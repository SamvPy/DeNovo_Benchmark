import tensorflow as tf
from tensorflow.keras import mixed_precision
import sys
import time
import numpy as np
from pyteomics import mgf as mgf
import math
import argparse
import os
from absl import logging

# Use the tensor cores
mixed_precision.set_global_policy('mixed_float16')


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
logging.set_verbosity(logging.ERROR)
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

MEMORY_LIMIT = 8000 # MiB
gpus = tf.config.list_physical_devices('GPU')

if len(gpus) > 1:
    gpu_i = 1
else:
    gpu_i = 0
tf.config.set_visible_devices(gpus[gpu_i], 'GPU')
tf.config.set_logical_device_configuration(
    gpus[gpu_i], 
    [tf.config.LogicalDeviceConfiguration(memory_limit=MEMORY_LIMIT)])

Proton = 1.007276035
H2O = 18.0105647
B_ION_OFFSET = Proton
Y_ION_OFFSET = H2O + Proton
RESOLUTION = 1000 # Max = 1000
MAXPEPTMASS = 5000

# ## Traning Data Generation
BATCH_SIZE = 128
BUFFER_SIZE = 3500000
SPECTRA_MAX_LENGTH = 1010
PEPTIDE_MAX_LENGTH = 50

amino2idx = {'<pad>': 0, '<bos>': 1, '<eos>': 2,
             'A': 3, 'C': 4, 'D': 5, 'E': 6, 'F': 7,
             'G': 8, 'H': 9, 'I': 10, 'K': 11, 'L': 12,
             'M': 13, 'N': 14, 'P': 15, 'Q': 16, 'R': 17,
             'S': 18, 'T': 19, 'V': 20, 'W': 21, 'Y': 22,
             'm': 23, 'n': 24, 'q': 25, 's': 26, 't': 27, 'y': 28}

idx2amino = {0: '<pad>',1: '<bos>', 2: '<eos>',
             3: 'A', 4: 'C', 5: 'D', 6: 'E', 7: 'F',
             8: 'G', 9: 'H', 10: 'I', 11: 'K', 12: 'L',
             13: 'M', 14: 'N', 15: 'P', 16: 'Q', 17: 'R',
             18: 'S', 19: 'T', 20: 'V', 21: 'W', 22: 'Y',
             23: 'm', 24: 'n', 25: 'q', 26: 's', 27: 't', 28: 'y'}


# High: D = n, E = q, L = I
# Low: K = Q (0.4), F = m (0.3)
idx2mass = {0: 0.00, 1: 0.00, 2: 0.00,
             3: 71.037, 4: 160.031, 5: 115.027, 6: 129.043, 7: 147.068,
             8: 57.021, 9: 137.059, 10: 113.084, 11: 128.095, 12: 113.084,
             13: 131.040, 14: 114.043, 15: 97.053, 16: 128.059, 17: 156.101,
             18: 87.032, 19: 101.048, 20: 99.068, 21: 186.079, 22: 163.063,
             23: 147.035, 24: 115.027, 25: 129.043, 26: 166.998, 27: 181.014, 28:243.030}

def tf_encode(pt):
    pt['inp'].set_shape([None, 2])
    return pt


def test_data_generator(test_file) :
    testfile = mgf.read(test_file)
    for spec in testfile :      
        peaklist = []
        
        charge = int(spec['params']['charge'][0])
        pmz_ = spec['params']['pepmass'][0]
        pmass = (pmz_ - Proton)*charge - H2O;
        
        inten = spec['intensity array']
        mz = spec['m/z array']
        
        scan_id = spec["params"]["scans"]

        if pmass > MAXPEPTMASS or len(mz) > 1000 :
            continue

        bion = pmass - H2O + B_ION_OFFSET;
        yion = pmass - H2O + Y_ION_OFFSET;
        
        pmass = round(pmass, 3)
        bion = round(bion, 3)
        yion = round(yion, 3)

        pmass_int = int(float(pmass)*RESOLUTION)
        peaklist.append([pmass_int, charge])
        
        #start
        peaklist.append([0, 0])
        
        peak_int = int(float(B_ION_OFFSET)*RESOLUTION)
        peaklist.append([peak_int, 101])
        
        peak_int = int(float(Y_ION_OFFSET)*RESOLUTION)
        peaklist.append([peak_int, 102])
        
        maxInten = 0.0
        for i in range(len(inten)) :
            inten[i] = math.sqrt(inten[i])
            maxInten = max(inten[i], maxInten)
        
        for i in range(len(mz)) :
            if mz[i] > pmass :
                continue
            peak_int = int(float(mz[i])*RESOLUTION)
            peaklist.append([peak_int, int(inten[i]/maxInten*100)])

        
        peak_int = int(float(bion)*RESOLUTION)
        peaklist.append([peak_int, 101])
        
        peak_int = int(float(yion)*RESOLUTION)
        peaklist.append([peak_int, 102])

        yield {'inp': peaklist, 'scan_id': scan_id}

def filter_max_length(x, sepc_max_length=SPECTRA_MAX_LENGTH):
    return tf.logical_and(tf.size(x['inp'][:, :1]) <= sepc_max_length, True)


# ## Optimizer
class CustomSchedule(tf.keras.optimizers.schedules.LearningRateSchedule):
    def __init__(self, d_model, warmup_steps):
        super(CustomSchedule, self).__init__()
    
        self.d_model = d_model
        self.warmup_steps = warmup_steps
    
    def __call__(self, step):
        d_model = tf.cast(self.d_model, tf.float32)
        arg1 = tf.math.rsqrt(step)
        arg2 = step * (self.warmup_steps ** -1.5)
        return tf.math.rsqrt(d_model) * tf.math.minimum(arg1, arg2)
    
    def get_config(self):
        config = {
        'd_model': self.d_model,
        'warmup_steps': self.warmup_steps,}
        return config

# ## Loss and metrics
loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True, reduction='none')

def loss_function(real, pred):
    mask = tf.math.logical_not(tf.math.equal(real, 0))
    loss_ = loss_object(real, pred)

    mask = tf.cast(mask, dtype=loss_.dtype)
    loss_ *= mask

    return tf.reduce_sum(loss_)/tf.reduce_sum(mask)


def accuracy_function(real, pred):
    accuracies = tf.equal(tf.cast(real, dtype=tf.int64), tf.argmax(pred, axis=2))

    mask = tf.math.logical_not(tf.math.equal(real, 0))
    accuracies = tf.math.logical_and(mask, accuracies)

    accuracies = tf.cast(accuracies, dtype=tf.float32)
    mask = tf.cast(mask, dtype=tf.float32)
    return tf.reduce_sum(accuracies)/tf.reduce_sum(mask)

def predict(spec, charge_, pmass, transformer, batch):
    encoder_input = tf.cast(spec, tf.int64)
    decoder_input = []
    
    check = []
    score = []
    smass = []
    count = 0
    
    smass.append([])
    smass.append([])
    check.append([])
    check.append([])
    score.append([])
    score.append([])
    
    for i in range(len(charge_)) : 
        pmass_int = int(pmass[i]*RESOLUTION)
    
        decoder_input.append([[1, pmass_int, 1, pmass_int]])
        check[0].append(False)
        check[1].append(False)
        smass[0].append(0)
        smass[1].append(0)
        score[0].append([])
        score[1].append([])
        
    decoder_input = tf.cast(decoder_input, tf.int64)
    inp = {'inp': encoder_input, 'dec_inp': decoder_input}
    for idx in range(PEPTIDE_MAX_LENGTH):
        predictions0, predictions1 = transformer.predict(inp, batch_size=batch)
        
        predictions0 = predictions0[: ,-1:, :]  # (batch_size, 1, vocab_size)
        predicted_id0 = tf.cast(tf.argmax(predictions0, axis=-1), tf.int32)
        softmax0 = tf.nn.softmax(predictions0)
        softmax0 = tf.math.reduce_max(softmax0, axis=-1).numpy()
        
        #reverse
        predictions1 = predictions1[: ,-1:, :]  # (batch_size, 1, vocab_size)
        predicted_id1 = tf.cast(tf.argmax(predictions1, axis=-1), tf.int32)
        softmax1 = tf.nn.softmax(predictions1)
        softmax1 = tf.math.reduce_max(softmax1, axis=-1).numpy()
        
        pred_ = []
        
        for i in range(len(predicted_id0)) :
            pred = []
            
            r = predicted_id0[i].numpy()[0]
            if r == 2 and not check[0][i] :
                count += 1
                check[0][i] = True
            
            smass[0][i] += idx2mass.get(r)
            rmass = pmass[i] - smass[0][i]
            
            mass_int = int(rmass*RESOLUTION)

            pred.append(r)
            pred.append(mass_int)

            score[0][i].append(softmax0[i][0])
            
            #reverse
            r = predicted_id1[i].numpy()[0]
            if r == 2 and not check[1][i] :
                count += 1
                check[1][i] = True
            
            smass[1][i] += idx2mass.get(r)
            rmass = pmass[i] - smass[1][i]
            
            mass_int = int(rmass*RESOLUTION)

            pred.append(r)
            pred.append(mass_int)
            
            pred_.append([pred])
            
            score[1][i].append(softmax1[i][0])
        
        if count == (len(predicted_id0)*2) :
            return inp['dec_inp'], score
        
        inp['dec_inp'] = tf.concat([inp['dec_inp'], pred_], axis=1)
              
        
    return inp['dec_inp'], score


def translate(result):
    result_ = result[:, 0:1]
    result_ = tf.squeeze(result_, axis=-1)
    seq = list(result_.numpy())
    predicted_sentence = ''
    for j in seq:
        if(j == 2) :
            break
        if(j == 0 or j == 1) :
            continue
        predicted_sentence += idx2amino.get(j)
        
    senten = ""
    for j in range(len(predicted_sentence)):
        senten += predicted_sentence[j]

    length = len(senten)
    senten = "['" + senten + "']"

    mass_int = result[:, 1:2]
    mass_int = tf.squeeze(mass_int, axis=-1)

    mass = float(mass_int.numpy()[length])/RESOLUTION
    
    #reverse
    result_ = result[:, 2:3]
    result_ = tf.squeeze(result_, axis=-1)
    seq = list(result_.numpy())
    predicted_sentence = ''
    for j in seq:
        if(j == 2) :
            break
        if(j == 0 or j == 1) :
            continue
        predicted_sentence += idx2amino.get(j)
        
    senten_R = ""
    for j in range(len(predicted_sentence)):
        senten_R += predicted_sentence[j]

    length = len(senten_R)
    senten_R = "['" + senten_R[::-1] + "']"

    mass_int = result[:, 3:4]
    mass_int = tf.squeeze(mass_int, axis=-1)

    mass_R = float(mass_int.numpy()[length])/RESOLUTION

    return senten, mass, senten_R, mass_R

def main(args) :
    BATCH_SIZE = args.batch_size

    if not args.no_multigpu :
        if args.no_nccl :
            strategy = tf.distribute.MirroredStrategy(cross_device_ops=tf.distribute.HierarchicalCopyAllReduce())
        else :
            strategy = tf.distribute.MirroredStrategy()
        BATCH_SIZE_PER_REPLICA = BATCH_SIZE
        GLOBAL_BATCH_SIZE = BATCH_SIZE_PER_REPLICA * strategy.num_replicas_in_sync
        options = tf.data.Options()
        options.experimental_distribute.auto_shard_policy = tf.data.experimental.AutoShardPolicy.DATA

    spec_data = tf.data.Dataset.from_generator(
        test_data_generator, 
        output_types={'inp': tf.int64, 'scan_id': tf.string},
        output_shapes={'inp': (None, 2), 'scan_id': ()}, 
        args=([args.spectrum_file])
    )

    test_data = spec_data.map(tf_encode, num_parallel_calls=tf.data.experimental.AUTOTUNE)
    test_data = test_data.filter(filter_max_length)
    test_data = test_data.cache()
    if not args.no_multigpu :
        test_data = test_data.with_options(options)
        test_data = test_data.padded_batch(GLOBAL_BATCH_SIZE)
    else :
        test_data = test_data.padded_batch(BATCH_SIZE)
    test_data = test_data.prefetch(tf.data.experimental.AUTOTUNE)

    custom_objects = {'loss_function' : loss_function, 'accuracy_function' : accuracy_function,
                 'CustomSchedule' : CustomSchedule}
    
    print("Loading Model... ", end='', flush = True)
    if not args.no_multigpu :
        with strategy.scope() :
            transformer = tf.keras.models.load_model(args.model_path, custom_objects=custom_objects)
            if not args.load_weights is None :
                print("Load Variables ... ", end='', flush = True)
                transformer.load_weights(args.load_weights)
    else :
        transformer = tf.keras.models.load_model(args.model_path, custom_objects=custom_objects)
        if not args.load_weights is None :
            print("Load Variables ... ", end='', flush = True)
            transformer.load_weights(args.load_weights)

    print("Done.", flush = True)

    total = 0
    for batch, (inp) in enumerate(test_data) :
        total += len(inp['inp'])

    resultfile = open(args.output_file, 'w')

    Mcount = 0
    index = 0
    end = 0.0
    for batch, (inp) in enumerate(test_data) :
        index = index + 1
        
        print("Loading Spectra... ", end='', flush = True)
        charge_= []
        pepmass= []
        for i in range(len(inp['inp'])) :
            charge_.append(inp['inp'][i][0][1].numpy())
            mass_int = inp['inp'][i][0][0]
            pepmass.append(float(mass_int.numpy())/RESOLUTION)
        spec = inp['inp']
        scan_ids = inp['scan_id']
        print("{}/{} Done.".format(Mcount+len(inp['inp']), total), flush = True)


        start = time.time()
        print('Procesing {}/{} ...'.format(Mcount+len(inp['inp']), total), flush = True)
        if not args.no_multigpu :
            with strategy.scope() :
                result, score = predict(spec, charge_, pepmass, transformer, GLOBAL_BATCH_SIZE)
        else :
            result, score = predict(spec, charge_, pepmass, transformer, BATCH_SIZE)
        end += (time.time()-start)
        print('Procesing {}/{} ... ({:.4f} Sec)'.format(Mcount+len(inp['inp']), total, end), flush = True)
        
        print('Writing Results ...', end='', flush = True)
        for idx in range(len(spec)) :
            Mcount = Mcount + 1
            senten, mass, senten_R, mass_R = translate(result[idx])
            prob = 1
            rprob = 1
            for i in score[0][idx] :
                prob *= i
            for i in score[1][idx] :
                rprob *= i
            with open(args.output_file, 'a') as resultfile :
                resultfile.write('{}\t{}\t{:.3f}\t{}\t{:.3f}\t{:.6f}\t{}\t{:.3f}\t{:.6f}\t{}\n'.format(Mcount, charge_[idx], pepmass[idx], senten, mass, prob, senten_R, mass_R, rprob, scan_ids[idx]))
        print("{}/{} Done.".format(Mcount, total), flush = True)
    print('Procesing {}/{} ... Done. ({:.4f} Sec)'.format(Mcount, total, end), flush = True)

if __name__ == '__main__':    
    parser = argparse.ArgumentParser(description="NovoB")
    parser.add_argument('-m', '--model_path', required=True, help='Saved Model Path (folder)')
    parser.add_argument('--load_weights', help='Load weights after loading model')
    parser.add_argument('-i', '--spectrum_file', required=True, help='spectrum file (.mgf)')
    parser.add_argument('-o', '--output_file', default='result.txt', help='output file (text file), default=result.txt')
    parser.add_argument('-b', '--batch_size', default=256, type=int, help='Batch size to predict peptides (on single GPU), default=256')
    parser.add_argument('-g', '--no_multigpu', action='store_true', help='Do not use multigpu')
    parser.add_argument('-n', '--no_nccl', action='store_true', help='Do not use nccl (when using multigpu)')
    args = parser.parse_args()
    main(args)
