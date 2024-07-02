import tensorflow as tf
import sys
import time
import numpy as np
from pyteomics import mgf as mgf
import math
import argparse
import os
from absl import logging

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
logging.set_verbosity(logging.ERROR)
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

Proton = 1.007276035
H2O = 18.0105647
B_ION_OFFSET = Proton
Y_ION_OFFSET = H2O + Proton
RESOLUTION = 1000 # Max = 1000
MAXPEPTMASS = 3000

# ## Traning Data Generation
BATCH_SIZE = 128
BUFFER_SIZE = 3000000
SPECTRA_MAX_LENGTH = 1010
PEPTIDE_MAX_LENGTH = 55

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

def preprocess(sequences, pmass):
    t_input = list(map(lambda sentence : list(''.join(sentence)), sequences))
    mass = 0.0
    r_mass = 0.0
    tm_input = []
    
    mass_int = int(pmass*RESOLUTION)
    tm_input.append([1, mass_int, 1, mass_int])
    for char in t_input :
        for f, r in zip(char, reversed(char)) :
            tmp = []
            mass += idx2mass.get(amino2idx.get(f))
            remass = pmass - mass
            mass_int = int(remass*RESOLUTION)
                
            tmp.append(amino2idx.get(f))
            tmp.append(mass_int)
    
            r_mass += idx2mass.get(amino2idx.get(r))
            remass = pmass - r_mass
            mass_int = int(remass*RESOLUTION)
          
            tmp.append(amino2idx.get(r))
            tmp.append(mass_int)
    
            tm_input.append(tmp)  
    
    remass = pmass - mass
    mass_int = int(remass*RESOLUTION)
  
    tmp = []
    tmp.append(2)
    tmp.append(mass_int)
    
    remass = pmass - r_mass
    mass_int = int(remass*RESOLUTION)
  
    tmp.append(2)
    tmp.append(mass_int)
    
    tm_input.append(tmp)
    return tm_input


def spec_data_generator(train_file) :
    trainfile = mgf.read(train_file)
    for spec in trainfile :      
        peaklist = []
        
        charge = int(spec['params']['charge'][0])
        pmz_ = spec['params']['pepmass'][0]
        pmass = (pmz_ - Proton)*charge - H2O;
        
        inten = spec['intensity array']
        mz = spec['m/z array']
        
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
        
        
        pept = []
        peptide = spec['params']['seq']
        peptide = peptide.replace("C(+57.02)", "C").replace("M(+15.99)", "m").replace("N(+.98)", "n").replace("Q(+.98)", "q")
        pept.append(peptide)

        pept = preprocess(pept, pmass)
        pept = np.reshape(pept, (-1, 4))
        
        pept_tar = pept[:-1, :]
        
        real = pept[:, :1]
        real1 = real[1:, :]
        real1 = np.squeeze(real1)
        real = pept[:, 2:3]
        real2 = real[1:, :]
        real2 = np.squeeze(real2)
        
        pept_tar = pept_tar.tolist()
        yield peaklist, pept_tar, real1, real2


def tf_encode(pt1, pt2, en1, en2):
    pt1.set_shape([None, 2])
    pt2.set_shape([None, 4])
    
    en1.set_shape([None])
    en2.set_shape([None])
    
    pt = {'inp' : pt1, 'dec_inp' : pt2}
    en = {'output_1' : en1, 'output_2' : en2}

    return pt, en

def test_data_generator(test_file) :   
    testfile = mgf.read(test_file)
    for spec in testfile :      
        peaklist = []
        
        charge = int(spec['params']['charge'][0])
        pmz_ = spec['params']['pepmass'][0]
        pmass = (pmz_ - Proton)*charge - H2O;
        
        inten = spec['intensity array']
        mz = spec['m/z array']
        
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
        
        
        pept = []
        peptide = spec['params']['seq']
        peptide = peptide.replace("C(+57.02)", "C").replace("M(+15.99)", "m").replace("N(+.98)", "n").replace("Q(+.98)", "q")
        pept.append(peptide)

        pept = preprocess(pept, pmass)
        pept = np.reshape(pept, (-1, 4))
        
        pept_tar = pept[:-1, :]
        
        real = pept[:, :1]
        real1 = real[1:, :]
        real1 = np.squeeze(real1)
        real = pept[:, 2:3]
        real2 = real[1:, :]
        real2 = np.squeeze(real2)
        
        pept_tar = pept_tar.tolist()
        yield peaklist, pept_tar, real1, real2

def filter_max_length(x, y, sepc_max_length=SPECTRA_MAX_LENGTH, pep_max_length=PEPTIDE_MAX_LENGTH):
    return tf.logical_and(tf.size(x['inp'][:, :1]) <= sepc_max_length,
                        tf.size(y['output_1'][:]) <= pep_max_length-1)

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
                spec_data_generator, output_types=(tf.int64, tf.int64, tf.int64, tf.int64), 
                output_shapes =( (None, 2), (None, 4), (None), (None) ),
                args=([args.learning_file]) )
    train_data = spec_data.map(tf_encode)
    train_data = train_data.filter(filter_max_length)
    train_data = train_data.cache()
    if not args.no_multigpu :
        train_data = train_data.with_options(options)
        train_data = train_data.shuffle(BUFFER_SIZE).padded_batch(GLOBAL_BATCH_SIZE)
    else :
        train_data = train_data.shuffle(sBUFFER_SIZE).padded_batch(BATCH_SIZE)
    train_data = train_data.prefetch(tf.data.experimental.AUTOTUNE)

    spec_data = tf.data.Dataset.from_generator(
                test_data_generator, output_types=(tf.int64, tf.int64, tf.int64, tf.int64),
                output_shapes =( (None, 2), (None, 4), (None), (None) ),
                args=([args.validation_file]) )
    test_data = spec_data.map(tf_encode)
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
    print("Loding Model... ", end='', flush = True)
    if not args.no_multigpu :
        with strategy.scope() :
            transformer = tf.keras.models.load_model(args.model_path, custom_objects=custom_objects)
    else :
        transformer = tf.keras.models.load_model(args.model_path, custom_objects=custom_objects)
    
    print("Done.", flush = True)

    print("Shuffling Spectra...", flush = True)
    if not args.no_multigpu :
        with strategy.scope() :
            transformer.fit(train_data, epochs=args.epochs, validation_data=test_data)
    else :
        transformer.fit(train_data, epochs=args.epochs, validation_data=test_data)

    print("Saving Model... ", end='', flush = True)
    if not args.save_weights is None :
        transformer.save_weights(args.save_weights)
    else :
        transformer.save(args.save_model)
    print("Done.", flush = True)



if __name__ == '__main__':    
    parser = argparse.ArgumentParser(description="NovoB")
    parser.add_argument('-m', '--model_path', required=True, help='Saved Model Path (folder)')
    parser.add_argument('--save_model', default='model', help='Path to save Model (folder), default = model')
    parser.add_argument('--save_weights', help='Path/file_name to save weights (path/file_name)')
    parser.add_argument('-l', '--learning_file', metavar='LEARNING_FILE.mgf', required=True, help='Learning File (.mgf)')
    parser.add_argument('-v', '--validation_file', required=True, metavar='VALIDATION_FILE.mgf', help='Validation File (.mgf)')
    parser.add_argument('-b', '--batch_size', default=128, type=int, help='Batch size for learning model (on single GPU), default=128')
    parser.add_argument('-e', '--epochs', default=30, type=int, help='the number of epochs, default=30')
    parser.add_argument('-g', '--no_multigpu', action='store_true', help='Do not use multigpu')
    parser.add_argument('-n', '--no_nccl', action='store_true', help='Do not use nccl (when using multigpu)')

    args = parser.parse_args()
    main(args)
