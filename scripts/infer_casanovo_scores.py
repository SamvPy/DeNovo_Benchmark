# Should be run with casanovo virtual environment
# denovo_utils should also be installed


'''Infer prediction scores from a casanovo model checkpoint given a psmlist.'''

import pandas as pd
import torch.nn as nn
import numpy as np
from psm_utils.io import read_file, write_file
from psm_utils import PSMList
import torch
from torch.utils.data import DataLoader, TensorDataset, Dataset
import json
import argparse

from pathlib import Path
from pyteomics import mgf
from pyteomics.mzml import PreIndexedMzML
from tqdm import tqdm
import os

from casanovo.casanovo import setup_logging, setup_model, ModelRunner
from casanovo.denovo.model import Spec2Pep
from denovo_utils.parsers.converters import DenovoEngineConverter
from denovo_utils.parsers.constants import MODIFICATION_MAPPING
import logging

logger = logging.getLogger(__name__)

def peptidoform_converter(
    peptidoform_in,
    peptidoform_out,
    peptidoform,
):
    proforma = peptidoform.proforma

    if isinstance(peptidoform_in, str):
        mod_in = MODIFICATION_MAPPING[peptidoform_in].values()
    elif isinstance(peptidoform_in, dict):
        mod_in = list(peptidoform_in.values())
        for mod_in_format, mod_in_reformat in peptidoform_in.items():
            proforma = proforma.replace(mod_in_format, mod_in_reformat)
    else:
        raise Exception('peptidoform_in should be mod_mapper or str')
    if isinstance(peptidoform_out, str):
        mod_out = MODIFICATION_MAPPING[peptidoform_out]
    elif isinstance(peptidoform_out, dict):
        mod_out = peptidoform_out
    else:
        raise Exception('peptidoform_out should be mod_mapper or str')

    if len(mod_in) > 0:
        unsupported_mods = [mod for mod in mod_in if mod not in mod_out.values()]
        for unsupported_mod in unsupported_mods:
            if unsupported_mod in proforma:
                return None
        
        for mod_out_reformat, mod_out_proforma in mod_out.items():
            proforma = proforma.replace(mod_out_proforma, mod_out_reformat)
    return proforma.split('/')[0]

def get_specid_index(path_peak):
    if path_peak.lower().endswith('mgf'):
        spectrum_ids = pd.DataFrame(
            pd.DataFrame(
                mgf.read(path_peak)
            )["params"].tolist()
        )["title"]
        scan_spectrum = {}

        # TODO: Check if this is correct
        for i, spectrum_id in enumerate(spectrum_ids.values):
            scan_spectrum[spectrum_id] = f"index={i}"
        return scan_spectrum
    
    elif path_peak.lower().endswith('mzml'):
        mzml_file = PreIndexedMzML(path_peak)
        scan_spectrum = {}
        for i, spectrum in tqdm(enumerate(mzml_file)):
            # The scan in casanovo equals
            # the scan in the id header of the mzml spectrum
            scan_spectrum[spectrum['id']] = spectrum['id'].split()[-1]
        return scan_spectrum
    
    else:
        raise Exception(f"path_peak must be mzML or MGF and filename should end accordingly. Path peak: {path_peak}")


def load_results(
    path_peak: str,
    psmlist=None,
    path_results=None,
    engine=None
) -> dict:
    
    if psmlist is None:
        parser = DenovoEngineConverter.select(engine)
        psmlist = parser.parse(
            result_path=path_results,
            mgf_path=path_peak
        )
    
    # Get mapping of casanovo based indexation and spectrum_id
    specid_index = get_specid_index(path_peak)

    id_peptide_dict = {
        specid_index[psm.spectrum_id]: peptidoform_converter(
            peptidoform_in=engine,
            peptidoform_out="casanovo",
            peptidoform=psm.peptidoform
        )
        for psm in psmlist
    }
    return id_peptide_dict, specid_index

class InterpretableCasanovo(nn.Module):
    def __init__(self, model: Spec2Pep):
        super().__init__()
        self.encoder = model.encoder
        self.decoder = model.decoder
        self.softmax = model.softmax
        self.max_length = model.max_length
        self.n_beams = model.n_beams
        self.model = model

    def forward(self, spectra, precursors, peptide=None, apply_softmax=False):
        
        memories, mem_masks = self.encoder(spectra)

        if peptide is None:
            pred, _ = self.decoder(None, precursors, memories, mem_masks)
        else:
            pred, _ = self.decoder(peptide, precursors, memories, mem_masks)

        if apply_softmax:
            pred = self.softmax(pred)
        return pred
    
def load_model_data(
    path_checkpoint=Path,
    path_config=Path,
    path_output=Path,
    path_peak=list[Path],
    indices=list
):
    output = setup_logging(
        output=path_output,
        verbosity="info"
    )
    config, model = setup_model(
        model=path_checkpoint,
        config=path_config,
        output=output,
        is_train=False
    )
    indices_out = []
    indices_all = []

    with ModelRunner(config, model) as runner:
        runner.initialize_trainer(train=False)
        runner.initialize_model(train=False)
        model_pt = runner.model
        test_index = runner._get_index(path_peak, False, "")
        runner.initialize_trainer(train=False)
        runner.initialize_data_module(test_index=test_index)
        runner.loaders.setup(stage="test", annotated=False)

        dl = runner.loaders._make_loader(
            runner.loaders.test_dataset,
            1
        )
        datapoints = []
        for entry in dl:
            index = entry[2][0][1]
            indices_all.append(index)
            if index in indices:
                indices_out.append(index)
                datapoints.append(entry)

    return InterpretableCasanovo(model_pt), datapoints, {
        'indices_out': indices_out, "indices_all": indices_all,
        'test_index': test_index, 'indices_in': indices
    }


def extract_scores(
    model,
    datapoint,
    peptide
):
    spec, prec, meta = datapoint

    scores = model.forward(
        spectra=spec,
        precursors=prec,
        peptide=peptide,
        apply_softmax=True
    ).detach().numpy()

    scores = scores[-1]

    tokenized = model.decoder.tokenize(peptide).detach().numpy()

    score_aa = []
    for i, aa_i in enumerate(tokenized):
        score_aa.append(scores[i][aa_i])

    peptide_score = np.mean(score_aa)
    aa_scores = (score_aa + peptide_score) / 2
    aa_scores = aa_scores[:-1]
    return peptide_score, aa_scores


def extract_scores_batch(
        model,
        spectra,
        precursors,
        peptides
    ):

    # Batch forward pass through the model
    scores_batch = model.forward(
        spectra=spectra,
        precursors=precursors,
        peptide=peptides,  # Assuming model can handle a batch of peptides
        apply_softmax=True
    ).detach().numpy()
   
    # Tokenize the entire batch of peptides at once (if supported)
    tokenized_batch = [model.decoder.tokenize(s).detach().numpy() for s in peptides]

    # Initialize lists to store peptide and amino acid scores
    peptide_scores = []
    aa_scores_list = []

    # Compute amino acid scores for all peptides in the batch using matrix operations
    for i in range(len(tokenized_batch)):
        tokenized = tokenized_batch[i]  # Tokenized amino acids for the i-th peptide
        scores = scores_batch[i]  # Corresponding score matrix for the i-th peptide

        # Retrieve scores for the tokenized amino acids
        score_aa = scores[np.arange(len(tokenized)), tokenized]  # Shape: (seq_len,)
        
        # Compute peptide score as the mean of amino acid scores
        peptide_score = np.mean(score_aa)
        peptide_scores.append(peptide_score)

        # Adjust amino acid scores
        aa_scores = (score_aa + peptide_score) / 2
        aa_scores = aa_scores[:-1]  # Remove the last token's score
        aa_scores_list.append(aa_scores)

    return peptide_scores, aa_scores_list


class CasanovoDataset(Dataset):
    def __init__(self, data, peptides):
        self.data = data
        self.peptides = peptides
        assert len(data) == len(peptides)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx], self.peptides[idx]

def prepare_batch(
    batch
):
    b, peptides = list(zip(*batch))
    spectra, precursors, spectrum_ids = list(zip(*b))
    spectra = [x.reshape(-1, 2) for x in spectra]
    precursors = [x.flatten() for x in precursors]
    spectra = torch.nn.utils.rnn.pad_sequence(spectra, batch_first=True)
    spectrum_ids = [x[0][1] for x in spectrum_ids]
    
    precursors = torch.stack(precursors)
    return spectra, precursors, peptides, spectrum_ids


if __name__ == "__main__":
    print("WARNING: For now only supports PSMLists with a single entry per spectrum.")
    parser = argparse.ArgumentParser(description="Rescore a psmlist with the casanovo model.")
    parser.add_argument('--mzml', required=True, help='Path to peak file.')
    parser.add_argument('--psmlist', required=True, help='Path to psmlist to rescore.')
    parser.add_argument('--path_model', required=True, help='Path to configuration file.')
    parser.add_argument('--path_config', required=True, help='Path to casanovo configuration file.')
    parser.add_argument('--modification_mapping', required=True, help='Provide either mapping dictionary or the engine name to map modifications. Alternatively, provide empty dictionary if mappings are already done.')
    parser.add_argument('--psmlist_filetype', required=False, default='infer', help='Reading setting for psmlist (see psm-utils filetype).')
    args = parser.parse_args()

    psm_list = read_file(args.psmlist, filetype=args.psmlist_filetype)
    with open(args.modification_mapping, 'r') as f:
        modification_mapping = json.load(f)

    # Load the results and create an index mapper
    id_peptide_dict, specid_index = load_results(
        path_peak=args.mzml,
        psmlist=psm_list,
        engine=modification_mapping
    )
    model, data, indices = load_model_data(
        path_checkpoint=args.path_model,
        path_config=args.path_config,
        path_output=Path('.'),
        path_peak=[args.mzml],
        indices=list(id_peptide_dict.keys())
    )

    index_specid = {v:k for k, v in specid_index.items()}

    # Extract the annotated peptides for each datapoint
    peptide_list = [
        id_peptide_dict[batch[-1][0][-1]] for batch in data
    ]

    dataset = CasanovoDataset(data, peptides=peptide_list)
    dl = DataLoader(
        dataset,
        batch_size=8,
        collate_fn=prepare_batch
    )

    result_dict_batched = {}

    for batch in tqdm(dl, desc="Inferring casanovo scores..."):
        pep_scores, aa_scores = extract_scores_batch(
            model=model,
            spectra=batch[0],
            precursors=batch[1],
            peptides=batch[2]
        )
        for pep_score, aa_score, scan in zip(pep_scores, aa_scores, batch[3]):
            result_dict_batched[index_specid[scan]] = {
                "peptide_score": pep_score,
                "aa_score": aa_score
            }
    
    for psm in psm_list:
        try:
            casanovo_score = result_dict_batched[psm['spectrum_id']]
        except KeyError:
            casanovo_score = np.nan
            logger.info(f"No score outputted for {psm['spectrum_id']}")

        psm['rescoring_features'].update({"Casanovo_score": casanovo_score['peptide_score']})
    
    write_file(
        psm_list=psm_list,
        filename="results.parquet",
        filetype='parquet',
        show_progressbar=True
    )