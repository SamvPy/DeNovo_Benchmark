"""
Get InstaNovo predictions, both from the original 5-beam search method and log probabilities
from the annotated peptides found in a sage search.
"""
from __future__ import annotations

import argparse
import logging
import os
import time
from pathlib import Path
from typing import Any

import pandas as pd
import numpy as np
import polars as pl
import torch
import yaml
from torch.utils.data import DataLoader
from tqdm import tqdm
import re
from psm_utils.io import read_file

from instanovo.constants import MASS_SCALE
from instanovo.inference.knapsack import Knapsack
from instanovo.inference.knapsack_beam_search import KnapsackBeamSearchDecoder
from instanovo.transformer.dataset import collate_batch
from instanovo.transformer.dataset import SpectrumDataset
from instanovo.transformer.model import InstaNovo
from instanovo.utils.metrics import Metrics
from instanovo.transformer.predict import _setup_knapsack, get_preds
from custom_scripts import PeptideSpecificDecoder
from instanovo.transformer.model import InstaNovo
import argparse

logging.basicConfig(filename='get_instanovo_predictions.log', level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Path to the input file containing filenames")
    parser.add_argument("ipc_dir")
    parser.add_argument("identifications_dir")
    parser.add_argument("output_dir")
    args = parser.parse_args()

    # Read input file
    with open(args.input_file, 'r') as f:
        filenames = f.read().splitlines()

    # Define paths
    ipc_dir = args.ipc_dir
    id_dir = args.identifications_dir
    output_dir = args.output_dir

    args = parser.parse_args()

    model, config = InstaNovo.load(
        "/home/sam/denovo_project/model_checkpoints/instanovo/base/instanovo_yeast.pt"
    )

    config["n_workers"] = 8
    config["subset"] = 1
    denovo = True
    knapsack_path = "/home/samva/Doctorate_project/Denovo_sideproject/checkpoints/instanovo_checkpoint/knapsack/"

    # Loop over filenames and run script for each pair of files
    for filename in filenames:
        # Construct paths for feather files, .sage.tsv files, and output
        ipc_path = os.path.join(ipc_dir, filename + ".ipc")
        id_path = os.path.join(id_dir, filename + ".sage.tsv")
        output_path = os.path.join(output_dir, filename + ".csv")
    
        get_matching_scores(model, config, ipc_path, id_path, output_path)


def get_matching_scores(model, config, ipc_path, id_path, output_path,
                        batch_size=128):
    preds = []
    targs = []
    probs = []
    prob_aa = []

    device = "cuda" if torch.cuda.is_available() else "cpu"
    device
    model = model.to(device).eval()
    decoder = PeptideSpecificDecoder(model)

    # Load and parse in the data
    # The identifications
    identifications = read_file(id_path, filetype="sage")
    identifications = identifications.get_rank1_psms()
    identifications = identifications.to_dataframe()

    identifications["scan_number"] = identifications.spectrum_id.apply(
        lambda x: int(re.search(r'(?<=scan=)\d+', x).group())
    )

    # The ipc dataset for instanovo modelling
    df = pl.read_ipc(ipc_path).to_pandas()
    index_cols = [
        "experiment_name",
        "scan_number",
        "id",
        "global_index",
        "spectrum_index",
        "file_index",
        "sample",
        "file",
        "index",
        "fileno",
    ]
    cols = [x for x in df.columns if x in index_cols]

    # Merge the datasets and only keep matching scans
    identifications_ipc_merge = pd.merge(
        identifications[
            ["peptidoform", "scan_number", "score", "qvalue", "is_decoy", "rescoring_features"]
        ], df, on="scan_number"
    )
    identifications_ipc_merge["modified_sequence"] = [str(x)[:-2] for x in identifications_ipc_merge.peptidoform.tolist()]
    identifications_ipc_merge["modified_sequence"] = identifications_ipc_merge.modified_sequence.apply(
        lambda x: x.replace(
            "[+57.0215]", "(+57.02)"
        ).replace(
            "[+15.9949]", "(+15.99)"
        ).replace(
            "[+42.010567]-", ""
        ).replace(
            "I", "L" # Add this in second run
        ).replace(
            "U", "C(+57.02)"
        )
    )

    pred_df = identifications_ipc_merge[cols].copy()

    logging.info(f"Shape identifications: {identifications.shape}; shape ipc: {df.shape}; Merged shape: {identifications_ipc_merge.shape}")

    s2i = model.decoder._aa2idx
    ds = SpectrumDataset(identifications_ipc_merge[df.columns], s2i, config["n_peaks"], return_str=False, eos_symbol="$")
    dl = DataLoader(ds, batch_size=batch_size, shuffle=False, collate_fn=collate_batch)

    start = time.time()
    for _, batch in tqdm(enumerate(dl), total=len(dl)):
        spectra, precursors, spectra_mask, peptides, _ = batch

        # Fix encoding issue
        peptides = peptides
    
        spectra = spectra.to(device)
        precursors = precursors.to(device)
        spectra_mask = spectra_mask.to(device)
        peptides = peptides.to(device)

        with torch.no_grad():
            p = decoder.decode(
                spectra=spectra,
                precursors=precursors,
                spectra_mask=spectra_mask,
                peptides_tgt=peptides,
            )

            preds += ["".join(model.decode(x))[1:] for x in p.sequence]
            probs += [x.detach().cpu() for x in p.log_probability]
            prob_aa += [x.detach().cpu() for x in p.log_probabilities]
            targs += ["".join(model.decode(x))[1:] for x in peptides]

    delta = time.time() - start
    logging.info(
        f"Average time per batch (bs={config['predict_batch_size']}): {delta/len(dl):.1f} seconds"
    )

    pred_df["preds"] = preds
    pred_df["log_probs"] = probs
    pred_df["aa_probs"] = prob_aa
    pred_df["targs"] = targs

    pred_df.to_csv(output_path[:-4]+"_targeted.csv", index=False)

if __name__ == "__main__":
    main()
