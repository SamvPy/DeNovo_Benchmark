"""
Get InstaNovo predictions, both from the original 5-beam search method and log probabilities
from the annotated peptides found in a sage search.
"""
from __future__ import annotations

import argparse
import logging
import os
import time
from typing import Any

import polars as pl
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm

from psm_utils.io import read_file

from instanovo.constants import MASS_SCALE
from instanovo.inference.knapsack import Knapsack
from instanovo.inference.knapsack_beam_search import KnapsackBeamSearchDecoder
from instanovo.transformer.dataset import collate_batch
from instanovo.transformer.dataset import SpectrumDataset
from instanovo.transformer.model import InstaNovo
from instanovo.transformer.predict import _setup_knapsack
from instanovo.transformer.model import InstaNovo
import argparse

logging.basicConfig(filename='run_denovo_instanovo.log', level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="Path to the input file containing filenames")
    parser.add_argument("ipc_dir")
    parser.add_argument("identifications_dir")
    parser.add_argument("output_dir")
    parser.add_argument("--n_workers", default=8)
    parser.add_argument("--batch_size", default=256)
    args = parser.parse_args()

    # Read input file
    with open(args.input_file, 'r') as f:
        filenames = f.read().splitlines()

    # Define paths
    ipc_dir = args.ipc_dir
    id_dir = args.identifications_dir
    output_dir = args.output_dir

    model, config = InstaNovo.load(
        "/home/sam/denovo_project/model_checkpoints/instanovo/base/instanovo_yeast.pt"
    )

    config["n_workers"] = int(args.n_workers)
    config["subset"] = 1
    config["predict_batch_size"] = int(args.batch_size)
    denovo =True
    knapsack_path = "/home/sam/denovo_project/model_checkpoints/instanovo/knapsack/"

    # Loop over filenames and run script for each pair of files
    for filename in filenames:
        # Construct paths for feather files, .sage.tsv files, and output
        ipc_path = os.path.join(ipc_dir, filename + ".ipc")
        id_path = os.path.join(id_dir, filename + ".sage.tsv")
        output_path = os.path.join(output_dir, filename + ".feather")
    
        get_top_preds(ipc_path, id_path, model, config, output_path, knapsack_path)

def get_top_preds(
    data_path: str,
    id_path: str,
    model: InstaNovo,
    config: dict[str, Any],
    output_path: str | None = None,
    knapsack_path: str | None = None,
    device: str = "cpu",
):
    logging.info(f"Loading data from {data_path}")
    df = pl.read_ipc(data_path)
    df = df.sample(fraction=config["subset"], seed=0)
    logging.info(
        f"Data loaded, evaluating {config['subset']*100:.1f}%, {df.shape[0]} samples in total."
    )

    df = df.to_pandas()
    identifications = read_file(id_path, filetype="sage")
    identifications = identifications.get_rank1_psms()
    identifications = identifications.to_dataframe()
    identifications["scan_number"] = identifications.spectrum_id.apply(lambda x: int(x[5:]))

    df = df[df.scan_number.isin(identifications.scan_number)]
    logging.info(f"Shape df: {df.shape}")

    vocab = list(config["residues"].keys())
    config["vocab"] = vocab
    s2i = {v: i for i, v in enumerate(vocab)}
    i2s = {i: v for i, v in enumerate(vocab)}

    ds = SpectrumDataset(df, s2i, config["n_peaks"], return_str=True)

    dl = DataLoader(
        ds,
        batch_size=config["predict_batch_size"],
        num_workers=config["n_workers"],
        shuffle=False,
        collate_fn=collate_batch,
    )

    model = model.to(device)
    model = model.eval()

    # setup decoder
    if knapsack_path is None or not os.path.exists(knapsack_path):
        logging.info("Knapsack path missing or not specified, generating...")
        knapsack = _setup_knapsack(model)
        decoder = KnapsackBeamSearchDecoder(model, knapsack)
        if knapsack_path is not None:
            logging.info(f"Saving knapsack to {knapsack_path}")
            knapsack.save(knapsack_path)
    else:
        logging.info("Knapsack path found. Loading...")
        decoder = KnapsackBeamSearchDecoder.from_file(model=model, path=knapsack_path)
    
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

    pred_df = df[cols].copy()

    preds = []
    targs = []
    probs = []

    start = time.time()
    for _, batch in tqdm(enumerate(dl), total=len(dl)):
        spectra, precursors, spectra_mask, peptides, _ = batch
        spectra = spectra.to(device)
        precursors = precursors.to(device)
        spectra_mask = spectra_mask.to(device)

        with torch.no_grad():
            p = decoder.decode(
                spectra=spectra,
                precursors=precursors,
                beam_size=3,
                max_length=config["max_length"],
                return_all_beams=True
            )

            preds += [["".join(i.sequence) for i in x] for x in p]
            probs += [[i.log_probability for i in x] for x in p]
            # targs += list(peptides)

    delta = time.time() - start
    logging.info(f"Time taken for {data_path} is {delta:.1f} seconds")
    logging.info(
        f"Average time per batch (bs={config['predict_batch_size']}): {delta/len(dl):.1f} seconds"
    )

    # pred_df["targets"] = targs
    pred_df["preds"] = preds
    pred_df["log_probs"] = probs

    pred_df.reset_index().to_feather(output_path)

def _setup_knapsack(model: InstaNovo) -> Knapsack:
    residue_masses = model.peptide_mass_calculator.masses
    residue_masses["$"] = 0
    residue_indices = model.decoder._aa2idx
    return Knapsack.construct_knapsack(
        residue_masses=residue_masses,
        residue_indices=residue_indices,
        max_mass=4000.00,
        mass_scale=MASS_SCALE,
    )

if __name__ == "__main__":
    main()