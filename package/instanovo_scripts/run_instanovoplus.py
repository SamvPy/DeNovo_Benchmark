# Environment: instanovo

import pandas as pd
import polars as pl
import argparse
import torch
import pyarrow as pa
import pyarrow.dataset as ds
import json
import os

from datasets import Dataset
from tqdm import tqdm
from torch.utils.data import DataLoader
from instanovo.diffusion.dataset import AnnotatedPolarsSpectrumDataset
from instanovo.diffusion.dataset import collate_batches
from instanovo.diffusion.multinomial_diffusion import MultinomialDiffusion
from instanovo.inference.diffusion import DiffusionDecoder


def pddf_to_dl(df: pd.DataFrame, model, decoder):

    peptides = df["modified_sequence"].tolist()
    diffusion_dataset = AnnotatedPolarsSpectrumDataset(
        pl.from_pandas(df), peptides=peptides
    )
     
    diffusion_dataloader = DataLoader(
        dataset=diffusion_dataset,
        batch_size=64,
        shuffle=False,
        collate_fn=collate_batches(
            residues=model.residues,
            max_length=model.config.max_length,
            time_steps=decoder.time_steps,
            annotated=True
        )
    )

    return diffusion_dataloader

def predict(batch, decoder, device):
    spectra, spectra_padding_mask, precursors, peptides, peptide_padding_mask, spectrum_id = batch
    spectra = spectra.to(device)
    spectra_padding_mask = spectra_padding_mask.to(device)
    precursors = precursors.to(device)
    peptides = peptides.to(device)
    peptide_padding_mask = peptide_padding_mask.to(device)

    with torch.no_grad():
        batch_predictions, batch_log_probs = decoder.decode(
            spectra=spectra,
            spectra_padding_mask=spectra_padding_mask,
            precursors=precursors,
            initial_sequence=peptides
        )
    return batch_predictions, batch_log_probs, spectrum_id

def main(args):
    device = args.device

    # Load models and checkpoints
    diffusion_model = MultinomialDiffusion.load(args.model_path)
    diffusion_model = diffusion_model.to(device).eval()
    diffusion_decoder = DiffusionDecoder(model=diffusion_model)

    # Load data
    diffusion_dataloader = pddf_to_dl(
        df=pd.read_feather(args.input_file),
        model=diffusion_model,
        decoder=diffusion_decoder
    )

    # Initialize storage objects for predictions
    predictions = []
    log_probs = []
    spectrum_index_list = []

    for batch in tqdm(diffusion_dataloader, total=len(diffusion_dataloader)):
        
        # Prediction step for each batch
        batch_predictions, batch_log_probs, spectrum_index = predict(
            batch=batch,
            decoder=diffusion_decoder,
            device=device
        )

        # Save predictions to list
        predictions.extend(["".join(sequence) for sequence in batch_predictions])
        log_probs.extend(batch_log_probs)
        spectrum_index_list.extend(spectrum_index)

    # Transform predictions to dataframe
    diffusion_predictions = pd.DataFrame({
        "spectrum_index": spectrum_index_list,
        "predictions": predictions,
        "log_probabilities": log_probs,
    })

    # Read mapping file
    with open(args.mapping_file, "r") as json_file:
        mapping_file = json.load(json_file)
    
    diffusion_predictions["spectrum_id"] = diffusion_predictions["spectrum_index"].apply(
        lambda x: mapping_file[str(x)]["title"]
    )

    if not os.path.exists(args.output_folder):
        os.mkdir(args.output_folder)

    diffusion_predictions.to_csv(
        os.path.join(args.output_folder, os.path.basename(args.input_file).split(".")[0]) + ".csv",
        index=False)
    diffusion_predictions.to_csv(
        os.path.basename(args.input_file).split(".")[0] + ".csv",
        index=False)

    print("Finished instanovo+ process.")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run Instanovo+ refinement using a diffusion model on a pickled input file.")
    parser.add_argument("-i", "--input_file", required=True, help="Pickled input file (https://github.com/instadeepai/InstaNovo)")
    parser.add_argument("-m", "--model_path", required=True, help="Path to model checkpoints of diffusion model.")
    parser.add_argument("-c", "--mapping_file", required=True, help="JSON-file mapping indices (int) to spectrum_ids (str).")
    parser.add_argument("-o", "--output_folder", default="", help="Folder to store the output in.")
    parser.add_argument("-d", "--device", default="cuda", help="Device to run the predictions on (cuda or cpu)")

    args = parser.parse_args()

    main(args)