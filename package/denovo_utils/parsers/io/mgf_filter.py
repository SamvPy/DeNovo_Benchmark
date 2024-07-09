# Run script with NovoB environment

from psm_utils import PSMList
from pyteomics import mgf
from glob import glob
from tqdm import tqdm
import argparse
import logging
import os

logger = logging.getLogger(__name__)
logging.basicConfig(filename="mgf_file_parsing.log", level=logging.INFO)

def filter_on_charge(mgf_file_path, output_dir, max_charge):

    mgf_file = mgf.read(mgf_file_path)
    mgf_new = []

    for i, spectrum in tqdm(enumerate(mgf_file), desc="Going through spectra..."):
        try:
            _ = spectrum["params"]["charge"][0]

            if _ <= max_charge:
                mgf_new.append(spectrum)

        except:
            continue
    
    mgf.write(
        mgf_new,
        output=output_dir + "/" + mgf_file_path.split("/")[-1],
        header=mgf_file.header,
    )

    logging.info(
        f"{mgf_file_path}: {i} -> {len(mgf_new)} spectra. (-{i+1-len(mgf_new)})"
    )

def main(args):
    if args.file_type == "file":
        filter_on_charge(args.mgf_folder, args.output_folder, args.max_charge)

    else:
        for file_path in tqdm(glob(args.mgf_folder + "/*.mgf"),
                              desc="Filtering mgf-files..."):
            filter_on_charge(file_path, args.output_folder, args.max_charge)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filter mgf files based")
    parser.add_argument('-f', '--mgf_folder', required=True, help="Folder to mgf files")
    parser.add_argument('-o', '--output_folder', required=True, help="Output folder to store mgf files")
    parser.add_argument('-t', '--file_type', default="dir", help="dir or file?")
    parser.add_argument('-c', '--max_charge', default=8, help="Filter out spectra in mgf with charges higher than -c")

    args = parser.parse_args()

    main(args)
