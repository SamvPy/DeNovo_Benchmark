from tqdm import tqdm
from glob import glob
from ..io import filter_on_charge
import argparse

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