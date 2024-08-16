import argparse
from ..io.writers import reformat_mgf




def main(args):

    reformat_mgf(
        path_mgf=args.input_file,
        output_dir=args.output_dir,
        min_peak_length=int(args.min_peak_length),
        max_charge=int(args.max_charge)
    )

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Reformat mgf, to make them compatible with expected de novo input."
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help="Path to the mgf-file which needs to be reformatted"
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        default="",
        help="Output folder to store reformatted mgf-file."
    )
    parser.add_argument(
        "-p",
        "--min_peak_length",
        default=15,
        help="Filter out spectra with less than min_peak_length peaks."
    )
    parser.add_argument(
        "-c",
        "--max_charge",
        default=8,
        help="Filter out spectra above max_charge."
    )
    
    args = parser.parse_args()
    
    main(args)