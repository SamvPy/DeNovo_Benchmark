import argparse
import os

import numpy as np

from ..constants import EXTENSIONS, MODIFICATION_MAPPING_TO_SPECTRALIS
from ..converters import DenovoEngineConverter
from ..io import MGFWriter, psmlist_to_mgf


def main(args):

    print('Processing for engine: {}'.format(args.denovo_engine))
    # # Parse the denovo_engines variable
    # if "all" in args.denovo_engines:
    #     engine_extension = EXTENSIONS
    # else:
    #     engine_extension = {}
    #     for engine in args.denovo_engines:
    #         try:
    #             engine_extension[engine] = EXTENSIONS[engine]
    #         except:
    #             raise KeyError(
    #                 f"Attribute '{engine}' was passed as argument, but only '{list(EXTENSIONS.keys())}' are allowed!"
    #             )

    # Set up variables
    filename = os.path.basename(args.mgf_path).split(".")[0]
    extension = '.'+os.path.basename(args.result_path).split('.')[-1]

    # If the results of a de novo engine are converted to a different format.
    # Can happen if the results are filtered before post-processing and stored in parquet format.
   
    if extension != EXTENSIONS[args.denovo_engine]:
        parser_format = extension[1:]
    else:
        parser_format = args.denovo_engine

    # Read in results first and split by psm candidate rank.
    if args.batches:
        
        parser = DenovoEngineConverter.select(parser_format)
        results = parser.parse(
            result_path=os.path.join(args.result_path),
            mgf_path=args.mgf_path
        )
        results_by_rank = {}

        ranks = results['rank']
        for rank in np.unique(ranks):
            results_by_rank[rank] = results[ranks==rank]

        for rank, psmlist in results_by_rank.items():
            # Init writer
            writer = MGFWriter(
                mgf_path=args.mgf_path,
                modification_mapping=MODIFICATION_MAPPING_TO_SPECTRALIS,
                output_folder=args.output_folder
            )

            # Load results
            writer.load_result(psmlist=psmlist, denovo_engine=args.denovo_engine)

            # Parse
            writer.merge()

            # Out
            out_filename = writer.filename + f'_{rank}.spectralis.mgf'
            writer.write(out_filename=out_filename)

    # Write all results in a single mgf file. MGF can become extremely big.
    else:
        # Set up writer
        writer = MGFWriter(
            mgf_path=args.mgf_path,
            modification_mapping=MODIFICATION_MAPPING_TO_SPECTRALIS,
            output_folder=args.output_folder,
        )

        # Read in results in the writer
        writer.read_result(
            result_path=args.result_path,
            denovo_engine=args.denovo_engine
        )

        writer.merge()

        out_filename = writer.filename + f'.spectralis.mgf'
        writer.write(out_filename=out_filename)

    # denovo_engine = DenovoEngineConverter.select(label=args.denovo_engine)
    # psmlist = denovo_engine.parse(
    #     result_path=args.result_path,
    #     mgf_path=args.mgf_path
    # )

    # psmlist_to_mgf(
    #     psmlist=psmlist,
    #     mgf_path=args.mgf_path,
    #     output_folder=args.output_folder,
    #     modification_mapping=MODIFICATION_MAPPING_TO_SPECTRALIS,
    #     exclusion_list=args.exclusion_list,
    #     annotated=True
    # )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Write mgf files, annotated by a given de novo search"
    )
    parser.add_argument(
        "-r",
        "--result_path",
        required=True,
        help="Path to the result file.",
    )
    parser.add_argument(
        "-m", "--mgf_path", required=True, help="Path to unannotated mgf-file"
    )
    parser.add_argument(
        "-d",
        "--denovo_engine",
        required=True,
        help="The denovo engine used to generate the files search result file.",
    )
    parser.add_argument(
        "-b",
        "--batches",
        default=True,
        help="If set to True, separate MGF-files will be created for each engine and psm-rank."
    )

    parser.add_argument(
        "-o",
        "--output_folder",
        default=f"{os.getcwd()}",
        help="Output folder to store annotated mgf-file",
    )
    parser.add_argument(
        "-x",
        "--exclusion_list",
        default=[],
        help="List containing tags in peptides that should be dropped (e. g. modification tags...)",
    )

    args = parser.parse_args()

    main(args)
