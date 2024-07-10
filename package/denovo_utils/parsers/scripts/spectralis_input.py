from ..converters import DenovoEngineConverter
from ..io import psmlist_to_mgf, MGFWriter
from ..constants import MODIFICATION_MAPPING_TO_SPECTRALIS, EXTENSIONS
import argparse
import os

def main(args):

    # Parse the denovo_engines variable
    if 'all' in args.denovo_engines:
        engine_extension = EXTENSIONS
    else:
        engine_extension = {}
        for engine in args.denovo_engines:
            try:
                engine_extension[engine] = EXTENSIONS[engine]
            except:
                raise KeyError(
                    f"Attribute '{engine}' was passed as argument, but only '{list(EXTENSIONS.keys())}' are allowed!"
                )

    # Set up variables
    filename = os.path.basename(args.mgf_path).split(".")[0]

    # Set up writer
    writer = MGFWriter(
        mgf_path=args.mgf_path,
        modification_mapping=MODIFICATION_MAPPING_TO_SPECTRALIS,
        output_folder=args.output_folder
    )

    # Read in results in the writer
    for denovo_engine, extension in engine_extension.items():
        writer.read_result(
            result_path=os.path.join(
                args.result_path, denovo_engine, filename+extension
            ),
            denovo_engine=denovo_engine
        )
    
    writer.merge()
    writer.write()

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

    parser = argparse.ArgumentParser(description="Write mgf files, annotated by a given de novo search")
    parser.add_argument('-r', '--result_path', required=True, help="Path to the root folder, containing all result files.")
    parser.add_argument('-m', '--mgf_path', required=True, help="Path to unannotated mgf-file")
    parser.add_argument('-d', '--denovo_engines', required=True, nargs='+', help="The denovo engine used to generate the files search result file.")
    parser.add_argument('-o', '--output_folder', default=f"{os.getcwd()}", help="Output folder to store annotated mgf-file")
    parser.add_argument('-x', '--exclusion_list', default=[], help="List containing tags in peptides that should be dropped (e. g. modification tags...)")

    args = parser.parse_args()

    main(args)