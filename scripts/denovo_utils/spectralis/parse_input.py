import pandas as pd
import os
import logging
import argparse

from pyteomics.mztab import MzTab
from pyteomics import mgf
from enum import Enum
from tqdm import tqdm
from psm_utils import PSM, PSMList, Peptidoform


logger = logging.getLogger(__name__)
logging.basicConfig(filename="denovo_output_parsing.log", level=logging.INFO)

tqdm.pandas()

# Define parser functions
MODIFICATION_MAPPING = {
    "casanovo": {
        "+57.021": "[UNIMOD:4]",
        "+15.995": "[UNIMOD:35]",
        "+0.984": "[UNIMOD:7]", # Check both D - N+0.984 and E - Q+0.984
        "+43.006-17.027": "[+25.980265]-", # [UNIMOD:Carbamylation][UNIMOD:Ammonia-loss]
        "+42.011": "[UNIMOD:1]-",
        "+43.006": "[UNIMOD:5]-",
        "-17.027": "[UNIMOD:385]-"
    },
    "instanovo": {
        "C(+57.02)": "[UNIMOD:4]",
        "M(ox)": "M[UNIMOD:35]",
        "M(+15.99)": "M[UNIMOD:35]",
        "N(+.98)": "N[UNIMOD:7]",
        "Q(+.98)": "Q[UNIMOD:7]"
    },
    "contranovo": {
        "C+57.021": "C[UNIMOD:4]",
        "+43.006-17.027": "[+25.980265]-",
        "-17.027+43.006": "[+25.980265]-",
        "+42.011": "[UNIMOD:1]-",  # Acetylation
        "+43.006": "[UNIMOD:5]-",  # 5
        "-17.027": "[UNIMOD:385]-",  # NH3 loss
        # AA mods:
        "M+15.995": "M[UNIMOD:35]",  # Met oxidation
        "N+0.984": "N[UNIMOD:7]", # Asn deamidation
        "Q+0.984": "Q[UNIMOD:7]"  # Gln deamidation
    },
    "pepnet": {
        "Z": "M[UNIMOD:35]",
        "C": "C[UNIMOD:4]"
    },
    "novob": {
        "C": "C[UNIMOD:4]",
        "m": "M[UNIMOD:35]",
        "n": "N[UNIMOD:7]",
        "q": "Q[UNIMOD:7]",
        "s": "S[UNIMOD:21]", 
        "t": "T[UNIMOD:21]",
        "y": "Y[UNIMOD:21]"
    },
    "general": {
        "acetylation": "[UNIMOD:1]"
    }
}

# Generate a comprehensive modification dictionary
def generate_modification_labels(mapping: dict):
    modification_dictionaries = list(mapping.values())
    all_modification_labels = {}

    for d in modification_dictionaries:
        for k, v in d.items():
            if k not in all_modification_labels:
                all_modification_labels[k] = v

    return all_modification_labels

ALL_MODIFICATION_LABELS = generate_modification_labels(MODIFICATION_MAPPING)

def generate_spectralis_mod_map() -> dict:
    spectralis_modifications = [
        "C[UNIMOD:4]", "M[UNIMOD:35]", 'Q[UNIMOD:7]', 'N[UNIMOD:7]'
    ]

    spectralis_mod_map = {
        'Q[UNIMOD:7]': 'E',
        'N[UNIMOD:7]': 'D'
    }

    for modification in ALL_MODIFICATION_LABELS.values():
        if modification not in spectralis_modifications:
            spectralis_mod_map[modification] = ""
    return spectralis_mod_map

MODIFICATION_MAPPING_TO_SPECTRALIS = generate_spectralis_mod_map()


def parse_peptidoform(peptide: str, mapping: dict, max_length=30):
    peptide_parsed = peptide
    for k, v in mapping.items():
        if ("-" in v) and (not peptide_parsed.startswith(k)):
            peptide_parsed = peptide_parsed.replace(k, v[:-1])
        else:
            peptide_parsed = peptide_parsed.replace(k, v)

    try:
        peptidoform = Peptidoform(peptide_parsed)
        if (len(peptidoform) > max_length) or (peptidoform.precursor_charge > 6) or (len(peptidoform) < 2):
            return None
        return peptidoform
    except:
        logging.info(f"Failed to parse: {peptide}")
        return None


def casanovo_parser(result_path: str, mgf_path: str, mapping: dict, max_length=30):

    # ASSUMPTION: 
    # The CasaNovo spectra_ref columns has prefix ms_run[1]:index= and the number is a count
    # of spectra, starting from 0 and going to n for a file with n spectra
    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    _ = mgf_file.pop("charge")
    result = pd.DataFrame(MzTab(result_path).spectrum_match_table).set_index("PSM_ID").reset_index().reset_index()
    run = os.path.basename(result_path)

    mgf_file = mgf_file.reset_index()
    result["index"] = result.spectra_ref.apply(lambda x: x.split("=")[-1])

    length_mgf = len(mgf_file)
    length_result = len(result)
    assert length_mgf == length_result

    # Fuse the metadata of the spectra with result file
    joined_file = pd.concat([mgf_file, result], axis=1)
    length_join = len(joined_file)

    # Sanity checks
    assert length_mgf == length_join
    assert length_result == length_join

    # Parse to psm utils type format
    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["sequence"]) + "/" + str(int(x["charge"])), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(
        lambda x: x[0]
    )

    # Parse peptidoforms with max_length == 30
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)
    
    psm_list = joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["title"],
                run=run,
                score=x["search_engine_score[1]"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source=x["search_engine"][0] + x["search_engine"][1],
                metadata={
                    "aa_scores": x["opt_ms_run[1]_aa_scores"],
                    "calc_mass_to_charge": x["calc_mass_to_charge"],
                    "spectra_ref": x["spectra_ref"],
                    "scans": x["scans"]
                }
            ),
            axis=1
        ).tolist()
    

    psmlist = PSMList(psm_list=psm_list)
    return psmlist

def instanovo_parser(result_path: str, mgf_path: str, mapping: dict, max_length=30):

    # ASSUMPTION: 
    # The output of Instanovo has the same length and order (in terms of spectra) as the mgf-file
    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    result = pd.read_csv(result_path)
    run = os.path.basename(result_path)

    length_mgf = len(mgf_file)
    length_result = len(result)
    assert length_mgf == length_result

    # Fuse the metadata of the spectra with result file
    joined_file = pd.concat([mgf_file, result], axis=1)
    length_join = len(joined_file)

    # Sanity checks
    assert length_mgf == length_join
    assert length_result == length_join

    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["preds"]) + "/" + str(int((x["charge"][0]))), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(
        lambda x: x[0]
    )
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    psmlist = PSMList(
        psm_list=joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["title"],
                run=run,
                score=x["log_probs"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source="InstaNovo",
                metadata={
                    "scans": x["scans"]
                }
            ),
            axis=1
        ).tolist()
    )
    return psmlist

def contranovo_parser(result_path: str, mgf_path: str, mapping: dict, max_length=30):

    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    _ = mgf_file.pop("charge")
    result = pd.DataFrame(MzTab(result_path).spectrum_match_table).set_index("PSM_ID").reset_index().reset_index()
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    joined_file = result.rename(
        columns={
            "spectra_ref": "title"
        }
    ).merge(mgf_file, on="title")

    # Sanity checks
    assert len(result) == len(joined_file)
    if len(mgf_file) < len(joined_file):
        logging.info(f"{result_path} for ContraNovo: dropped {len(mgf_file)-len(joined_file)} spectra.")

    # Parse to psm utils type format
    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["sequence"]) + "/" + str(int((x["charge"]))), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(
        lambda x: x[0]
    )
    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    psmlist = PSMList(
        psm_list=joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["title"],
                run=run,
                score=x["search_engine_score[1]"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source="ContraNovo",
                metadata={
                    "aa_scores": x["opt_ms_run[1]_aa_scores"],
                    "calc_mass_to_charge": x["calc_mass_to_charge"],
                    "scans": x["scans"]
                }
            ),
            axis=1
        ).tolist()
    )
    return psmlist

def novob_parser(result_path: str, mgf_path: str, mapping: dict, max_length=30):

    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    _ = mgf_file.pop("charge")
    result = pd.read_csv(result_path, sep="\t", header=None).rename(
        columns={
            0: "Mcount",
            1: "charge",
            2: "peptide_mass",
            3: "sequence_forward",
            4: "mass_forward",
            5: "probability_forward",
            6: "sequence_reverse",
            7: "mass_reverse",
            8: "probability_reverse",
            9: "scans"
        }
    )
    result["scans"] = result.apply(lambda x: x["scans"][2:-1], axis=1) # Remove b' ... '
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    joined_file = result.merge(mgf_file, on="scans")

    # Sanity checks
    assert len(result) == len(joined_file)
    if len(mgf_file) < len(joined_file):
        logging.info(f"{result_path} for ContraNovo: dropped {len(mgf_file)-len(joined_file)} spectra.")

    # Parse to psm utils type format
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(
        lambda x: x[0]
    )

    def select_top_PSM(x, max_length):
        if x["probability_forward"] > x["probability_reverse"]:
            peptidoform = eval(x["sequence_forward"])[0] + "/" + str(int(x["charge"]))
            mass_error = x["mass_forward"]
            proba = x["probability_forward"]
        else:
            peptidoform = eval(x["sequence_reverse"])[0] + "/" + str(int(x["charge"]))
            mass_error = x["mass_reverse"]
            proba = x["probability_reverse"]

        peptide = parse_peptidoform(peptidoform, mapping, max_length)
        if peptide:
            return None
        else:
            return PSM(
                peptidoform=parse_peptidoform(peptidoform, mapping),
                spectrum_id=x["title"],
                run=run,
                score=proba,
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source="NovoB",
                metadata={
                    "ppm_error": mass_error,
                    "scans": x["scans"]
                }
            )

    joined_file = joined_file.progress_apply(
        lambda x: select_top_PSM(x, max_length),
        axis=1
    )
    joined_file = joined_file.dropna()

    psmlist = PSMList(
        psm_list=joined_file.tolist()
    )
    return psmlist

def pepnet_parser(result_path: str, mgf_path: str, mapping: dict, max_length=30):
    
    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    result = pd.read_csv(result_path, sep="\t").rename(
        columns={"TITLE": "title"}
    )
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    joined_file = result.merge(mgf_file, on="title")

    # Sanity checks
    assert len(result) == len(joined_file)

    joined_file["peptidoform"] = joined_file.apply(
        lambda x: str(x["DENOVO"]) + "/" + str(int((x["charge"][0]))), axis=1
    )
    joined_file["precursor_mz"] = joined_file["pepmass"].apply(
        lambda x: x[0]
    )

    joined_file["peptidoform"] = joined_file["peptidoform"].apply(
        lambda x: parse_peptidoform(x, mapping, max_length)
    )
    joined_file = joined_file.dropna(subset=["peptidoform"]).reset_index(drop=True)

    psmlist = PSMList(
        psm_list=joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=x["peptidoform"],
                spectrum_id=x["title"],
                run=run,
                score=x["Score"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source="PepNet",
                metadata={
                    "positional_scores": x["Positional Score"],
                    "ppm_error": x["PPM Difference"]
                }
            ),
            axis=1
        ).tolist()
    )
    return psmlist


class DenovoEngineNotSupported(ValueError):
    def __init__(self, parameter, allowed_values):
        self.parameter = parameter
        self.allowed_values = allowed_values
        super().__init__(f"Invalid parameter value: {parameter}. Allowed values are: {allowed_values}")


# Define supported parsers for de novo search engines as an Enum with associated parser functions
class DenovoEngine(Enum):
    CASANOVO = ("casanovo", casanovo_parser)
    INSTANOVO = ("instanovo", instanovo_parser)
    CONTRANOVO = ("contranovo", contranovo_parser)
    NOVOB = ("novob", novob_parser)
    PEPNET = ("pepnet", pepnet_parser)

    def __init__(self, label, parser_func):
        self.label = label
        self.parser_func = parser_func

    def parse(self, result_path: str, mgf_path: str):
        return self.parser_func(result_path, mgf_path, MODIFICATION_MAPPING[self.label])

    @classmethod
    def from_string(cls, label: str):
        for engine in cls:
            if engine.label == label:
                return engine
        raise DenovoEngineNotSupported(label, [e.label for e in cls])



def denovo_to_psmlist(result_path: str, mgf_path: str, denovo: DenovoEngine) -> PSMList:
    """
    Parse results to PSM utils using the specified de novo engine.

    Args:
        results (str): The results to be parsed.
        denovo (str): The de novo engine to use. Must be one of the supported engines.

    Raises:
        DenovoEngineNotSupported: If the specified de novo engine is not supported.
    """

    denovo_engine = DenovoEngine.from_string(denovo)
    psmlist = denovo_engine.parse(result_path, mgf_path)
    return psmlist


def annotate_psms_to_mgf(psmlist: PSMList, mgf_path: str, output_folder: str, modification_mapping: dict, exclusion_list = []):
    filename = os.path.basename(mgf_path).split(".")[0]
    mgf_file = mgf.read(mgf_path)
    annotation_dict = psmlist.to_dataframe().loc[
        :, ["peptidoform", "spectrum_id"]
        ].set_index("spectrum_id").to_dict("index")
    
    mgf_annotated = []

    for spectrum in mgf_file:
        
        # When the peptide prediction didnt fit the requirements, skip the spectrum.
        try:
            annotated_peptide = annotation_dict[
                spectrum["params"]["title"]
            ]["peptidoform"].proforma
        except:
            continue

        for k, v in modification_mapping.items():
            annotated_peptide = annotated_peptide.replace(k, v)
        
        for exclude in exclusion_list:
            if exclude in annotated_peptide:
                continue
        
        spectrum["params"]["seq"] = annotated_peptide.split("/")[0]
        mgf_annotated.append(spectrum)

    mgf.write(
        mgf_annotated,
        output= os.path.join(output_folder, filename) + "_annotated.mgf",
        header=mgf_file.header
    )


def main(args):
    
    psmlist = denovo_to_psmlist(
        result_path=args.search_result,
        mgf_path=args.mgf_file,
        denovo=args.denovo_engine
    )   
    annotate_psms_to_mgf(
        psmlist=psmlist,
        mgf_path=args.mgf_file,
        output_folder=args.output_folder,
        modification_mapping=MODIFICATION_MAPPING_TO_SPECTRALIS,
        exclusion_list=args.exclusion_list
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Write mgf files, annotated by a given de novo search")
    parser.add_argument('-r', '--search_result', required=True, help="Path to result file generated by a de novo search engine.")
    parser.add_argument('-m', '--mgf_file', required=True, help="Path to unannotated mgf-file")
    parser.add_argument('-d', '--denovo_engine', required=True, help="The denovo engine used to generate the files search result file.")
    parser.add_argument('-o', '--output_folder', default="", help="Output folder to store annotated mgf-file")
    parser.add_argument('-x', '--exclusion_list', default=[], help="List containing tags in peptides that should be dropped (e. g. modification tags...)")

    args = parser.parse_args()

    main(args)