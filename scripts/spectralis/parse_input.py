from enum import Enum
import pandas as pd
from pyteomics.mztab import MzTab
from pyteomics import mgf
import os
import logging
from tqdm import tqdm
from psm_utils import PSM, PSMList


logger = logging.getLogger(__name__)
logging.basicConfig(filename="denovo_output_parsing.log", level=logging.INFO)

tqdm.pandas()

# Define parser functions
MODIFICATION_MAPPING = {
    "casanova": {
        "+57.021": "[U:Carbamidomethyl]",
        "+15.995": "[U:Oxidation]",
        "+0.984": "[U:Deamidation]", # Check both D - N+0.984 and E - Q+0.984
        "+42.011": "[U:Acetylation]-",
        "+43.006": "[U:Carbamylation]-",
        "-17.027": "[U:Ammonia-loss]-",
        "+43.006-17.027": "[+25.980265]-" # [U:Carbamylation][U:Ammonia-loss]
    },
    "instanovo": {
        "C(+57.02)": "[U:Carbamidomethyl]",
        "M(ox)": "M[U:Oxidation]",
        "M(+15.99)": "M[U:Oxidation]",
        "N(+.98)": "N[U:Deamidation]",
        "Q(+.98)": "Q[U:Deamidation]"
    },
    "contranovo": {
        "C[+57.021]": "C[U:Carbamidomethyl]",
        "+42.011": "[U:Acetylation]-",  # Acetylation
        "+43.006": "[U:Carbamylation]-",  # Carbamylation
        "-17.027": "[U:Ammonia-loss]-",  # NH3 loss
        "+43.006-17.027": "[+25.980265]-",
        # AA mods:
        "M+15.995": "M[U:Oxidation]",  # Met Oxidation
        "N+0.984": "N[U:Deamidation]", # Asn Deamidation
        "Q+0.984": "Q[U:Deamidation]"  # Gln Deamidation
    },
    "pepnet": {
        "Z": "M[U:Oxidation]",
        "C": "C[U:Carbamidomethyl]"
    },
    "novob": {
        "C": "C[U:Carbamidomethyl]",
        "m": "M[U:Oxidation]",
        "n": "N[U:Deamidation]",
        "q": "Q[U:Deamidation]",
        "s": "S[U:Phosphorylation]", 
        "t": "T[U:Phosphorylation]",
        "y": "Y[U:Phosphorylation]"
    }
}


def parse_peptidoform(peptide, mapping):
    peptide_parsed = peptide
    for k, v in mapping.items():
        peptide_parsed = peptide_parsed.replace(k, v)
    return peptide_parsed


def casanovo_parser(result_path: str, mgf_path: str, mapping: dict):

    # ASSUMPTION: 
    # The output of CasaNovo has the same length and order (in terms of spectra) as the mgf-file
    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    result = pd.DataFrame(MzTab(result_path).spectrum_match_table).set_index("PSM_ID").reset_index().reset_index()
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

    # Parse to psm utils type format
    result["peptidoform"] = result.apply(
        lambda x: str(x["sequence"]) + "/" + str(x["charge"]), axis=1
    )
    result["precursor_mz"] = result["pepmass"].apply(
        lambda x: x[0]
    )

    psmlist = PSMList(
        result.progress_apply(
            lambda x: PSM(
                peptidoform=parse_peptidoform(x["peptidoform"], mapping),
                spectrum_id=x["title"],
                run=run,
                score=x["search_engine_score[1]"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source=x["search_engine"],
                metadata={
                    "aa_scores": x["opt_ms_run[1]_aa_scores"],
                    "calc_mass_to_charge": x["calc_mass_to_charge"],
                    "spectra_ref": x["spectra_ref"],
                    "scans": x["scans"]
                }
            ),
            axis=1
        ).tolist()
    )
    return psmlist

def instanovo_parser(result_path: str, mgf_path: str, mapping: dict):

    # ASSUMPTION: 
    # The output of Instanovo has the same length and order (in terms of spectra) as the mgf-file
    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    result = pd.read_csv(result_path).reset_index()
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

    result["peptidoform"] = result.apply(
        lambda x: str(x["preds"]) + "/" + str(int((x["charge"][0])))
    )
    result["precursor_mz"] = result["pepmass"].apply(
        lambda x: x[0]
    )

    psmlist = PSMList(
        result.progress_apply(
            lambda x: PSM(
                peptidoform=parse_peptidoform(x["peptidoform"], mapping),
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

def contranovo_parser(result_path: str, mgf_path: str, mapping: dict):

    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    result = pd.DataFrame(MzTab(result_path).spectrum_match_table).set_index("PSM_ID").reset_index().reset_index()
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    joined_file = result.rename(
        {
            "spectra_ref": "title"
        }
    ).merge(mgf_file, on="title")

    # Sanity checks
    assert len(result) == len(joined_file)
    if len(mgf_file) < len(joined_file):
        logging.info(f"{result_path} for ContraNovo: dropped {len(mgf_file)-len(joined_file)} spectra.")

    # Parse to psm utils type format
    result["peptidoform"] = result.apply(
        lambda x: str(x["preds"]) + "/" + str(int((x["charge"][0])))
    )
    result["precursor_mz"] = result["pepmass"].apply(
        lambda x: x[0]
    )

    psmlist = PSMList(
        joined_file.progress_apply(
            lambda x: PSM(
                peptidoform=parse_peptidoform(x["peptidoform"], mapping),
                spectrum_id=x["title"],
                run=run,
                score=x["search_engine_score[1]"],
                precursor_mz=x["precursor_mz"],
                retention_time=x["rtinseconds"],
                source=x["search_engine"],
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

def novob_parser(result_path: str, mgf_path: str, mapping: dict):

    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    result = pd.read_csv(result_path, sep="\t", header=None).rename(
        columns={
            0: "Mcount",
            1: "charge",
            2: "pepmass",
            3: "sequence_forward",
            4: "mass_forward",
            5: "probability_forward",
            6: "sequence_reverse",
            7: "mass_reverse",
            8: "probability_reverse",
            9: "scans"
        }
    )
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    joined_file = result.merge(mgf_file, on="scans")

    # Sanity checks
    assert len(result) == len(joined_file)
    if len(mgf_file) < len(joined_file):
        logging.info(f"{result_path} for ContraNovo: dropped {len(mgf_file)-len(joined_file)} spectra.")

    # Parse to psm utils type format
    result["precursor_mz"] = result["pepmass"].apply(
        lambda x: x[0]
    )

    def select_top_PSM(x):
        if x["probability_forward"] > x["probability_reverse"]:
            peptidoform = eval(x["sequence_forward"])[0] + "/" + str(x["charge"])
            mass_error = x["mass_forward"]
            proba = x["probability_forward"]
        else:
            peptidoform = eval(x["sequence_reverse"])[0] + "/" + str(x["charge"])
            mass_error = x["mass_reverse"]
            proba = x["probability_reverse"]
        
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
                "calc_mass_to_charge": x["calc_mass_to_charge"],
                "scans": x["scans"]
            }
        )

    psmlist = PSMList(
        joined_file.progress_apply(
            lambda x: select_top_PSM(x),
            axis=1
        ).tolist()
    )
    return psmlist

def pepnet_parser(result_path: str, mgf_path: str, mapping: dict):
    
    mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())
    result = pd.read_csv(result_path, sep="\t").rename(
        columns={"TITLE", "title"}
    )
    run = os.path.basename(result_path)

    # Fuse the metadata of the spectra with result file
    joined_file = result.merge(mgf_file, on="title")

    # Sanity checks
    assert len(result) == len(joined_file)

    result["peptidoform"] = result.apply(
        lambda x: x["DENOVO"] + "/" + str(int((x["charge"][0])))
    )
    result["precursor_mz"] = result["pepmass"].apply(
        lambda x: x[0]
    )

    psmlist = PSMList(
        result.progress_apply(
            lambda x: PSM(
                peptidoform=parse_peptidoform(x["peptidoform"], mapping),
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
        self.parser_func(result_path, mgf_path, MODIFICATION_MAPPING[self.label])

    @classmethod
    def from_string(cls, label: str):
        for engine in cls:
            if engine.label == label:
                return engine
        raise DenovoEngineNotSupported(label, [e.label for e in cls])



def denovo_to_psmlist(result_path: str, mgf_path: str, denovo: DenovoEngine):
    """
    Parse results to PSM utils using the specified de novo engine.

    Args:
        results (str): The results to be parsed.
        denovo (str): The de novo engine to use. Must be one of the supported engines.

    Raises:
        DenovoEngineNotSupported: If the specified de novo engine is not supported.
    """

    denovo_engine = DenovoEngine.from_string(denovo)
    denovo_engine.parse(result_path, mgf_path)