
from psm_utils import PSMList, Peptidoform
import pandas as pd

def spectralis_parser(psm_list: PSMList, result_path):

    spectralis_df = pd.read_csv(result_path)
    spectralis_df["same_sequence"] = spectralis_df.apply(
        lambda x: x["peptide_spectralis-ea"]==x["peptide_init"], axis=1
    )
    spectralis_df["spectrum_id"] = spectralis_df["scans"].apply(lambda x: x.split("||")[0])
    spectralis_dict = spectralis_df.set_index("spectrum_id").to_dict("index")

    for psm in psm_list:

        spectralis_entry = spectralis_dict[psm.spectrum_id]
        peptide_init = Peptidoform(
            spectralis_entry["peptide_init"].replace("I", "L").replace("UNLMOD", "UNIMOD")
            +"/"+str(psm.peptidoform.precursor_charge)
        )
        unchanged = peptide_init==psm.peptidoform

        psm.provenance_data.update({"spectralis": {
            "spectralis_peptide_init": peptide_init,
            "spectralis_peptide_init_unchanged": unchanged,
            "spectralis_refined": spectralis_entry["same_sequence"],
            "spectralis_peptide_ea": Peptidoform(
                spectralis_entry["peptide_spectralis-ea"]+"/"+str(psm.peptidoform.precursor_charge)
            ),
            "spectralis_score_init": spectralis_entry["score_init"],
            "spectralis_score_ea": spectralis_entry["score_spectralis-ea"]
        }})

def instanovoplus_parser(psm_list: PSMList, result_path):
    pass

REFINEMENT_PARSERS = {
    "spectralis": spectralis_parser,
    "instanovoplus": instanovoplus_parser
}

def add_refinement(psm_list: PSMList, result_path, refiner):
    """Refinement results are stored in nested dict in provenance data."""
    if refiner not in ["instanovoplus", "spectralis"]:
        raise NotImplementedError(
            f"Refiners only implemented for {list(REFINEMENT_PARSERS.keys())}"
        )

    REFINEMENT_PARSERS[refiner](
        psm_list=psm_list,
        result_path=result_path
    )
    