import pyopenms as oms
from pyteomics import mgf
import os
import pandas as pd
from pyteomics.mass import calculate_mass
from psm_utils import PSM, PSMList

from spectrum_utils.proforma import Proteoform, Modification
from spectrum_utils.fragment_annotation import get_theoretical_fragments
from psm_utils import Peptidoform
from spectrum_utils import proforma

from denovo_utils.analysis import calculate_hyperscore
from denovo_utils.parsers import proforma_to_theoretical_spectrum
from denovo_utils.parsers.converters import SpectralisParser, DenovoEngineConverter
from ms2rescore.feature_generators import BasicFeatureGenerator, MS2PIPFeatureGenerator, DeepLCFeatureGenerator

import deeplc
from deeplc.plot import scatter
from ms2rescore.report.charts import (
    calculate_feature_qvalues,
    feature_ecdf_auc_bar,
    fdr_plot,
    ms2pip_correlation,
)

from tqdm import tqdm

import spectrum_utils.plot as sup

from denovo_utils.utils.pandas import get_psm_type, get_spectralis_score
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np


def peptidoform_has_modification(peptidoform, allowed_modificiations=["[UNIMOD:4]", "[UNIMOD:35]"]):
    peptidoform_str = peptidoform.proforma
    for allowed_modification in allowed_modificiations:
        peptidoform_str = peptidoform_str.replace(allowed_modification, "")
    
    return "[" in peptidoform_str

def hyperscore_difference(row, reference):
    try:
        reference_hyperscore = reference[row["spectrum_id"]]
        return float(row["hyperscore"])-float(reference_hyperscore)
    except:
        return None
    
def evaluate_prediction_isobaricity(
        row, ground_truth_peptide, ground_truth_hyperscore
):
    try:
        sequence_match = ground_truth_peptide[
            row["spectrum_id"]
        ] == row["peptide"]
        if sequence_match:
            return "Match"
        
        ref_hyperscore = ground_truth_hyperscore[row["spectrum_id"]]
        if row["hyperscore"] == ref_hyperscore:
            return "Isobaric"
        
        elif row["hyperscore"] > ref_hyperscore:
            return "Better"

        elif row["hyperscore"] < ref_hyperscore:
            return "Worse"

        else:
            return "Error?"

    except:
        return "Unpredicted"
    
def main():
    filenames = [
        # 'F01_Fraction2',
        # 'F01_Fraction4',
        # 'S14_Rep2',
        # 'S14_Rep1',
        # 'S08',
        # 'S03',
        # 'S14_Rep3',
        # 'F07_Fraction4',
        # 'S11_Fraction3',
        # 'S11_Fraction1',
        # 'S07',
        # 'F07_Fraction3',
        # 'F08_Rep2',
        # 'F07_Fraction2',
        # 'F07_Fraction1',
        # 'S11_Fraction2',
        # 'S11_Fraction4',
        # 'F08_Rep1',
        # 'F01_Fraction1',
        # 'F01_Fraction3',
        # 'F08_Rep3',
        'F06',
        #'S05'
    ]

    keep_cols = [
        "proforma",
        "sequence",
        "spectrum_id",
        "run",
        "engine",
        "score",
        "qvalue",
        "spectralis_score",
        "hyperscore",
        "is_decoy",
        "has_modification",
        "rescoring_features"
    ]

    for filename in filenames:
        print(filename)
        root_data="/home/samva/Doctorate/data_directory/denovo_project"
        mgf_path=os.path.join(root_data, "mgf_filtered", filename + ".mgf")
        results_dir=os.path.join(root_data, "denovo_results")

        parser_spectralis = SpectralisParser(
            mgf_path=mgf_path,
            results_dir=results_dir
        )

        # Casanovo, instanovo, pepnet, contranovo ran together
        parser_spectralis.parse(
            path_spectralis=os.path.join(
                results_dir,
                "refinement/spectralis/pt1", filename + "_annotated_rescoring.csv"
            )
        )

        # NovoB, Novor, PepNovo+ ran together
        parser_spectralis.parse(
            path_spectralis=os.path.join(
                results_dir,
                "refinement/spectralis/pt2", filename + "_annotated_rescoring.csv"
            )
        )

        # Sage results ran separately
        parser_spectralis.parse(
            path_spectralis=os.path.join(
                results_dir,
                "refinement/spectralis/pt3", filename + "_annotated_rescoring.csv"
            )
        )

        psmlist = parser_spectralis.psmlist
        psmlist["run"] = [filename]*len(psmlist)
        decoy_status = psmlist["is_decoy"] 
        decoy_status = np.where(decoy_status == None, False, decoy_status)
        psmlist["is_decoy"] = decoy_status 

        psmlist["qvalue"] = [1 if x is None else x for x in psmlist["qvalue"]]
        spectrum_ids_to_keep = psmlist[
            (psmlist["source"]=="sage") &
            (psmlist["qvalue"]<.01)
        ]["spectrum_id"]

        mgf_file = mgf.read(mgf_path)

        basic_fgen = BasicFeatureGenerator()

        basic_fgen.add_features(psmlist)
        print("Added basic features.")
        # ms2pip_fgen.add_features(psmlist)
        # print("Added MS2PIP features.")
        # deeplc_fgen.add_features(psmlist)
        # print("Added DeepLC features.")

        for psm in tqdm(psmlist):
            hyperscore = calculate_hyperscore(
                psm=psm,
                mgf_file=mgf_file,
                engine="pyopenms"
            )
            psm["rescoring_features"].update(
                {"hyperscore": hyperscore}
            )
        print("Added hyperscore.")

        df = psmlist.to_dataframe()
        df = df[df.spectrum_id.isin(spectrum_ids_to_keep)].reset_index(drop=True)
        df["spectralis_score"] = df.apply(get_spectralis_score, axis=1)
        df["hyperscore"] = pd.DataFrame(df["rescoring_features"].tolist())["hyperscore"]
        df["psm_type"] = df.apply(get_psm_type, axis=1)
        df["proforma"] = df.peptidoform.apply(lambda x: x.proforma)
        df["peptide"] = df.peptidoform.apply(lambda x: x.sequence.replace("I", "L"))
        df["has_modification"] = df.peptidoform.apply(peptidoform_has_modification)

        ground_truth_peptide = df.loc[df.source=="sage", ["spectrum_id", "peptide"]].set_index("spectrum_id").to_dict()["peptide"]
        ground_truth_hyperscore = df.loc[df.source=="sage", ["spectrum_id", "hyperscore"]].set_index("spectrum_id").to_dict()["hyperscore"]

        df["hyperscore_diff"] = df.apply(lambda x: hyperscore_difference(x, reference=ground_truth_hyperscore), axis=1)
        df["match_type"] = df.progress_apply(
            lambda x: evaluate_prediction_isobaricity(
                x, 
                ground_truth_peptide=ground_truth_peptide, 
                ground_truth_hyperscore=ground_truth_hyperscore
            ), axis=1
        )

        df.to_pickle(os.path.join("/home/samva/Doctorate/DeNovo_Benchmark/notebooks/analysis/filtered_results", filename+".pkl"))

if __name__ == "__main__":
    main()