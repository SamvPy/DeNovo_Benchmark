# import pyopenms as oms
# from pyteomics import mgf
# import os
# import pandas as pd
# from pyteomics.mass import calculate_mass
# from psm_utils import PSM, PSMList

# from spectrum_utils.proforma import Proteoform, Modification
# from spectrum_utils.fragment_annotation import get_theoretical_fragments
# from psm_utils import Peptidoform
# from spectrum_utils import proforma

# from denovo_utils.analysis import calculate_hyperscore
# from denovo_utils.parsers import proforma_to_theoretical_spectrum
# from denovo_utils.parsers.converters import SpectralisParser, DenovoEngineConverter
# from ms2rescore.feature_generators import BasicFeatureGenerator, MS2PIPFeatureGenerator, DeepLCFeatureGenerator

# import deeplc
# from deeplc.plot import scatter
# from ms2rescore.report.charts import (
#     calculate_feature_qvalues,
#     feature_ecdf_auc_bar,
#     fdr_plot,
#     ms2pip_correlation,
# )

# from tqdm import tqdm

# import spectrum_utils.plot as sup

# from denovo_utils.utils.pandas import get_psm_type, get_spectralis_score
# import seaborn as sns
# from matplotlib import pyplot as plt
# import numpy as np


# def peptidoform_has_modification(peptidoform, allowed_modificiations=["[UNIMOD:4]", "[UNIMOD:35]"]):
#     peptidoform_str = peptidoform.proforma
#     for allowed_modification in allowed_modificiations:
#         peptidoform_str = peptidoform_str.replace(allowed_modification, "")
    
#     return "[" in peptidoform_str

# def hyperscore_difference(row, reference):
#     try:
#         reference_hyperscore = reference[row["spectrum_id"]]
#         return float(row["hyperscore"])-float(reference_hyperscore)
#     except:
#         return None
    
# def evaluate_prediction_isobaricity(
#         row, ground_truth_peptide, ground_truth_hyperscore
# ):
#     try:
#         sequence_match = ground_truth_peptide[
#             row["spectrum_id"]
#         ] == row["peptide"]
#         if sequence_match:
#             return "Match"
        
#         ref_hyperscore = ground_truth_hyperscore[row["spectrum_id"]]
#         if row["hyperscore"] == ref_hyperscore:
#             return "Isobaric"
        
#         elif row["hyperscore"] > ref_hyperscore:
#             return "Better"

#         elif row["hyperscore"] < ref_hyperscore:
#             return "Worse"

#         else:
#             return "Error?"

#     except:
#         return "Unpredicted"
    
# def main():
#     filenames = [
#         # 'F01_Fraction2',
#         # 'F01_Fraction4',
#         # 'S14_Rep2',
#         # 'S14_Rep1',
#         # 'S08',
#         # 'S03',
#         # 'S14_Rep3',
#         # 'F07_Fraction4',
#         # 'S11_Fraction3',
#         # 'S11_Fraction1',
#         # 'S07',
#         # 'F07_Fraction3',
#         # 'F08_Rep2',
#         # 'F07_Fraction2',
#         # 'F07_Fraction1',
#         # 'S11_Fraction2',
#         # 'S11_Fraction4',
#         # 'F08_Rep1',
#         # 'F01_Fraction1',
#         # 'F01_Fraction3',
#         # 'F08_Rep3',
#         'F06',
#         #'S05'
#     ]

#     keep_cols = [
#         "proforma",
#         "sequence",
#         "spectrum_id",
#         "run",
#         "engine",
#         "score",
#         "qvalue",
#         "spectralis_score",
#         "hyperscore",
#         "is_decoy",
#         "has_modification",
#         "rescoring_features"
#     ]

#     for filename in filenames:
#         print(filename)
#         root_data="/home/samva/Doctorate/data_directory/denovo_project"
#         mgf_path=os.path.join(root_data, "mgf_filtered", filename + ".mgf")
#         results_dir=os.path.join(root_data, "denovo_results")

#         parser_spectralis = SpectralisParser(
#             mgf_path=mgf_path,
#             results_dir=results_dir
#         )

#         # Casanovo, instanovo, pepnet, contranovo ran together
#         parser_spectralis.parse(
#             path_spectralis=os.path.join(
#                 results_dir,
#                 "refinement/spectralis/pt1", filename + "_annotated_rescoring.csv"
#             )
#         )

#         # NovoB, Novor, PepNovo+ ran together
#         parser_spectralis.parse(
#             path_spectralis=os.path.join(
#                 results_dir,
#                 "refinement/spectralis/pt2", filename + "_annotated_rescoring.csv"
#             )
#         )

#         # Sage results ran separately
#         parser_spectralis.parse(
#             path_spectralis=os.path.join(
#                 results_dir,
#                 "refinement/spectralis/pt3", filename + "_annotated_rescoring.csv"
#             )
#         )

#         psmlist = parser_spectralis.psmlist
#         psmlist["run"] = [filename]*len(psmlist)
#         decoy_status = psmlist["is_decoy"] 
#         decoy_status = np.where(decoy_status == None, False, decoy_status)
#         psmlist["is_decoy"] = decoy_status 

#         psmlist["qvalue"] = [1 if x is None else x for x in psmlist["qvalue"]]
#         spectrum_ids_to_keep = psmlist[
#             (psmlist["source"]=="sage") &
#             (psmlist["qvalue"]<.01)
#         ]["spectrum_id"]

#         mgf_file = mgf.read(mgf_path)

#         basic_fgen = BasicFeatureGenerator()

#         basic_fgen.add_features(psmlist)
#         print("Added basic features.")
#         # ms2pip_fgen.add_features(psmlist)
#         # print("Added MS2PIP features.")
#         # deeplc_fgen.add_features(psmlist)
#         # print("Added DeepLC features.")

#         for psm in tqdm(psmlist):
#             hyperscore = calculate_hyperscore(
#                 psm=psm,
#                 mgf_file=mgf_file,
#                 engine="pyopenms"
#             )
#             psm["rescoring_features"].update(
#                 {"hyperscore": hyperscore}
#             )
#         print("Added hyperscore.")

#         df = psmlist.to_dataframe()
#         df = df[df.spectrum_id.isin(spectrum_ids_to_keep)].reset_index(drop=True)
#         df["spectralis_score"] = df.apply(get_spectralis_score, axis=1)
#         df["hyperscore"] = pd.DataFrame(df["rescoring_features"].tolist())["hyperscore"]
#         df["psm_type"] = df.apply(get_psm_type, axis=1)
#         df["proforma"] = df.peptidoform.apply(lambda x: x.proforma)
#         df["peptide"] = df.peptidoform.apply(lambda x: x.sequence.replace("I", "L"))
#         df["has_modification"] = df.peptidoform.apply(peptidoform_has_modification)

#         ground_truth_peptide = df.loc[df.source=="sage", ["spectrum_id", "peptide"]].set_index("spectrum_id").to_dict()["peptide"]
#         ground_truth_hyperscore = df.loc[df.source=="sage", ["spectrum_id", "hyperscore"]].set_index("spectrum_id").to_dict()["hyperscore"]

#         df["hyperscore_diff"] = df.apply(lambda x: hyperscore_difference(x, reference=ground_truth_hyperscore), axis=1)
#         df["match_type"] = df.progress_apply(
#             lambda x: evaluate_prediction_isobaricity(
#                 x, 
#                 ground_truth_peptide=ground_truth_peptide, 
#                 ground_truth_hyperscore=ground_truth_hyperscore
#             ), axis=1
#         )

#         df.to_pickle(os.path.join("/home/samva/Doctorate/DeNovo_Benchmark/notebooks/analysis/filtered_results", filename+".pkl"))

# if __name__ == "__main__":
#     main()

from denovo_utils.parsers import DenovoEngineConverter
from denovo_utils.utils.pandas import collapse_casanovo_score
from denovo_utils.parsers.constants import EXTENSIONS
from denovo_utils.analysis.feature_generation import calculate_hyperscore
from denovo_utils.analysis.missing_fragmentations import MissingFragmentationSiteFGen
from spectrum_utils.utils import mass_diff
from pyteomics import mgf
import pandas as pd
import os
from glob import glob
from tqdm import tqdm
import numpy as np

from denovo_utils.analysis.missing_fragmentations import get_annotated_spectrum
from denovo_utils.utils.annotation import AnnotatedSpectrum

def psm_to_annotated_spectrum(psm, mgf_file):
    annotated_spec = get_annotated_spectrum(
        mgf_file,
        psm,
        ions="byp",
        neutral_losses=True
    )
    return AnnotatedSpectrum(
        annotated_spec,
        peplen=len(psm["peptidoform"]),
        process=True,
        add_neutral_losses=["-H2O", "-NH3"]
    )


def psmlist_to_annotated_spectra(psmlist, mgf_file, out_path, save=True):
    spec_dict = {}
    for i, psm in enumerate(psmlist):
        spec_dict[i] = {
            "spectrum_id": psm["spectrum_id"],
            "annotated_spectrum": psm_to_annotated_spectrum(
                psm=psm,
                mgf_file=mgf_file
            ),
            "run": psm["run"],
            "source": psm["source"]
        }
    df = pd.DataFrame(spec_dict)

    if save:
        df.to_feather(out_path)
    return df


def calculate_ppms(ion_mz_delta):
    y_ppm = []
    for y in ion_mz_delta["y1"]:
        if y != 0:
            y_ppm.append(y)
    for y in ion_mz_delta["y2"]:
        if y != 0:
            y_ppm.append(y)

    b_ppm = []
    for b in ion_mz_delta["b1"]:
        if b != 0:
            b_ppm.append(b)
    for b in ion_mz_delta["b2"]:
        if b != 0:
            b_ppm.append(b)

    if y_ppm:
        y_ppm_mean = np.mean(y_ppm)
    else:
        y_ppm_mean = 0.
    
    if b_ppm:
        b_ppm_mean = np.mean(b_ppm)
    else:
        b_ppm_mean = 0.
    return b_ppm_mean, y_ppm_mean

def calculate_pct_intensities(tic, isotope_matrix):
    y_intensity = np.sum(isotope_matrix["y1"]) + np.sum(isotope_matrix["y2"])
    b_intensity = np.sum(isotope_matrix["b1"]) + np.sum(isotope_matrix["b2"])
    p_intensity = np.sum(isotope_matrix["p"])

    return {
        "y_pct": y_intensity/tic,
        "b_pct": b_intensity/tic,
        "p_pct": p_intensity/tic,
        "fragmented_pct": (
            y_intensity + b_intensity
        ) / (y_intensity + b_intensity + p_intensity)
    }

def calculate_peak_features(annotated_spectrum):
    n_peaks = len(annotated_spectrum.spectrum.mz)

    annotated_peaks = 0
    for ion_type in ["b1", "b2", "y1", "y2", "p"]:
        annotated_peaks += np.sum(annotated_spectrum.isotope_matrix[ion_type]>0)
    
    return {
        "n_peaks": n_peaks,
        "pct_peaks_annotated": annotated_peaks/n_peaks
    }
    
def add_features(psm, missfrag_fgen):
    metadata, features = missfrag_fgen._calculate_features(psm=psm, ions="byp")
    features["hyperscore"] = calculate_hyperscore(
        psm=psm,
        mgf_file=missfrag_fgen.mgf_file,
        engine="pyopenms"
    )
    features["tic"] = metadata["annotated_spectrum"].tic
    features["precursor_err_ppm"] = mass_diff(
        psm.peptidoform.theoretical_mz,
        missfrag_fgen.mgf_file.get_by_id(psm["spectrum_id"])["params"]["pepmass"][0],
        mode_is_da=False
    )
    features["precursor_err_da"] = mass_diff(
        psm.peptidoform.theoretical_mz,
        missfrag_fgen.mgf_file.get_by_id(psm["spectrum_id"])["params"]["pepmass"][0],
        mode_is_da=True
    )
    features["ppm_b_mean"], features["ppm_y_mean"] = calculate_ppms(metadata["annotated_spectrum"].ion_mz_delta)
    
    intensity_features = calculate_pct_intensities(
        metadata["annotated_spectrum"].tic, 
        metadata["annotated_spectrum"].isotope_matrix
    )
    features.update(intensity_features)

    features["precursor_peak_present"] = np.sum(metadata["annotated_spectrum"].isotope_matrix["p"]) > 0

    peak_features = calculate_peak_features(metadata["annotated_spectrum"])
    features.update(peak_features)

    psm["rescoring_features"].update(features)

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
        ] == row["sequence"]
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

def merge_results(
        root_ground_truth,
        root_denovo_results,
        root_mgf,
        root_out,
        engines,
        ground_truth_filetype="sage",
        only_identified_spectra=True
    ):
    for mgf_path in glob(os.path.join(root_mgf, "*.mgf")):
        print(mgf_path)
        missfrag_fgen = MissingFragmentationSiteFGen(spectrum_path=mgf_path, only_by_ions=False)

        df_list = []

        filename = os.path.basename(mgf_path).split(".")[0]

        parser = DenovoEngineConverter.select(ground_truth_filetype)
        ground_truth_psmlist = parser.parse(
            os.path.join(root_ground_truth, filename + EXTENSIONS[ground_truth_filetype]),
            mgf_path=mgf_path
        )

        if only_identified_spectra:
            ground_truth_psmlist = ground_truth_psmlist[
                (ground_truth_psmlist["qvalue"] < .01) &
                (~ground_truth_psmlist["is_decoy"])
            ]
            spectrum_ids_identified = ground_truth_psmlist["spectrum_id"].tolist()

        for psm in tqdm(ground_truth_psmlist):
            add_features(psm, missfrag_fgen)

        ground_truth_df = ground_truth_psmlist.to_dataframe()
        ground_truth_df["hyperscore"] = ground_truth_df["rescoring_features"].apply(lambda x: x["hyperscore"])
        ground_truth_df["proforma"] = ground_truth_df.peptidoform.apply(lambda x: x.proforma
                                                              .replace("I", "L")
                                                              .replace("N[UNIMOD:7]", "D")
                                                              .replace("Q[UNIMOD:7]", "E")
                                                            )
        ground_truth_df["sequence"] = ground_truth_df.peptidoform.apply(lambda x: x.sequence.replace("I", "L")
                                                                                  .replace("N[UNIMOD:7]", "D")
                                                                                  .replace("Q[UNIMOD:7]", "E"))

        ground_truth_sequence = ground_truth_df.loc[
            ground_truth_df.source=="sage",
            ["spectrum_id", "sequence"]
        ].set_index("spectrum_id").to_dict()["sequence"]

        ground_truth_hyperscore = ground_truth_df.loc[
            ground_truth_df.source=="sage",
            ["spectrum_id", "hyperscore"]
        ].set_index("spectrum_id").to_dict()["hyperscore"]

        
        df_list.append(ground_truth_df)

        for engine in engines:
            print(f"\t{engine}")
            parser = DenovoEngineConverter.select(engine)
            psmlist = parser.parse(
                os.path.join(root_denovo_results, engine, f"{filename}.{engine}.csv"),
                mgf_path=mgf_path
            )

            for psm in tqdm(psmlist):
                add_features(psm, missfrag_fgen)

            df = psmlist.to_dataframe()
    
            if only_identified_spectra:
                df = df[df.spectrum_id.isin(spectrum_ids_identified)]
            
            df_list.append(df)
        
        df_results = pd.concat(df_list, ignore_index=True)
        df_results["proforma"] = df_results.peptidoform.apply(lambda x: x.proforma
                                                              .replace("I", "L")
                                                              .replace("N[UNIMOD:7]", "D")
                                                              .replace("Q[UNIMOD:7]", "E")
                                                            )
        df_results["sequence"] = df_results.peptidoform.apply(lambda x: x.sequence.replace("I", "L")
                                                                                  .replace("N[UNIMOD:7]", "D")
                                                                                  .replace("Q[UNIMOD:7]", "E"))
        df_results["run"] = filename
        df_results["score"] = df_results.apply(collapse_casanovo_score, axis=1)
        
        df_results["hyperscore_diff"] = df_results.apply(
            lambda x: hyperscore_difference(x, ground_truth_hyperscore), axis=1
        )
        df_results["match_type"] = df_results.apply(
            lambda x: evaluate_prediction_isobaricity(
                x,
                ground_truth_peptide=ground_truth_sequence,
                ground_truth_hyperscore=ground_truth_hyperscore
            ), axis=1
        )


        df_results = df_results.loc[
            :,
            [
                "spectrum_id",
                "proforma",
                "sequence",
                "match_type",
                "run",
                "source",
                "score",
                "hyperscore",
                "hyperscore_diff",
                "qvalue",
                "is_decoy",
                "rescoring_features",
                "precursor_mz",
                "retention_time",
                "protein_list"
            ]
        ]
        df_results.to_pickle(os.path.join(root_out, filename+".pkl"))


if __name__ == "__main__":
    root_ground_truth = "/home/samva/Doctorate/data_directory/PXD028735/search_results/identification"
    root_denovo_results = "/home/samva/Doctorate/data_directory/PXD028735/denovo_results"
    root_mgf = "/home/samva/Doctorate/data_directory/PXD028735/mgf/Orbitrap_QE/reformatted"
    root_out = "/home/samva/Doctorate/data_directory/PXD028735/denovo_results/merged2"
    engines = [
        "casanovo",
        "instanovo",
        "contranovo",
        "pepnet",
        "novob"
    ]
    psmlist = merge_results(
        root_ground_truth=root_ground_truth,
        root_denovo_results=root_denovo_results,
        root_mgf=root_mgf,
        root_out=root_out,
        engines=engines,
        only_identified_spectra=False
    )