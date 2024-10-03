import logging
import os
import re
from copy import deepcopy
from typing import Dict

import pandas as pd
from psm_utils import PSMList
from pyteomics import mgf
from tqdm import tqdm


from ..converters import DenovoEngineConverter
from ..exceptions import NoResultsToMergeException, SeparatorCharacterInTitle

logger = logging.getLogger(__name__)
logging.basicConfig(filename="denovo_writer.log", level=logging.INFO)


# TODO: Split up the MGFWriter into a reader and a writer to separate functionalities more
class MGFWriter:
    def __init__(
        self,
        mgf_path: str,
        modification_mapping: dict = {},
        exclusion_list: list = [],
        separator: str = "|",
        output_folder: str = "",
    ) -> None:
        self.mgf_path = mgf_path
        self.filename = os.path.basename(mgf_path).split(".")[0]
        self.results: Dict[str, PSMList] = {}
        self.parsed = False
        self.separator = separator
        self.modification_mapping = modification_mapping
        self.exclusion_list = exclusion_list
        self.output_folder = output_folder
        assert os.path.isdir(self.output_folder)

    def _parse_results_merged(self, df: pd.DataFrame) -> dict:
        if not hasattr(self, "results_merged"):
            AttributeError(
                "No attribute 'results_merged' initialized. First read in de novo results and merge them with the appropriate methods!"
            )

        df["peptidoform"] = df["peptidoform"].apply(lambda x: x.proforma)
        df = (
            df[["spectrum_id", "peptidoform", "source"]]
            .groupby(by=["spectrum_id", "peptidoform"])
            .apply(lambda x: list(x["source"]))
            .reset_index()
            .rename(columns={0: "search_engine"})
        )

        if len(df[df["spectrum_id"].apply(lambda x: "|" in x)]) != 0:
            titles_with_char = df[df["spectrum_id"].apply(lambda x: "|" in x)][
                "spectrum_id"
            ]
            raise SeparatorCharacterInTitle(
                char=self.separator, title=titles_with_char.iloc[0]
            )

        df["title"] = df.apply(
            lambda x: x["spectrum_id"]
            + 2 * self.separator
            + self.separator.join(x["search_engine"]),
            axis=1,
        )

        return df.loc[:, ["spectrum_id", "title", "peptidoform"]]

    def _index_spectra(self, df: pd.DataFrame, col: str):
        """
        Given a dataframe with duplicate entries in 'col',
        create a mapping from 'col' to indexes in the dataframe.

        During this process, the dataframe (df) is reindexed and should be overwritten.

        Parameters
        ----------
        df: pd.DataFrame
            The dataframe used to create an index dictionary
        col: str
            The column in the dataframe to map indices with

        Returns
        -------
        tuple
            dict:
                A mapping of entries in 'col' with the dataframe indices
            pd.DataFrame:
                The re-indexed dataframe (df)
        """
        df = df.reset_index(drop=True)
        mapping = {}

        for k, v in df[col].to_dict().items():
            if v not in mapping.keys():
                mapping[v] = []
            mapping[v].append(k)
        return mapping, df

    def read_result(self, result_path: str, denovo_engine: str) -> None:
        parser = DenovoEngineConverter.select(denovo_engine)
        psmlist = parser.parse(result_path=result_path, mgf_path=self.mgf_path)
        self.results[denovo_engine] = psmlist

    def merge(self, parse=True):
        """
        Merge the results of the de novo search engines and parse for writing if required.

        Parsing will remove all non-required information.
        """

        if len(self.results) == 0:
            raise NoResultsToMergeException()

        results_merged = pd.concat(
            [psm_df.to_dataframe() for psm_df in self.results.values()],
            ignore_index=True,
        )

        if not parse:
            self.results_merged = results_merged
            self.parsed = False

        results_merged = self._parse_results_merged(df=results_merged)
        self.mapping, self.results_merged = self._index_spectra(
            df=results_merged, col="spectrum_id"
        )
        self.parsed = True

    @property
    def mgf_file(self):
        return mgf.read(self.mgf_path)

    def write(self):

        if not hasattr(self, "results") and not self.parsed:
            raise AttributeError(
                "First read in results and merge/parse them with appropriate methods!"
            )

        mgf_annotated = []
        for spectrum in self.mgf_file:
            try:
                indices = self.mapping[spectrum["params"]["title"]]
            except KeyError:
                continue

            for i in indices:

                spectrum_new = deepcopy(spectrum)

                entry = self.results_merged.loc[i, :]
                spectrum_new = proforma_to_spectralis(
                    proforma=entry["peptidoform"],
                    spectrum=spectrum_new,
                    modification_mapping=self.modification_mapping,
                    exclusion_list=self.exclusion_list,
                )

                # Spectralis uses the scans as output identifier, so this should be unique.
                spectrum_new["params"]["scans"] = entry["title"]
                mgf_annotated.append(spectrum_new)

        mgf.write(
            spectra=mgf_annotated,
            output=os.path.join(self.output_folder, self.filename + "_annotated.mgf"),
            header=self.mgf_file.header,
        )


def parse_title(spectrum: dict):
    title_pattern = r"controllerType=\d+ controllerNumber=\d+ scan=(\d+)"
    try:
        title = spectrum["params"]["title"]
    except KeyError:
        return False, False
    
    match = re.search(title_pattern, title)
    if match:
        title_id = match.group(0)
        scan_number = match.group(1)
        return title_id, scan_number
    else:
        return False, False


def parse_charge(spectrum: dict):

    try:
        charge = spectrum["params"]["charge"]
    except KeyError:
        return False
    
    if len(charge)>1:
        logging.warning(f"{spectrum['params']['title']}: Has multiple charges ({charge})")
    return int(charge[0])


def parse_pepmass(spectrum: dict, silence_warning=True):

    try:
        pepmass = spectrum["params"]["pepmass"]
        if len(pepmass)>1 and not silence_warning:
            logging.warning(f"{spectrum['params']['title']}: Pepmass field holds multiple values. ({pepmass})")
        return pepmass
    except KeyError:
        raise Exception("No pepmass field in the precursor info of MGF spectrum.")


def check_required_fields(spectrum: dict):

    start_new_index = False

    # Check title field
    title_id, _ = parse_title(spectrum)
    if title_id == False:
        start_new_index = True
    
    # Check charge field
    _ = parse_charge(spectrum)

    # Check pepmass
    _ = parse_pepmass(spectrum, silence_warning=False)

    return start_new_index


def reformat_mgf(
        path_mgf,
        output_dir="",
        min_peak_length=15,
        max_charge=8
    ):

    mgf_file = mgf.read(path_mgf)
    filename = os.path.basename(path_mgf)
    filename_out = filename.replace(".mgf", "") + "_reformatted.mgf"
    logging.info(f"\n\nParsing {filename}")

    total_spectra = len(mgf_file)
    start_new_index = check_required_fields(spectrum=mgf_file[0])

    spectrum_list = []

    skipped_missing_charge = 0
    skipped_above_max_charge = 0
    skipped_by_peak_length = 0

    for i, spectrum in enumerate(mgf_file):
        
        if len(spectrum["m/z array"]) < min_peak_length:
            skipped_by_peak_length += 1
            continue

        if start_new_index:
            title, scan_id = f"controllerType=0 controllerNumber=1 scan={i}", i
        else:
            title, scan_id = parse_title(spectrum=spectrum)
        
        charge = parse_charge(spectrum=spectrum)
        if charge==False:
            skipped_missing_charge += 1
            continue

        if charge > max_charge:
            skipped_above_max_charge += 1
            continue

        pepmass = parse_pepmass(spectrum=spectrum)
        rtinseconds = spectrum["params"]["rtinseconds"]
            
        params = {
            "title": title,
            "pepmass": pepmass,
            "rtinseconds": rtinseconds,
            "charge": str(charge)+"+",
            "scans": scan_id,
        }
        spectrum["params"] = params
        spectrum_list.append(spectrum)

    if spectrum_list == []:
        logging.warning(f"{filename}: No spectra left after precursor info validation and filtering.")

    logging.info(
        f"\nTotal spectra: {total_spectra}; After filtering: {len(spectrum_list)}\n\t"
        f"Filtered spectra: {total_spectra-len(spectrum_list)}\n\t" \
        f"Peaks below{min_peak_length}: {skipped_by_peak_length}\n\t" \
        f"No charge field: {skipped_missing_charge}\n\t" \
        f"Above max charge ({max_charge}): {skipped_above_max_charge}"
    )

    if (skipped_missing_charge+skipped_by_peak_length) == total_spectra:
        raise Exception("Use a charge predictor for annotation! All MGFs have missing charge field.")

    mgf.write(
        spectra=spectrum_list,
        output=os.path.join(output_dir, filename_out),
        header=mgf_file.header
    )


def proforma_to_spectralis(
    proforma: str, spectrum: dict, modification_mapping: dict = [], exclusion_list=[]
) -> dict:

    annotated_peptide = proforma
    for k, v in modification_mapping.items():
        annotated_peptide = annotated_peptide.replace(k, v)

    for exclude in exclusion_list:
        if exclude in annotated_peptide:
            continue

    spectrum["params"]["seq"] = annotated_peptide.split("/")[0]
    return spectrum


def psmlist_to_mgf(
    psmlist: PSMList,
    mgf_path: str,
    output_folder: str,
    modification_mapping: dict = [],
    exclusion_list: dict = [],
    annotated=True,
):
    """
    From a psm-utils PSMList with annotated sequences, write results to the accomponaying MGF-file.

    The generated MGF-file is compatible as Spectralis input.

    Parameters
    ----------
    psmlist: PSMList
        The PSMs that need to be written to an MGF file.
    mgf_path: str
        The mgf-path to the corresponding raw data.
    output_folder: str
        A path to a folder where the MGF file will be stored.
    modification_mapping: dict
        A modification mapper defined in denovo_utils.constants. Maps modifications to the desired format.
    exclusion_list: dict
        A list of strings which, if in the peptide sequence, prompts the writer to drop the PSM.
    annotated:
        Whether to write the peptide sequences in the MGF-file as
    """

    filename = os.path.basename(mgf_path).split(".")[0]
    mgf_file = mgf.read(mgf_path)
    annotation_dict = (
        psmlist.to_dataframe()
        .loc[:, ["peptidoform", "spectrum_id"]]
        .set_index("spectrum_id")
        .to_dict("index")
    )

    mgf_annotated = []

    for spectrum in mgf_file:

        # When the peptide prediction didnt fit the requirements, skip the spectrum.
        try:
            annotated_peptide = annotation_dict[spectrum["params"]["title"]][
                "peptidoform"
            ].proforma
        except:
            continue

        spectrum = proforma_to_spectralis(
            proforma=annotated_peptide,
            spectrum=spectrum,
            modification_mapping=modification_mapping,
            exclusion_list=exclusion_list,
        )
        mgf_annotated.append(spectrum)

    mgf.write(
        mgf_annotated,
        output=os.path.join(output_folder, filename) + "_annotated.mgf",
        header=mgf_file.header,
    )


def filter_on_charge(mgf_file_path, output_dir, max_charge):

    filename = os.path.basename(mgf_file_path)
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
        output=os.path.join(output_dir, filename),
        header=mgf_file.header,
    )

    logging.info(
        f"{mgf_file_path}: {i} -> {len(mgf_new)} spectra. (-{i+1-len(mgf_new)})"
    )
