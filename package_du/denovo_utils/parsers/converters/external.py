"""psm_utils-based search result parser."""

from psm_utils import PSMList
from psm_utils.io import read_file

import pandas as pd
from pyteomics import mgf

from .utils import mzml_reader

def psmutils_parser(
    result_path: str,
    mgf_path: str,
    mapping: dict,
    max_length=30,
    label: str = "infer",
    im: bool = False,
    **kwargs
) -> PSMList:
    """
    Return a `PSMList` from any search result supported within the psm_utils package.

    Parameters
    ----------
    result_path: str
        Path to the search results.
    mgf_path: str (unused)
        Path to the MGF-file that was used to generate search results
    mapping: dict
        Mapping dictionary for converting engine-specific to
        unified modification and amino acid labels.
    max_length: int
        Any peptide sequence longer than this value will be ignored.
    label: str
        The engine used within the psm_utils.io.read_file() function.

    Return
    ------
    psm_utils.PSMList
        The PSMList, representing the search results.
    """
    psmlist = read_file(result_path, filetype=label)
    psmlist.rename_modifications(mapping)

    # The precursor_mz could be different
    if mgf_path.lower().endswith('.mzml'):
        mgf_file = mzml_reader(mgf_path, im=im)
        # In mzml format, the indexation with spectra_ref is done with the id (title) part
        mgf_file['index'] = mgf_file['title'].apply(lambda x: int(x.split('scan=')[-1]))
        mgf_file = mgf_file.set_index('index')
    else:
        mgf_file = pd.DataFrame(pd.DataFrame(mgf.read(mgf_path))["params"].tolist())

    mgf_file['precursor_mz'] = mgf_file['pepmass'].apply(lambda x: x[0])
    mgf_file = mgf_file[['title', 'precursor_mz']]
    mgf_file = mgf_file.set_index('title')
    mgf_dict = mgf_file.to_dict('index')

    new_psm_list = []
    for psm in psmlist:
        if len(psm["peptidoform"]) <= max_length:
            psm['precursor_mz'] = mgf_dict[psm['spectrum_id']]['precursor_mz']
            new_psm_list.append(psm)

    return PSMList(psm_list=new_psm_list)
