"""psm_utils-based search result parser."""

from psm_utils import PSMList
from psm_utils.io import read_file


def psmutils_parser(
    result_path: str,
    mgf_path: str,
    mapping: dict,
    max_length=30,
    label: str = "infer",
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

    new_psm_list = []
    for psm in psmlist:
        if len(psm["peptidoform"]) <= max_length:
            new_psm_list.append(psm)

    return PSMList(psm_list=new_psm_list)
