"""Adanovo parser function."""

from .casanovo import casanovo_parser
import logging
from psm_utils import PSMList

logger = logging.getLogger(__name__)
logging.basicConfig(filename="denovo_output_parsing.log", level=logging.INFO)


def adanovo_parser(
    result_path: str, mgf_path: str, mapping: dict, max_length: int = 30, im: bool=False, **kwargs
) -> PSMList:
    """
    Return a `PSMList` from a casanovo search result and its associated spectral file.

    Parameters
    ----------
    result_path: str
        Path to the search results.
    mgf_path: str
        Path to the MGF-file that was used to generate search results
    mapping: dict
        Mapping dictionary for converting engine-specific to
        unified modification and amino acid labels.
    max_length: int
        Any peptide sequence longer than this value will be ignored.

    Return
    ------
    psm_utils.PSMList
        The PSMList, representing the search results.
    """
    psm_list = casanovo_parser(
        result_path=result_path,
        mgf_path=mgf_path,
        mapping=mapping,
        max_length=max_length,
        im=im
    )
    if len(psm_list) == 0:
        return PSMList(psm_list=[])
    psm_list['source'] = ['AdaNovo'] * len(psm_list)
    return psm_list