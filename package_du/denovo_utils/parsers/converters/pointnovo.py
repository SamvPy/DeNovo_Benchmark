"""PointNovo parser function."""

from psm_utils import PSMList


def pointnovo_parser(
    result_path: str, mgf_path: str, mapping: dict, max_length=30, **kwargs
) -> PSMList:
    """
    Return a `PSMList` from a ContraNovo search result and its MGF spectral file.

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
    raise NotImplementedError("pointnovo_parser not implemented yet!")
