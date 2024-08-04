"""Class which collects all search result parsers."""

from enum import Enum

from psm_utils import PSMList
from psm_utils.io import FILETYPES

from ..constants import MODIFICATION_MAPPING
from ..exceptions import DenovoEngineNotSupported
from .casanovo import casanovo_parser
from .contranovo import contranovo_parser
from .external import psmutils_parser
from .instanovo import instanovo_parser
from .novob import novob_parser
from .novor import novor_parser
from .pepnet import pepnet_parser
from .pepnovo import pepnovo_parser
from .pointnovo import pointnovo_parser


# Define supported parsers for de novo search engines as
# an Enum with associated parser functions
class DenovoEngineConverter(Enum):
    """
    Enum holding all search result parsers.

    Each parser accepts as input the search results and mgf-file and
    outputs a PSMList. This class centralizes all these disparate parsers.

    When using this `Enum`, always use the `select` classmethod to access a parser.

    Example
    -------
    >>> parser = DenovoEngine.select('casanovo')
    >>> psmlist = parser.parse('result.mztab', 'file.mgf')
    """

    CASANOVO = ("casanovo", casanovo_parser)
    INSTANOVO = ("instanovo", instanovo_parser)
    CONTRANOVO = ("contranovo", contranovo_parser)
    NOVOB = ("novob", novob_parser)
    PEPNET = ("pepnet", pepnet_parser)
    NOVOR = ("novor", novor_parser)
    PEPNOVO = ("pepnovo", pepnovo_parser)
    POINTNOVO = ("pointnovo", pointnovo_parser)
    PMSUTILS = ("psm-utils", psmutils_parser)

    def __init__(self, label: str, parser_func: callable) -> None:
        """
        Initialize the enum. Should never be done by the user.

        Parameters
        ----------
        label: str
            The name of the parser.
        parser_func: callable
            The parser function.

        Attributes
        ----------
        label: str
            The name of the parser.
        parser_func: callable
            The parser function.
        psm_utils_parser: str | None
            Defines the engine to use within the psm-utils package for result parsing.
        """
        self.label = label
        self.parser_func = parser_func
        self.psm_utils_parser = None

    def parse(self, result_path: str, mgf_path: str, max_length: int = 30) -> PSMList:
        """
        Parse the results from a specified de novo search engine.

        Parameters
        ----------
        result_path : str
            Path to the result file (filetype dependent on search engine).
        mgf_path : str
            Path to the mgf file.

        Returns
        -------
        list
            A PSMList
            https://psm-utils.readthedocs.io/en/latest/api/psm_utils/#psm_utils.PSMList
        """
        if isinstance(self.psm_utils_parser, str):
            inference_label = self.psm_utils_parser
        else:
            inference_label = self.label

        if inference_label not in MODIFICATION_MAPPING.keys():
            print("No modification mapping defined for '{}'.".format(inference_label))
            mapping = {}
        else:
            mapping = MODIFICATION_MAPPING[inference_label]

        return self.parser_func(
            result_path=result_path,
            mgf_path=mgf_path,
            mapping=mapping,
            max_length=max_length,
            label=self.psm_utils_parser,
        )

    @classmethod
    def select(cls, label: str):
        """
        Select a de novo search engine result parser by its label.

        Each converter accepts as input a search result and spectral file (mgf)
        and outputs a `PSMList`.

        Parameters
        ----------
        label : str
            The label of the de novo search engine. Options include:
            - 'casanovo'
            - 'instanovo'
            - 'contranovo'
            - 'novob'
            - 'pepnet'
            - 'novor'
            - 'pepnovo'
            - 'pointnovo'
            - 'psm-utils' or any engine supported within psm-utils, e.g. sage.

        Returns
        -------
        DenovoEngine
            The corresponding `DenovoEngineConverter` enum member.

        Raises
        ------
        DenovoEngineNotSupported
            If the provided label does not match any known de novo search engine.

        Examples
        --------
        >>> parser = DenovoEngineConverter.select('casanovo')
        >>> psmlist = parser.parse('result.mztab', 'file.mgf')
        """
        # The label should be different according to support psm-utils reader
        if label in FILETYPES.keys():
            for engine in cls:
                if engine.label == "psm-utils":
                    engine.psm_utils_parser = label
                    return engine
            raise Exception("Unexpected error.")

        for engine in cls:
            if engine.label == label:
                return engine

        supported_filetypes = [e.label for e in cls] + list(FILETYPES.keys())
        raise DenovoEngineNotSupported(label, supported_filetypes)
