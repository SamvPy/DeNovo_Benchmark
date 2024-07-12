from .casanovo import casanovo_parser
from .instanovo import instanovo_parser
from .contranovo import contranovo_parser
from .novob import novob_parser
from .pepnet import pepnet_parser
from .novor import novor_parser
from .pepnovo import pepnovo_parser
from .pointnovo import pointnovo_parser

from enum import Enum
from psm_utils import PSMList

from ..constants import MODIFICATION_MAPPING
from ..exceptions import DenovoEngineNotSupported

# Define supported parsers for de novo search engines as an Enum with associated parser functions
class DenovoEngineConverter(Enum):
    """
    Examples
    --------
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

    def __init__(self, label, parser_func):
        self.label = label
        self.parser_func = parser_func

    def parse(self, result_path: str, mgf_path: str) -> PSMList:
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
            A PSMList (https://psm-utils.readthedocs.io/en/latest/api/psm_utils/#psm_utils.PSMList)
        """
        return self.parser_func(result_path, mgf_path, MODIFICATION_MAPPING[self.label])

    @classmethod
    def select(cls, label: str):
        """
        Select a de novo search engine result parser by its label.

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
        
        Returns
        -------
        DenovoEngine
            The corresponding `DenovoEngine` enum member.
        
        Raises
        ------
        DenovoEngineNotSupported
            If the provided label does not match any known de novo search engine.
        
        Examples
        --------
        >>> parser = DenovoEngineConverter.select('casanovo')
        >>> psmlist = parser.parse('result.mztab', 'file.mgf')
        """
        for engine in cls:
            if engine.label == label:
                return engine
        raise DenovoEngineNotSupported(label, [e.label for e in cls])