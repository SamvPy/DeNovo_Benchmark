from psm_utils.io import read_file
from psm_utils import PSMList

def psmutils_parser(result_path: str, mgf_path:str, mapping: dict, max_length=30, label:str="infer", **kwargs):
    psmlist = read_file(
        result_path,
        filetype=label
    )
    psmlist.rename_modifications(mapping)

    new_psm_list = []
    for psm in psmlist:
        if len(psm["peptidoform"]) <= max_length:
            new_psm_list.append(psm)

    return PSMList(psm_list=new_psm_list)