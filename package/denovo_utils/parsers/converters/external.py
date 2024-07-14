from psm_utils.io import read_file

def psmutils_parser(result_path: str, mgf_path:str, mapping=None, max_length=None):
    return read_file(
        result_path,
        filetype="infer"
    )