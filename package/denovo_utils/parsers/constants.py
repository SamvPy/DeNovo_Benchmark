
### GENERAL STUFF

MODIFICATION_MAPPING = {
    "casanovo": {
        "+57.021": "[UNIMOD:4]",
        "+15.995": "[UNIMOD:35]",
        "+0.984": "[UNIMOD:7]", # Check both D - N+0.984 and E - Q+0.984
        "+43.006-17.027": "[+25.980265]-", # [UNIMOD:Carbamylation][UNIMOD:Ammonia-loss]
        "+42.011": "[UNIMOD:1]-",
        "+43.006": "[UNIMOD:5]-",
        "-17.027": "[UNIMOD:385]-"
    },
    "instanovo": {
        "C": "C[UNIMOD:4]",
        "M(ox)": "M[UNIMOD:35]",
        "M(+15.99)": "M[UNIMOD:35]",
        "N(+.98)": "N[UNIMOD:7]",
        "Q(+.98)": "Q[UNIMOD:7]"
    },
    "contranovo": {
        "C+57.021": "C[UNIMOD:4]",
        "+43.006-17.027": "[+25.980265]-",
        "-17.027+43.006": "[+25.980265]-",
        "+42.011": "[UNIMOD:1]-",  # Acetylation
        "+43.006": "[UNIMOD:5]-",  # 5
        "-17.027": "[UNIMOD:385]-",  # NH3 loss
        # AA mods:
        "M+15.995": "M[UNIMOD:35]",  # Met oxidation
        "N+0.984": "N[UNIMOD:7]", # Asn deamidation
        "Q+0.984": "Q[UNIMOD:7]"  # Gln deamidation
    },
    "pepnet": {
        "Z": "M[UNIMOD:35]",
        "C": "C[UNIMOD:4]"
    },
    "novob": {
        "C": "C[UNIMOD:4]",
        "m": "M[UNIMOD:35]",
        "n": "N[UNIMOD:7]",
        "q": "Q[UNIMOD:7]",
        "s": "S[UNIMOD:21]", 
        "t": "T[UNIMOD:21]",
        "y": "Y[UNIMOD:21]"
    },
    "novor": {
        "C(1)": "C[UNIMOD:4]",
        "M(0)": "M[UNIMOD:35]"
    },
    "pepnovo": {
        "M+16": "M[UNIMOD:35]",
        "C": "C[UNIMOD:4]"
    },
    "general": {
        "acetylation": "[UNIMOD:1]"
    }
}

# Generate a comprehensive modification dictionary
def generate_all_modification_labels(mapping: dict) -> dict:
    """
    Generate a comprehensive dictionary of all modification labels.

    Parameters
    ----------
    mapping : dict
        The modification mappings for different engines.

    Returns
    -------
    dict
        A dictionary of all unique modification labels.
    """
    modification_dictionaries = list(mapping.values())
    all_modification_labels = {}

    for d in modification_dictionaries:
        for k, v in d.items():
            if k not in all_modification_labels:
                all_modification_labels[k] = v

    return all_modification_labels

ALL_MODIFICATION_LABELS = generate_all_modification_labels(MODIFICATION_MAPPING)


EXTENSIONS = {
    "casanovo": ".mztab",
    "contranovo": ".mztab",
    "instanovo": ".csv",
    "novob": ".tsv",
    "pepnet": ".tsv",
    "novor": ".csv",
    "pepnovo": ".mgf.out",
    "directag": ".tags"
}


# def generate_spectralis_mod_map(all_labels: dict) -> dict:
#     """
#     Generate a mapping of modifications to Spectralis format.

#     Parameters
#     ----------
#     all_labels : dict
#         A dictionary of all unique modification labels.

#     Returns
#     -------
#     dict
#         A dictionary mapping modifications to Spectralis format.
#     """
#     spectralis_modifications = [
#         "C[UNIMOD:4]", "M[UNIMOD:35]", 'Q[UNIMOD:7]', 'N[UNIMOD:7]'
#     ]

#     spectralis_mod_map = {
#         'Q[UNIMOD:7]': 'E',
#         'N[UNIMOD:7]': 'D'
#     }

#     for modification in all_labels.values():
#         if modification not in spectralis_modifications:
#             spectralis_mod_map[modification] = ""
#     return spectralis_mod_map

# MODIFICATION_MAPPING_TO_SPECTRALIS = generate_spectralis_mod_map(ALL_MODIFICATION_LABELS)


### SPECTRALIS STUFF
MODIFICATION_MAPPING_TO_SPECTRALIS = {
    'Q[UNIMOD:7]': 'E',
    'N[UNIMOD:7]': 'D',
    'C[UNIMOD:4]': 'C', # L257 in spectralis_master of spectralis codebase
    '[+25.980265]-': '',
    '[UNIMOD:1]-': '',
    '[UNIMOD:5]-': '',
    '[UNIMOD:385]-': '',
    '[UNIMOD:35]': '',
    'S[UNIMOD:21]': 'S',
    'T[UNIMOD:21]': 'T',
    'Y[UNIMOD:21]': 'Y',
    '[UNIMOD:1]': '',
    '[UNIMOD:5]': '',
    '[+25.980265]': '',
    '[UNIMOD:385]': '',
    '-': ''
}
ENGINES = [
    "Casanovo4.2.0",
    "InstaNovo",
    "PepNet",
    "ContraNovo"
]
ENGINES_MAPPING = {
    "Casanovo4.2.0": "casanovo",
    "InstaNovo": "instanovo",
    "PepNet": "pepnet",
    "ContraNovo": "contranovo"
}

### PEPNOVO STUFF
PEPNOVO_COLUMN_MAPPING = {
    "#Index": "index",
    "RnkScr": "score",
    "PnvScr": "score_pepnovo",
    "N-Gap": "n_shift",
    "C-Gap": "c_shift",
    "[M+H]": "mass",
    "Charge": "charge",
    "Sequence": "peptide"
}
PEPNOVO_COLUMNS = [
    'index',
    'score',
    'score_pepnovo',
    'n_shift',
    'c_shift',
    'mass',
    'charge',
    'peptide',
    'index',
    'scans',
    'title',
    'sqs'
]
no_solutions = "# No solutions found."
spectrum_error_1 = "# Could not process spectrum..."
spectrum_error_2 = "#Problem reading spectrum..."
too_few_peaks = "# too few peaks..."

spectrum_errors = {
    no_solutions: "no_solution",
    spectrum_error_1: "process_error",
    spectrum_error_2: "reading_error",
    too_few_peaks: "sparse_spectrum"
}