
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
        "C(+57.02)": "[UNIMOD:4]",
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

def generate_spectralis_mod_map(all_labels: dict) -> dict:
    """
    Generate a mapping of modifications to Spectralis format.

    Parameters
    ----------
    all_labels : dict
        A dictionary of all unique modification labels.

    Returns
    -------
    dict
        A dictionary mapping modifications to Spectralis format.
    """
    spectralis_modifications = [
        "C[UNIMOD:4]", "M[UNIMOD:35]", 'Q[UNIMOD:7]', 'N[UNIMOD:7]'
    ]

    spectralis_mod_map = {
        'Q[UNIMOD:7]': 'E',
        'N[UNIMOD:7]': 'D'
    }

    for modification in all_labels.values():
        if modification not in spectralis_modifications:
            spectralis_mod_map[modification] = ""
    return spectralis_mod_map

MODIFICATION_MAPPING_TO_SPECTRALIS = generate_spectralis_mod_map(ALL_MODIFICATION_LABELS)