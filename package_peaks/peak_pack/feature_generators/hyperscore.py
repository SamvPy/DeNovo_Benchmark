from tqdm import tqdm

from ms2rescore.feature_generators.base import FeatureGeneratorBase
from typing import List
from psm_utils import PSMList

import numpy as np

def factorial(n):
    """
    Compute factorial of n using a loop.
    Parameters:
        n (int): Non-negative integer.
    Returns:
        int: Factorial of n.
    """
    result = 1
    for i in range(1, n + 1):
        result *= i
    return result
 
def calculate_hyperscore(n_y, n_b, y_ion_intensities, b_ion_intensities):
    """
    Calculate the hyperscore for a peptide-spectrum match.
    Parameters:
        n_y (int): Number of matched y-ions.
        n_b (int): Number of matched b-ions.
        y_ion_intensities (list): Intensities of matched y-ions.
        b_ion_intensities (list): Intensities of matched b-ions.
    Returns:
        float: Calculated hyperscore.
    """
    # Calculate the product of y-ion and b-ion intensities
    product_y = np.sum(y_ion_intensities) if y_ion_intensities else 1
    product_b = np.sum(b_ion_intensities) if b_ion_intensities else 1
 
    # Calculate factorial using custom function
    factorial_y = factorial(n_y)
    factorial_b = factorial(n_b)
 
    # Compute hyperscore
    hyperscore = np.log(factorial_y * factorial_b * (product_y+product_b))
    return hyperscore

def infer_fragment_identity(frag, allow_ion_types=['b', 'y']):
    ion = frag.ion

    is_allowed = False
    for allowed_ion_type in allow_ion_types:
        if allowed_ion_type in ion:
            is_allowed=True
            break
    
    if not is_allowed:
        return False
    # Radicals
    if "·" in ion:
        return False
    if frag.neutral_loss is not None:
        return False
    if frag.charge > 2:
        return False
    
    return ion[0]

def get_by_fragments(annotated_spectrum):
    b_intensities = []
    y_intensities = []
    for peak in annotated_spectrum:
        
        for fragment in peak.annotation:

            ion_type = infer_fragment_identity(fragment)
            
            if ion_type == 'b':
                b_intensities.append(peak.intensity)
            if ion_type == 'y':
                y_intensities.append(peak.intensity)
    return b_intensities, y_intensities


class HyperscoreGenerator(FeatureGeneratorBase):
    """MS2Rescore type feature generator for ppm-errors"""
    def __init__(self, *args, **kwargs):
        """Initialize feature generator class."""
        self.config = kwargs
        super().__init__(*args, **kwargs)

    @property
    def feature_names(self) -> List[str]:
        """Names of features added to rescoring_features dict."""
        return [
            "hyperscore"
        ]

    @property
    def name(self) -> str:
        return "hyperscore"
    @property
    def input_type(self) -> str:
        return "spectrum_vector"
    
    def add_features(self, psm_list: PSMList, spectra: List) -> None:
        """Compute and add rescoring features to psmlist."""

        for psm, spectrum in tqdm(zip(psm_list, spectra)):
            
            b, y = spectrum.by
            hs = calculate_hyperscore(
                n_y=len(y),
                n_b=len(b),
                y_ion_intensities=y,
                b_ion_intensities=b
            )
            
            psm.rescoring_features.update({"hyperscore": hs})
