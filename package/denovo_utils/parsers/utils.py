from psm_utils import Peptidoform
from psm_utils.peptidoform import PeptidoformException
from spectrum_utils import proforma
from spectrum_utils.fragment_annotation import get_theoretical_fragments

import re

import pyopenms as oms
import numpy as np

# Add missing modification
modification = oms.ResidueModification()
modification.setAverageMass(-17.0305)
modification.setFullName("Loss of ammonia")
modification.setMonoMass(-17.026549)
modification.setName("N-oxobutanoic")
modification.setId("UniMod:385")

db = oms.ModificationsDB()
db.addModification(modification)

def proforma_to_oms(peptide: str):
    """
    Parses a peptide sequence in proforma format to pyOpenMS compatible format

    Parameter
    ---------
    peptide: str
        Peptide string in proforma format
    
    Returns
    -------
    AASequence (pyOpenMS):
        A peptide sequence in pyOpenMS format
    int:
        charge of the peptide
    """
    # Check if peptide is in proforma format
    try:
        peptidoform = Peptidoform(peptide)
    except:
        raise PeptidoformException(f"({peptide}) not in proforma format.")

    # Error with UNIMOD:385
    peptide=peptide.replace("UNIMOD:385", "-17.027")

    # Reformat unimod modifications
    pattern_unimod = r'\[UNIMOD:(\d+)\]'
    def replace_unimod(match):
        return f"(UniMod:{match.group(1)})"
    
    peptide_oms_str = re.sub(
        pattern=pattern_unimod,
        repl=replace_unimod,
        string=peptide
    )

    # Parse N-terminal modifications
    if ")-" in peptide_oms_str:
        nterm_modification, peptide_oms_str = peptide_oms_str.split(")-")
        nterm_modification += ")"
        peptide_oms_str = "." + nterm_modification + peptide_oms_str + "."
    elif "]-" in peptide_oms_str:
        nterm_modification, peptide_oms_str = peptide_oms_str.split("]-")
        nterm_modification += "]"
        peptide_oms_str = "." + nterm_modification + peptide_oms_str + "."

    # Split the charge from the peptide string
    if "/" in peptide_oms_str:
        peptide_oms_str, charge = peptide_oms_str.split("/")
    else:
        charge = None

    peptide_oms = oms.AASequence.fromString(peptide_oms_str)

    return peptide_oms, charge

def proforma_to_theoretical_spectrum(peptide: str, engine="spectrum-utils"):
    """
    Create a theoretical spectrum from a peptide sequence.

    Parameter
    ---------
    peptide: str
        Peptide sequence in proforma format
    engine: str
        The engine to use to create theoretical spectrum.
        Can only be 'pyopenms' or 'spectrum-utils' (default)
        
    Return
    ------
    MSSpectrum
        Spectrum object in pyOpenMS format
    """
    if engine=="spectrum-utils":
        proforma_str = proforma.parse(peptide)[0]
        theoretical_fragments = get_theoretical_fragments(
            proteoform=proforma_str
        )
        mz_array = [peak[1] for peak in theoretical_fragments]
        intensity_array = [1.0]*len(mz_array)

        spectrum = oms.MSSpectrum()
        spectrum.set_peaks(
            [
                mz_array,
                intensity_array
            ]
        )
        return spectrum

    elif engine=="pyopenms":
        # Reformat peptide sequence in pyOpenMS format
        peptide_oms, charge = proforma_to_oms(peptide=peptide)
        
        # Initialize the required objects to create the spectrum
        spectrum = oms.MSSpectrum()
        tsg = oms.TheoreticalSpectrumGenerator()
        p = oms.Param()

        p.setValue("add_b_ions", "true")
        p.setValue("add_metainfo", "true")
        tsg.setParameters(param=p)

        # Create the theoretical spectrum
        tsg.getSpectrum(
            spec=spectrum,
            peptide=peptide_oms,
            min_charge=1,
            max_charge=2
        )
        return spectrum
    
    else:
        raise Exception("Engine not suported for theoretical spectrum generation.")
