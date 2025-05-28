"""Utility functions for proforma peptide sequence parsing."""

import logging
import re
from typing import Literal, Optional

import pyopenms as oms
from psm_utils import Peptidoform
from psm_utils.peptidoform import PeptidoformException
from pyteomics.proforma import ProFormaError
from spectrum_utils import proforma
from spectrum_utils.fragment_annotation import get_theoretical_fragments

# Add missing modification
modification_385 = oms.ResidueModification()
modification_385.setAverageMass(-17.0305)
modification_385.setFullName("Loss of ammonia")
modification_385.setMonoMass(-17.026549)
modification_385.setName("N-oxobutanoic")
modification_385.setId("UniMod:385")

modification_5 = oms.ResidueModification()
modification_5.setAverageMass(43.0247)
modification_5.setFullName("Carbamylation")
modification_5.setMonoMass(43.005814)
modification_5.setName("Carbamyl")
modification_5.setId("UniMod:5")

db = oms.ModificationsDB()
db.addModification(modification_385)
db.addModification(modification_5)


def proforma_to_oms(peptide: str) -> tuple[oms.AASequence, Optional[int]]:
    """
    Parse a peptide sequence in proforma format to pyOpenMS compatible format.

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
        _ = Peptidoform(peptide)
    except ProFormaError:
        raise PeptidoformException(f"({peptide}) not in proforma format.")

    # Error with UNIMOD:385
    peptide = peptide.replace("UNIMOD:385", "-17.027")
    peptide = peptide.replace("UNIMOD:5", "43.005814")

    # Reformat unimod modifications
    pattern_unimod = r"\[UNIMOD:(\d+)\]"

    def replace_unimod(match):
        return f"(UniMod:{match.group(1)})"

    peptide_oms_str = re.sub(
        pattern=pattern_unimod, repl=replace_unimod, string=peptide
    )

    # Parse N-terminal modifications
    if ")-" in peptide_oms_str:
        peptide_oms_list = peptide_oms_str.split(")-")
        nterm_modification, peptide_oms_str = peptide_oms_list[-2], peptide_oms_list[-1]
        nterm_modification += ")"
        peptide_oms_str = "." + nterm_modification + peptide_oms_str + "."
    elif "]-" in peptide_oms_str:
        peptide_oms_list = peptide_oms_str.split("]-")
        nterm_modification, peptide_oms_str = peptide_oms_list[-2], peptide_oms_list[-1]
        nterm_modification += "]"
        peptide_oms_str = "." + nterm_modification + peptide_oms_str + "."

    # Split the charge from the peptide string
    if "/" in peptide_oms_str:
        peptide_oms_str, charge = peptide_oms_str.split("/")
    else:
        charge = None

    peptide_oms = oms.AASequence.fromString(peptide_oms_str)

    return peptide_oms, charge


def proforma_to_theoretical_spectrum(
    peptide: str, engine: Literal["spectrum-utils", "pyopenms"] = "spectrum-utils"
) -> oms.MSSpectrum:
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
    if engine == "spectrum-utils":
        proforma_str = proforma.parse(peptide)[0]
        theoretical_fragments = get_theoretical_fragments(proteoform=proforma_str)
        mz_array = [peak[1] for peak in theoretical_fragments]
        intensity_array = [1.0] * len(mz_array)

        spectrum = oms.MSSpectrum()
        spectrum.set_peaks([mz_array, intensity_array])
        return spectrum

    elif engine == "pyopenms":
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
        tsg.getSpectrum(spec=spectrum, peptide=peptide_oms, min_charge=1, max_charge=2)
        return spectrum

    else:
        raise Exception("Engine not suported for theoretical spectrum generation.")


def parse_peptidoform(
    peptide: str, mapping: dict, max_length=30
) -> Optional[Peptidoform]:
    """
    Parse a peptide string into a psm-utils peptidoform.

    Also parsed modifications with the mapping dictionary.

    Parameters
    ----------
    peptide: str
        Peptide string in proforma format
    mapping: dict
        Mapping of modifications or residues.
    max_length: int
        Ignores peptide sequences with longer residues.
    """
    peptide_parsed = peptide

    # If missing prediction for particalur spectrum
    if peptide.startswith('nan/'):
        return None

    for k, v in mapping.items():
        # Make sure N-terminal mods occur N-terminally!
        if ("-" in v):
            if peptide_parsed.startswith(k):
                peptide_parsed = peptide_parsed.replace(k, v, 1)
            if k in peptide_parsed[1:]: # Ensures that this mod only occurs N-terminally
                logging.warning(f"Mapping string is not consistent for {peptide} ({v})")
                return None
        else:
            peptide_parsed = peptide_parsed.replace(k, v)

    try:
        peptidoform = Peptidoform(peptide_parsed)
        if (
            (len(peptidoform) > max_length)
            or (peptidoform.precursor_charge > 6)
            or (len(peptidoform) < 2)
        ):
            return None
        return peptidoform
    except PeptidoformException:
        logging.warning(f"Failed to parse: {peptide}")
        return None
