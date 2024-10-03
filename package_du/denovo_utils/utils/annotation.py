from numba import types
from numba.typed import Dict
from alphapept.constants import isotopes
from alphapept.chem import dict_to_dist
from pyteomics import mass as ptmass
import spectrum_utils.spectrum as sus
import numpy as np
from matplotlib import pyplot as plt

NEUTRAL_LOSSES = [
    "-HCOOH",
    "-CH4OS",
    "-CO",
    "-CO2",
    "-H2O",
    "-NH3",
    "-H3PO4",
    "-C2H4O2S",
    "-H",
    "-HCONH2",
    "-HPO3",
    "-SO3"
]

def select_lowest_error_annotation(annotation_list):
    """Return the index of the annotation with lowest Da error."""
    if len(annotation_list)==1:
        return 0
    
    lowest = np.inf
    best_annotation_index = None
    for i, annotation in enumerate(annotation_list):
        mz_delta = annotation.mz_delta[0]
        if np.abs(mz_delta) < lowest:
            lowest = mz_delta
            best_annotation_index = i
    
    return best_annotation_index

def select_backbone_first(annotation_list):
    """Give preference to backbone ions before looking at neutral losses."""
    if len(annotation_list)==1:
        return 0
    
    backbone_i = []
    backbone_annot = []

    nl_i = []
    nl_annot = []

    for i, annotation in enumerate(annotation_list):
        if annotation.neutral_loss:
            nl_i.append(i)
            nl_annot.append(annotation)
        else:
            backbone_i.append(i)
            backbone_annot.append(annotation)
    
    if len(backbone_i)==1:
        return backbone_i[0]
    elif len(backbone_i) > 1:
        i = select_lowest_error_annotation(backbone_annot)
        return backbone_i[i]
    else:
        i = select_lowest_error_annotation(nl_annot)
        return nl_i[i]

class AnnotatedSpectrum():
    def __init__(
            self,
            spectrum: sus.MsmsSpectrum=None,
            peplen: int=None,
            tolerance: float = .02,
            process=False,
            fill=None,
            verbose=False,
            add_neutral_losses = NEUTRAL_LOSSES,
            mode_multiple_annotations=select_backbone_first
        ):
        """
        The AnnotatedSpectrum object stores and allows easy retrieval of annotated ions.

        These include the following:
        - abcxyz: both singly and doubly charged
        - p: 1-5 charges
        
        They are stored in a dictionary according to peptide length
        """
        if isinstance(fill, AnnotatedSpectrum):
            self.fill_spectrum(fill)

        else:
            self.spectrum = spectrum
            self.ion_length = peplen-1
            self.tolerance = tolerance
            self.verbose = verbose
            self.multiple_annotations = []
            self.tic = np.sum(self.spectrum.intensity)
            self.annotated = False
            self.neutral_losses = add_neutral_losses
            self.mode_multiple_annotations = mode_multiple_annotations
        self._initialize_ion_vectors()

        if process:
            self.annotate_ions()
            self.add_isotopes()
    
    def fill_spectrum(self, spectrum):
        self.spectrum = spectrum.spectrum
        self.ion_length = spectrum.ion_length
        self.tolerance = spectrum.tolerance
        self.verbose = spectrum.verbose
        self.multiple_annotations = []
        self.tic = np.sum(self.spectrum.intensity)
        self.annotated = False
        self.mode_multiple_annotations = self.mode_multiple_annotations
        self._initialize_ion_vectors()

    def _get_ion_vector(self):
        return {
            "b1": np.zeros(self.ion_length),
            "b2": np.zeros(self.ion_length),
            "y1": np.zeros(self.ion_length),
            "y2": np.zeros(self.ion_length),
            "a1": np.zeros(self.ion_length),
            "a2": np.zeros(self.ion_length),
            "c1": np.zeros(self.ion_length),
            "c2": np.zeros(self.ion_length),
            "z1": np.zeros(self.ion_length),
            "z2": np.zeros(self.ion_length),
            "x1": np.zeros(self.ion_length),
            "x2": np.zeros(self.ion_length),
            "p": np.zeros(self.spectrum.precursor_charge),
            "I": [],
        }
    
    def _get_isotope_matrix(self, n_isotopes):
        return {
            "b1": np.full((self.ion_length, n_isotopes), .0),
            "b2": np.full((self.ion_length, n_isotopes), .0),
            "y1": np.full((self.ion_length, n_isotopes), .0),
            "y2": np.full((self.ion_length, n_isotopes), .0),
            "a1": np.full((self.ion_length, n_isotopes), .0),
            "a2": np.full((self.ion_length, n_isotopes), .0),
            "c1": np.full((self.ion_length, n_isotopes), .0),
            "c2": np.full((self.ion_length, n_isotopes), .0),
            "z1": np.full((self.ion_length, n_isotopes), .0),
            "z2": np.full((self.ion_length, n_isotopes), .0),
            "x1": np.full((self.ion_length, n_isotopes), .0),
            "x2": np.full((self.ion_length, n_isotopes), .0),
            "p": np.full((self.spectrum.precursor_charge, n_isotopes), .0),
            "I": [],
        }

    def _initialize_ion_vectors(self):
        self.ion_vectors = self._get_ion_vector()
        self.ion_mz = self._get_ion_vector()
        self.ion_mz_delta = self._get_ion_vector()

        # Neutral loss tables
        self.nl_ion_vectors = {
            neutral_loss: self._get_ion_vector() for neutral_loss in NEUTRAL_LOSSES
        }
        self.nl_ion_mz = {
            neutral_loss: self._get_ion_vector() for neutral_loss in NEUTRAL_LOSSES
        }
        self.nl_ion_mz_delta = {
            neutral_loss: self._get_ion_vector() for neutral_loss in NEUTRAL_LOSSES
        }

    def handle_multiple_annotations(self, dict_spectrum):
        """
        Make sure every entry only holds 1 annotation.

        Multiple annotations are handled with a function passed to the mode parameter.
        By default, it selects the one with the lowest absolute Dalton error.

        Returns
        -------
        dict
            mz_list, intensity_list, annotation_list
        """
        mz_list = []
        intensity_list = []
        annotation_list = []

        for peak_i in range(len(dict_spectrum['mz'])):

            annotations = dict_spectrum['annotation'][peak_i].fragment_annotations
            mz = dict_spectrum['mz'][peak_i]
            intensity = dict_spectrum['intensity'][peak_i]


            if len(annotations)>0:
                if len(annotations)>1:
                    self.multiple_annotations.append(
                        {
                            'annotations': annotations,
                            'mz': mz,
                            'intensity': intensity,
                        }
                    )
                
                best_annotation_index = self.mode_multiple_annotations(annotations)
                mz_list.append(mz)
                intensity_list.append(intensity)
                annotation_list.append(annotations[best_annotation_index])

        return {
            "intensity": np.array(intensity_list),
            "mz": np.array(mz_list),
            "annotation": np.array(annotation_list),
        }
    
    def add_ion_type_annotation(self, annotation_list):
        """
        Return a list of corresponding ion types.

        - abc xyz ions: ion_type(charge)
        - p ions: p
        - I ions: I
        """
        ion_type_annotation = []
        ion_type_annotation_nl = []
        nl_type = []
        for a in annotation_list:
            ion_type = a.ion_type
            nl = a.neutral_loss
            
            # Whether to add in neutral loss table
            if nl is not None:
                if ion_type[0] not in ["p", "I"]:
                    ion_type = ion_type[0]+str(a.charge)
                ion_type_annotation_nl.append(ion_type)
                nl_type.append(nl)
                ion_type_annotation.append("-")
                continue

            ion_type_annotation_nl.append("-")
            nl_type.append("-")
            if ion_type[0] == "p":
                if a.charge <= self.spectrum.precursor_charge:
                    ion_type_annotation.append("p")
                else:
                    ion_type_annotation.append("-")
            elif ion_type[0] == "I":
                ion_type_annotation.append("I")
            else:
                ion_type_annotation.append(ion_type[0]+str(a.charge))
        return np.array(ion_type_annotation), np.array(ion_type_annotation_nl), np.array(nl_type)

    def add_sorting_index(self, annotation_list):
        """
        Sort the abc xyz ions by ion_type number. The p ions by charge. The I ions by mz.
        """
        sorted_index = []
        encountered_I = 0
        for a in annotation_list:
            ion_type = a.ion_type
            if ion_type[0] == "p":
                if a.charge <= self.spectrum.precursor_charge:
                    sorted_index.append(a.charge)
                else:
                    if self.verbose:
                        print(f"Warning: Higher charged precursor ion found '{a}'")
                    sorted_index.append(0)
            elif ion_type[0] == "I":
                encountered_I += 1
                sorted_index.append(encountered_I)
            else:
                sorted_index.append(int(ion_type[1:]))


        # -1 for correct indexing starting from 0
        return np.array(sorted_index)-1
    
    def add_mz_delta(self, annotation_list):

        mz_deltas = []
        for a in annotation_list:
            ion_type = a.ion_type
            if ion_type[0] == "p":
                if a.charge <= self.spectrum.precursor_charge:
                    mz_deltas.append(np.abs(a.mz_delta[0]))
                else:
                    mz_deltas.append("-")
            else:
                mz_deltas.append(np.abs(a.mz_delta[0]))
        return np.array(mz_deltas)

    def fill_nl_ion_vectors(self, spec_dict, ion_type):
        
        for nl_type in self.neutral_losses:
            mask = (spec_dict["ion_type_nl"]==ion_type) & (spec_dict["nl_type"]==nl_type)
            sorting_index = spec_dict["ion_n"][mask]
            
            self.nl_ion_vectors[nl_type][ion_type][sorting_index] = spec_dict["intensity"][mask]
            self.nl_ion_mz_delta[nl_type][ion_type][sorting_index] = spec_dict["mz_delta"][mask]
            self.nl_ion_mz[nl_type][ion_type][sorting_index] = spec_dict["mz"][mask]

    def annotate_ions(self, ions="a1 a2 b1 b2 c1 c2 x1 x2 y1 y2 z1 z2 p I".split()):
        spec_dict = self.handle_multiple_annotations(
            {
                "intensity": self.spectrum.intensity,
                "mz": self.spectrum.mz,
                "annotation": self.spectrum.annotation
            }
        )
        spec_dict["ion_type"], spec_dict["ion_type_nl"], spec_dict["nl_type"] = self.add_ion_type_annotation(spec_dict["annotation"])
        spec_dict["ion_n"] = self.add_sorting_index(spec_dict["annotation"])
        spec_dict["mz_delta"] = self.add_mz_delta(spec_dict["annotation"])
        
        if self.verbose:
            print(spec_dict)
        # Fill the singly charged a b c x y z ion vectors

        masks = []
        for ion_type in ions:

            if np.isin(np.array([ion_type]), spec_dict["ion_type_nl"]).any():
                self.fill_nl_ion_vectors(spec_dict, ion_type)

            if not np.isin(np.array([ion_type]), spec_dict["ion_type"]).any():
                continue
    
            # TODO: Handle Immonium ion annotations
            if ion_type == "I":
                continue

            # Store in array of 5 sorted by charge
            mask = spec_dict["ion_type"]==ion_type 
            sorting_index = spec_dict["ion_n"][mask]
            self.ion_vectors[ion_type][sorting_index] = spec_dict["intensity"][mask]
            self.ion_mz_delta[ion_type][sorting_index] = spec_dict["mz_delta"][mask]
            self.ion_mz[ion_type][sorting_index] = spec_dict["mz"][mask]
            masks.append(mask)

        ignored_mask = sum(masks)==0
        self.ignored_ions = {
            "intensity": spec_dict["intensity"][ignored_mask],
            "mz": spec_dict["mz"][ignored_mask],
            "annotation": spec_dict["annotation"][ignored_mask]
        }
        self.annotated = True
    
    def return_ion_type(self, ion_type):
        if ion_type not in self.ion_vectors.keys():
            raise Exception(f"Ion type {ion_type} not supported.")
        
        return np.array([i for i in self.ion_vectors[ion_type]])
    
    def count_ions(self, ion_type, pct=False):
        if ion_type not in self.ion_vectors.keys():
            raise Exception(f"Ion type {ion_type} not supported.")
        
        n = sum(self.ion_vectors[ion_type]!=0)
        if pct:
            return n/self.ion_length
        else:
            return n
        
    def add_isotopes(self, isom=1.003355, n_isotopes=5):
        self.isotope_matrix = self._get_isotope_matrix(n_isotopes=n_isotopes)
        self.isotope_matrix_mz = self._get_isotope_matrix(n_isotopes=n_isotopes)

        # Singly charged fragment ions
        for ion_type in self.ion_vectors.keys():
            
            # Find isotope peaks for special ion type: precursor ion
            if ion_type == "p":
                for ion_n, charge_state in enumerate(range(1, self.spectrum.precursor_charge+1)):
                    if self.ion_vectors[ion_type][ion_n] == 0:
                        continue
                    
                    # Get monoisotopic intensity for precursor ion
                    monoisotopic_intensity = self.ion_vectors[ion_type][ion_n]
                    monoisotopic_mass = self.ion_mz[ion_type][ion_n]
                    self.isotope_matrix[ion_type][ion_n][0] = monoisotopic_intensity

                    # Register other isotopes
                    for i in range(1, n_isotopes):

                        # Get the intensity of next isotope and store if found
                        next_isotope = monoisotopic_mass + i*(isom/charge_state)
                        isotope_intensity, isotope_mz = self.extract_isotope(
                            mz=next_isotope,
                            tolerance=self.tolerance
                        )
                        if isotope_intensity:
                            self.isotope_matrix[ion_type][ion_n][i] = isotope_intensity
                            self.isotope_matrix_mz[ion_type][ion_n][i] = isotope_mz
                        
                        # Stop searching for other isotopes if previous not found
                        else:
                            break
                continue

            # Skip "I" ions for isotope peak searching
            elif not ion_type.endswith("1") and not ion_type.endswith("2"):
                continue

            # Find isotope peaks for singly and doubly charged ion type series
            # Similar code as for precursor ions
            charge_state = int(ion_type[-1])
            for ion_n in range(self.ion_length):

                # If monoisotopic mass is not found, do not search for other isotopes
                if self.ion_vectors[ion_type][ion_n] == 0:
                    continue

                # Register monoisotopic mass
                monoisotopic_intensity = self.ion_vectors[ion_type][ion_n]
                monoisotopic_mass = self.ion_mz[ion_type][ion_n]
                self.isotope_matrix[ion_type][ion_n][0] = monoisotopic_intensity

                # Register other isotopes
                for i in range(1, n_isotopes):
                    next_isotope = monoisotopic_mass+ i*(isom/charge_state)
                    isotope_intensity, isotope_mz = self.extract_isotope(
                        mz=next_isotope, 
                        tolerance=self.tolerance
                    )
                    if isotope_intensity:
                        self.isotope_matrix[ion_type][ion_n][i] = isotope_intensity
                        self.isotope_matrix_mz[ion_type][ion_n][i] = isotope_mz
                    else:
                        break

    def extract_isotope(self, mz, tolerance, handle_multiple=max):
        candidate_intensities = self.spectrum.intensity[
            (self.spectrum.mz < mz+tolerance) & 
            (self.spectrum.mz > mz-tolerance)
        ]
        candidate_mzs = self.spectrum.mz[
            (self.spectrum.mz < mz+tolerance) & 
            (self.spectrum.mz > mz-tolerance)
        ]
        if candidate_intensities.any():
            registered_intensity = handle_multiple(candidate_intensities)
            max_i = np.argmax(candidate_intensities)
            registered_mz = candidate_mzs[max_i]
            return registered_intensity, registered_mz
        else:
            return False, "_"


def aaseq_to_isotope_dist(peptide, isotope_constants):

    # Generate atom dictionary
    atom_dict = ptmass.BasicComposition()
    for aa in peptide:
        if aa=="U":
            aa="A"
        atom_dict += ptmass.std_aa_comp[aa]
    
    # Convert to numba format
    numba_typed_dict = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.int64
    )

    for key, value in atom_dict.items():
        numba_typed_dict[key] = value

    return dict_to_dist(numba_typed_dict, isotope_constants)


def get_annotated_spectrum(spectrum, ions="abcxyzIp"):
    """Annotate the peaks of a spectrum"""
    return sus.MsmsSpectrum(
            identifier=spectrum["spectrum_id"],
            precursor_mz=spectrum["precursor_mz"],
            precursor_charge=spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["charge state"],
            mz=spectrum["m/z array"],
            intensity=spectrum["intensity array"]
        ).annotate_proforma(str(spectrum["peptidoform"]),
                            .02, "Da", ions, neutral_losses=False, max_ion_charge=5)


def plot_isotope_mse(spec, isotopes, ion_type, ion_n, log_transform=False):
    peptide = spec.spectrum.proforma[:-2]
    if ion_type[0] in "x y z".split():
        peptide_part = peptide[-ion_n:]
    elif ion_type == "p":
        peptide_part = peptide
    else:
        peptide_part = peptide[:ion_n]
    
    print(peptide_part)
    observed = spec.isotope_matrix[ion_type][ion_n]
    observed = observed / max(observed)
    expected = aaseq_to_isotope_dist(peptide_part, isotopes).intensities[:5]

    if log_transform:
        observed = np.log2(observed + .001)
        expected = np.log2(expected + .001)

    plt.stem(observed, markerfmt = "g")
    plt.stem(expected, markerfmt = "r")