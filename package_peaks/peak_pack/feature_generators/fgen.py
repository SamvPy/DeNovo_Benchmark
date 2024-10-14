from ms2rescore.feature_generators.base import FeatureGeneratorBase
from typing import List
from psm_utils import PSM, PSMList
import logging
import numpy as np
from pyteomics import mgf

from .explained_intensity import ExplainedIntensityFeatures
from .hyperscore import HyperscoreGenerator
from .missing_fragmentation import MissingFragmentationFeatures
from .peak_features import PeakFeatures
from .ppm_error import PPMFeatures

from ..annotation.evidence import PeptideEvidence
from ..annotation.process import parallel_annotate_psms
from ..utils import ion_dict_to_matrix

class PeakFeatureGenerator(FeatureGeneratorBase):
    """MS2Rescore type feature generator for peak-type features."""

    def __init__(self, *args, **kwargs):
        """
        Initialize feature generator class.
        
        Specify which feature types should be included.
        """
        self.config = kwargs
        self.fgens = [
            PeakFeatures(),
            MissingFragmentationFeatures(),
            HyperscoreGenerator(**kwargs),
            PPMFeatures(),
            ExplainedIntensityFeatures()
        ]
        
        self.annotate = not ((len(self.fgens)==1) and (self.fgens[0].name=="general"))
        self.args = args

    @property
    def feature_names(self) -> List[str]:
        """Names of features added to rescoring_features dict."""
        return [
            # Type 1
            "tic",
            "peak_count"
            # Type 2
            "missing_frag_count",
            "missing_frag_pct",
            "missing_frag_sites",
            "missing_frag_longest_sequence",
            "missing_frag_from_n",
            "missing_frag_from_c",
            # Type 3
            "explained_y_pct",
            "explained_b_pct",
            "explained_by_pct",
            "explained_nl_pct",
            "explained_all_pct",
            # Type 4
            "hyperscore"
        ]

    def add_features(self, psm_list: PSMList) -> None:
        """Compute and add rescoring features to psmlist."""

        logging.info("Retrieving spectra...")
        mgf_file = mgf.read(self.config["spectrum_path"])
        spectra = mgf_file.get_by_ids([psm.spectrum_id for psm in psm_list])
        logging.info("Done.")

        if self.annotate:            
            # Get annoated spectrum vectors
            pes = []
            svs = parallel_annotate_psms(
                psm_list=psm_list,
                mgf_spectra=spectra,
                config=self.config,
            )
        
            # Get list of peptide evidences (missing fragmentation sites)
            for psm, sv in zip(psm_list, svs):
                matrix = np.concatenate(
                    [
                        ion_dict_to_matrix(
                            ion_dict=m,
                            ion_types=self.config["ion_types"],
                            n=sv.n
                        ) 
                        for m in sv.experimental_intensity.values()
                    ],
                    axis=0
                )

                pe = PeptideEvidence(
                    peptidoform=sv.peptidoform,
                    ion_matrix=matrix,
                    evidence_labels=[
                        ion_type + nl 
                        for ion_type in self.config["ion_types_evidence"] 
                        for nl in self.config["neutral_losses_evidence"] 
                    ]
                )
                psm.metadata.update({"peptide_evidence": pe})
                pes.append(pe)
            
        # Use the feature generators to generate a specific feature set on the annotated spectra and peptide evidence objects.
        # These features are stored inside the psms in the psm_list passed as parameters
        for fgen in self.fgens:
            if fgen.input_type=="peptide_evidence":
                fgen.add_features(psm_list, pes)
            elif fgen.input_type=="spectrum_vector":
                fgen.add_features(psm_list, svs)
            elif fgen.input_type=="psm_list":
                fgen.add_features(psm_list, spectra)
            else:
                raise Exception(f"Unrecognized input type: {fgen.input_type}")