from typing import Optional, List
from peak_pack.annotation import PeptideEvidence
from .psm import PSM
from .spectrum import Spectrum
from ..parsers import DenovoEngineConverter


class Run:
    def __init__(self, run_id):
        self.run_id = run_id
        self.spectra: dict[str, Spectrum] = {}  # Dictionary to map spectrum_id to Spectrum object
    
    def __repr__(self):
        engines = []
        for spectra in self.spectra.values():
            for psm in spectra.psm_candidates:
                if psm.engine_name not in engines:
                    engines.append(psm.engine_name)
        return f"{len(self.spectra)} spectra loaded from {len(engines)} engines ({engines})."

    def load_data(self, psmlist, is_ground_truth=False, filter_by_gt=True):

        for psm_ in psmlist.get_rank1_psms():
            # Only ground truth psms can add a Spectrum object to the run.
            if is_ground_truth:
                if psm_['is_decoy'] or psm_['qvalue']>.01:
                    continue
                spectrum = Spectrum(psm_['spectrum_id'])
                self.add_spectrum(spectrum=spectrum)

            spectrum = self.get_spectrum(psm_["spectrum_id"])
            
            if not isinstance(spectrum, Spectrum):
                if filter_by_gt:
                    continue
                # Add spectra if psms without a gt still needs to be added.
                else:
                    spectrum = Spectrum(psm_["spectrum_id"])
                    self.add_spectrum(spectrum)
            
            psm = PSM(
                peptidoform=psm_['peptidoform'],
                score=psm_['score'],
                engine_name=psm_['source'],
                aa_score=extract_metadata(psm_, 'aa_scores', infer_datatype=True),
                peptide_evidence=PeptideEvidence.load(
                    extract_metadata(psm_, 'peptide_evidence')
                )
            )
            
            # Add ms2rescore score if rescored.
            if 'score_ms2rescore' in psm_['provenance_data'].keys():
                psm.scores.add_score(
                    score=psm_['provenance_data']['score_ms2rescore'],
                    metadata="ms2rescore",
                    score_type="peptide"
                )

            spectrum.add_psm(psm, is_ground_truth=is_ground_truth)


    def load_refinement(self, psmlist, filter_by_gt=True):
        raise NotImplementedError()

    def add_spectrum(self, spectrum: Spectrum):
        if spectrum.spectrum_id in self.spectra.keys():
            raise Exception(f"Cannot add an existing spectrum. {spectrum.spectrum_id}")
        self.spectra[spectrum.spectrum_id] = spectrum
    
    def get_spectrum(self, spectrum_id) -> Spectrum:
        return self.spectra.get(spectrum_id, None)
    
    def get_common_spectra(self, engines, include_gt=True) -> 'Run':
        spectra = {}
        for spectrum_id, spectrum in self.spectra.items():
            psms = []
            common = True
            for engine in engines:
                psm = spectrum.get_psms_by_engine(engine_name=engine)
                if len(psm) == 0:
                    common = False
                    break
                psms += psm

            if common:
                if include_gt:
                    psms += [spectrum.psm_gt]
                spectra[spectrum_id] = spectrum
        
        new_run = Run(self.run_id)
        new_run.spectra = spectra
        return new_run


def extract_metadata(psm, key, infer_datatype=False):
    if key in psm['metadata'].keys():
        el = psm['metadata'][key]
        if infer_datatype:
            return eval(el)
        return el
    return None