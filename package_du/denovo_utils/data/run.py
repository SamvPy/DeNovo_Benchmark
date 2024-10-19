from .spectrum import Spectrum
from ..parsers import DenovoEngineConverter


class Run:
    def __init__(self, run_id):
        self.run_id = run_id
        self.spectra = {}  # Dictionary to map spectrum_id to Spectrum object
    
    def load_data(self, path_mgf, path_results, engine):
        parser = DenovoEngineConverter.select(engine)
        psm_list = parser.parse(
            result_path=path_results,
            mgf_path=path_mgf
        )


    def add_spectrum(self, spectrum: Spectrum):
        self.spectra[spectrum.spectrum_id] = spectrum
    
    def get_spectrum(self, spectrum_id):
        return self.spectra.get(spectrum_id, None)