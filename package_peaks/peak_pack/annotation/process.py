from .annotation import SpectrumVector

import multiprocessing
from tqdm import tqdm


class PSMProcessor:
    def __init__(self, config):
        self.config = config
    
    def process_psm(self, psm_mgf):
        sv = SpectrumVector(
            ion_types = self.config["ion_types"],
            neutral_losses = self.config["neutral_losses"]
        )
        sv.parse(psm_mgf[0], psm_mgf[1])

        return sv


# Use multiprocessing.Pool for parallel execution
def parallel_annotate_psms(psm_list, mgf_spectra, config):
    """
    Parse a bunch of spectra in parallel.
    
    psmlist and mgf_spectra should be in the same order regarding spectrum_ids.
    """
    results = []

    psm_processor = PSMProcessor(config=config)

    with multiprocessing.Pool(processes=config["max_workers"]) as pool, tqdm(total=len(psm_list)) as pbar:
        # Use pool.imap to process the function with multiple inputs in parallel
        for result in pool.imap(psm_processor.process_psm, zip(psm_list, mgf_spectra)):
            pbar.update()
            pbar.refresh()

            results.append(
                result
            )
    return results