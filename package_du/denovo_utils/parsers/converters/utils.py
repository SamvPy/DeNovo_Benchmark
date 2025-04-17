from pyteomics import mzml
import pandas as pd
from tqdm import tqdm

def mzml_reader(path, im):
    mzml_file = mzml.read(path)
    ion_mobility = None
    entries = []
    for spectrum in tqdm(mzml_file):

        title = spectrum['id']
        selected_ion = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
        pepmass = selected_ion['selected ion m/z']
        rtinseconds = spectrum['scanList']['scan'][0]['scan start time'] * 60
        charge = [selected_ion['charge state']]
        scans = title.split('scan=')[-1]

        # If contains ion mobility, extract it
        if im:
            ion_mobility = spectrum['scanList']['scan'][0]['inverse reduced ion mobility']

        entries.append(
            {
                'title': title,
                'pepmass': (pepmass, None),
                'rtinseconds': rtinseconds,
                'charge': charge,
                'scans': scans,
                'ion_mobility': ion_mobility
            }
        )
    return pd.DataFrame(entries)

def infer_rank(spectrum_id, spectrum_id_dict):
    'Rank is believed to be in order.'
    if spectrum_id not in spectrum_id_dict.keys():
        spectrum_id_dict[spectrum_id] = 0
    
    spectrum_id_dict[spectrum_id] += 1
    return spectrum_id_dict[spectrum_id]