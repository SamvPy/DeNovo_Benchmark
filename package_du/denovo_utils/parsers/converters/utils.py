from pyteomics import mzml
import pandas as pd
from tqdm import tqdm

def mzml_reader(path):
    mzml_file = mzml.read(path)
    
    entries = []
    for spectrum in tqdm(mzml_file):

        title = spectrum['id']
        selected_ion = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
        pepmass = selected_ion['selected ion m/z'] * selected_ion['charge state']
        rtinseconds = spectrum['scanList']['scan'][0]['scan start time'] * 60
        charge = [selected_ion['charge state']]
        scans = title.split('scan=')[-1]

        entries.append(
            {
                'title': title,
                'pepmass': (pepmass, None),
                'rtinseconds': rtinseconds,
                'charge': charge,
                'scans': scans
            }
        )
    return pd.DataFrame(entries)