### Data index for Eawag data

# imports
import pandas as pd


def get_data_index(lake, dtype):
    """
    Get data index of available processed data.

    Parameters
    ----------
    lake : str
        Lake (i.e., Zug, Lucerne).
    dtype : str
        Data type (i.e., ctd, microstructure, mooring).

    Returns
    -------
    data_index : pd.DataFrame
        Index of processed data.
    """
    di_path = f'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{dtype}.json'

    if dtype == 'ctd':
        data_index = pd.read_json(di_path).sort_values(by=['date', 'time'], ascending=True).reset_index(drop=True)
    elif dtype == 'mooring':
        data_index = pd.read_json(di_path).sort_values(by=['date', 'location', 'depth'], ascending=True).reset_index(drop=True)
    else:
        raise ValueError('Only ctd and mooring data indices implemented.')
    
    return data_index