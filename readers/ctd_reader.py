### Class for reading CTD data

# imports
import pandas as pd
import xarray as xr
import pyrsktools as rsk
import os



class CTDReader():
    SID_BRAND_MAP = {'1807': 'Sea&Sun', '2023': 'Sea&Sun', '66131': 'RBR'}
    COLS_DROP_RBR = ['temperature1']
    COLS_MAP_RBR = {}

    def __init__(self, serial_id, lake, location, year, date, bathy_file, datalakes=False):
        """
        Initialize CTDReader object.

        Parameters
        ----------
        serial_id : str
            Serial number of CTD.
        lake : str
            Lake where CTD profiles.
        year : str
            Year of CTD profiles. 
        date : str
            Date (YYYYMMDD) of CTD profiles.
        bathy_file : str
            File path to bathymetry data.
        datalakes : bool
            Toggle whether to read from Eawag drive or DataLakes.
        """
        self.serial_id = serial_id
        self.lake = lake
        self.location = location
        self.year = year
        self.date = date
        self.bathy_file = bathy_file
        self.datalakes = datalakes

    
    def separate_profiles_RBR_L0(self):
        """
        Load raw (L0) CTD from RBR CTD, separate profiles, and write to individual files.

        Returns
        -------
        profiles : list
            List of output file paths and xarray Datasets for each profile.
        """
        brand = self.SID_BRAND_MAP[self.serial_id]
        if brand != 'RBR':
            raise ValueError("Only RBR CTD L0 data requires profile separation.")
        
        froot, ext = os.path.splitext(self.fpath)
        if ext != '.rsk':
            raise ValueError("Only RBR CTD L0 .rsk files require profile separation.")
        
        with rsk.RSK(self.fpath) as f:
            f.readdata()
            data = f.data

            # find indices of each profile
            f.computeprofiles()
            profiles_idx = f.getprofilesindices(direction='both')

        # separate profiles and write out
        profiles = []
        for i, p_idx in enumerate(profiles_idx):
            p_data = pd.DataFrame(data[p_idx])
            p_data = p_data.set_index('timestamp')
            ds = xr.Dataset.from_dataframe(p_data)
            out_path = froot + f'_{i+1}.nc'
            profiles.append((out_path, ds))

            ds.to_netcdf(out_path, mode='w', format='NETCDF4')

        return profiles

        


    def load_from_L0(self):
        """
        Load raw (L0) CTD data into xarray Dataset.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by CTD.
        """
        brand = self.SID_BRAND_MAP[self.serial_id]

        if brand == 'Sea&Sun':
            raise NotImplementedError
        elif brand == 'RBR':
            pass
        else:
            raise NotImplementedError("Only Sea&Sun and RBR CTDs are handled.")