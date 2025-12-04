### Class for reading CTD data

# imports
import pandas as pd
import numpy as np
import xarray as xr
import pyrsktools as rsk
import os

from .profile_reader import ProfileReader


class CTDReader(ProfileReader):
    SID_BRAND_MAP = {'1807': 'Sea&Sun', '2023': 'Sea&Sun', '66131': 'RBR'}
    COLS_DROP_RBR = ['temperature1', 'speed_of_sound']
    COLS_MAP_RBR = {}
    COLS_DROP_SEASUN = ['CAP25', 'DO_ml', 'pH', 'Redox', 'Vbatt', 'IntD', 'IntT']
    COLS_MAP_SEASUN = {
        'Press': 'press',
        'Temp': 'temp',
        'Cond': 'cond',
        'Turb': 'turb',
        'sat': 'd_oxygen_sat',
        'DO_mg': 'd_oxygen_conc',
        'pH_Tc': 'pH',
    }
    VAR_ATTRS = {
        'press': {'units': 'dbar', 'long_name': 'Pressure'},
        'temp': {'units': '°C', 'long_name': 'Temperature'},
        'cond': {'units': 'mS/cm', 'long_name': 'Conductivity'},
        'turb': {'units': 'FTU', 'long_name': 'Turbidity'},              # FTU = NTU
        'd_oxygen_sat': {'units': '%', 'long_name': 'Dissolved Oxygen Saturation'},
        'd_oxygen_conc': {'units': 'mg/l', 'long_name': 'Dissolved Oxygen Concentration'},
        'pH': {'units': 'pH', 'long_name': 'Potential of Hydrogen'},
        'depth': {'units': 'm', 'long_name': 'Depth'},
        'salin': {'units': 'PSU', 'long_name': 'Salinity'}              # PSU = ppt
    }


    def __init__(self, serial_id, lake, year, date, bathy_file, datalakes=False):
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
        super.__init__(lake, year, date, bathy_file, datalakes)
        self.serial_id = serial_id


    # ---------- Utility ----------

    def assign_attributes(self, ds):
        """
        Assign metadata to xarray Dataset attributes.
        """
        ds.assign_attrs(
            serial_id = self.serial_id
        )

        return ds

    
    # ---------- Organization ----------

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
        air_idx_start = 0
        for i, p_idx in enumerate(profiles_idx):
            # calculate air pressure prior to downcast
            air_pressure = np.mean(data[air_idx_start: min(p_idx)]['pressure'])
            air_idx_start = max(p_idx) + 1

            p_data = pd.DataFrame(data[p_idx])
            p_data = p_data.set_index('timestamp')
            ds = xr.Dataset.from_dataframe(p_data)
            ds = ds.assign_coords(air_pressure=air_pressure)
            out_path = froot + f'_{i+1}.nc'
            profiles.append((out_path, ds))

            ds.to_netcdf(out_path, mode='w', format='NETCDF4')

        return profiles


    # ---------- Processing ----------

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
            data = self.parse_sea_and_sun_L0()
        elif brand == 'RBR':
            pass
        else:
            raise NotImplementedError("Only Sea&Sun and RBR CTDs are handled.")
        
        data = data.set_index('time')
        ds = xr.Dataset.from_dataframe(data)

        # assign attributes to data variables
        for var, attrs in self.VAR_ATTRS.items():
            if var in ds:
                ds[var].attrs.update(attrs)

        return ds
        

    def parse_sea_and_sun_L0(self):
        """
        Parse raw (L0) data from Sea & Sun CTD.

        Returns
        -------
        data : pd.DataFrame
            Data from Sea & Sun CTD.
        """
        with open(self.fpath, encoding='latin1', errors='ignore') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if 'Lines' in line:
                break
        else:
            raise ValueError("Start of data not found within file.")
        
        cols = lines[i+2].replace(";", "").split()[1:]
        #units = lines[i+3].replace(";", "").replace("[", "").replace("]", "").split()

        data = pd.read_csv(self.fpath, sep=r'\s+', header=None, skiprows=i+5, names=cols, engine='python', encoding='cp1252')
        data['time'] = pd.to_datetime(data['IntD'] + " " + data['IntT'], format='%d.%m.%Y %H:%M:%S.%f')
        data['time'] = data['time'].dt.tz_localize('UTC').dt.tz_convert('Europe/Zurich')

        data = data.drop(self.COLS_DROP_SEASUN, axis=1)
        data = data.rename(columns=self.COLS_MAP_SEASUN)

        return data
    

    def parse_seabird_L0(self):
        """
        Parse raw (L0) data from Seabird CTD.
        """
        raise NotImplementedError
    

    def parse_RBR_L0(self):
        """
        Parse raw (L0) data from RBR CTD.
        """
        raise NotImplementedError