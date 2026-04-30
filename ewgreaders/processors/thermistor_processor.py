### Class for processing Thermistor data

# imports
import json
import os
from glob import glob
import pyrsktools as rsk
import pandas as pd
import xarray as xr
import warnings


class ThermistorProcessor:
    MD_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_md.json'
    DPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/'
    DIPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/mooring.json'
    THERMISTORS = ['rbr_temp', 'rbr_duet']
    COLS_MAP = {'timestamp': 'time'}
    VARS_MAP = {'temperature': 'temp', 'pressure': 'press'}
    VAR_ATTRS = {
        'time': {'long_name': 'Coordinated Universal Time (UTC)'},
        'temp': {'units': '°C', 'long_name': 'Temperature'},
        'press': {'units': 'dbar', 'long_name': 'Pressure'}
    }

    
    def __init__(self, lake, year, date, location, serial_id):
        """
        Initialize ThermistorProcessor object.

        Parameters
        ----------
        lake : str
            Lake where thermistor is deployed.
        year : str
            Year of thermistor retrieval. 
        date : str
            Date (YYYYMMDD) of thermistor retrieval.
        location : str
            Location code within lake of thermistor deployment.
        serial_id : str
            Serial number of thermistor.
        """
        self.lake = lake
        self.year = year
        self.date = date
        self.location = location
        self.serial_id = serial_id

        self.md_file = self.locate_md_file()
        self.sensor = self.get_sensor_type()
        self.depth = self.get_depth()
        self.dpath_L0, self.dpath_L1, self.dpath_L2 = self.locate_data_dirs()
        

    
    # ---------- Metadata ----------
    
    def locate_md_file(self):
        """
        Locate metadata file.

        Returns
        -------
        md_path : str
            File path to metadata JSON file.
        """
        return self.MD_PATH.format(lake=self.lake, year=self.year, date=self.date, location=self.location)
    

    def open_md_file(self):
        """
        Open metadata file.

        Returns
        -------
        md : dict
            Mooring metadata.
        """
        with open(self.md_file, 'r') as f:
            md = json.load(f)

        return md
    

    def get_sensor_type(self):
        """
        Parse metadata file for sensor type.

        Returns
        -------
        sensor : str
            Type of sensor.
        """
        md = self.open_md_file()
        for i in md['instruments']:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.THERMISTORS:
                return i['instrument']
            
        raise ValueError(f'{self.serial_id} sensor not found')
    
    
    def get_mab(self):
        """
        Parse metadata file for sensor type.

        Returns
        -------
        sensor : str
            Type of sensor.
        """
        md = self.open_md_file()
        for i in md['instruments']:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.THERMISTORS:
                return i['mab']
            
        raise ValueError(f'{self.serial_id} sensor not found')
    
    
    def get_total_depth(self):
        """
        Parse metadata file for lake depth at mooring location.
        
        Returns
        -------
        total_depth : float
            Lake depth at mooring location.
        """
        md = self.open_md_file()

        return md['lake_depth']
    
    
    def get_depth(self):
        """
        Calculate depth from total depth and mab metadata.

        Returns
        -------
        depth : float
            Depth [m] of sensor.
        """
        mab = self.get_mab()
        total_depth = self.get_total_depth()

        return total_depth - mab
    

    # ---------- Navigation ----------

    def locate_data_dirs(self):
        """
        Locate data directories for L0, L1, and L2 data.

        Returns
        -------
        dpath_L0 : str
            Path to L0 data directory.
        dpath_L1 : str
            Path to L1 data directory.
        dpath_L2 : str
            Path to L2 data directory.
        """
        dpath = self.DPATH.format(lake=self.lake, location=self.location, year=self.year, date=self.date)

        return os.path.join(dpath, 'L0'), os.path.join(dpath, 'L1'), os.path.join(dpath, 'L2')


    def locate_file_L0(self):
        """
        Locate file with raw (L0) thermistor data.

        Returns
        -------
        fpath_L0 : str
            Path to L0 data file.
        """
        fpaths = glob(f'{self.dpath_L0}/*{self.serial_id}*.rsk')
        
        if len(fpaths) != 1:
            raise IndexError(f'Could not find single data file for {self.serial_id}.')
        
        return fpaths[0]
    

    # ---------- L0 to L1 ----------

    def parse_L0(self):
        """
        Load raw (L0) thermistor data into xarray Dataset.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by thermistor.
        """
        fpath_L0 = self.locate_file_L0()

        if self.sensor in  ['rbr_temp', 'rbr_duet']:
            with rsk.RSK(fpath_L0) as f:
                f.readdata()
                data = pd.DataFrame(f.data)

            data = data.rename(columns=self.COLS_MAP)
            data = data.set_index('time')
            ds = xr.Dataset.from_dataframe(data)
        else:
            raise NotImplementedError("Only rbr_temp and rbr_duet sensors are handled.")
        
        return ds
    
    
    # ---------- L1 to L2 ----------

    def organize_data_vars(self, ds):
        """
        Rename data variables.

        Parameters
        ----------
        ds : xr.Dataset
            Thermistor data.
        
        Returns
        -------
        ds : xr.Dataset
            Thermistor data with desired data variables.
        """
        if self.sensor in ['rbr_temp', 'rbr_duet']:
            vars_map = {k: v for k, v in self.VARS_MAP.items() if k in ds.data_vars}
        else:
            raise NotImplementedError('Only rbr_temp and rbr_duet sensors are handled.')

        return ds.rename_vars(vars_map)
    

    @staticmethod
    def calculate_depth(press, p_atm=10.1325):
        """
        Calculate depth of thermistor from pressure date.
        Approximate depth = pressure - air pressure

        Parameters
        ----------
        press : xr.DataArray
            Pressure [dbar].
        p_atm : float
            Atmospheric pressure [dbar].

        Returns
        -------
        depth : xr.DataArray
            Depth below water surface [m].
        """
        air_pressure = press.where(press <= p_atm).min().item()
        depth = (press - air_pressure).rename('depth')
        
        return depth

    
    
    def assign_attributes(self, ds):
        """
        Assign attributes to data variables and to dataset.
        Add depth and serial id coordinates.

        Parameters
        ----------
        ds : xr.Dataset
            Thermistor data.

        Returns
        -------
        ds : xr.Dataset
            Thermistor data with attributes.
        """
        # data variables
        for var, attrs in self.VAR_ATTRS.items():
            if var in ds:
                ds[var].attrs.update(attrs)

        # dataset
        md = self.open_md_file()
        md_xr = {
            'location': md['mooring'],
            'xsc': md['xsc'],
            'ysc': md['ysc'],
            'lake_depth': md['lake_depth'],
            'deployment': pd.to_datetime(md['deployment']).date(),
            'retrieval': pd.to_datetime(md['retrieval']).date(),
            'sensor': self.sensor,
            'serial_id': self.serial_id,
            'depth': self.depth
        }
        ds = ds.assign_attrs(md_xr)

        # add depth and serial id coordinates
        ds = ds.assign_coords(depth=self.depth, serial_id=self.serial_id)

        return ds


    def quality_assurance(self):
        """
        Run quality assurance on L1 thermistor data.

        Parameters
        ----------
        ds : xr.Dataset
            L1 thermistor data.

        Returns
        -------
        ds : xr.Dataset
            Processed (L2) thermistor data.
        """
        ds = self.organize_data_vars(ds)
        ds['depth'] = self.calculate_depth(ds['press'])
        ds = self.assign_attributes(ds)
        
        return ds


    # ---------- Writing ----------

    def write_to_nc(self, ds, level, overwrite=True):
        """
        Write xarray Dataset to .nc file.

        Parameters
        ----------
        ds : xr.Dataset
            Thermistor data.
        level : str
            L1 or L2.
        overwrite : bool
            If True, overwrite existing L1 data.

        Returns
        -------
        fpath : str
            File path to written data.
        """
        if level == 'L1':
            fpath = os.path.join(self.dpath_L1, f'{self.sensor}_{self.serial_id}_L1.nc')
        elif level == 'L2':
            fpath = os.path.join(self.dpath_L2, f'{self.sensor}_{self.serial_id}_L2.nc')
        else:
            raise ValueError('Writing level must be L1 or L2.')

        if os.path.exists(fpath) and not overwrite:
            warnings.warn(f'{fpath} already exists and overwrite = False.')
        else:
            ds.to_netcdf(fpath)

        return fpath
    

    # ---------- Data Index ----------

    def update_data_index(self):
        """
        Update data index after processing new thermistor.

        Returns
        -------
        data_index : pd.DataFrame
            Mooring data index.
        di_path : str
            File path to mooring data index.
        """
        raise NotImplementedError
    

    # ---------- Pipeline ----------

    def process(self, update=False):
        """
        Process raw (L0) thermistor data.  Convert to xarray and write to .nc (L1).
        Run quality assurance and write to .nc (L2).

        Parameters
        ----------
        update : bool
            If True, update data index with newly processed data.
        """
        ds = self.parse_L0()
        fpath_L1 = self.write_to_nc(ds, 'L1')
        ds_qa = self.quality_assurance(ds)
        fpath_L2 = self.write_to_nc(ds_qa, 'L2')
        if update:
            data_index, di_path = self.update_data_index()