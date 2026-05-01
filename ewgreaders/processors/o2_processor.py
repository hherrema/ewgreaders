### Class for processing Oxygen Logger data

# imports
import json
import os
from glob import glob
import pandas as pd
import pyrsktools as rsk
import xarray as xr
import warnings


class O2Processor:
    MD_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_md.json'
    DT_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_dt.csv'
    DPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/'
    DIPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/mooring.json'
    OXYGEN_LOGGERS = ['minidot', 'rbr_do']
    COLS_DT_MINIDOT = ['UTC_Date_&_Time', 'Coordinated Universal Time']
    COLS_INT_MINIDOT = ['Unix Timestamp']
    COLS_FLOAT_MINIDOT = ['Battery', 'Temperature', 'Dissolved Oxygen', 'Dissolved Oxygen Saturation', 'Q']
    COLS_MAP_MINIDOT = {'UTC_Date_&_Time': 'time'}
    VARS_DROP_MINIDOT = ['Unix Timestamp', 'Coordinated Universal Time', 'Battery', 'Q']
    VARS_MAP_MINIDOT = {
        'Temperature': 'temp', 
        'Dissolved Oxygen': 'do2_conc', 
        'Dissolved Oxygen Saturation': 'do2_sat'
    }
    COLS_MAP_RBR_DO = {'timestamp': 'time', 'dissolved_o2_saturation': 'do2_sat'}
    VAR_ATTRS = {
        'time': {'long_name': 'Coodinated Universal Time (UTC)'},
        'do2_conc': {'units': 'mg/l', 'long_name': 'Dissolved Oxygen Concentration'},
        'do2_sat': {'units': '%', 'long_name': 'Dissolved Oxygen Saturation'},
    }


    def __init__(self, lake, year, date, location, serial_id):
        """
        Initialize O2Processor object.

        Parameters
        ----------
        lake : str
            Lake where oxygen logger is deployed.
        year : str
            Year of oxygen logger retrieval. 
        date : str
            Date (YYYYMMDD) of oxygen logger retrieval.
        location : str
            Location code within lake of oxygen logger deployment.
        serial_id : str
            Serial number of oxygen logger.
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
            if i['serial_id'] == self.serial_id and i['instrument'] in self.OXYGEN_LOGGERS:
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
            if i['serial_id'] == self.serial_id and i['instrument'] in self.OXYGEN_LOGGERS:
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
    
    
    def get_depth(self, dt=True):
        """
        Parse depth table for instrument depth.

        Parameters
        ----------
        dt : bool
            If False, calculate depth from total depth and mab metadata.

        Returns
        -------
        depth : float
            Depth [m] of sensor.
        """
        if dt:
            dt_path = self.DT_PATH.format(lake=self.lake, year=self.year, date=self.date, location=self.location)
            depth_table = pd.read_csv(dt_path, dtype={'serial_id': str})
            depth = depth_table[depth_table['serial_id'] == self.serial_id].iloc[0].depth
        else:
            mab = self.get_mab()
            total_depth = self.get_total_depth()
            depth = total_depth - mab

        return depth
    

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
        Locate file with raw (L0) oxygen logger data.

        Returns
        -------
        fpath_L0 : str
            Path to L0 data file.
        """
        if self.sensor == 'minidot':
            fpath_L0 = f'{self.dpath_L0}/7450-{self.serial_id}/Cat.txt'
        elif self.sensor == 'rbr_do':
            fpaths = glob(f'{self.dpath_L0}/*{self.serial_id}*.rsk')
            if len(fpaths) != 1:
                raise IndexError(f'Could not find single data file for {self.serial_id}.')
            fpath_L0 = fpaths[0]
        else:
            raise NotImplementedError("Only minidot and rbr_do sensors are handled.")
        
        return fpath_L0
    

    # ---------- L0 to L1 ----------
    
    def parse_minidot_L0(self, fpath_L0):
        """
        Parse raw (L0) data from Minidot oxygen logger.

        Parameters
        ----------
        fpath_L0 : str
            File path to raw (L0) Minidot oxygen logger data.

        Returns
        -------
        data : pd.DataFrame
            Data from Minidot oxygen logger.
        """
        with open(fpath_L0, 'r') as f:
            lines = [x[:-1] for x in f if len(x.split(',')) > 1]

        # extract column names
        cols = [x.lstrip(' ') for x in lines[0].split(',')]

        data = []
        for line in lines[2:]:
            data.append([x.lstrip(' ') for x in line.split(',')])
        data = pd.DataFrame(data, columns=cols)

        # cast to proper datatypes, map to time dimension
        data[self.COLS_DT_MINIDOT] = data[self.COLS_DT_MINIDOT].apply(pd.to_datetime)
        data[self.COLS_INT_MINIDOT] = data[self.COLS_INT_MINIDOT].astype(int)
        data[self.COLS_FLOAT_MINIDOT] = data[self.COLS_FLOAT_MINIDOT].astype(float)
        data = data.rename(columns=self.COLS_MAP_MINIDOT)

        return data
    
    
    def parse_RBR_DO_L0(self, fpath_L0):
        """
        Parse raw (L0) data from RBR_DO oxygen logger.

        Parameters
        ----------
        fpath_L0 : str
            File path to raw (L0) RBR_DO oxygen logger data.

        Returns
        -------
        data : pd.DataFrame
            Data from RBR_DO oxygen logger.
        """
        raise NotImplementedError
        with rsk.RSK(fpath_L0) as f:
            f.readdata()
            data = pd.DataFrame(f.data)

        data = data.rename(columns=self.COLS_MAP_RBR_DO)

        return data
    

    def parse_L0(self):
        """
        Load raw (L0) oxygen logger data into xarray Dataset.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by oxygen logger.
        fpath_L0 : str
            File path to raw (L0) oxygen logger data.
        """
        fpath_L0 = self.locate_file_L0()

        if self.sensor == 'minidot':
            data = self.parse_minidot_L0(fpath_L0)
        elif self.sensor == 'rbr_do':
            data = self.parse_RBR_DO_L0(fpath_L0)
        else:
            raise NotImplementedError("Only minidot and rbr_do sensors are handled.")
        
        data = data.set_index('time')
        ds = xr.Dataset.from_dataframe(data)

        return ds, fpath_L0

    
    # ---------- L1 to L2 ----------

    def organize_data_vars(self, ds):
        """
        Drop and rename data variables.

        Parameters
        ----------
        ds : xr.Dataset
            Oxygen logger data.
        
        Returns
        -------
        ds : xr.Dataset
            Oxygen logger data with desired data variables.
        """
        if self.sensor == 'minidot':
            ds = ds.drop_vars(self.VARS_DROP_MINIDOT)
            vars_map = {k: v for k, v in self.VARS_MAP_MINIDOT.items() if k in ds.data_vars}
        else:
            raise NotImplementedError('Only minidot sensors are handled.')

        return ds.rename_vars(vars_map)
    
    
    def assign_attributes(self, ds):
        """
        Assign attributes to data variables and to dataset.
        Add depth and serial id coordinates.

        Parameters
        ----------
        ds : xr.Dataset
            Oxygen logger data.

        Returns
        -------
        ds : xr.Dataset
            Oxygen logger data with attributes.
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
    

    def quality_assurance(self, ds):
        """
        Run quality assurance on L1 oxygen logger data.

        Parameters
        ----------
        ds : xr.Dataset
            L1 oxygen logger data.

        Returns
        -------
        ds : xr.Dataset
            Processed (L2) oxygen logger data.
        """
        ds = self.organize_data_vars(ds)
        ds = self.assign_attributes(ds)
        
        return ds
    

    # ---------- Writing ----------

    def write_to_nc(self, ds, level, overwrite=True):
        """
        Write xarray Dataset to .nc file.

        Parameters
        ----------
        ds : xr.Dataset
            Oxygen logger data.
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
        Update data index after processing new oxygen logger.

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
        Process raw (L0) oxygen logger data.  Convert to xarray and write to .nc (L1).
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