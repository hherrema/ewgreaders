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
    DT_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_dt.csv'
    BATHY_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/bathymetry.nc'
    DPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/'
    DIPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/mooring.json'
    THERMISTORS = ['rbr_temp', 'rbr_duet']
    COLS_MAP = {'timestamp': 'time'}
    VARS_MAP = {'temperature': 'temp', 'pressure': 'press'}
    VAR_ATTRS = {
        'time': {'long_name': 'Coordinated Universal Time (UTC)'},
        'temp': {'units': '°C', 'long_name': 'Temperature'},
        'press': {'units': 'dbar', 'long_name': 'Pressure'},
        'depth': {'units': 'm', 'long_name': 'Depth'},
        'serial_id': {'long_name': 'Serial ID'}
    }

    
    def __init__(self, lake, year, date, location, serial_id, t_offset=None):
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
        t_offset : str
            Time offset from sensor clock to correct time (e.g., +/-HH:MM:SS).
        """
        self.lake = lake
        self.year = year
        self.date = date
        self.location = location
        self.serial_id = serial_id
        self.t_offset = t_offset

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
    

    def get_swiss_coords(self, oom=True):
        """
        Parse metadata file for Swiss coordinates of mooring location.

        Parameters
        ----------
        oom : bool
            Toggle to add order of magnitude (2, 1) to (x, y) coordinates.

        Returns
        -------
        xsc : int
            Longitude coordinate.
        ysc : int
            Latitude coordinate.
        """
        md = self.open_md_file()

        xsc = md['xsc']
        ysc = md['ysc']

        if oom:
            xsc = int(xsc + 2e6)
            ysc = int(ysc + 1e6)

        return xsc, ysc
    

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
    
    
    def get_total_depth(self, from_bathy=False):
        """
        Parse metadata file for lake depth at mooring location.

        Parameters
        ----------
        from_bathy : bool
            If True, get total depth from bathymetry file.
        
        Returns
        -------
        total_depth : float
            Lake depth at mooring location.
        """
        if from_bathy:
            bathy = xr.open_dataset(self.BATHY_PATH.format(lake=self.lake))
            xsc, ysc = self.get_swiss_coords()
            total_depth = bathy.sel(xsc=xsc, ysc=ysc).depth.item()
        else:
            md = self.open_md_file()
            total_depth = md['lake_depth']

        return total_depth
    
    
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
            total_depth = self.get_total_depth(from_bathy=True)
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
        # add depth and serial id coordinates
        ds = ds.assign_coords(depth=self.depth, serial_id=self.serial_id)

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
            'bathy_depth': self.get_total_depth(from_bathy=True),
            'deployment': md['deployment'],
            'retrieval': md['retrieval'],
            'sensor': self.sensor,
            'serial_id': self.serial_id,
            'depth': self.depth,
            't_offset': str(self.t_offset)
        }
        ds = ds.assign_attrs(md_xr)

        return ds
    

    def derive_vars(self, ds):
        """
        Process L1 thermistor data to derive depth and assign attributes.
        
        Parameters
        ----------
        ds : xr.Dataset
            Thermistor data.

        Returns
        -------
        ds : xr.Dataset
            Processed thermistor data.
        """
        ds = self.organize_data_vars(ds)
        if self.sensor == 'rbr_duet':
            ds['depth'] = self.calculate_depth(ds['press'])
        ds = self.assign_attributes(ds)

        return ds


    def correct_clock_offset(self, ds):
        """
        Apply correction to sensor clock offset.

        Parameters
        ----------
        ds : xr.Dataset
            Thermistor data.

        Returns
        -------
        ds : xr.Dataset
            Thermistor data with corrected time dimension.
        """
        ds['time'] = ds['time'] + pd.to_timedelta(self.t_offset)

        return ds
    

    def quality_assurance(self, ds):
        """
        Run quality assurance on L1 thermistor data.

        Parameters
        ----------
        ds : xr.Dataset
            L1 thermistor data.

        Returns
        -------
        ds : xr.Dataset
            Quality assured (L2) thermistor data.
        """
        if self.t_offset:
            ds = self.correct_clock_offset(ds)
        
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
        root = f'Q:/Messdaten/Aphys_Hypothesis_data/{self.lake}/'
        years = [d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))]

        di = []
        for yr in years:
            root_yr = f'Q:/Messdaten/Aphys_Hypothesis_data/{self.lake}/{yr}/Mooring/'

            dates = os.listdir(root_yr)
            for dt in dates:
                root_date = os.path.join(root_yr, dt)

                md_files = glob(f'{root_date}/*_md.json')
                for md_file in md_files:
                    loc = os.path.basename(md_file).split('_md')[0]
                    dt_path = self.DT_PATH.format(lake=self.lake, year=yr, date=dt, location=loc)
                    
                    # REMOVE ONCE DEPTH TABLES MADE FOR ALL DEPLOYMENTS
                    if not os.path.exists(dt_path):
                        continue

                    depth_table = pd.read_csv(dt_path, dtype={'serial_id': str})
                    dp_L2 = os.path.join(self.DPATH.format(lake=self.lake, year=yr, date=dt, location=loc), 'L2')
                    
                    for fp in os.listdir(dp_L2):
                        serial_id = fp.split('_L2')[0].split('_')[-1]
                        dt_sel = depth_table[depth_table['serial_id'] == serial_id].iloc[0]
                        with open(md_file, 'r') as f:
                            md = json.load(f)

                        di.append({
                            'lake': self.lake,
                            'date': pd.to_datetime(dt),
                            'location': loc,
                            'xsc': md['xsc'],
                            'ysc': md['ysc'],
                            'deploy': pd.to_datetime(md['deployment']),
                            'retrieve': pd.to_datetime(md['retrieval']),
                            'sensor': dt_sel['instrument'],
                            'serial_id': serial_id,
                            'depth': dt_sel['depth']
                        })

        di_path = self.DIPATH.format(lake=self.lake)
        with open(di_path, 'w') as f:
            json.dump(di, f, indent=2)

        return pd.DataFrame(di).sort_values(by=['date', 'location', 'depth'], ascending=True).reset_index(drop=True), di_path
    

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
        # L0 to L1
        ds = self.parse_L0()
        fpath_L1 = self.write_to_nc(ds, 'L1')

        # L1 to L2
        ds = self.derive_vars(ds)
        ds = self.quality_assurance(ds)
        fpath_L2 = self.write_to_nc(ds, 'L2')

        # data index
        if update:
            data_index, di_path = self.update_data_index()