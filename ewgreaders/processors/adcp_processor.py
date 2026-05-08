### Class for processing ADCP data

# imports
import json
import os
from glob import glob
import dolfyn as dlfn
import pandas as pd
import xarray as xr
import numpy as np
import warnings


class ADCPProcessor:
    MD_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_md.json'
    DT_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_dt.csv'
    BATHY_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/bathymetry.nc'
    DPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/'
    DIPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/mooring.json'
    ADCPS = ['adcp']
    VARS_DROP = ['number']
    VARS_MAP = {
        'salinity': 'salin',
        'pressure': 'press',
        'pressure_std': 'press_std',
        'beam2inst_orientmat': 'rotmat'
    }
    VAR_ATTRS = {
        'time': {'long_name': 'Coodinated Universal Time (UTC)'},
        'dir': {'long_name': 'Reference frame (Earth)'},
        'range': {'units': 'm', 'long_name': 'Distance to center of range bin'},
        'beam': {'long_name': 'Beam reference frame'},
        'earth': {'long_name': 'Earth reference frame'},
        'inst': {'long_name': 'Instrument reference frame'},
        'depth': {'units': 'm', 'long_name': 'Depth'},
        'builtin_test_fail': {'long_name': 'Built-in test failure flag'},
        'c_sound': {'units': 'm/s', 'long_name': 'Speed of sound in water'},
        'z': {'units': 'm', 'long_name': 'Depth'},
        'pitch': {'units': 'degrees', 'long_name': 'Pitch angle'},
        'roll': {'units': 'degrees', 'long_name': 'Roll angle'},
        'heading': {'units': 'degrees', 'long_name': 'Heading angle'},
        'temp': {'units': '°C', 'long_name': 'Temperature'},
        'salin': {'units': 'PSU', 'long_name': 'Salinity'},
        'min_preping_wait': {'units': 's', 'long_name': 'Minimum pre-ping wait time (i.e., time between measurements)'},
        'heading_std': {'units': 'degrees', 'long_name': 'Heading angle standard deviation'},
        'pitch_std': {'units': 'degrees', 'long_name': 'Pitch angle standard deviation'},
        'roll_std': {'units': 'degrees', 'long_name': 'Roll angle standard deviation'},
        'press': {'units': 'dbar', 'long_name': 'Pressure'},
        'press_std': {'units': 'dbar', 'long_name': 'Pressure standard deviation'},
        'vel': {'units': 'm/s', 'long_name': 'Current velocity'},
        'amp': {'units': 'counts', 'long_name': 'Acoustic signal amplitude (intensity)'},
        'corr': {'units': 'counts', 'long_name': 'Acoustic signal correlation (beam consistency)'},
        'prcnt_gd': {'units': '%', 'long_name': 'Proportion of acceptable singal returns'}, 
        'rotmat': {'long_name': 'Rotation matrix'},
        'orientmat': {'long_name': 'Orientation matrix'},
        'serial_id': {'long_name': 'Serial ID'}
    }


    def __init__(self, lake, year, date, location, serial_id):
        """
        Initialize ADCPProcessor object.

        Parameters
        ----------
        lake : str
            Lake where ADCP is deployed.
        year : str
            Year of ADCP retrieval. 
        date : str
            Date (YYYYMMDD) of ADCP retrieval.
        location : str
            Location code within lake of ADCP deployment.
        serial_id : str
            Serial number of ADCP.
        """
        self.lake = lake
        self.year = year
        self.date = date
        self.location = location
        self.serial_id = serial_id

        self.md_file = self.locate_md_file()
        self.sensor = self.get_sensor_type()
        self.depth = self.get_depth()
        self.orientation = self.get_orientation()
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
        return self.MD_PATH.format(lake=self.lake, location=self.location, year=self.year, date=self.date)
    

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
            if i['serial_id'] == self.serial_id and i['instrument'] in self.ADCPS:
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
            if i['serial_id'] == self.serial_id and i['instrument'] in self.ADCPS:
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
            bathy = xr.open_dataset(self.BATHY_PATH.formate(lake=self.lake))
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
    

    def get_orientation(self):
        """
        Parse metadata file for ADCP orientation.
        
        Returns
        -------
        orientation : str
            Orientation {'up', 'down'} of ADCP.
        """
        md = self.open_md_file()
        for i in md['instruments']:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.ADCPS:
                if 'up' in i['comments'].lower():
                    return 'up'
                elif 'down' in i['comments'].lower():
                    return 'down'
            
        raise ValueError(f'{self.serial_id} orientation not found')
    

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
        Locate file with raw (L0) ADCP data.

        Returns
        -------
        fpath_L0 : str
            Path to L0 data file.
        """
        fpaths = glob(f'{self.dpath_L0}/*{self.serial_id}*.000')
        
        if len(fpaths) != 1:
            raise IndexError(f'Could not find single data file for {self.serial_id}.')
        
        return fpaths[0]
    

    # ---------- L0 to L1 ----------    

    def parse_L0(self):
        """
        Load raw (L0) ADCP data into xarray Dataset.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by ADCP.
        """
        fpath_L0 = self.locate_file_L0()

        ds = dlfn.read(fpath_L0)
        del ds.time.attrs['units']  # remove time units to avoid writing error

        return ds

    
    # ---------- L1 to L2 ----------        

    def calculate_depth(self, ds):
        """
        Calculate depth from ADCP depth and range bins.  
        Change depth data variable to z, assign depth coordinate, swap range dimension for depth.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with depth dimension.
        """
        ds = ds.rename({'depth': 'z'})
        
        if self.orientation == 'up':
            ds = ds.assign_coords(depth=self.depth - ds['range'])
        elif self.orientation == 'down':
            ds = ds.assign_coords(depth=self.depth + ds['range'])
        else:
            raise ValueError('ADCP orientation must be up or down.')
        
        ds = ds.swap_dims({'range': 'depth'})

        # check if ADCP range reaches lake surface or bottom
        total_depth = self.get_total_depth(from_bathy=True)
        if self.orientation == 'up':
            self.surfbot = ds['depth'].min().item() <= 0
        elif self.orientation == 'down':
            self.surfbot = ds['depth'].max().item() >= total_depth

        return ds
        
    
    def organize_data_vars(self, ds):
        """
        Drop and rename data variables.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        
        Returns
        -------
        ds : xr.Dataset
            ADCP data with desired data variables.
        """
        ds = ds.drop_vars(self.VARS_DROP)
        vars_map = {k: v for k, v in self.VARS_MAP.items() if k in ds.data_vars}

        return ds.rename_vars(vars_map)
    
    
    def assign_attributes(self, ds):
        """
        Assign attributes to data variables and to dataset.
        Add depth and serial id coordinates.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with attributes.
        """
        # add serial id coordinate
        ds = ds.assign_coords(serial_id=self.serial_id)

        # data variables
        for var, attrs in self.VAR_ATTRS.items():
            if var in ds:
                ds[var].attrs = attrs

        # dataset
        md = self.open_md_file()
        md_xr = {
            'location': md['mooring'],
            'xsc': md['xsc'],
            'ysc': md['ysc'],
            'lake_depth': md['lake_depth'],
            'bathy_depth': self.get_total_depth(from_bathy=True),
            'deployment': md['deployment'],
            'retrieval':md['retrieval'],
            'sensor': self.sensor,
            'serial_id': self.serial_id,
            'depth': self.depth,
            'orientation': self.orientation,
        }

        # include attributes already in dataset
        ds_attrs = ds.attrs
        ds_attrs.pop('serialnum', None)
        ds_attrs.pop('orientation', None)

        md_xr.update(ds_attrs)

        ds = ds.assign_attrs(md_xr)

        return ds
        

    def derive_vars(self, ds):
        """
        Process L1 ADCP data to derive depth, assign attributes.
        
        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        ds : xr.Dataset
            Processed ADCP data.
        """
        ds = self.calculate_depth(ds)
        ds = self.organize_data_vars(ds)
        ds = self.assign_attributes(ds)

        return ds

    # ---------- Quality Assurance ----------

    def qa_interface_surface(self, ds):
        """
        Filter data impacted by lake surface.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with surface interface filtered.
        """
        dist_sidelobe = self.depth * (1 - np.cos(ds.attrs['beam_angle'] * np.pi / 180))

        return ds.where(ds.depth >= dist_sidelobe, drop=True)
    
    
    def qa_interface_bottom(self, ds):
        """
        Filter data impacted by lake bottom.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with bottom interface filtered.
        """
        total_depth = self.get_total_depth(from_bathy=True)
        dist_sidelobe = (total_depth - self.depth) * (1 - np.cos(ds.attrs['beam_angle'] * np.pi / 180))
        
        return ds.where(ds.depth <= total_depth - dist_sidelobe, drop=True)
    

    def qa_min_corr(self, ds, corr_thresh=64):
        """
        Filter data with at least one beam correlation < 64.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        corr_thresh : int
            Threshold for beam correlations.

        Returns
        -------
        ds : xr.Dataset
            ADCP with minimum correlations filtered.
        """
        corr1 = ds.corr.sel(beam=1) >= corr_thresh
        corr2 = ds.corr.sel(beam=2) >= corr_thresh
        corr3 = ds.corr.sel(beam=3) >= corr_thresh
        corr4 = ds.corr.sel(beam=4) >= corr_thresh

        return ds.where(corr1 & corr2 & corr3 & corr4)
    
    
    def qa_pg14(self, ds, pg14_thresh=25):
        """
        Filter data with PG1 + PG4 < 25%.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        pg14_thresh : int
            Threshold for minimum good data percentage.

        Returns
        -------
        ds : xr.Dataset
            ADCP with good data filtered.
        """
        pg14 = ds.prcnt_gd.sel(beam=1) + ds.prcnt_gd.sel(beam=4)

        return ds.where(pg14 >= pg14_thresh)
    

    def qa_pg3(self, ds, pg3_thresh=25):
        """
        Filter data with PG3 > 25 %.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        pg3_thresh : int
            Threshold for maximum bad data percentage.

        Returns
        -------
        ds : xr.Dataset
            ADCP with bad data filtered.
        """
        pg3 = ds.prcnt_gd.sel(beam=3)

        return ds.where(pg3 <= pg3_thresh)
    

    def qa_vel_error(self, ds, velerr_thresh=0.05):
        """
        Filter data with velocity error > 0.05 m/s.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        velerr_thresh : float
            Threshold for maximum velocity error.

        Returns
        -------
        ds : xr.Dataset
            ADCP with velocity error filtered.
        """
        velerr = abs(ds.vel.sel(dir='err'))

        return ds.where(velerr <= velerr_thresh)
    

    def qa_tilt(self, ds, tilt_thresh=15):
        """
        Filter data with pitch or roll angle > 15°.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        tilt_thresh : int
            Threshold for maximum tilt.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with pitch and roll angle filtered.
        """
        pitch = abs(ds['pitch']) <= tilt_thresh
        roll = abs(ds['roll']) <= tilt_thresh

        return ds.where(pitch & roll)
    

    def qa_echo_amp_diff(self, ds, surfbot_toggle=False, ead_thresh=30):
        """
        Filter data with at least one beam with vertical echo difference between 
        consecutive bins > 30.

        Only required if beam reaches surface or bottom.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        surfbot_toggle : bool
            Toggle whether ADCP range reaches lake surface or bottom.
        ead_thresh : int
            Threshold for minimum vertical echo difference.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with vertical echo difference filtered.
        """
        if not self.surfbot and not surfbot_toggle:
            return ds
        
        echo_amp_diff = ds.amp.diff(dim='depth')
        ead1 = echo_amp_diff.sel(beam=1) <= ead_thresh
        ead2 = echo_amp_diff.sel(beam=2) <= ead_thresh
        ead3 = echo_amp_diff.sel(beam=3) <= ead_thresh
        ead4 = echo_amp_diff.sel(beam=4) <= ead_thresh

        return ds.where(ead1 & ead2 & ead3 & ead4)


    def qa_corr_stdev(self, ds, stdev_thresh=0.01, scale=100):  # threshold may be too harsh (test not from manual)
        """
        Filter data with standard deviation of 4 beams' correlations > 0.01.
        
        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        stdev_thresh : float
            Threshold for maximum standard deviation.
        scale : int
            Scale correlations to [0, 1]

        Returns
        -------
        ds : xr.Dataset
            ADCP data with correlation standard deviation filtered.
        """
        corr1 = ds.corr.sel(beam=1) / scale
        corr2 = ds.corr.sel(beam=2) / scale
        corr3 = ds.corr.sel(beam=3) / scale
        corr4 = ds.corr.sel(beam=4) / scale
        corr = xr.concat([corr1, corr2, corr3, corr4], dim='beam')
        corr_stdev = corr.std(dim='beam')

        return ds.where(corr_stdev <= stdev_thresh)
    

    def quality_assurance(self, ds):
        """
        Run quality assurance on L1 ADCP data.

        Parameters
        ----------
        ds : xr.Dataset
            L1 ADCP data.

        Returns
        -------
        ds : xr.Dataset
            Quality assured (L2) ADCP data.
        """
        ds = ds.where(~ds['builtin_test_fail'], drop=True)   # built-in qa test
        ds = self.qa_interface_surface(ds)
        ds = self.qa_interface_bottom(ds)
        ds = self.qa_min_corr(ds)
        ds = self.qa_pg14(ds)
        ds = self.qa_pg3(ds)
        ds = self.qa_vel_error(ds)
        ds = self.qa_tilt(ds)
        ds = self.qa_echo_amp_diff(ds)
        #ds = self.qa_corr_stdev(ds)

        return ds
    


    # ---------- Writing ----------

    def write_to_nc(self, ds, level, overwrite=True):
        """
        Write xarray Dataset to .nc file.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
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
        Update data index after processing new ADCP.

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
        Process raw (L0) ADCP data.  Convert to xarray and write to .nc (L1).
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