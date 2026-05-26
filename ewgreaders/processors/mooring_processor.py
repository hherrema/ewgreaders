### Base class for processing mooring data

# imports
import json
from datetime import datetime
import os
from glob import glob
import xarray as xr
import pandas as pd
import numpy as np
import scipy
import pyrsktools as rsk
import dolfyn as dlfn
import warnings


class MooringProcessor:
    MD_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_md.json'
    DT_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_dt.csv'
    BATHY_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/bathymetry.nc'
    DPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/'    
    ADCPS = ['adcp']
    THERMISTORS = ['rbr_temp', 'rbr_duet']
    OXYGEN_LOGGERS = ['minidot', 'rbr_do']
    COLS_MAP_RBR_DUET = {'timestamp': 'time'}


    def __init__(self, lake, year, date, location):
        """
        Initialize MooringProcessor object.
        
        Parameters
        ----------
        lake : str
            Lake where mooring is deployed.
        year : str
            Year of mooring retrieval.
        date : str
            Date (YYYYMMDD) of mooring retrieval.
        location : str
            Location code within lake of mooring deployment.
        """
        self.lake = lake
        self.year = year
        self.date = date
        self.location = location

        self.md_file = self.locate_md_file()
        self.dpath_L0, self.dpath_L1, self.dpath_L2 = self.locate_data_dirs()
        self.deploy, self.retrieve = self.get_deploy_retrieve_dates()
        self.xsc, self.ysc = self.get_swiss_coords()


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
    

    def get_deploy_retrieve_dates(self):
        """
        Parse metadata file for depolyment and retrieval dates.

        Returns
        -------
        deploy : datetime
            Date of mooring deployment.
        retrieve : datetime
            Date of mooring retrieval.
        """
        md = self.open_md_file()

        deploy = datetime.strptime(md['deployment'], '%d.%m.%Y').date()
        retrieve = datetime.strptime(md['retrieval'], '%d.%m.%Y').date()

        return deploy, retrieve
    

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
            total_depth = bathy.sel(xsc=self.xsc, ysc=self.ysc).depth.item()
        else:
            md = self.open_md_file()
            total_depth = md['lake_depth']

        return total_depth


    def get_instruments(self, pandas=False):
        """
        Parse metadata file for all instruments.

        Parameters
        ----------
        pandas : bool
            If True, return pandas DataFrame.

        Returns
        -------
        instruments : list
            Metadata dictionaries for all instruments on mooring.
        """
        md = self.open_md_file()

        if pandas:
            return pd.DataFrame(md['instruments'])
        else:
            return md['instruments']
        

    def get_adcps(self):
        """
        Parse metadata file for ADCPs.

        Returns
        -------
        adcps : list
            Metadata dictionaries for all ADCPs on mooring.
        """
        md = self.open_md_file()

        return [i for i in md['instruments'] if i['instrument'] == 'adcp']
    

    def get_thermistors(self):
        """
        Parse metadata file for thermistors.

        Returns
        -------
        thermistors : list
            Metadata dictionaries for all thermistors on mooring.
        """
        md = self.open_md_file()

        return [i for i in md['instruments'] if i['instrument'] in self.THERMISTORS]
    

    def get_oxygen_loggers(self):
        """
        Parse metadata file for oxygen loggers.

        Returns
        -------
        oxygen_loggers : list
            Metadata dictionaries for all oxygen loggers on mooring.
        """
        md = self.open_md_file()

        return [i for i in md['instruments'] if i['instrument'] in self.OXYGEN_LOGGERS]
    

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
        dpath = self.DPATH.format(lake=self.lake, year=self.year, date=self.date, location=self.location)

        return os.path.join(dpath, 'L0'), os.path.join(dpath, 'L1'), os.path.join(dpath, 'L2')
    

    def locate_file_L0_rbr_duet(self, serial_id):
        """
        Locate file with raw (L0) rbr_duet thermistor data.

        Parameters
        ----------
        serial_id : str
            Thermistor serial ID.

        Returns
        -------
        fpath_L0 : str
            Path to L0 data file.
        """
        fpaths = glob(f'{self.dpath_L0}/*{serial_id}*.rsk')
        
        if len(fpaths) != 1:
            raise FileNotFoundError(f'Could not find single data file for {serial_id}.')
        
        return fpaths[0]
    

    def locate_file_L0_adcp(self, serial_id):
        """
        Locate file with raw (L0) ADCP data.

        Parameters
        ----------
        serial_id : str
            ADCP serial ID.

        Returns
        -------
        fpath_L0 : str
            Path to L0 data file.
        """
        fpaths = glob(f'{self.dpath_L0}/*{serial_id}*.000')
        
        if len(fpaths) != 1:
            raise FileNotFoundError(f'Could not find single data file for {serial_id}.')
        
        return fpaths[0]

    

    # ---------- L0 ----------

    def parse_L0_rbr_duet(self, serial_id):
        """
        Load raw (L0) rbr_duet thermistor data into xarray Dataset.

        Parameters
        ----------
        serial_id : str
            Thermistor serial ID.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by thermistor.
        """
        fpath_L0 = self.locate_file_L0_rbr_duet(serial_id)

        with rsk.RSK(fpath_L0) as f:
            f.readdata()
            data = pd.DataFrame(f.data)

        data = data.rename(columns=self.COLS_MAP_RBR_DUET)
        data = data.set_index('time')
        
        return xr.Dataset.from_dataframe(data)
    

    def parse_L0_adcp(self, serial_id):
        """
        Load raw (L0) ADCP data into xarray Dataset.

        Parameters
        ----------
        serial_id : str
            ADCP serial ID.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by ADCP.
        """
        fpath_L0 = self.locate_file_L0_adcp(serial_id)

        return dlfn.read(fpath_L0)
    

    # ---------- Depth Regression ----------

    @staticmethod
    def calculate_depth_rbr_duet(press, p_atm=10.1325):
        """
        Calculate depth of thermistor from pressure date.
        Approximate depth = pressure - air pressure.
        Median depth approximates given thermistor records primarily in water.

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

        return depth.median().item()


    @staticmethod
    def calculate_depth_adcp(ds):
        """
        Calculate depth from ADCP readings.  
        Median depth approximates given ADCP records primarily in water.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        depth : float
            Depth below water surface [m].
        """
        # check if ADCP reads pressure (i.e., depth)
        if 'pressure' in ds.data_vars and (ds['pressure'] != 0).any() and ds['depth'].std() > 0:
            return ds['depth'].median().item()
        else:
            raise KeyError('ADCP does not measure depth.')
        

    def create_depth_table(self):
        """
        Extract depths from instruments on moorings with pressure sensors (RBR duet, ADCP).
        Run linear regression to determine depths of other sensors.

        Returns
        -------
        depth_table : pd.DataFrame
            Depths [m] of mooring instruments.
        """
        total_depth = self.get_total_depth()
        instruments = self.get_instruments()

        # extract sensor depths
        depth_table = []
        for i in instruments:
            if i['instrument'] == 'rbr_duet':
                try:
                    ds = self.parse_L0_rbr_duet(i['serial_id'])
                    depth_sensor = self.calculate_depth_rbr_duet(ds['pressure'])
                except FileNotFoundError:
                    depth_sensor = np.nan
            elif i['instrument'] == 'adcp':
                try:
                    ds = self.parse_L0_adcp(i['serial_id'])
                    depth_sensor = self.calculate_depth_adcp(ds)
                except (FileNotFoundError, KeyError):
                    depth_sensor = np.nan
            else:
                depth_sensor = np.nan

            depth_table.append({
                'instrument': i['instrument'],
                'serial_id': i['serial_id'],
                'depth_md': total_depth - i['mab'],
                'depth_sensor': depth_sensor
            })

        depth_table = pd.DataFrame(depth_table)
        #print(depth_table)

        # use sensor values for instruments at same depth (mean is really taking the non-NaN entry)
        depth_table['depth_sensor'] = depth_table.groupby('depth_md')['depth_sensor'].transform('mean')
        
        # vertical translation
        if depth_table['depth_sensor'].nunique() == 1:
            warnings.warn("Only 1 sensor depth, applying vertical translation.")
            b = (depth_table['depth_sensor'] - depth_table['depth_md']).mean()
            depth_table['depth_t'] = depth_table['depth_md'] + b
            depth_table['depth'] = depth_table['depth_sensor'].fillna(depth_table['depth_t']).round(1)

        # linear regression
        elif depth_table['depth_sensor'].nunique() > 1:
            m, b, _, _, _ = scipy.stats.linregress(depth_table.depth_md, depth_table.depth_sensor, nan_policy='omit')
            depth_table['depth_lr'] = m * depth_table['depth_md'] + b
            depth_table['depth'] = depth_table['depth_sensor'].fillna(depth_table['depth_lr']).round(1)

        else:
            raise ValueError('No sensor depths, regression not calculated.')

        return depth_table
    

    def write_depth_table(self, depth_table):
        """
        Write table with instrument depths to .csv.

        Parameters
        ----------
        depth_table : pd.DataFrame
            Depths [m] of mooring instruments.
        """
        dt_path = self.DT_PATH.format(lake=self.lake, year=self.year, date=self.date, location=self.location)

        depth_table.to_csv(dt_path, index=False)