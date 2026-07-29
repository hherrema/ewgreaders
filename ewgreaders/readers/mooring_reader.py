### Base class for reading mooring data

# imports
import json
from datetime import datetime
import xarray as xr
import pandas as pd
import numpy as np


class MooringReader:
    MD_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_md.json'
    DT_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_dt.csv'
    BATHY_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/bathymetry.nc'
    DPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/'
    ADCPS = ['adcp']
    THERMISTORS = ['rbr_temp', 'rbr_duet']
    OXYGEN_LOGGERS = ['minidot', 'rbr_do']

    def __init__(self, lake, date, location):
        """
        Initialize MooringReader object.

        Parameters
        ----------
        lake : str
            Lake where mooring is deployed.
        date : str
            Date (YYYY-MM-DD) of mooring retrieval.
        location : str
            Location code within lake of mooring deployment.
        """
        self.lake = lake
        self.year = str(date.year)
        self.date = date.strftime('%Y%m%d')
        self.location = location

        self.md_file = self.locate_md_file()
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
    

    # ---------- Reading ----------

    @staticmethod
    def create_instrument_chain(datasets):
        """
        Concatenate individual instrument data into single Dataset with all instruments.
        Works for thermistors and oxygen loggers.

        Parameters
        ----------
        datasets : list
            List of xarray Datasets from individual instruments.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by all instruments on mooring.
        """
        # align time samples
        ds_aligned = xr.align(*datasets, join='inner')

        return xr.concat(ds_aligned, dim='depth')
    

    @staticmethod
    def create_adcp_chain(datasets, gap=None, drop=None):
        """
        Concatenate individual ADCP data into a single Dataset with all ADCPs.

        Parameters
        ----------
        datasets : list
            List of xarray Datasets from individual ADCPs.
        gap : np.array
            Depths of gaps between ADCP ranges.
        drop : np.array
            Depth bins with known bad data to remove.

        Returns
        -------
        ds_adcp : xr.Dataset
            Dataset of data recorded by all ADCPs on mooring.
        """
        # reindex to shared time axis
        t0 = max(ds.time.values[0] for ds in datasets)
        tf = min(ds.time.values[-1] for ds in datasets)
        time_shared = datasets[0].time.sel(time=slice(t0, tf)).values
        all_adcp_aligned = [ds.reindex(time=time_shared, method="nearest") for ds in datasets]

        # concatenate datasets and average common range bins
        ds_adcp = xr.concat(all_adcp_aligned, dim='depth')
        ds_adcp = ds_adcp.sortby('depth')
        ds_adcp = ds_adcp.groupby('depth').mean()

        # don't interplate over gap from opposite looking ADCPs in double frame
        if isinstance(gap, np.ndarray):
            full_depth = np.sort(np.concatenate([ds_adcp.depth.values, gap]))
            ds_adcp = ds_adcp.reindex(depth=full_depth)

        # remove known bad depths
        if isinstance(drop, np.ndarray):
            ds_adcp = ds_adcp.drop_sel(depth=drop)
        
        return ds_adcp
    
