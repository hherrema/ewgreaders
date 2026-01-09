### Base class for reading mooring data

# imports
import os
import json
from datetime import datetime
import xarray as xr


class MooringReader:
    MD_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_md.json'
    BATHY_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/bathymetry.nc'
    DPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/'
    THERMISTORS = ['rbr_temp', 'rbr_duet']
    OXYGEN_LOGGERS = ['minidot', 'rbr_do']

    def __init__(self, lake, location, year, date, datalakes=False):
        """
        Initialize MooringReader object.

        Parameters
        ----------
        lake : str
            Lake where mooring is deployed.
        location : str
            Location code within lake of mooring deployment.
        year : str
            Year of mooring retrieval.
        date : str
            Date (YYYYMMDD) of mooring retrieval.
        datalakes : bool
            Toggle whether to read from Eawag drive or DataLakes.
        """
        self.lake = lake
        self.location = location
        self.year = year
        self.date = date
        self.datalakes = datalakes

        self.md_file = self.locate_md_file()
        self.md = self.open_md_file()
        self.total_depth = self.get_total_depth()

        self.dpath_L0, self.dpath_L1, self.dpath_L2 = self.locate_data_dirs()


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
        xsc = self.md['longitude']
        ysc = self.md['latitude']

        if oom:
            xsc = int(xsc + 2e6)
            ysc = int(ysc + 1e6)

        return xsc, ysc
    
    
    def get_total_depth(self):
        """
        Parse metadata file for total depth at mooring location.

        Returns
        -------
        total_depth : int
            Depth [m] of lake at mooring locaton.
        """
        return self.md['depth']
    

    def set_total_depth(self, total_depth=None):
        """
        Set depth of lake at position of mooring.

        Parameters
        ----------
        total_depth : float
            Depth of lake at mooring position, calculated manually.
        """
        if not total_depth:
            bathy = xr.open_dataset(self.bathy_file)
            total_depth = bathy.sel(xsc=self.xsc, ysc=self.ysc).depth.item()

        return total_depth
    

    def get_rope_extension_prcnt(self):
        """
        Parse metadata file for rope extension percentage.  Rope extneds under weight of mooring.

        Returns
        -------
        rope_extension_prcnt : int
            Precentage to increase every 5 m interval by.
        """
        return self.md['rope_extension_prcnt']
    

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
        deploy = datetime.strptime(self.md['deployment'], '%d.%m.%Y').date()
        retrieve = datetime.strptime(self.md['retrieval'], '%d.%m.%Y').date()

        return deploy, retrieve
    
    
    def get_adcps(self):
        """
        Parse metadata file for ADCPs.

        Returns
        -------
        adcps : list
            Metadata dictionaries for all ADCPs on mooring.
        """
        return [i for i in self.md['instruments'] if i['instrument'] == 'adcp']
    

    def get_thermistors(self):
        """
        Parse metadata file for thermistors.

        Returns
        -------
        thermistors : list
            Metadata dictionaries for all thermistors on mooring.
        """
        return [i for i in self.md['instruments'] if i['instrument'] in self.THERMISTORS]
    

    def get_oxygen_loggers(self):
        """
        Parse metadata file for oxygen loggers.

        Returns
        -------
        oxygen_loggers : list
            Metadata dictionaries for all oxygen loggers on mooring.
        """
        return [i for i in self.md['instruments'] if i['instrument'] in self.OXYGEN_LOGGERS]


    def get_mab(self):
        """
        Parse metadata file for meters above bottom of instrument.
        """
        raise NotImplementedError("Subclasses must implement method.")
    

    def set_depth(self):
        """"
        Set depth of instrument.
        """
        raise NotImplementedError("Subclasses must implement method.")
    

    def create_instrument_chain(self, datasets):
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
    
    
