### Base class for reading mooring data

# imports
import json
from datetime import datetime
import xarray as xr


class MooringReader:
    MD_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_md.json'

    def __init__(self, lake, location, year, date, bathy_file, datalakes=False):
        """
        Initialized MooringReader object.

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
        bathy_file : str
            File path to bathymetry data.
        datalakes : bool
            Toggle whether to read from Eawag drive or DataLakes.
        """
        self.lake = lake
        self.location = location
        self.year = year
        self.date = date
        self.bathy_file = bathy_file
        self.datalakes = datalakes

        self.md_file = self.locate_md_file()


    def locate_md_file(self):
        """
        Locate metadata file.

        Returns
        -------
        md_path : str
            File path to metadata JSON file.
        """
        return self.MD_PATH.format(lake=self.lake, location=self.location, year=self.year, date=self.date)


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
        with open(self.md_file, 'r') as f:
            md = json.load(f)

        xsc = md['longitude']
        ysc = md['latitude']

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
        with open(self.md_file, 'r') as f:
            md = json.load(f)

        return md['depth']
    

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
        with open(self.md_file, 'r') as f:
            md = json.load(f)

        return md['rope_extension_prcnt']
    

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
        with open(self.md_file, 'r') as f:
            md = json.load(f)

        deploy = datetime.strptime(md['deployment'], '%d.%m.%Y').date()
        retrieve = datetime.strptime(md['retrieval'], '%d.%m.%Y').date()

        return deploy, retrieve
    

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
    
    
