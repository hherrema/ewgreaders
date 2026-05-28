### Class for reading Thermistor data

# imports
from glob import glob
import xarray as xr


class ThermistorReader():
    DPATH_L2 = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/L2/'

    def __init__(self, lake, date, location, serial_id):
        """
        Initialize ThermistorReader object.

        Parameters
        ----------
        lake : str
            Lake where thermistor is deployed.
        date : str
            Date (YYYY-MM-DD) of thermistor retrieval.
        location : str
            Location code within lake of thermistor deployment.
        serial_id : str
            Serial number of thermistor.
        """
        self.lake = lake
        self.year = str(date.year)
        self.date = date.strftime('%Y%m%d')
        self.location = location
        self.serial_id = serial_id


    # ---------- Navigation ----------

    def locate_file_L2(self):
        """
        Locate file with processed (L2) thermistor data.
        
        Returns
        -------
        fpath_L2 : str
            Path to L2 data file.
        """
        dpath_L2 = self.DPATH_L2.format(lake=self.lake, year=self.year, date=self.date, location=self.location)
        fpaths = glob(f'{dpath_L2}/*{self.serial_id}_L2.nc')

        if len(fpaths) != 1:
            raise FileNotFoundError('Could not locate single L2 file.')
        
        return fpaths[0]
    

    # ---------- Reading ----------
    def load(self):
        """
        Load processed (L2) thermistor data.

        Returns
        -------
        ds : xr.Dataset
            Thermistor data.
        """
        fpath_L2 = self.locate_file_L2()

        ds = xr.open_dataset(fpath_L2)

        return ds.load()