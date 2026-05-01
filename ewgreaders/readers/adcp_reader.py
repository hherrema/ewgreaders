### Class for reading ADCP data

# imports
import xarray as xr
from glob import glob

class ADCPReader():
    DPATH_L2 = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/L2/'

    def __init__(self, lake, year, date, location, serial_id):
        """
        Initialize ADCPReader object.

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


    # ---------- Navigation ----------

    def locate_file_L2(self):
        """
        Locate file with processed (L2) ADCP data.
        
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
        Load processed (L2) ADCP data.

        Returns
        -------
        ds : xr.Dataset
            ADCP data.
        """
        fpath_L2 = self.locate_file_L2()

        return xr.open_dataset(fpath_L2)
