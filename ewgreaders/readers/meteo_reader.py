### Class for reading meteorological data

# imports
import xarray as xr


class MeteoReader:
    METEO_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/Meteo/meteo_{station}.nc'

    def __init__(self, lake, station):
        """
        Initialize MeteoReader object.

        Parameters
        ----------
        lake : str
            Lake name.
        station : str
            MeteoSwiss recording station.
        """
        self.lake = lake
        self.station = station


    # ---------- Reading ----------
    def load(self):
        """
        Load meteorological data.

        Returns
        -------
        ds = xr.Dataset
            Meteorological data.
        """
        fpath = self.METEO_PATH.format(lake=self.lake, station=self.station)

        ds = xr.open_dataset(fpath)

        return ds.load()