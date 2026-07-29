### Class for reading bathymetry data

# imports
import xarray as xr


class BathyReader:
    BATHY_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/bathymetry.nc'

    def __init__(self, lake):
        """
        Initialize BathyReader object.

        Parameters
        ----------
        lake : str
            Lake name.
        """
        self.lake = lake


    # ---------- Reading ----------

    def load(self):
        """
        Load bathymetry data.
        
        Returns
        -------
        bathy : xr.Dataset
            Bathymetry data.
        """
        bathy_path = self.BATHY_PATH.format(lake=self.lake)

        bathy = xr.open_dataset(bathy_path)

        return bathy.load()