### Class for reading microstructure data

# imports
import os
import xarray as xr


class uReader():
    DPATH_L2 = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Microstructure/{date}/L2/'

    def __init__(self, lake, date, fname, profile_num):
        """
        Initialize uReader object.

        Parameters
        ----------
        lake : str
            Lake of microstructure profile.
        date : str
            Date (YYYY-MM-DD) of microstructure profile.
        fname : str
            File name of profile.
        profile_num : str
            Profile number within file.
        """
        self.lake = lake
        self.year = str(date.year)
        self.date = date.strftime('%Y%m%d')
        self.fname = fname
        self.profile_num = profile_num


    # ---------- Navigation ----------

    def locate_file_L2(self):
        """
        Locate file with processed (L2) microstructure data.

        Returns
        -------
        fpath_L2 : str
            Path to L2 data file.
        """
        dpath_L2 = self.DPATH_L2.format(lake=self.lake, year=self.year, date=self.date)
        
        return os.path.join(dpath_L2, f'L2_{self.fname}_down_prof{self.profile_num}.nc')
    

    # ---------- Reading ----------

    def load(self):
        """
        Load processed (L2) microstructure data.

        Returns
        -------
        ds : xr.Dataset
            Microstructure data.
        """
        fpath_L2 = self.locate_file_L2()

        ds = xr.open_dataset(fpath_L2)

        return ds.load()