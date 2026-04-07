### Class for reading CTD data

# imports
import pandas as pd
import numpy as np
import xarray as xr
import pyrsktools as rsk
import os
from glob import glob

from .profile_reader import ProfileReader


class CTDReader():
    DPATH_L2 = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/CTD/{date}/L2/'

    def __init__(self, lake, date, fname):
        """
        Initialize CTDReader object.

        Parameters
        ----------
        lake : str
            Lake where CTD profiles.
        date : str
            Date (YYYY-MM-DD) of CTD profile.
        fanme : str
            File name of profile.
        """
        self.lake = lake
        self.year = str(date.year)
        self.date = date.strftime('%Y%m%d')
        self.fname = fname


    # ---------- Navigation ----------

    def locate_file_L2(self):
        """
        Locate file with processed (L2) CTD data.

        Returns
        -------
        fpath_L2 : str
            Path to L2 data file.
        """
        dpath_L2 = self.DPATH_L2.format(lake=self.lake, year=self.year, date=self.date)
        fpaths = glob(f'{dpath_L2}/*{self.fname}*.nc')

        if len(fpaths) != 1:
            raise FileNotFoundError('Could not locate single L2 file.')
        
        return fpaths[0]

    
    # ---------- Reading ----------
    
    def load(self):
        """
        Load processed (L2) CTD data.

        Returns
        -------
        ds : xr.Dataset
            CTD data.
        """
        fpath_L2 = self.locate_file_L2()

        return xr.open_dataset(fpath_L2)