### Class for reading Thermistor data

# imports
import json
import numpy as np
import pandas as pd
import xarray as xr
import pyrsktools as rsk

from .mooring_reader import MooringReader


class ThermistorReader(MooringReader):
    THERMISTORS = ['rbr_temp', 'rbr_duet']
    COLS_MAP = {'timestamp': 'time', 'temperature': 'temp'}

    def __init__(self, serial_id, lake, location, year, date, bathy_file, datalakes=False):
        """
        Initialize ThermistorReader object.

        Parameters
        ----------
        serial_id : str
            Serial number of thermistor.
        lake : str
            Lake where thermistor is deployed.
        location : str
            Location code within lake of thermistor deployment.
        year : str
            Year of thermistor retrieval. 
        date : str
            Date (YYYYMMDD) of thermistor retrieval.
        bathy_file : str
            File path to bathymetry data.
        datalakes : bool
            Toggle whether to read from Eawag drive or DataLakes.
        """
        super.__init__(lake, location, year, date, bathy_file, datalakes)
        self.serial_id = serial_id


    def get_mab(self):
        """
        Parse metadata file for meters above bottom of thermistor.

        Returns
        -------
        mab : float
            Meters above bottom for thermistor.
        """
        with open(self.md_file, 'r') as f:
            md = json.load(f)

        instruments = md['instruments']
        for i in instruments:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.THERMISTORS:
                return i['mab']
            
        return np.nan
    

    def set_depth(self):
        """
        Set depth of thermistor.

        Returns
        -------
        depth : float
            Depth below water surface of thermistor.
        """
        return self.total_depth - self.mab
    

    def get_sensor_type(self):
        """
        Parse metadata file for sensor type.

        Returns
        -------
        sensor : str
            Type of thermistor.
        """
        with open(self.md_file, 'r') as f:
            md = json.load(f)

        instruments = md['instruments']
        for i in instruments:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.THERMISTORS:
                return i['instrument']
            
        return None

    

    def load_from_L0(self):
        """
        Load raw (L0) thermistor data into xarray Dataset.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by thermistor.
        """
        if self.sensor in  ['rbr_temp', 'rbr_duet']:
            with rsk.RSK(self.fpath) as f:
                f.readdata()
                data = pd.DataFrame(f.data)

            data = data.rename(columns=self.COLS_MAP)
            data = data.set_index('time')

            ds = xr.Dataset.from_dataframe(data)
            ds = ds.assign_coords(depth=self.depth, serial_id=self.serial_id)
        else:
            raise NotImplementedError("Only rbr_temp and rbr_duet sensors are handled.")
        
        return ds