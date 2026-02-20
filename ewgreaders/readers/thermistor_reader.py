### Class for reading Thermistor data

# imports
import json
from glob import glob
import numpy as np
import pandas as pd
import xarray as xr
import pyrsktools as rsk

from .mooring_reader import MooringReader


class ThermistorReader(MooringReader):
    COLS_MAP = {'timestamp': 'time', 'temperature': 'temp'}

    def __init__(self, serial_id, lake, location, year, date, datalakes=False):
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
        datalakes : bool
            Toggle whether to read from Eawag drive or DataLakes.
        """
        super().__init__(lake, location, year, date, datalakes)
        self.serial_id = serial_id

        self.fpath_L0 = self.locate_data_file('L0')
        self.sensor = self.get_sensor_type()
        self.mab = self.get_mab()
        self.depth = self.set_depth()
    

    def get_sensor_type(self):
        """
        Parse metadata file for sensor type.

        Returns
        -------
        sensor : str
            Type of thermistor.
        """
        instruments = self.md['instruments']
        for i in instruments:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.THERMISTORS:
                return i['instrument']
            
        raise ValueError(f'{self.serial_id} sensor not found')
    

    def get_mab(self):
        """
        Parse metadata file for meters above bottom of thermistor.

        Returns
        -------
        mab : float
            Meters above bottom for thermistor.
        """
        instruments = self.md['instruments']
        for i in instruments:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.THERMISTORS:
                return i['mab']
            
        raise ValueError(f'{self.serial_id} meters above bottom not found.')
    

    def set_depth(self):
        """
        Set depth of thermistor.

        Returns
        -------
        depth : float
            Depth below water surface of thermistor.
        """
        return self.total_depth - self.mab
    

    def calculate_rbr_duet_depth(self, ds, p_atm=10.1325):
        """
        Calculate depth of thermistor with pressure data.
        Approximate depth as pressure - air pressure.

        Parameters
        ----------
        ds : xr.Dataset
            Thermistor data from RBR duet.
        p_atm : float
            Atmospheric pressure [dbar].

        Returns
        -------
        depth : float
            Sensor depth in lake [m].
        pressure : float
            Sensor pressure while in water [dbar].
        air_pressure : float
            Sensor pressure while in air [dbar].
        """
        # mean pressure while in air
        air_pressure = ds.where(ds.pressure <= p_atm).pressure.mean().item()
        if np.isnan(air_pressure):
            air_pressure = p_atm

        # median pressure while in water
        pressure = ds.where(ds.pressure > p_atm).pressure.median().item()

        # approximate depth as gauge pressure
        depth = pressure - air_pressure

        return depth, pressure, air_pressure


    def locate_data_file(self, level):
        """
        Locate file with thermistor data.

        Parameters
        ----------
        level : str
            Level of data.

        Returns
        -------
        fpath : str
            Path to data file.
        """
        if level == 'L0':
            fpaths = glob(f'{self.dpath_L0}/*{self.serial_id}*.rsk')
        elif level == 'L1':
            raise NotImplementedError
        elif level == 'L2':
            raise NotImplementedError
        else:
            raise ValueError("Data level must be L0, L1, or L2.")
        
        if len(fpaths) != 1:
            raise IndexError(f'Could not find single data file for {self.serial_id}.')
        
        return fpaths[0]
    

    def load_from_L0(self):
        """
        Load raw (L0) thermistor data into xarray Dataset.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by thermistor.
        """
        if self.sensor in  ['rbr_temp', 'rbr_duet']:
            with rsk.RSK(self.fpath_L0) as f:
                f.readdata()
                data = pd.DataFrame(f.data)

            data = data.rename(columns=self.COLS_MAP)
            data = data.set_index('time')

            ds = xr.Dataset.from_dataframe(data)
            ds = ds.assign_coords(depth=self.depth, serial_id=self.serial_id)
        else:
            raise NotImplementedError("Only rbr_temp and rbr_duet sensors are handled.")
        
        return ds