### Class for reading Oxygen logger data

# imports
import xarray as xr
import numpy as np
import pandas as pd
import pyrsktools as rsk
import json

from .mooring_reader import MooringReader


class O2Reader(MooringReader):
    OXYGEN_LOGGERS = ['minidot', 'rbr_do']
    COLS_DROP_MINIDOT = ['Unix Timestamp', 'Coordinated Universal Time', 'Battery', 'Q']
    COLS_MAP_MINIDOT = {
        'UTC_Date_&_Time': 'time', 
        'Temperature': 'temp', 
        'Dissolved Oxygen': 'd_oxygen_conc', 
        'Dissolved Oxygen Saturation': 'd_oxygen_sat'
    }
    COLS_MAP_RBR_DO = {'timestamp': 'time', 'dissolved_o2_saturation': 'd_oxygen_sat'}
    
    def __init__(self, serial_id, lake, location, year, date, bathy_file, datalakes=False):
        """
        Initialize O2Reader object.

        Parameters
        ----------
        serial_id : str
            Serial number of oxygen logger.
        lake : str
            Lake where oxygen logger is deployed.
        location : str
            Location code within lake of oxygen logger deployment.
        year : str
            Year of oxygen logger retrieval. 
        date : str
            Date (YYYYMMDD) of oxygen logger retrieval.
        bathy_file : str
            File path to bathymetry data.
        datalakes : bool
            Toggle whether to read from Eawag drive or DataLakes.
        """
        super.__init__(lake, location, year, date, bathy_file, datalakes)
        self.serial_id = serial_id


    def get_mab(self):
        """
        Parse metadata file for meters above bottom of oxygen logger.

        Returns
        -------
        mab : float
            Meters above bottom for oxygen logger.
        """
        with open(self.md_file, 'r') as f:
            md = json.load(f)

        instruments = md['instruments']
        for i in instruments:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.OXYGEN_LOGGERS:
                return i['mab']
            
        return np.nan
    

    def set_depth(self):
        """
        Set depth of oxygen logger.

        Returns
        -------
        depth : float
            Depth below water surface of oxygen logger.
        """
        return self.total_depth - self.mab
    

    def get_sensor_type(self):
        """
        Parse metadata file for sensor type.

        Returns
        -------
        sensor : str
            Type of oxygen logger.
        """
        with open(self.md_file, 'r') as f:
            md = json.load(f)

        instruments = md['instruments']
        for i in instruments:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.OXYGEN_LOGGERS:
                return i['instrument']
            
        return None
    

    def load_from_L0(self):
        """
        Load raw (L0) oxygen logger data into xarray Dataset.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by oxygen logger.
        """
        if self.sensor == 'minidot':
            data = self.parse_minidot_L0()
        elif self.sensor == 'rbr_do':
            data = self.parse_RBR_DO_L0()
        else:
            raise NotImplementedError("Only minidot and rbr_do sensors are handled.")
        
        data = data.set_index('time')
        ds = xr.Dataset.from_dataframe(data)
        ds = ds.assign_coords(depth=self.depth, serial_id=self.serial_id)

        return ds
    
    
    def parse_minidot_L0(self):
        """
        Parse raw (L0) data from Minidot oxygen logger.

        Returns
        -------
        data : pd.DataFrame
            Data from Minidot oxygen logger.
        """
        with open(self.fpath, 'r') as f:
                lines = [x[:-1] for x in f if len(x.split(',')) > 1]

        # extract colum names
        cols = [x.lstrip(' ') for x in lines[0].split(',')]

        data = []
        for line in lines[2:]:
            data.append([x.lstrip(' ') for x in line.split(',')])
        data = pd.DataFrame(data, columns=cols)

        data = data.drop(self.COLS_DROP_MINIDOT, axis=1)
        data = data.rename(columns=self.COLS_MAP_MINIDOT)
        data['time'] = pd.to_datetime(data['time'])
        data['temp'] = data['temp'].astype(float)
        data['d_oxygen'] = data['d_oxygen'].astype(float)
        data['d_oxygen_sat'] = data['d_oxygen_sat'].astype(float)

        return data
    

    def parse_RBR_DO_L0(self):
        """
        Parse raw (L0) data from RBR_DO oxygen logger.

        Returns
        -------
        data : pd.DataFrame
            Data from RBR_DO oxygen logger.
        """
        with rsk.RSK(self.fpath) as f:
            f.readdata()
            data = pd.DataFrame(f.data)

        data = data.rename(columns=self.COLS_MAP_RBR_DO)

        return data