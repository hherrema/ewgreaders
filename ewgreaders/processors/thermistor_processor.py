### Class for processing Thermistor data

# imports
import pyrsktools as rsk
import xarray as xr
import pandas as pd
from glob import glob
import os
import warnings

from .mooring_processor import MooringProcessor


class ThermistorProcessor(MooringProcessor):
    COLS_MAP = {'timestamp': 'time', 'temperature': 'temp'}

    
    def __init__(self, serial_id, lake, location, year, date):
        """
        Initialize ThermistorProcessor object.

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
        """
        super().__init__(lake, location, year, date)
        self.serial_id = serial_id

        self.fpath_L0 = self.locate_data_file()
        self.sensor = self.get_sensor_type()


    def locate_data_file(self):
        """
        Locate file with thermistor data.

        Returns
        -------
        fpath : str
            Path to data file.
        """
        fpaths = glob(f'{self.dpath_L0}/*{self.serial_id}*.rsk')
        
        if len(fpaths) != 1:
            raise IndexError(f'Could not find single data file for {self.serial_id}.')
        
        return fpaths[0]
    

    def parse_L0(self):
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
    

    def write_to_L1(self, ds, overwrite=True):
        """
        Write xarray Dataset to L1 .nc file.

        Parameters
        ----------
        ds : xr.Dataset
            Thermistor data.
        overwrite : bool
            If True, overwrite existing L1 data.

        Returns
        -------
        fpath_L1 : str
            File path to L1 data.
        """
        fpath_L1 = os.path.join(self.dpath_L1, f'{self.sensor_type}_{self.serial_id}_L1.nc')

        if os.path.exists(fpath_L1) and not overwrite:
            warnings.warn(f'{fpath_L1} already exists and overwrite = False.')
        else:
            ds.to_netcdf(fpath_L1)

        return fpath_L1
    
    def write_to_L2(self, ds, overwrite=True):
        """
        Write xarray Dataset to L2 .nc file.

        Parameters
        ----------
        ds : xr.Dataset
            Thermistor data.
        overwrite : bool
            If True, overwrite existing L1 data.

        Returns
        -------
        fpath_L2 : str
            File path to L2 data.
        """
        fpath_L2 = os.path.join(self.dpath_L2, f'{self.sensor_type}_{self.serial_id}_L2.nc')

        if os.path.exists(fpath_L2) and not overwrite:
            warnings.warn(f'{fpath_L2} already exists and overwrite = False.')
        else:
            ds.to_netcdf(fpath_L2)

        return fpath_L2

    
    def process(self):
        """
        Process raw (L0) thermistor data.
        """