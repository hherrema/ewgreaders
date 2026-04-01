### Class for processing ADCP data

# imports
import json
import os
from glob import glob
import dolfyn as dlfn
import warnings


class ADCPProcessor:
    MD_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_md.json'
    DPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/'
    ADCPS = ['adcp']


    def __init__(self, lake, location, year, date, serial_id):
        """
        Initialize ADCPProcessor object.

        Parameters
        ----------
        lake : str
            Lake where ADCP is deployed.
        location : str
            Location code within lake of ADCP deployment.
        year : str
            Year of ADCP retrieval. 
        date : str
            Date (YYYYMMDD) of ADCP retrieval.
        serial_id : str
            Serial number of ADCP.
        """
        self.lake = lake
        self.location = location
        self.year = year
        self.date = date
        self.serial_id = serial_id

        self.sensor = self.get_sensor_type()
        self.depth = self.get_depth()
        self.dpath_L0, self.dpath_L1, self.dpath_L2 = self.locate_data_dirs()
        

    
    # ---------- Metadata ----------
    
    def locate_md_file(self):
        """
        Locate metadata file.

        Returns
        -------
        md_path : str
            File path to metadata JSON file.
        """
        return self.MD_PATH.format(lake=self.lake, location=self.location, year=self.year, date=self.date)
    

    def open_md_file(self):
        """
        Open metadata file.

        Returns
        -------
        md : dict
            Mooring metadata.
        """
        with open(self.md_file, 'r') as f:
            md = json.load(f)

        return md
    

    def get_sensor_type(self):
        """
        Parse metadata file for sensor type.

        Returns
        -------
        sensor : str
            Type of sensor.
        """
        md = self.open_md_file()
        for i in md['instruments']:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.OXYGEN_LOGGERS:
                return i['instrument']
            
        raise ValueError(f'{self.serial_id} sensor not found')
    

    def get_mab(self):
        """
        Parse metadata file for sensor type.

        Returns
        -------
        sensor : str
            Type of sensor.
        """
        md = self.open_md_file()
        for i in md['instruments']:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.OXYGEN_LOGGERS:
                return i['mab']
            
        raise ValueError(f'{self.serial_id} sensor not found')
    
    
    def get_total_depth(self):
        """
        Parse metadata file for lake depth at mooring location.
        
        Returns
        -------
        total_depth : float
            Lake depth at mooring location.
        """
        md = self.open_md_file()

        return md['depth']
    
    
    def get_depth(self):
        """
        Calculate depth from total depth and mab metadata.

        Returns
        -------
        depth : float
            Depth [m] of sensor.
        """
        mab = self.get_mab()
        total_depth = self.get_total_depth()

        return total_depth - mab
    

    # ---------- Navigation ----------

    def locate_data_dirs(self):
        """
        Locate data directories for L0, L1, and L2 data.

        Returns
        -------
        dpath_L0 : str
            Path to L0 data directory.
        dpath_L1 : str
            Path to L1 data directory.
        dpath_L2 : str
            Path to L2 data directory.
        """
        dpath = self.DPATH.format(lake=self.lake, location=self.location, year=self.year, date=self.date)

        return os.path.join(dpath, 'L0'), os.path.join(dpath, 'L1'), os.path.join(dpath, 'L2')


    def locate_file_L0(self):
        """
        Locate file with raw (L0) ADCP data.

        Returns
        -------
        fpath_L0 : str
            Path to L0 data file.
        """
        fpaths = glob(f'{self.dpath_L0}/*{self.serial_id}*.000')
        
        if len(fpaths) != 1:
            raise IndexError(f'Could not find single data file for {self.serial_id}.')
        
        return fpaths[0]
    

    # ---------- L0 to L1 ----------    

    def parse_L0(self):
        """
        Load raw (L0) ADCP data into xarray Dataset.

        Returns
        -------
        ds : xr.Dataset
            Dataset of data recorded by ADCP.
        """
        fpath_L0 = self.locate_file_L0()

        return dlfn.read(fpath_L0)

    
    # ---------- L1 to L2 ----------

    def quality_assurance(self):
        """
        Run quality assurance on L1 ADCP data.
        """
        raise NotImplementedError
    


    # ---------- Writing ----------

    def write_to_nc(self, ds, level, overwrite=True):
        """
        Write xarray Dataset to .nc file.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        level : str
            L1 or L2.
        overwrite : bool
            If True, overwrite existing L1 data.

        Returns
        -------
        fpath : str
            File path to written data.
        """
        if level == 'L1':
            fpath = os.path.join(self.dpath_L1, f'{self.sensor}_{self.serial_id}_L1.nc')
        elif level == 'L2':
            fpath = os.path.join(self.dpath_L2, f'{self.sensor}_{self.serial_id}_L2.nc')
        else:
            raise ValueError('Writing level must be L1 or L2.')

        if os.path.exists(fpath) and not overwrite:
            warnings.warn(f'{fpath} already exists and overwrite = False.')
        else:
            ds.to_netcdf(fpath)

        return fpath
    

    # ---------- Pipeline ----------

    def process(self):
        """
        Process raw (L0) ADCP data.  Convert to xarray and write to .nc (L1).
        Run quality assurance and write to .nc (L2).
        """
        ds = self.parse_L0()
        fpath_L1 = self.write_to_nc(ds, 'L1')
        ds_qa = self.quality_assurance(ds)
        fpath_L2 = self.write_to_nc(ds_qa, 'L2')