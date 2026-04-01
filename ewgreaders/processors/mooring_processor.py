### Base class for processing mooring data

# imports
import json
import os


class MooringProcessor:
    MD_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}_md.json'
    BATHY_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/bathymetry.nc'
    DPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/Mooring/{date}/{location}/'
    ADCPS = ['adcp']
    THERMISTORS = ['rbr_temp', 'rbr_duet']
    OXYGEN_LOGGERS = ['minidot', 'rbr_do']


    def __init__(self, lake, location, year, date):
        """
        Initialize MooringProcessor object.
        
        Parameters
        ----------
        lake : str
            Lake where mooring is deployed.
        location : str
            Location code within lake of mooring deployment.
        year : str
            Year of mooring retrieval.
        date : str
            Date (YYYYMMDD) of mooring retrieval.
        """
        self.lake = lake
        self.location = location
        self.year = year
        self.date = date


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


    def get_adcps(self):
        """
        Parse metadata file for ADCPs.

        Returns
        -------
        adcps : list
            Metadata dictionaries for all ADCPs on mooring.
        """
        md = self.open_md_file()
        return [i for i in md['instruments'] if i['instrument'] == 'adcp']
    

    def get_thermistors(self):
        """
        Parse metadata file for thermistors.

        Returns
        -------
        thermistors : list
            Metadata dictionaries for all thermistors on mooring.
        """
        md = self.open_md_file()
        return [i for i in md['instruments'] if i['instrument'] in self.THERMISTORS]
    

    def get_oxygen_loggers(self):
        """
        Parse metadata file for oxygen loggers.

        Returns
        -------
        oxygen_loggers : list
            Metadata dictionaries for all oxygen loggers on mooring.
        """
        md = self.open_md_file()
        return [i for i in md['instruments'] if i['instrument'] in self.OXYGEN_LOGGERS]
    

    def get_sensor_type(self, serial_id):
        """
        Parse metadata file for sensor type.

        Returns
        -------
        sensor : str
            Type of sensor.
        """
        md = self.open_md_file()
        for i in md['instruments']:
            if i['serial_id'] == self.serial_id and i['instrument'] in self.ADCPS + self.THERMISTORS + self.OXYGEN_LOGGERS:
                return i['instrument']
            
        raise ValueError(f'{self.serial_id} sensor not found')
    

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