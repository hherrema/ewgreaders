### Class for processing CTD data

# imports
import json
import os
import pandas as pd
import xarray as xr
import warnings
import gsw as sw
from glob import glob
import math


class CTDProcessor:
    MD_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/CTD/{date}/{fname}_md.json'
    DPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/{year}/CTD/{date}/'
    DIPATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/ctd.json'
    #SID_BRAND_MAP = {'1807': 'Sea&Sun', '2023': 'Sea&Sun', '66131': 'RBR'}
    #COLS_DROP_RBR = ['temperature1', 'speed_of_sound']
    #COLS_MAP_RBR = {}
    VARS_DROP_SEASUN = ['DO_ml', 'pH', 'Vbatt', 'IntD', 'IntT']
    VARS_MAP_SEASUN = {
        'Press': 'press',
        'Temp': 'temp',
        'Cond': 'cond',
        'CAP25': 'cond25',
        'DO_mg': 'do2_conc',
        'sat': 'do2_sat',
        'Turb': 'turb',
        'pH_Tc': 'pH',
        'Chl_A': 'chl_a',
        'Redox': 'redox'
    }
    VAR_ATTRS = {
        'time': {'long_name': 'Coordinated Universal Time (UTC)'},
        'press': {'units': 'dbar', 'long_name': 'Pressure'},
        'temp': {'units': '°C', 'long_name': 'Temperature'},
        'cond': {'units': 'mS/cm', 'long_name': 'Conductivity'},
        'cond25': {'units': 'mS/cm', 'long_name': 'Conductivity corrected to 25°C'},
        'do2_conc': {'units': 'mg/l', 'long_name': 'Dissolved Oxygen Concentration'},
        'do2_sat': {'units': '%', 'long_name': 'Dissolved Oxygen Saturation'},
        'turb': {'units': 'FTU', 'long_name': 'Turbidity'},              # FTU = NTU
        'pH': {'units': 'pH', 'long_name': 'Potential of Hydrogen'},
        'chl_a': {'units': 'µg/L', 'long_name': 'Chlorophyll A Concentration'},
        'redox': {'units': 'mV', 'long_name': 'Oxidation-Reduction Potential'},
        'depth': {'units': 'm', 'long_name': 'Depth'},
        'salin': {'units': 'PSU', 'long_name': 'Salinity'},              # PSU = ppt
        'cond20': {'units': 'mS/cm', 'long_name': 'Conductivity corrected to 20°C'},
        'rho': {'units': 'kg/m3', 'long_name': 'Density'},
        'ptemp': {'units': '°C', 'long_name': 'Potential Temperature'},
        'prho': {'units': 'kg/m3', 'long_name': 'Potential Density'},
        'fall_speed': {'units': 'm/s', 'long_name': 'Falling Speed'}
    }

    def __init__(self, lake, year, date, fname):
        """
        Initialize CTDProcessor object.

        Parameters
        ----------
        lake : str
            Lake where CTD profiles.
        year : str
            Year of CTD profile. 
        date : str
            Date (YYYYMMDD) of CTD profile.
        fname : str
            File name of raw (L0) data.
        """
        self.lake = lake
        self.year = year
        self.date = date
        self.fname = fname

        self.md_file = self.locate_md_file()
        self.sensor = self.get_sensor_type()
        self.serial_id = self.get_serial_id()
        self.lat, self.lon = self.ch1903_to_latlng()
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
        return self.MD_PATH.format(lake=self.lake, year=self.year, date=self.date, fname=self.fname)
    
    
    def open_md_file(self):
        """
        Open metadata file.

        Returns
        -------
        md : dict
            Profile metadata.
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

        return md['sensor']

    
    def get_serial_id(self):
        """
        Parse metadata file for sensor serial number.

        Returns
        -------
        serial_id : str
            Sensor serial number.
        """
        md = self.open_md_file()

        return md['serial_id']
    

    def ch1903_to_latlng(self):
        """
        Convert CH1903 Swiss coordinates to latitude and longitude.

        Returns
        -------
        lat : float
            Latitude [°].
        lon : float
            Longitude [°].
        """
        md = self.open_md_file()
        xsc = md['xsc']
        ysc = md['ysc']

        # remove leading 2 (xsc) and 1 (ysc)
        if xsc > 2e6:
            xsc = xsc - int(2e6)
        if ysc > 1e6:
            ysc = ysc - int(1e6)

        x = (xsc - 600000) / 1000000
        y = (ysc - 200000) / 1000000

        lat = (16.9023892 
            + 3.238272 * y 
            - 0.270978 * x ** 2 
            - 0.002528 * y ** 2 
            - 0.0447 * x ** 2 * y 
            - 0.014 * y ** 3)
        
        lng = (2.6779094 
            + 4.728982 * x 
            + 0.791484 * x * y 
            + 0.1306 * x * y ** 2 
            - 0.0436 * x ** 3)
        
        lat = (lat * 100) / 36
        lng = (lng * 100) / 36

        return lat, lng
    

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
        dpath = self.DPATH.format(lake=self.lake, year=self.year, date=self.date)

        return os.path.join(dpath, 'L0'), os.path.join(dpath, 'L1'), os.path.join(dpath, 'L2')
    

    def locate_file_L0(self):
        """
        Locate file with raw (L0) CTD data.

        Returns
        -------
        fpath_L0 : str
            Path to L0 data file.
        """
        if self.sensor == 'sea&sun':
            fpath_L0 = f'{self.dpath_L0}/{self.fname}.TOB'
        else:
            raise NotImplementedError('Only sea&sun sensors are handled.')
        
        return fpath_L0
    

    # ---------- L0 to L1 ----------
    """
    def separate_profiles_RBR_L0(self):
        Load raw (L0) CTD from RBR CTD, separate profiles, and write to individual files.

        Returns
        -------
        profiles : list
            List of output file paths and xarray Datasets for each profile.

        brand = self.SID_BRAND_MAP[self.serial_id]
        if brand != 'RBR':
            raise ValueError("Only RBR CTD L0 data requires profile separation.")
        
        froot, ext = os.path.splitext(self.fpath)
        if ext != '.rsk':
            raise ValueError("Only RBR CTD L0 .rsk files require profile separation.")
        
        with rsk.RSK(self.fpath) as f:
            f.readdata()
            data = f.data

            # find indices of each profile
            f.computeprofiles()
            profiles_idx = f.getprofilesindices(direction='both')

        # separate profiles and write out
        profiles = []
        air_idx_start = 0
        for i, p_idx in enumerate(profiles_idx):
            # calculate air pressure prior to downcast
            air_pressure = np.mean(data[air_idx_start: min(p_idx)]['pressure'])
            air_idx_start = max(p_idx) + 1

            p_data = pd.DataFrame(data[p_idx])
            p_data = p_data.set_index('timestamp')
            ds = xr.Dataset.from_dataframe(p_data)
            ds = ds.assign_coords(air_pressure=air_pressure)
            out_path = froot + f'_{i+1}.nc'
            profiles.append((out_path, ds))

            ds.to_netcdf(out_path, mode='w', format='NETCDF4')

        return profiles
    """

    def parse_sea_and_sun_L0(self, fpath_L0):
        """
        Parse raw (L0) data from Sea&Sun CTD.

        Parameters
        ----------
        fpath_L0 : str
            File path to raw (L0) Sea&Sun CTD data.

        Returns
        -------
        data : pd.DataFrame
            Data from Sea&Sun CTD.
        """
        with open(fpath_L0, encoding='latin1', errors='ignore') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if 'Lines' in line:
                break
        else:
            raise ValueError("Start of data not found within file.")
        
        cols = lines[i+2].replace(";", "").split()[1:]
        #units = lines[i+3].replace(";", "").replace("[", "").replace("]", "").split()

        data = pd.read_csv(fpath_L0, sep=r'\s+', header=None, skiprows=i+5, names=cols, engine='python', encoding='cp1252')
        data['time'] = pd.to_datetime(data['IntD'] + " " + data['IntT'], format='%d.%m.%Y %H:%M:%S.%f')
        #data['time_local'] = data['time'].dt.tz_localize('UTC').dt.tz_convert('Europe/Zurich')

        return data
    
    
    def parse_seabird_L0(self, fpath_L0):
        """
        Parse raw (L0) data from Seabird CTD.

        Parameters
        ----------
        fpath_L0 : str
            File path to raw (L0) Seabird CTD data.
        """
        raise NotImplementedError
    

    def parse_RBR_L0(self, fpath_L0):
        """
        Parse raw (L0) data from RBR CTD.

        Parameters
        ----------
        fpath_L0 : str
            File path to raw (L0) RBR CTD data.
        """
        raise NotImplementedError
    

    def parse_L0(self):
        """
        Load raw (L0) CTD data into xarray Dataset.

        Returns
        -------
        ds : xr.Dataset
            CTD data.
        """
        fpath_L0 = self.locate_file_L0()

        if self.sensor == 'sea&sun':
            data = self.parse_sea_and_sun_L0(fpath_L0)
        else:
            raise NotImplementedError("Only sea&sun sensors are handled.")
        
        data = data.set_index('time')
        ds = xr.Dataset.from_dataframe(data)

        return ds, fpath_L0
    

    # ---------- L1 to L2 ----------

    def extract_downcast(self, ds):
        """
        Parse profile for downcast. Exclude 0.5 dbar at surface and bottom.

        Parameters
        ----------
        ds : xr.Dataset
            CTD profile.

        Returns
        -------
        ds_downcast : xr.Dataset
            CTD downcast.
        air_pressure : float
            Air pressure [dbar].
        """
        # surface index from air pressure
        air_pressure = ds['Press'].min().item()
        idx_surface = (ds['Press'] > air_pressure + 0.5).argmax(dim='time').item()

        # bottom index from max pressure
        max_pressure = ds['Press'].max().item()
        idx_bottom = (ds['Press'] < max_pressure - 0.5).argmin(dim='time').item()

        return ds.isel(time=slice(idx_surface, idx_bottom)), air_pressure
    
    
    def organize_data_vars(self, ds):
        """
        Drop and rename data variables.

        Parameters
        ----------
        ds : xr.Dataset
            CTD downcast.
        
        Returns
        -------
        ds : xr.Dataset
            CTD data with desired data variables.
        """
        if self.sensor == 'sea&sun':
            ds = ds.drop_vars(self.VARS_DROP_SEASUN)
            vars_map = {k: v for k, v in self.VARS_MAP_SEASUN.items() if k in ds.data_vars}
        else:
            raise NotImplementedError('Only sea&sun sensors are handled.')

        return ds.rename_vars(vars_map)
    

    @staticmethod
    def salinity(temp, cond, cond_coef=0.874e-3):
        """
        Calculate salinity from conductivity.
        
        Parameters
        ----------
        temp : array_like
            Water temperature [°C]
        cond : array_like
            Conductivity [mS/cm].               # CHECK UNITS FOR EACH CTD
        cond_coef : float
            Coefficient to calculate salinity from conductivity at 20°C.

        Returns
        -------
        salin : array_like
            Salinity [PSU].
        cond20 : array_like
            Conductivity corrected to 20°C.
        """
        temp_correction = (1.8626 
                           - 0.052908 * temp 
                           + 0.00093057 * temp ** 2 
                           - 6.78e-6 * temp ** 3)
        
        cond20 = temp_correction * cond * 1e3    # conductivity at 20 °C
        salin =  cond_coef * cond20

        return salin, cond20
        

    @staticmethod
    def density(temp, salin):
        """
        Calculate water density.

        Parameters
        ----------
        temp : array_like
            Water temperature [°C].
        salin : array_like
            Water salinity [PSU].

        Returns
        -------
        rho : array_like
            Water density [kg/m3].
        """
        rho = 1e3 * (0.9998395 
                    + 6.7914e-5 * temp 
                    - 9.0894e-6 * temp ** 2 
                    + 1.0171e-7 * temp ** 3 
                    - 1.2846e-9 * temp ** 4 
                    + 1.1592e-11 * temp ** 5 
                    - 5.0125e-14 * temp ** 6 
                    + salin * (8.181e-4 
                                - 3.85e-6 * temp 
                                + 4.96e-8 * temp ** 2))
        
        return rho
    

    @staticmethod
    def calculate_depth(press, air_pressure, rho, lat):
        """
        Calculate depth below water surface.

        Parameters
        ----------
        press : array_like
            Pressure [dbar].
        air_pressure : float
            Air pressure [dbar].
        rho : array_like
            Water density [kg/m3]
        lat : float
            Latitude [°].
        
        Returns
        -------
        depth : array_like
            Depth below water surface [m].
        """
        pressure_adjusted = press - air_pressure

        if math.isnan(lat):
            g = 9.81
        else:
            g = sw.grav(lat, pressure_adjusted)

        depth = 1e4 * pressure_adjusted / (rho * g)

        return depth
    

    @staticmethod
    def potential_temperature(temp, salin, press, air_pressure, p_ref=0):
        """
        Calculate potential temperature.

        Parameters
        ----------
        temp : array_like
            Water temperature [°C].
        salin : array_like
            Water salinity [PSU].
        press : array_like
            Pressure [dbar].
        air_pressure : float
            Air pressure [dbar].
        p_ref : float
            Reference pressure [dbar].
        
        Returns
        -------
        pt : array_like
            Potential temperature [°C].
        """
        pressure_adjusted = press - air_pressure
        
        return sw.pt_from_t(salin, temp, pressure_adjusted, p_ref)
    
    
    def assign_attributes(self, ds):
        """
        Assign attributes to data variables and to dataset.

        Parameters
        ----------
        ds : xr.Dataset
            CTD data.

        Returns
        -------
        ds : xr.Dataset
            CTD data with attributes.
        """
        # data variables
        for var, attrs in self.VAR_ATTRS.items():
            if var in ds:
                ds[var].attrs.update(attrs)

        # dataset
        md = self.open_md_file()
        md_xr = {k: (str(v) if isinstance(v, bool) else v) for k, v in md.items()}
        ds = ds.assign_attrs(md_xr)

        return ds

    
    def quality_assurance(self, ds):
        """
        Run quality assurance on L1 CTD data.

        Parameters
        ----------
        ds : xr.Dataset
            L1 CTD data.

        Returns
        -------
        ds : xr.Dataset
            Processed (L2) CTD data.
        """
        ds, air_pressure = self.extract_downcast(ds)
        ds = self.organize_data_vars(ds)
        ds[['salin', 'cond20']] = self.salinity(ds['temp'], ds['cond'])
        ds['rho'] = self.density(ds['temp'], ds['salin'])
        ds['depth'] = self.calculate_depth(ds['press'], air_pressure, ds['rho'], self.lat)
        ds['ptemp'] = self.potential_temperature(ds['temp'], ds['salin'], ds['press'], air_pressure)
        ds['prho'] = self.density(ds['ptemp'], ds['salin'])
        ds['fall_speed'] = ds['depth'].differentiate("time") * 1.e9    # 1e9 to convert fron ns to s
        ds = ds.swap_dims({'time': 'depth'})
        ds = self.assign_attributes(ds)
        
        return ds
    

    # ---------- Writing ----------

    def write_to_nc(self, ds, level, overwrite=True):
        """
        Write xarray Dataset to .nc file.

        Parameters
        ----------
        ds : xr.Dataset
            CTD data.
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
            fpath = os.path.join(self.dpath_L1, f'{self.sensor}_{self.serial_id}_{self.fname}_L1.nc')
        elif level == 'L2':
            fpath = os.path.join(self.dpath_L2, f'{self.sensor}_{self.serial_id}_{self.fname}_L2.nc')
        else:
            raise ValueError('Writing level must be L1 or L2.')

        if os.path.exists(fpath) and not overwrite:
            warnings.warn(f'{fpath} already exists and overwrite = False.')
        else:
            ds.to_netcdf(fpath, mode='w', format='NETCDF4')

        return fpath
    

    # ---------- Data Index ----------

    def update_data_index(self):
        """
        Update data index after processing new CTD profile.

        Returns
        -------
        data_index : pd.DataFrame
            CTD data index.
        di_path : str
            File path to CTD data index.
        """
        root = f'Q:/Messdaten/Aphys_Hypothesis_data/{self.lake}/'
        years = [d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))]

        di = []
        for yr in years:
            root_yr = f'Q:/Messdaten/Aphys_Hypothesis_data/{self.lake}/{yr}/CTD/'

            dates = os.listdir(root_yr)
            for dt in dates:
                root_date = os.path.join(root_yr, dt)
                dp_L2 = os.path.join(self.DPATH.format(lake=self.lake, year=yr, date=dt), 'L2')

                md_files = glob(f'{root_date}/*_md.json')
                for md_file in md_files:
                    fn = os.path.basename(md_file).split('_md')[0]
                    fp_L2 = glob(f'{dp_L2}/*{fn}*')
                    if len(fp_L2) == 1:
                        with open(md_file, 'r') as f:
                            md = json.load(f)

                        di.append({
                            'lake': md['lake'],
                            'date': md['date'],
                            'time': md['time'],
                            'profile_loc': md['profile_loc'],
                            'xsc': md['xsc'],
                            'ysc': md['ysc'],
                            'sensor': md['sensor'],
                            'serial_id': md['serial_id'],
                            'fname': fn
                        })

        di_path = self.DIPATH.format(lake=self.lake)
        with open(di_path, 'w') as f:
            json.dump(di, f, indent=2)

        return pd.DataFrame(di).sort_values(by=['date', 'time'], ascending=True).reset_index(drop=True), di_path


    # ---------- Pipeline ----------

    def process(self, update=False):
        """
        Process raw (L0) CTD data.  Convert to xarray and write to .nc (L1).
        Run quality assurance and write to .nc (L2).

        Parameters
        ----------
        update : bool
            If True, update data index with newly processed profile.
        """
        ds, self.fpath_L0 = self.parse_L0()
        self.fpath_L1 = self.write_to_nc(ds, 'L1')
        ds_qa = self.quality_assurance(ds)
        self.fpath_L2 = self.write_to_nc(ds_qa, 'L2')
        if update:
            data_index, di_path = self.update_data_index()