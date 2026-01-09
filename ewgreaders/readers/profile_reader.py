### Base class for reading profile data

# imports
import gsw as sw
import xarray as xr


class ProfileReader:
    BATHY_PATH = 'Q:/Messdaten/Aphys_Hypothesis_data/{lake}/bathymetry.nc'

    def __init__(self, lake, year, date, fpath, datalakes=False):
        """
        Initialize ProfileReader object.

        Parameters
        ----------
        lake : str
            Lake where profile occurs.
        year : str
            Year of profile. 
        date : str
            Date (YYYYMMDD) of profile.
        fpath : str
            File path to profile.
        datalakes : bool
            Toggle whether to read from Eawag drive or DataLakes.
        """
        self.lake = lake
        self.year = year
        self.date = date
        self.fpath = fpath
        self.datalakes = datalakes

    
    # ---------- Utility ----------
    def load_bathymetry(self):
        """
        Load lake bathymetry data.

        Returns
        -------
        bathy : xr.Dataset
            Lake bathymetry.
        """
        return xr.open_dataset(self.BATHY_PATH.format(lake=self.lake))



    def calculate_depth(pressure, air_pressure, rho, lat):
        """
        Calculate depth below water surface.

        Parameters
        ----------
        pressure : array_like
            Pressure [dbar].
        air_pressure : float
            Air pressure [dbar].
        rho : array_like
            Water density [kg/m^3]
        lat : float
            Latitude [°].
        
        Returns
        -------
        depth : array_like
            Depth below water surface [m].
        """
        pressure_adjusted = pressure - air_pressure
        depth = 1e4 * pressure_adjusted / (rho * sw.g(lat))

        return depth


    def calculate_salinity(temperature, conductivity, cond_coef=0.874e-3):
        """
        Calculate salinity from conductivity.
        
        Parameters
        ----------
        temperature : array_like
            Water temperature [°C]
        conductivity : array_like
            Conductivity [mS/cm].               # CHECK UNITS FOR EACH CTD
        cond_coef : float
            Coeffiction to calculate salinity from conductivity at 20 °C.
        """
        temp_correction = (1.8626 
                           - 0.052908 * temperature 
                           + 0.00093057 * temperature ** 2 
                           - 6.78e-6 * temperature ** 3)
        
        conductivity20 = temp_correction * conductivity * 1e3    # conductivity at 20 °C
        salinity =  cond_coef * conductivity20

        return salinity
        


    def calculate_density(temperature, salinity):
        """
        Calculate water density.

        Parameters
        ----------
        temperature : array_like
            Water temperature [°C].
        salinity : array_like
            Water salinity.
        """
        rho = 1e3 * (0.9998395 
                    + 6.7914e-5 * temperature 
                    - 9.0894e-6 * temperature ** 2 
                    + 1.0171e-7 * temperature ** 3 
                    - 1.2846e-9 * temperature ** 4 
                    + 1.1592e-11 * temperature ** 5 
                    - 5.0125e-14 * temperature ** 6 
                    + salinity * (8.181e-4 
                                - 3.85e-6 * temperature 
                                + 4.96e-8 * temperature ** 2))
        
        return rho
