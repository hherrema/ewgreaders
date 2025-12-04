### Utility methods from limnology data

# imports
import seawater as sw

def ch1903_to_latlng(xsc, ysc):
    """
    Convert CH1903 Swiss coordinates to latitude and longitude.
    
    Parameters
    ----------
    xsc : int
        x Swiss coordinate.
    ysc : int
        y Swiss coordinate.

    Returns
    -------
    lat : float
        Latitude [°].
    lon : float
        Longitude [°].
    """
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

