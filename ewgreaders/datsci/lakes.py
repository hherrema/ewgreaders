### Methods for analyzing lake data

# imports
import xarray as xr


def thorpe_displacement(depth, temp):
    """
    Calculate Thorpe displacement.
    
    Parameters
    ----------
    depth : array_like
        Vertical profile depth values.
    temp : array_like
        Vertical profile temperature values.
    """
    # ensure depth increases monotonically
    raise NotImplementedError