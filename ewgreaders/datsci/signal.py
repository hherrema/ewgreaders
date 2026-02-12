### Methods for signal processing of lake data

# imports
import xarray as xr
import numpy as np
import scipy


def rolling_average(data, depth, ra_window):
    """
    Compute rolling average.

    Parameters
    ----------
    data : array_like
        Data variable to compute rolling average of.
    depth : array_like
        Depth below water surface.
    ra_window : float
        Depth window for rolling average.

    Returns
    -------
    data_ra : np.array
        Rolling average of data.
    """
    data_ra = []
    for d in depth:
        mask = (depth >= d - (ra_window/2)) & (depth <= d + (ra_window/2))
        data_ra.append(np.mean(data[mask]))

    return np.array(data_ra)


def savitzky_golay(arr):
    """
    Apply Savitzky-Golay filter to smooth array.

    Parameters
    ----------
    arr : array_like
        Array of values to smooth.

    Returns
    -------
    arr_sg : array_like
        Smoothed array.
    """
    window = int(np.ceil(len(arr)/10) // 2 * 2 + 1)
    polyorder = min(3, window)
    
    return scipy.signal.savgol_filter(arr, window, polyorder, mode='nearest')


def order_profile(var, depth, descending):
    """
    Order profile.  Stable so repeat values maintain original order.

    Parameters
    ----------
    var : xr.DataArray
        Variable profile to sort.
    depth : xr.DataArray
        Depth below water surface, increasing.
    descending : bool
        True if var decreases with depth, False if var increases with depth.
    """
    if descending:
        s = -1
    else:
        s = 1

    var_ascending = var*s
    idx = var_ascending.argsort(kind='mergesort').values
    var_sorted = var_ascending.isel(time=idx)
    depth_sorted = depth.isel(time=idx)

    return var_sorted*s, depth_sorted


def valid_depths(ds, thresh):
    """
    Filter depths below threshold of non-nan values.

    Parameters
    ----------
    ds : xr.DataArray
        Data with a depth dimension.
    thresh : float
        Theshold for percentage of data with non-nan values.

    Returns
    -------
    ds : xr.DataArray
        Filtered data with only valid depths.
    """
    valid_depths = ds.notnull().mean(dim='time')

    # ADCP data has range coordinate
    if 'range' in ds.coords:
        ds = ds.sel(range=valid_depths >= thresh)

    # thermistor and oxygen logger data has depth coordinate
    else:
        ds = ds.sel(depth=valid_depths >= thresh)

    return ds