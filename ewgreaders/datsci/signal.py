### Methods for signal processing of lake data

# imports
import xarray as xr
import numpy as np
import scipy


def rolling_average_z(da, ra_window):
    """
    Compute rolling average along depth dimension.

    Parameters
    ----------
    da : xr.DataArray
        Data to compute rolling average of.
    ra_window : float
        Depth window for rolling average [m].

    Returns
    -------
    da_ra : xr.DataArray
        Rolling average of data.
    """
    da_ra = [da.sel(depth=slice(d - (ra_window/2), d + (ra_window/2))).mean() for d in da.depth]

    return xr.DataArray(da_ra, dims=da.dims, coords=da.coords, name=da.name)


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


def order_profile(da, surfmax):
    """
    Order profile.  Stable so repeat values maintain original order.

    Parameters
    ----------
    da : xr.DataArray
        Profile to sort.
    surfmax : bool
        True if var is max at surface, False if var is max at bottom.
    """
    if surfmax:
        s = -1
    else:
        s = 1

    da_ascending = da*s
    idx = da_ascending.argsort(kind='mergesort').values
    da_sorted = da_ascending.isel(depth=idx)

    return da_sorted*s


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
    
    return ds.sel(depth=valid_depths >= thresh)