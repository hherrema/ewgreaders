### Methods for signal processing of lake data

# imports
import xarray as xr
import numpy as np
import scipy


# ---------- Time Series ---------- #

def rolling_average_t():
    raise NotImplementedError


# ---------- Depth Profiles ---------- #

def rolling_average_z(da, ra_window=1):
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


def binned_average_z(da, bin_size=1):
    """
    Compute binned average along depth dimension.

    Parameters
    ----------
    da : xr.DataArray
        Data to compute binned average of.
    bin_size : float
        Depth bin size [m].

    Returns
    -------
    da_ba : xr.DataArray
        Binned average of data.
    """
    bins = np.arange(0, da.depth.max() + bin_size, bin_size)

    da_ba = da.groupby_bins('depth', bins, labels=bins[:-1] + bin_size/2).mean()
    
    return da_ba.rename({'depth_bins': 'depth'})


def order_profile(da, surfmax, reassign_depth=False):
    """
    Order profile.  Stable so repeat values maintain original order.

    Parameters
    ----------
    da : xr.DataArray
        Profile to sort.
    surfmax : bool
        True if var is max at surface, False if var is max at bottom.
    reassign_depth: bool
        If True, reassign original depth dimension.

    Returns
    -------
    da_sorted : xr.DataArray
        Ordered profile.
    """
    if surfmax:
        s = -1
    else:
        s = 1

    # extract depth before sorting
    depth = da['depth']

    da_ascending = da * s
    idx = da_ascending.argsort(kind='mergesort').values
    da_sorted = da_ascending.isel(depth=idx) * s

    if reassign_depth:
        da_sorted = da_sorted.assign_coords(depth=depth)

    return da_sorted


# ---------- Depth-Binned Data ---------- #

def valid_depths(ds, thresh):
    """
    Filter depths below threshold of non-nan values.
    Currently implemented for ADCP data.

    Parameters
    ----------
    ds : xr.Dataset
        Data with a depth dimension.
    thresh : float
        Theshold for percentage of data with non-nan values.

    Returns
    -------
    ds : xr.Dataset
        Filtered data with only valid depths.
    """
    prop_valid = ds['vel'].notnull().mean(dim=['time' ,'dir'])

    return ds.sel(depth=prop_valid >= thresh)


# ---------- Misc ---------- #

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





