### Fluid dynamics for lake data

# imports
import xarray as xr
import numpy as np


# ---------- Utility ---------- #


# ---------- CTD ---------- #
def calculate_brunt_vaisala_frequency(ds):
    """
    Calculate Brunt-Väisälä (buoyancy) frequency.
    N^2 = (g / ρ) * (dρ/dz)

    Parameters
    ----------
    ds : xr.Dataset
        CTD data.

    Returns
    -------
    N : xr.DataArray
    """
    raise NotImplementedError


# ---------- ADCP ---------- #

def calculate_TKE(ds):
    """
    Calculate turbulent kinetic energy (TKE).
    TKE = (u'^2 + v'^2 + w'^2) / 2

    Parameters
    ----------
    ds : xr.Dataset
        ADCP data.

    Returns
    -------
    tke : xr.DataArray
        Turbulent kinetic energy timeseries as function of depth.
    tke : 
    """
    u = ds.vel.sel(dir='E')
    v = ds.vel.sel(dir='N')
    w = ds.vel.sel(dir='U')

    # calculate velocity perturbations
    u_prime = u - u.mean(dim='time')
    v_prime = v - v.mean(dim='time')
    w_prime = w - w.mean(dim='time')

    tke = 0.5 * (u_prime**2 + v_prime**2 + w_prime**2)

    return tke.drop_vars('dir').rename('tke')


def calculate_froude_number(ds, densimetric=False):
    """
    Calculate Froude number.
    Fr = u / sqrt(g*L)

    Parameters
    ----------
    ds : xr.Dataset
        ADCP data.
    densimetric : bool, optional
        Toggle densimetric Froude number (reduced gravity).

    Returns
    -------
    Fr : xr.DataArray
        Froude number timeseries.
    """
    g = 9.81
    if densimetric:
        raise NotImplementedError
    
    L = (ds.range.max().item() - ds.range.min().item()) + ds.attrs['cell_size']   # half cell size at top and bottom

    # horizontal speed
    u = ds.vel.sel(dir='E')
    v = ds.vel.sel(dir='N')
    speed = xr.ufuncs.sqrt(u**2 + v**2)

    Fr = (speed / np.sqrt(g * L)).mean(dim='range')    # average over depth

    return Fr.rename('Fr')