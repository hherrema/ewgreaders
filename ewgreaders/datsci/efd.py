### Fluid dynamics for lake data

# imports
import xarray as xr
import numpy as np
from .signal import order_profile


# ---------- Utility ---------- #


# ---------- CTD ---------- #

def brunt_vaisala_frequency(ds, sort=True):
    """
    Calculate Brunt-Väisälä (buoyancy) frequency.
    N^2 = (g / ρ) * (dρ/dz)
    z represents depth, with z = 0 at the surface.

    Parameters
    ----------
    ds : xr.Dataset
        CTD data.
    sort : bool
        If True, sort profile to monotonic order.

    Returns
    -------
    N2 : xr.DataArray
    """
    g = 9.81

    # mask for valid depths, quality checked temperature and conductivity (salinity)
    mask = (ds['depth'].notnull()) & (ds['Temp_qual'] == 0) & (ds['Cond_qual'] == 0)
    depth = ds['depth'][mask]
    rho = ds['rho'][mask]

    # sort density profile, use original depth
    if sort:
        rho, _ = order_profile(rho, depth, False)

    drhodz = np.gradient(rho, depth)
    N2 = (g/rho) * drhodz

    # return N2.rename('N2').assign_coords(depth=depth)
    return xr.DataArray(N2.values, dims=["time"], coords={"time": depth.time, "depth": ("time", depth.values)}, name="N2")


# ---------- ADCP ---------- #

def TKE(ds):
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


def froude_number(ds, densimetric=False):
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


def vertical_shear(ds):
    """
    Calculate vertical shear.
    S2 = (du/dz^2) + (dv/dz)^2

    Parameters
    ----------
    ds : xr.Dataset
        ADCP data.
    
    Returns
    -------
    S2 : xr.DataArray
        Vertical shear timeseries as a function of depth.
    """
    u = ds.vel.sel(dir='E')
    v = ds.vel.sel(dir='N')

    # vertical derivatives of horizontal velocity
    dudz = u.differentiate('range')
    dvdz = v.differentiate('range')

    S2 = dudz**2 + dvdz**2

    return S2.rename('shear')


# ---------- CTD + ADCP ----------
def richardson_number(ds_adcp):
    """
    Calculate Richardson number.
    Ri = N^2 / S^2
    """
    raise NotImplementedError