### Methods for analyzing lake data

# imports
import xarray as xr
import numpy as np
import scipy
import math


# ---------- Utility ---------- #
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


# ---------- CTD ---------- #


def locate_thermocline(ds, ra_window=1):
    """
    Locate thermocline, defined as maximum vertical temperature gradient.

    Parameters
    ----------
    ds : xr.Dataset 
        Single profile CTD data.
    ra_window : float
        Depth window for rolling average.

    Returns
    -------
    thermocline_depth : float
        Depth of thermocline.
    """
    # only keep samples with valid depths, quality checked temperature
    mask = ds['depth'].notnull() & (ds['Temp_qual'] == 0)
    depth = ds['depth'][mask].values
    temp = ds['Temp'][mask].values

    # rolling average to smooth temperature profile
    temp_ra = rolling_average(temp, depth, ra_window)

    # locate maximum negative temperature gradient
    dTdz = np.gradient(temp_ra, depth)
    thermocline_idx = np.nanargmin(dTdz)
    thermocline_depth = depth[thermocline_idx]

    return thermocline_depth


def locate_epi_meta_hypolimnion(ds):
    """
    Locate epilimnion, metalimnion, and hypolimnion regions, based on temperature profile.

    Parameters
    ----------
    ds : xr.Dataset
        Single profile CTD data.

    Returns
    -------
    epilimnion : tuple
        Start and end depths of epilimnion.
    metalimnion : tuple
        Start and end depths of metalimnion.
    hypolimnion : tuple
        Start and end depths of hypolimnion.
    """
    raise NotImplementedError


def calculate_transect_min_dox(datasets):
    """
    Calculate minimum dissolved oxygen concentration across profiles from same transect.

    Parameters
    ----------
    datasets : list
        Datasets for profiles taken during same transect.

    Returns
    -------
    min_dox : float
        Minimum dissolved oxygen concentration during transect.
    min_dox_depth: float
        Depth of minimum dissolved oxygen concentraton.
    max_transect_depth : float
        Maximum depth during transect.
    """
    min_dox, min_dox_depth, max_transect_depth = np.inf, np.inf, -np.inf
    for ds in datasets:
        # only keep samples with valid depths, quality checked dissolved oxygen
        mask = ds['depth'].notnull() & (ds['DO_mg_qual'] == 0)
        depth = ds['depth'][mask]
        dox_conc = ds['DO_mg'][mask]

        # calculate profile values
        min_dox_p = dox_conc.min().item()
        min_dox_idx = np.where(dox_conc == min_dox_p)[0]
        min_dox_depth_p = depth[min_dox_idx].min().item()
        max_depth_p = depth.max().item()

        # update transect minimum dissolved oxygen concentration
        if min_dox_p < min_dox:
            min_dox = min_dox_p
            min_dox_depth = min_dox_depth_p

        # update transect maximum depth
        if max_depth_p > max_transect_depth:
            max_transect_depth = max_depth_p

    return min_dox, min_dox_depth, max_transect_depth


def locate_anoxia(ds, epsilon=0.05, transect_min_dox=None):
    """
    Locate start of anoxic zone.

    Parameters
    ----------
    ds : xr.Dataset
        Single profile CTD data.
    epsilon : float
        Allowed deviation from minimum dissolved oxygen concentration.  SHOULD BE TAKEN AS MEASUREMENT ERROR OF CTD.
    transect_min_dox : float
        Minimum dissolved oxygen during transect.

    Returns
    -------
    anoxic_depth : float
        Shallowest depth of anoxic conditions.
    anoxic_dox : float
        Minimum dissolved oxygen concentration.
    anoxic_sem : float
        Standard error of dissolved oxygen measurements below anoxic depth.
    """
    # only keep samples with valid depths, quality checked dissolved oxygen
    mask = ds['depth'].notnull() & (ds['DO_mg_qual'] == 0)
    depth = ds['depth'][mask]
    dox_conc = ds['DO_mg'][mask]

    # compare with transect minimum
    if transect_min_dox:
        anoxic_idx = np.where(dox_conc <= transect_min_dox + epsilon)[0]
    # compare within profile
    else:        
        anoxic_idx = np.where(dox_conc <= min(dox_conc.min().item(), epsilon))[0]

    if len(anoxic_idx) == 0:
        raise IndexError('No anoxic depths located.')
    
    anoxic_depth = depth[anoxic_idx].min().item()
    anoxic_dox = dox_conc[anoxic_idx.min():]

    return anoxic_depth, np.min(anoxic_dox), scipy.stats.sem(anoxic_dox)


def thorpe_displacement(ds):
    """
    Calculate Thorpe displacement.
    
    Parameters
    ----------
    ds : xr.Dataset
        Single profile CTD data.
    """
    # ensure depth increases monotonically
    raise NotImplementedError


# ---------- ADCP ----------

def locate_interface_bounds(depth, mxsc, mysc, angle_rad, bathy):
    """
    Locate land boundaries of flux interface.

    TO DO: edge case where d < depth is not true land boundary (i.e., ridge).

    Parameters
    ----------
    depth : float
        Depth of ADCP bin.
    mxsc : int
        x Swiss coordinate of mooring location.
    mysc : int
        y Swiss coordinate of mooring location.
    angle_rad : float
        Angle of interface in radians.
    bathy : xr.Dataset
        Lake bathymetry.

    Returns
    -------
    xsc1 : float
        x Swiss coordinate of boundary 1.
    ysc1 : float
        y Swiss coordiate of boundary 1.
    xsc2 : float
        x Swiss coordinate of boundary 2.
    ysc2 : float
        y Swiss coordiate of boundary 2.
    distance : float
        Distance between boundary points.

    """
    # positive direction
    xsc1, ysc1 = mxsc, mysc
    while xsc1 < bathy.xsc.max().item() and ysc1 < bathy.ysc.max().item():
        d = bathy.sel(xsc=xsc1, ysc=ysc1, method='nearest').depth.item()

        # boundary reached
        if d < depth or math.isnan(d):
            break

        xsc1 = xsc1 + np.cos(angle_rad)
        ysc1 = ysc1 + np.sin(angle_rad)

    # negative direction
    xsc2, ysc2 = mxsc, mysc
    while xsc2 > bathy.xsc.min().item() and ysc2 > bathy.ysc.min().item():
        d = bathy.sel(xsc=xsc2, ysc=ysc2, method='nearest').depth.item()

        # boundary reached
        if d < depth or math.isnan(d):
            break

        xsc2 = xsc2 - np.cos(angle_rad)
        ysc2 = ysc2 - np.sin(angle_rad)

    # distance between boundary points
    distance = np.sqrt((xsc1 - xsc2)**2 + (ysc1 - ysc2)**2)

    return xsc1, ysc1, xsc2, ysc2, distance


def calculate_flux(ds, depth, mxsc, mysc, angle_deg, bathy):
    """
    Calculate flux from ADCP velocity data.

    Parameters
    ----------
    ds : xr.Dataset
        ADCP data.
    depth : float
        Depth of ADCP bin.
    mxsc : int
        x Swiss coordinate of mooring location.
    mysc : int
        y Swiss coordinate of mooring location.
    angle_deg : float
        Angle of interface in degrees.  0° indicates E-W interface.  90° indicates N-S interface.
    bathy : xr.Dataset
        Lake bathymetry.

    Returns
    -------
    flux : xr.DataArray
        Flux timeseries across interface [m^3/s].
    """
    if angle_deg < -90 or angle_deg >= 90:
        raise ValueError("Interface angle must be in range [-90°, 90°)")
    
    angle_rad = np.deg2rad(angle_deg)
    xsc1, ysc1, xsc2, ysc2, distance = locate_interface_bounds(depth, mxsc, mysc, angle_rad, bathy)
    h = ds.attrs['cell_size']
    da = distance * h               # area of cross-section rectangle [m^2]

    u = ds.vel.sel(dir='E').sel(range=depth, method='nearest')
    v = ds.vel.sel(dir='N').sel(range=depth, method='nearest')

    u_orth = u * (-1) * np.sin(angle_rad)
    v_orth = v * np.cos(angle_rad)

    # flux = velocity * cross-section area
    flux = (u_orth + v_orth) * da
    
    return flux.rename('flux')


def calculate_flow_angle(ds):
    """
    Calculate angle of flow from ADCP velocity data.
    Map {'E', 'N', 'W', 'S'} to angles {0°, 90°, 180°, 270°}.

    Parameters
    ----------
    ds : xr.DataArray
        ADCP data.

    Returns
    -------
    flow_angle : xr.DataArray
        Flow angle timeseries [°].
    """
    u = ds.vel.sel(dir='E')
    v = ds.vel.sel(dir='N')

    flow_angle = np.rad2deg(np.arctan2(v, u)) % 360

    return flow_angle.rename('angle')


def calculate_excursion_length(ds):
    """
    Calculate distance a water particle moves during a contiguous period of unidirectional motion.
    
    Reference: LawrEtal1997

    Parameters
    ----------
    ds : xr.Dataset
        ADCP data.

    Returns
    -------
    """
    raise NotImplementedError