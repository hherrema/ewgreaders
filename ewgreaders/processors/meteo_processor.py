### Class for processing meteo station data

# imports
import pandas as pd
import xarray as xr


class MeteoProcessor:
    COLS_MAP = {
    'reference_timestamp': 'time',
    'tre200s0': 'air_temp',
    'fkl010z0': 'wind_speed',
    'dkl010z0': 'wind_direction'
    }
    VAR_ATTRS = {
        'time': {'long_name': 'Coordinated Universal Time (UTC)'},
        'air_temp': {'units': '°C', 'long_name': 'Air Temperature (2 meters above ground)'},
        'wind_speed': {'units': 'm/s', 'long_name': 'Wind Speed (10 minute average)'},
        'wind_direction': {'units': '°', 'long_name': 'Wind Direction (10 minute average)'}
    }

    df = pd.read_csv(meteo_path, sep=';')
    df = df.dropna(axis='columns', how='all')
    df = df[cols_map.keys()].rename(columns=cols_map)
    df['time'] = pd.to_datetime(df['time'], format='%d.%m.%Y %H:%M')
    df = df.set_index('time')
    ds = xr.Dataset(df)

    for var, attrs in var_attrs.items():
    if var in ds:
        ds[var].attrs.update(attrs)

    ds = ds.assign_attrs({'station': 'Cham'})