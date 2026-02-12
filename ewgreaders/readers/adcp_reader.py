### Class for reading ADCP data

# imports
import dolfyn as dlfn
import xarray as xr
import numpy as np
import json
from glob import glob

from .mooring_reader import MooringReader


class ADCPReader(MooringReader):

    def __init__(self, serial_id, lake, location, year, date, datalakes=False):
        """
        Initialize ADCPReader object.

        Parameters
        ----------
        serial_id : str
            Serial number of ADCP.
        lake : str
            Lake where ADCP is deployed.
        location : str
            Location code within lake of ADCP deployment.
        year : str
            Year of ADCP retrieval.
        date : str
            Date (YYYYMMDD) of ADCP retrieval.
        datalakes : bool
            Toggle whether to read from Eawag drive or DataLakes.
        """
        super().__init__(lake, location, year, date, datalakes)
        self.serial_id = serial_id

        self.mab = self.get_mab()
        self.depth = self.get_depth()
        self.orientation = self.get_orientation()


    # ---------- Metadata ---------- #
    
    def get_mab(self):
        """
        Parse metadata file for meters above bottom of ADCP.

        Returns
        -------
        mab : float
            Meters above bottom for ADCP.
        """
        instruments = self.md['instruments']
        for i in instruments:
            if i['serial_id'] == self.serial_id and i['instrument'] == 'adcp':
                return i['mab']
            
        raise ValueError(f'{self.serial_id} meters above bottom not found.')
    

    def get_depth(self, ds=None):
        """
        Parse metadata file for depth of ADCP.
        Use ADCP data if ADCP measures pressure.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        
        Returns
        -------
        depth : float
            Depth below water surface of ADCP.
        """
        if ds and 'pressure' in ds.data_vars:
            return ds.depth.where(ds.depth != 0, drop=True).mean().item()

        return self.total_depth - self.mab
    

    def get_orientation(self, ds=None):
        """
        Parse metadata file for orientation of ADCP.
        Use data attribute if ADCP data provided.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        orientation : str
            Orientation {up, down} of ADCP.
        """
        if ds and 'orientation' in ds.attrs.keys():
            pass
            #return ds.attrs['orientation']          # orientation not always correct in data file
        
        instruments = self.md['instruments']
        for i in instruments:
            if i['serial_id'] == self.serial_id and i['instrument'] == 'adcp':
                if "up" in i['comments']:
                    return "up"
                elif "down" in i['comments']:
                    return "down"
        
        raise ValueError(f'{self.serial_id} orientation not found.')


    def range_to_depth(self, ds):
        """
        Convert ADCP range values to depths.  Set attribute if surface or bottom is reached.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with range values converted to depths.
        """

        if self.orientation == 'up':
            ds['range'] = np.round(self.depth - ds.range.values, 1)
        elif self.orientation == 'down':
            ds['range'] = np.round(self.depth + ds.range.values, 1)

        self.surfbot = (ds['range'].max().item() >= self.total_depth) | (ds['range'].min().item() <= 0)

        return ds


    # ---------- Reading ---------- #

    def locate_data_file(self, level):
        """
        Locate file with ADCP data.

        Parameters
        ----------
        level : str
            Level of data.

        Returns
        -------
        fpath : list
            Path(s) to data file(s).
        """
        if level == 'L0':
            fpath = glob(f'{self.dpath_L0}/*{self.serial_id}*.000')
        elif level == 'L1':
            raise NotImplementedError
        elif level == 'L2':
            raise NotImplementedError
        else:
            raise ValueError("Data level must be L0, L1, or L2.")
        
        if len(fpath) == 0:
            raise IndexError(f'Could not find data file for {self.serial_id}.')
        
        return fpath
    

    def load_from_L0(self, fpath):
        """
        Load raw (L0) ADCP data into xarray Dataset.

        Parameters
        ----------
        fpath : list
            Path(s) to data file(s).

        Returns
        -------
        ds : xr.Dataset
            Dataset of all data recorded by ADCP.
        """
        if len(fpath) == 1:
            ds = dlfn.read(fpath[0])
        else:
            raise NotImplementedError
        
        # set depth and orientation attributes
        self.depth = self.get_depth(ds)
        self.orientation = self.get_orientation(ds)

        return ds



    # ---------- Quality Assurance ---------- #

    def run_qa(self, ds):
        """
        Apply all quality assurance filters to ADCP data.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        ds_qa : xr.Dataset
            Quality-assured ADCP data.
        """
        ds = self.qa_interface_surface(ds)
        ds = self.qa_interface_bottom(ds)
        ds = self.qa_min_corr(ds)
        ds = self.qa_pg14(ds)
        ds = self.qa_pg3(ds)
        ds = self.qa_vel_error(ds)
        ds = self.qa_tilt(ds)
        ds = self.qa_echo_amp_diff(ds)
        #ds = self.qa_corr_stdev(ds)

        return ds



    def qa_interface_surface(self, ds):
        """
        Filter data impacted by lake surface.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with surface interface filtered.
        """
        dist_sidelobe = self.depth * (1 - np.cos(ds.attrs['beam_angle'] * np.pi / 180))

        return ds.where(ds.range >= dist_sidelobe, drop=True)
    
    
    def qa_interface_bottom(self, ds):
        """
        Filter data impacted by lake bottom.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with bottom interface filtered.
        """
        dist_sidelobe = (self.total_depth - self.depth) * (1 - np.cos(ds.attrs['beam_angle'] * np.pi / 180))
        
        return ds.where(ds.range <= self.total_depth - dist_sidelobe, drop=True)
    

    def qa_min_corr(self, ds, corr_thresh=64):
        """
        Filter data with at least one beam correlation < 64.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        corr_thresh : int
            Threshold for beam correlations.

        Returns
        -------
        ds : xr.Dataset
            ADCP with minimum correlations filtered.
        """
        corr1 = ds.corr.sel(beam=1) >= corr_thresh
        corr2 = ds.corr.sel(beam=2) >= corr_thresh
        corr3 = ds.corr.sel(beam=3) >= corr_thresh
        corr4 = ds.corr.sel(beam=4) >= corr_thresh

        return ds.where(corr1 & corr2 & corr3 & corr4)
    
    
    def qa_pg14(self, ds, pg14_thresh=25):
        """
        Filter data with PG1 + PG4 < 25%.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        pg14_thresh : int
            Threshold for minimum good data percentage.

        Returns
        -------
        ds : xr.Dataset
            ADCP with good data filtered.
        """
        pg14 = ds.prcnt_gd.sel(beam=1) + ds.prcnt_gd.sel(beam=4)

        return ds.where(pg14 >= pg14_thresh)
    

    def qa_pg3(self, ds, pg3_thresh=25):
        """
        Filter data with PG3 > 25 %.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        pg3_thresh : int
            Threshold for maximum bad data percentage.

        Returns
        -------
        ds : xr.Dataset
            ADCP with bad data filtered.
        """
        pg3 = ds.prcnt_gd.sel(beam=3)

        return ds.where(pg3 <= pg3_thresh)
    

    def qa_vel_error(self, ds, velerr_thresh=0.05):
        """
        Filter data with velocity error > 0.05 m/s.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        velerr_thresh : float
            Threshold for maximum velocity error.

        Returns
        -------
        ds : xr.Dataset
            ADCP with velocity error filtered.
        """
        velerr = abs(ds.vel.sel(dir='err'))

        return ds.where(velerr <= velerr_thresh)
    

    def qa_tilt(self, ds, tilt_thresh=15):
        """
        Filter data with pitch or roll angle > 15°.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        tilt_thresh : int
            Threshold for maximum tilt.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with pitch and roll angle filtered.
        """
        pitch = abs(ds['pitch']) <= tilt_thresh
        roll = abs(ds['roll']) <= tilt_thresh

        return ds.where(pitch & roll)
    

    def qa_echo_amp_diff(self, ds, surfbot_toggle=False, ead_thresh=30):
        """
        Filter data with at least one beam with vertical echo difference between 
        consecutive bins > 30.

        Only required if beam reaches surface or bottom.

        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        surfbot_toggle : bool
            Toggle whether ADCP range reaches lake surface or bottom.
        ead_thresh : int
            Threshold for minimum vertical echo difference.

        Returns
        -------
        ds : xr.Dataset
            ADCP data with vertical echo difference filtered.
        """
        if not self.surfbot and not surfbot_toggle:
            return ds
        
        echo_amp_diff = ds.amp.diff(dim='range')
        ead1 = echo_amp_diff.sel(beam=1) <= ead_thresh
        ead2 = echo_amp_diff.sel(beam=2) <= ead_thresh
        ead3 = echo_amp_diff.sel(beam=3) <= ead_thresh
        ead4 = echo_amp_diff.sel(beam=4) <= ead_thresh

        return ds.where(ead1 & ead2 & ead3 & ead4)


    def qa_corr_stdev(self, ds, stdev_thresh=0.01, scale=100):  # threshold may be too harsh (test not from manual)
        """
        Filter data with standard deviation of 4 beams' correlations > 0.01.
        
        Parameters
        ----------
        ds : xr.Dataset
            ADCP data.
        stdev_thresh : float
            Threshold for maximum standard deviation.
        scale : int
            Scale correlations to [0, 1]

        Returns
        -------
        ds : xr.Dataset
            ADCP data with correlation standard deviation filtered.
        """
        corr1 = ds.corr.sel(beam=1) / scale
        corr2 = ds.corr.sel(beam=2) / scale
        corr3 = ds.corr.sel(beam=3) / scale
        corr4 = ds.corr.sel(beam=4) / scale
        corr = xr.concat([corr1, corr2, corr3, corr4], dim='beam')
        corr_stdev = corr.std(dim='beam')

        return ds.where(corr_stdev <= stdev_thresh)