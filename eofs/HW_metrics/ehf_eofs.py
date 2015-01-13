'''Calculate the EOFs of EHF heatwave metrics.

EHF heat wave metrics have been calculated from AWAP Tmax and Tmin data using a
90th percentile with respect to the entire analysis period. Metrics are annual
heat wave characteristics. They are HWF (frequency), HWD (duration), HWA 
(amplitude), HWM (magnitude), HWN (number) and HWT (timing).
'''

def load_heat_waves(filename):
    '''Load heat wave metrics from a netcdf file.

    Arguments 
    filename -- name of file containing heat wave metrics.
    maskname -- name of file containing AWAP land-sea mask.

    Returns
    hwf -- frequency
    hwn -- number
    hwd -- duration
    hwa -- amplitude
    hwm -- magnitude
    hwt -- timing
    '''
    from netCDF4 import Dataset
    ncin = Dataset(filename, 'r')
    hwf = ncin.variables['HWF_EHF'][:]
    hwn = ncin.variables['HWN_EHF'][:]
    hwd = ncin.variables['HWD_EHF'][:]
    hwa = ncin.variables['HWA_EHF'][:]
    hwm = ncin.variables['HWM_EHF'][:]
    hwt = ncin.variables['HWT_EHF'][:]
    lat = ncin.variables['lat'][:]
    lon = ncin.variables['lon'][:]
    times = ncin.variables['Times'][:]
    return hwf, hwn, hwd, hwa, hwm, hwt, lat, lon, times

if __name__ == "__main__":
    from numpy import arange
    # Load the heat wave metrics.
    fname = '/home/nfs/z5032520/heatwaves-svn/indices/EHF_metrics_1911-2013.nc' 
    hwf, hwn, hwd, hwa, hwm, hwt, lat, lon, times = load_heat_waves(fname)
    # Years
    years = arange(1911,2013)

