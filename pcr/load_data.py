"""Load_data contains functions to load datafiles for a heat wave PCA.
"""

def load_heat_waves(filename):
    """Load heat wave metrics from a netcdf file.

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
    """
    from netCDF4 import Dataset
    import numpy as np
    ncin = Dataset(filename, 'r')
    maskfile = ('/srv/ccrc/data35/z5032520/AWAP/mask/varmask.nc')
    mask1file = ('/srv/ccrc/data35/z5032520/AWAP/mask/AWAP_Land-Sea-Mask_0.5deg.nc')
    varmasknc = Dataset(maskfile,'r')
    varmask = varmasknc.variables['mask'][:]
    hwf = ncin.variables['HWF_EHF'][:]
    mask2 = np.empty(hwf.shape)
    for n in range(hwf.shape[0]):
        mask2[n, :, :] = varmask
    mask1nc = Dataset(mask1file,'r')
    mask = abs(mask1nc.variables['LSM'][:]-1)
    mask1 = np.empty(hwf.shape)
    for n in range(hwf.shape[0]):
        mask1[n, :, :] = mask
    hwf = np.ma.array(ncin.variables['HWF_EHF'][:], mask=mask1)
    hwn = np.ma.array(ncin.variables['HWN_EHF'][:], mask=mask1)
    hwf[hwn==0] = 0
    hwd = np.ma.array(ncin.variables['HWD_EHF'][:], mask=mask1)
    hwd[hwn==0] = 0
    hwa = np.ma.array(ncin.variables['HWA_EHF'][:], mask=mask1)
    hwa[hwn==0] = 0
    hwm = np.ma.array(ncin.variables['HWM_EHF'][:], mask=mask1)
    hwm[hwn==0] = 0
    hwt = np.ma.array(ncin.variables['HWT_EHF'][:], mask=mask1)
    hwt[hwn==0] = 0
    lat = ncin.variables['lat'][:]
    lon = ncin.variables['lon'][:]
    times = ncin.variables['time'][:]
    hwf = np.ma.array(hwf, mask=mask2)
    hwn = np.ma.array(hwn, mask=mask2)
    hwd = np.ma.array(hwd, mask=mask2)
    hwt = np.ma.array(hwt, mask=mask2)
    hwa = np.ma.array(hwa, mask=mask2)
    hwm = np.ma.array(hwm, mask=mask2)
    hwa = np.ma.array(hwa, mask=mask2)
    hwm = np.ma.array(hwm, mask=mask2)
    return hwf, hwn, hwd, hwa, hwm, hwt, lat, lon, times

def load_index(fname, standardize=False):
    """Load an index from a file in a single column. 
    
    Arguments
    fname -- name of csv file containing monthly nino3.4 data.
    standardize -- if True the series is standardized.

    Returns
    series -- pandas time series of index.
    """
    from numpy import genfromtxt, nan
    from pandas import date_range, Series
    data = genfromtxt(fname)
    years = data[:,0].astype(int)
    dates = date_range(str(years[0]), str(years[-1]+1), freq='M')
    series = Series(data[:,-1], index=dates)
    series = series.replace(-99.9, nan)
    series = series.replace(-999, nan)
    if standardize == True:
        mu = series.mean()
        sigma = series.std()
        series = (series - mu)/sigma
    return series

def load_index2(fname, standardize=False):
    """Load an index with timestamps from a text file with multiple
    columns by month.

    Arguments
    fname -- name of whitespace delimited file.
    standardize -- if Ture the series is standardized.

    Returns
    series -- pandas time series of index.
    """
    from numpy import genfromtxt, nan
    from pandas import date_range, Series
    data = genfromtxt(fname)
    years = data[:,0].astype(int)
    index = data[:,1:]
    dates = date_range(str(years[0]), str(years[-1]+1), freq='M')
    series = index.reshape(index.shape[0]*index.shape[1])
    series = Series(series, index=dates)
    series = series.replace(-99.9, nan)
    series = series.replace(-999, nan)
    if standardize == True:
        mu = series.mean()
        sigma = series.std()
        series = (series - mu)/sigma
    return series
