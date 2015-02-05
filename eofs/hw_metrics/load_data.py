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
    from numpy import ma, empty
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
    masknc = Dataset('../mask/varmask.nc','r')
    mask = masknc.variables['mask'][:]
    mask2 = empty(hwf.shape)
    for n in range(hwf.shape[0]):
        mask2[n, :, :] = mask
    mask1 = hwf.mask
    hwf = ma.array(hwf, mask=mask2)
    hwn = ma.array(hwn, mask=mask2)
    hwd = ma.array(hwd, mask=mask2)
    hwa = ma.array(hwa, mask=mask2)
    hwm = ma.array(hwm, mask=mask2)
    hwt = ma.array(hwt, mask=mask2)
    for itime in range(times.size):
        for ilon in range(lon.size):
            for ilat in range(lat.size):
                if not mask1[itime,ilat,ilon]:
                    if hwa.mask[itime,ilat,ilon]:
                        hwa.mask[itime,ilat,ilon] = False
                        hwa[itime,ilat,ilon] = 0
                    if hwm.mask[itime,ilat,ilon]:
                        hwm.mask[itime,ilat,ilon] = False
                        hwm[itime,ilat,ilon] = 0
    return hwf, hwn, hwd, hwa, hwm, hwt, lat, lon, times

def load_index(fname, standardize=True):
    """Load an index from a CSV file. 
    
    The nino3.4 data is normalised by its mean and standard deviation after
    it is loaded.

    Arguments
    fname -- name of csv file containing monthly nino3.4 data

    Returns
    series -- pandas time series of index.
    """
    from numpy import genfromtxt
    from pandas import date_range, Series
    data = genfromtxt(fname, delimiter=',')
    series = data.reshape(data.shape[0]*data.shape[1])
    dates = date_range('1869-12', \
            '2013-12', freq='M').shift(15, freq='D')
    series = Series(series, index=dates)
    if standardize == True:
        mu = series.mean()
        sigma = series.std()
        series = (series - mu)/sigma
    return series


def load_index2(fname, standardize=True):
    """Load an index with timestamps from a text file.

    Arguments
    fname -- name of whitespace delimited file.

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
