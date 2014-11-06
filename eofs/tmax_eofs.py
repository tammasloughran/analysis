'''Calculate the eofs of tmax from AWAP data.

This is a crappy test script.
'''


def deseason(data, time, n=7):
    '''Remove seasonality from a data matrix.

    deseason removes seasonality from a numpy matrix of data with a
    moving window that is +-n timesteps wide. At the moment it 
    converts slices of a numpy matrix into a pandas series and then
    performs the moving average. A numpy matrix of deseasonalized
    data is returned.

    Arguments
    data -- 3D masked numpy array containing data to be deseasoned.
    time -- time axis data from .nc file.
    n -- size of the window for the moving average.

    Returns
    data - the deseasoned numpy data matrix.
    '''
    from pandas import Series, rolling_mean
    from numpy import isnan
    # Get the starting and end dates of the time axis.
    dates = get_dates(time, frequency='D')
    # Main deseasoning loop.
    tdim, ydim, xdim = data.shape
    mask = data.mask
    for x in range(0, xdim, 1):
        for y in range(0, ydim, 1):
            if not data.mask[0, y, x]:
                # Select point x,y.
                a_series = data[:, y, x]
                a_series = Series(a_series, index=dates)
                # Find the rolling mean.
                base = a_series['1960-12':'1991-01']
                roll_mean = rolling_mean(base, window=2*n+1, center=True)
                # Select base period. '61 to '90 is the standard BoM period.
                base = roll_mean['1961':'1990']
                for month in range(1, 13, 1):
                    for day in range(1, 32, 1):
                        # Find average for each day of year over base period.
                        average = base[((base.index.month==month)&\
                                        (base.index.day==day))].mean()
                        if not isnan(average):
                            # If average is NaN, then it must be a day
                            # that doesn't exist. eg 30 Feb.
                            # Subtract the average from each respective day of
                            # year of the series.
                            a_series[((a_series.index.month==month)&\
                                      (a_series.index.day==day))] -= average
                # Return the series to the data matrix.
                data[:, y, x] = a_series
    return data


def do_rotation(pcs, eofs, space='State'):
    '''Prepare data and perform varimax rotation.
    
    do_rotation reshapes the EOFs for rotation, then applies the 
    rotation and reshapes back into the original form.
    Rotation can be done on the PCs in sample space or on the EOFs in
    state space. In state space the masked values in the EOFs are removed
    before rotation and returned afterwards.

    Arguments
    pcs -- matrix of PCs allong the second axis.
    eofs -- matrix of EOFs allong the first axis.
    space -- type of varimax rotaion. Either 'state' or 'sample'.

    Returns
    pcs -- rotated PCs.
    eofs -- rotated EOFs.
    '''
    from numpy import dot, ma, where, ones, isnan, NaN
    import rotate
    if space == 'sample':
        pcs, R = rotate.varimax(pcs)
        nmaps, ny, nx = eofs.shape
        ngridpoints = nx*ny
        eofs2d = eofs.reshape([nmaps, ngridpoints])
        rot_eofs = dot(R, eofs2d)
        eofs = rot_eofs.reshape([nmaps, ny, nx])
    elif space == 'state':
        nmaps, ny, nx = eofs.shape
        ngridpoints = nx * ny
        eofs2d = eofs.reshape([nmaps, ngridpoints])
        nonMissingIndex = where(isnan(eofs2d.data[0]) == False)[0]
        dataNoMissing = eofs2d.data[:, nonMissingIndex]
        rot_eofs_nomiss, R = rotate.varimax(dataNoMissing.T, normalize=True)
        rotated_eofs = ones([nmaps, ngridpoints]) * NaN
        rotated_eofs = rotated_eofs.astype(eofs2d.dtype)
        rotated_eofs[:, nonMissingIndex] = rot_eofs_nomiss.T
        rotated_eofs = rotated_eofs.reshape([nmaps, ny, nx])
        eofs = ma.masked_array(rotated_eofs, eofs.mask)
        pcs = dot(pcs, R)
    return pcs, eofs


def get_dates(time, frequency='M'):
    '''Create a date_range object from time axis array.

    Arguments
    time -- the array containing the time axis.
    frequency -- the required frequency of the date_range object.

    Returns
    dates -- a date_range object.
    '''
    from pandas import date_range
    string = str(int(time[0]))
    start_date = string[6:] + '/' + string[4:6] + '/' + string[0:4]
    string = str(time[-1])
    end_date = string[6:] + '/' + string[4:6] + '/' + string[0:4]
    dates = date_range(start_date, end_date, freq=frequency)
    return dates


def load_data(filename, maskname):
    '''Load the the tmax data and apply a mask.

    Arguments
    filename -- name of the file containing the data.
    maskname -- name of the file containing the land sea mask.

    Returns
    temp2 -- masked t_max data.
    lons -- longitudes.
    lats -- latitudes.
    nctime -- times.
    '''
    from netCDF4 import Dataset
    from numpy import arange, empty, ma
    # Load the data from file.
    ncin = Dataset(filename, 'r')
    temp = ncin.variables['tmax'][:]
    lons = ncin.variables['lon'][:]
    lats = ncin.variables['lat'][:]
    nctime = ncin.variables['time'][:]
    ncin.close()
    # Load the land sea mask.
    maskfile = Dataset(maskname, 'r')
    mask = maskfile.variables['LSM'][:]
    maskfile.close()
    # The mask needs to be inverted. In numpy 1=no data and 0 = data.
    # I know it's a bit silly but that's how it is.
    mask = abs(mask - 1)
    dims = temp.shape
    mask2 = empty(dims)
    for n in range(dims[0]):
        mask2[n, :, :] = mask
    mask = None
    del mask
    # Finally apply the mask.
    temp2 = ma.masked_array(temp, mask2)
    return temp2, lons, lats, nctime


def load_index(fname):
    '''Load nino 3.4 index. 
    
    The nino3.4 data is normalised by its standard deviation after
    it is loaded.

    Arguments
    fname -- name of csv file containing monthly nino3.4 data

    Returns
    ninonorm -- normalised nino3.4 data.
    '''
    from numpy import genfromtxt
    from pandas import date_range, Series
    ninodata = genfromtxt(fname, delimiter=',')
    ninodata = ninodata.reshape(ninodata.shape[0]*ninodata.shape[1])
    ninodates = date_range('1869-12', \
            '2013-12', freq='M').shift(15, freq='D')
    nino34 = Series(ninodata,index=ninodates)
    ninostd = nino34.std()
    ninonorm = nino34/ninostd
    return ninonorm


def plot_eigenvalues(eigens, errors):
    '''Make a scree plot of the first 20 eigenvalues.
    
    It also plots the error bars from the Nort test. 
    PC/EOF pairs whos error pars do not overlap are significant.

    Arguments
    eigens -- array of eigenvalues.
    errors -- the corresponding array of errors from a North test.
    '''
    import matplotlib.pyplot as plt
    plt.figure()
    neig = range(len(eigens))
    plt.errorbar(range(1,20), eigens[0:19], yerr=errors[0:19], fmt='.')
    plt.xlabel('PC')
    plt.ylabel('Eigenvalue')
    plt.title('Scree Plot')
    plt.savefig('scree.v0.6.eps')


def plot_eofs(eofs, lon, lat, name):
    '''Plots the eof patterns and saves them in .eps format.

    Arguments
    eofs -- matrix of EOFs.
    lon -- array of longitudes from the netcdf file.
    lat -- array of latitudes from the netcdf rile.
    name -- output plot file name.
    '''
    from numpy import arange, meshgrid
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    # Use an equidistant cylyndrical map projection.
    levs = arange(-1.8, 1.9, 0.2)
    parallels = arange(-40., -9., 10.)
    meridians = arange(120., 160., 10.,)
    string = 'EOF '
    plt.figure()
    fig, axes = plt.subplots(nrows=2, ncols=2)
    count = 0
    for ax in axes.flat:
        map_axes = Basemap(ax=ax, projection='cyl',
            llcrnrlat=-44, urcrnrlat=-10,
            llcrnrlon=112, urcrnrlon=156)
        x, y = map_axes(*meshgrid(lon, lat))
        cs = map_axes.contourf(x, y, eofs[count, :, :].squeeze(),
            levels=levs, cmap=plt.cm.RdBu_r)
        map_axes.drawparallels(parallels, labels=[True,False,False,False])
        map_axes.drawmeridians(meridians, labels=[False,False,False,True])
        map_axes.drawcoastlines()
        ax.set_title(string + str(count+1))
        count = count + 1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.12, 0.05, 0.76])
    cb = plt.colorbar(cs, cax=cbar_ax, orientation='vertical')
    plotfilename = name + 'v0.6.eps'
    plt.savefig(plotfilename, format='eps')


def plot_pcs(pcs, time):
    '''Plot the 1st and 2nd principle component time series.

    Arguments
    pcs -- matrix of PCs.
    time -- array of time axis.
    '''
    import matplotlib.pyplot as plt
    from pandas import date_range, Series
    dates = get_dates(time, frequency='M')
    pc1 = Series(pcs[:,0], index=dates)
    pc2 = Series(pcs[:,1], index=dates)
    plt.figure()
    fig, axes = plt.subplots(nrows=2, sharex=True, squeeze=True)
    pc1.plot(ax=axes[0], style='k', title='PC 1', lw=0.5)
    pc2.plot(ax=axes[1], style='k', title='PC 2', lw=0.5)
    fig.text(0.06, 0.5, 'Normalized Score', ha='center', 
            va='center', rotation='vertical')
    plt.xlabel('Date')
    plt.savefig('pcs_0.6.eps',format='eps')


def resample(data, time):
    '''Calculate the montly means allong the first axis of the data matrix.

    Arguments
    data -- the data matrix to resample into monthly means.
    time -- the original time axis.

    Returns
    resampledata -- the resampled data matrix as monthly means.
    '''
    from pandas import Series
    from numpy import array, ma, zeros
    dates = get_dates(time, frequency='D')
    tdim, ydim, xdim = data.shape
    nothing = Series(zeros(tdim), index=dates)
    nothing = nothing.resample('M', how='mean')
    newdim = len(nothing)
    del nothing
    resampledata = zeros([newdim, ydim, xdim])
    newmask = data.mask[0:newdim, :, :]
    for x in range(0, xdim, 1):
        for y in range(0, ydim, 1):
            a_series = data[:, y, x]
            a_series = Series(a_series, index=dates)
            resampledata[:, y, x] = array(a_series.resample('M', how='mean'))
    resampledata = ma.masked_array(resampledata, newmask)
    return resampledata


def theil_detrend(data):
    '''Detrend data with a Theil slope esimator.

    The Theil slope estimator uses up lots of memory so
    only use it if the dataset is small. Otherwise use 
    scipy.signal.detrend. The Theil slope estimator is also known
    as the Sen's Kendall or a Theil-Sen slope estimator.

    Arguments
    data -- the data to be detrended.
    '''
    tdim, ydim, xdim = data.shape
    for x in range(0, xdim, 1):
        for y in range(0, ydim, 1):
            if not data.mask[0, y, x]:
                a_series = data[:, y, x]
                medslope, medintercept = theilslopes(a_series, \
                        arange(len(a_series)))
                trend = medslope*(arange(len(a_series))) + medintercept
                a_series -= trend
                data[:, y, x] = a_series
    return data


if __name__ == "__main__":
    from eofs.standard import Eof
    from numpy import dot, load, arange, cos, sqrt, deg2rad, newaxis
    import numpy as np
    import sys
    from scipy.signal import detrend
    import scipy.stats as stats

    if sys.argv[1] == '--reload':
        # Load the t_max data.
        print 'Loading data'
        filename = 'AWAP_TX_1911-2011_0.5deg.nc'
        maskname = 'AWAP_Land-Sea-Mask_0.5deg.nc'
        t_max, lon, lat, time = load_data(filename, maskname)
        # Detrend
        t_max = detrend(t_max, axis=0, type='linear')
        # We are not interested in the seasonal cycle so this will be
        # removed with a +-7 (15 day) moving window average.
        print 'Deseasoning'
        t_max = deseason(t_max, time,n=7)
        # Resample the data to monthly means.
        t_max = resample(t_max)
        # Save the masked and deseasoned data to file.
        t_max.dump('masked_deseasoned_tmax_0.5d')
        lon.dump('lons_0.5d')
        lat.dump('lats_0.5d')
        time.dump('times_0.5d')
    elif sys.argv[1] == '--continue':
        print 'Loading masked and deseasoned data.'
        t_max = load('masked_deseasoned_tmax_0.5d')
        lon = load('lons_0.5d')
        lat = load('lats_0.5d')
        time = load('times_0.5d')
    else:
        print 'Specify --load or --continue'
        sys.exit()

    # Calculate weightings.
    coslat = cos(deg2rad(lat)).clip(0., 1.)
    wgts = sqrt(coslat)[..., newaxis]

    # Set up the EOF solver.
    print 'Setting up solver.'
    #solver = Eof(t_max)
    solver = Eof(t_max, weights=wgts)
    # PCs.
    print 'Calculating PCs'
    pcs = solver.pcs(pcscaling=1, npcs=4)
    # Explained variance.
    print 'Calculating explained variance'
    explained_variance = solver.varianceFraction()
    # North test errors.
    print 'Calculating North test errors'
    errors = solver.northTest(vfscaled=True)
    # Eigenvalues.
    print 'Calculating eigenvalues'
    eigens = solver.eigenvalues()
    # EOFs and EOFs expressed as covariance and correlation.
    print 'Calculating EOFs...'
    eofs = solver.eofs(eofscaling=2, neofs=4)
    print 'as variance...'
    eofs_covariance = solver.eofsAsCovariance(pcscaling=1, neofs=4)
    print 'as correlation'
    eofs_correlation = solver.eofsAsCorrelation(neofs=4)

    # Apply rotation to PCs and EOFs.
    print 'Rotating PCs and EOFs'
    pcs, eofs = do_rotation(pcs, eofs, space='sample')

    # load nino index and correlate to pc
    nino34 = load_index('NINO_3.4_monthly_index.csv')
    dmi = load_index('DMI_monthly_index.csv')
    ninoslice = nino34['1911-01':'2011-12']
    dmislice = dmi['1911-01':'2011-12']
    print 'Correlations are'
    ninopc1rho, p1 = stats.spearmanr(ninoslice, pcs[:,0])
    dmipc1rho, p2 = stats.spearmanr(dmislice,pcs[:,0])
    ninopc2rho, p3 = stats.spearmanr(ninoslice, pcs[:,1])
    dmipc2rho, p4 = stats.spearmanr(dmislice,pcs[:,1])
    print '     nino3.4         dmi'
    print 'PC1:', ninopc1rho, dmipc1rho
    print 'PC2:', ninopc2rho, dmipc2rho

    # Plotting.
    print 'Plotting'
    plot_eigenvalues(explained_variance, errors)
    plot_eofs(eofs, lon, lat, 'EOFs')
    plot_eofs(eofs_covariance, lon, lat, 'EOFs_Covariance')
    plot_eofs(eofs_correlation, lon, lat, 'EOFs_Correlation')
    plot_pcs(pcs, time)
    print 'Done!'
