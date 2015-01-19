'''Calculate the EOFs of EHF heatwave metrics.

EHF heat wave metrics have been calculated from AWAP Tmax and Tmin data using a
90th percentile with respect to the entire analysis period. Metrics are annual
heat wave characteristics. They are HWF (frequency), HWD (duration), HWA 
(amplitude), HWM (magnitude), HWN (number) and HWT (timing).
'''

def get_dates(time, frequency='M'):
    '''Create a date_range object from time axis array.

    Arguments
    time -- the array containing the time axis.
    frequency -- the required frequency of the date_range object.

    Returns
    dates -- a date_range object.
    '''
    from pandas import date_range
    #string = str(int(time[0]))
    #start_date = string[6:] + '/' + string[4:6] + '/' + string[0:4]
    #string = str(time[-1])
    #end_date = string[6:] + '/' + string[4:6] + '/' + string[0:4]
    start_date = time[0]
    end_date = time[-1]
    dates = date_range(start_date, end_date, freq=frequency)
    return dates

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
    levs = arange(-6, 6.5, 0.5)
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
    dates = get_dates(time, frequency='A')
    pc1 = Series(pcs[:,0], index=dates)
    #pc1 = pc1.resample('A', how='mean')
    pc2 = Series(pcs[:,1], index=dates)
    #pc2 = pc2.resample('A', how='mean')
    plt.figure()
    fig, axes = plt.subplots(nrows=2, sharex=True, squeeze=True)
    pc1.plot(ax=axes[0], style='k', title='PC 1', lw=0.5)
    pc2.plot(ax=axes[1], style='k', title='PC 2', lw=0.5)
    fig.text(0.06, 0.5, 'Normalized Score', ha='center', 
            va='center', rotation='vertical')
    plt.xlabel('Date')
    plt.savefig('pcs_0.6.eps',format='eps')

if __name__ == "__main__":
    from numpy import arange, cos, sqrt, deg2rad, newaxis, dot
    from eofs.standard import Eof
    import rotate
    start_year = 1911
    end_year = 2012
    # Load the heat wave metrics.
    fname = '/srv/ccrc/data35/z5032520/AWAP/yearly/ehfhw/CCRC_NARCliM_1911-20\
    14_EHFheatwaves_summer_AWAP0.5deg.nc' 
    hwf, hwn, hwd, hwa, hwm, hwt, lat, lon, times = load_heat_waves(fname)
    # Years
    years = arange(start_year,end_year+1)
    # Calculate weightings.
    coslat = cos(deg2rad(lat)).clip(0.,1.)
    wgts = sqrt(coslat)[..., newaxis]
    # Set up solver
    retain = 4
    solver = Eof(hwf, weights=wgts)
    pcs = solver.pcs(pcscaling=1, npcs=retain)
    explained_variance = solver.varianceFraction()
    errors = solver.northTest(vfscaled=True)
    eigens = solver.eigenvalues()
    eofs = solver.eofs(eofscaling=2, neofs=retain)
    eofs_covariance = solver.eofsAsCovariance(pcscaling=1, neofs=retain)
    eofs_correlation = solver.eofsAsCorrelation(neofs=retain)
    # Apply rotation to PCs and EOFs.
    eofs2 = eofs
    pcs, eofs = rotate.do_rotation(pcs, eofs, space='state')
    # Plotting.
    plot_eigenvalues(explained_variance, errors)
    plot_eofs(eofs, lon, lat, 'Rotated_EOFs')
    plot_eofs(eofs2, lon, lat, 'EOFs')
    plot_eofs(eofs_covariance, lon, lat, 'EOFs_Covariance')
    plot_eofs(eofs_correlation, lon, lat, 'EOFs_Correlation')
    plot_pcs(pcs, times)
