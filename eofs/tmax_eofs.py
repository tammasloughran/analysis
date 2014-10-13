"""
Calculate the eofs of tmax from AWAP data.
This is a crappy test script.
"""

def deseason(data, time, n=7):
    """
    deseason(data, time, n) deseasonalizes a numpy matrix of data with a
    moving window that is +-n timesteps wide. At the moment it 
    converts slices of a numpy matrix into a pandas series and then
    performs the moving average. A numpy matrix of deseasonalized
    data is returned. This is probably not the fastest way to do 
    this but, whatever man...

    Arguments
    data - 3D masked numpy array containing data to be deseasoned.
    time - time axis data from .nc file.
    n - the size of the window for the moving average. Defaults to 7 (i.e. 15 days).

    Returns
    data - the deseasoned data matrix.
    """
    from pandas import date_range, Series, rolling_mean
    from numpy import array, isnan
    # Get the starting and end dates of the time axis.
    string = str(int(time[0]))
    start_date = string[6:]+"/"+string[4:6]+"/"+string[0:4]
    string = str(time[-1])
    end_date = string[6:]+"/"+string[4:6]+"/"+string[0:4]
    # Define a range of dates for the data.
    dates = date_range(start_date, end_date, freq='D')
    # Main deseasoning loop.
    tdim, ydim, xdim = data.shape
    for x in range(0,xdim,1):
        for y in range(0,ydim,1):
            if not data.mask[0,y,x]:
                # Find the rolling mean at x,y.
                a_series = Series(data[:,y,x], index=dates)
                base = a_series['1960-12':'1991-01']
                roll_mean = rolling_mean(base, window=2*n+1, center=True)
                # Select the base period. '61 to '90 is the standard BoM base period.
                base = roll_mean['1961':'1990']
                for month in range(1,13,1):
                    for day in range(1,32,1):
                        # Find the average for each day of the year over the base period.
                        average = base[((base.index.month==month)&(base.index.day==day))].mean()
                        if not isnan(average):
                            # If average is NaN, then it must be a day that doesn't exist. eg 30 Feb.
                            # Subtract the average from each respective day of year of the series.
                            a_series[((a_series.index.month==month)&(a_series.index.day==day))] -= average
                # Return the series to the data matrix.
                data[:,y,x] = array(a_series)
    return data

def load_data(filename,maskname):
    """
    load_data(filename, maskname) loads t_max and dimension data from the 
    file in filename and applies a land sea mask.

    Arguments
    filename - name of the file containing the data.
    maskname - name of the file containing the land sea mask.

    Returns
    temp2 - masked t_max data.
    lons - longitudes.
    lats - latitudes.
    nctime - times.
    """
    from netCDF4 import Dataset
    from numpy import arange, empty, ma
    # Load the data from file.
    ncin = Dataset(filename, 'r')
    temp = ncin.variables['tmax'][:]
    lons = ncin.variables['lon'][:]
    lats = ncin.variables['lat'][:]
    nctime = ncin.variables['time'][:]
    # The time axis in the nc file is all screwy so 
    # I'll just use a fake time axis for now.
    # nctime = arange(nctime.size)
    ncin.close()
    # Load the land sea mask.
    maskfile = Dataset(maskname, 'r')
    mask = maskfile.variables['LSM'][:]
    maskfile.close()
    # The mask needs to be inverted. In numpy 1=no data and 0 = data.
    # I know it's a bit silly but that's how it is.
    mask = abs(mask - 1)
    # I dont know of any convenient, memory eficient way of broadcasting 
    # the 2D mask onto 3D data. Using numpy operators to do this slows 
    # the computer to a halt. So I define a 3D mask, fill it in and 
    # release the old mask.
    dims = temp.shape
    mask2 = empty(dims)
    for n in range(dims[0]):
        mask2[n,:,:] = mask
    mask = None
    del mask
    # Finally apply the mask.
    temp2 = ma.masked_array(temp, mask2)
    temp = None
    del temp
    mask2 = None
    del mask2
    return temp2, lons, lats, nctime

def plot_pcs(pcs,time):
    """
    plot_pcs(pcs,time) plots the 1st and 2nd principle component time series.
    """
    import matplotlib.pyplot as plt
    from pandas import date_range, Series
    string = str(int(time[0]))
    start_date = string[6:]+"/"+string[4:6]+"/"+string[0:4]
    string = str(time[-1])
    end_date = string[6:]+"/"+string[4:6]+"/"+string[0:4]
    # Define a range of dates for the data.
    dates = date_range(start_date, end_date, freq='D')
    pc1 = Series(pcs[:,0], index=dates)
    pc1monthly = pc1.resample('M', how='mean')
    pc2 = Series(pcs[:,1], index=dates)
    pc2monthly = pc2.resample('M', how='mean')
    plt.figure()
    pc1monthly.plot(style='b')
    pc2monthly.plot(style='g')
    plt.title("PC1 (blue) & PC2 (green)")
    plt.xlabel("Date")
    plt.ylabel("Normalized Score")
    plt.savefig('pcs_0.5d.v0.3.eps', format='eps')

def plot_eofs(eofs,lon,lat):
    """
    plot_eofs(eofs,lon,lat) plots the eof patterns and saves them
    in .eps format.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from numpy import arange, meshgrid

    # Use an equidistant cylyndrical map projection.
    levs = arange(-4.,4.5,0.5)
    parallels = arange(-40., -9., 10.)
    meridians = arange(120., 160., 10.,)
    string = "EOF "

    plt.figure()
    fig, axes = plt.subplots(nrows=2, ncols=2)
    count = 0
    for ax in axes.flat:
        map_axes = Basemap(ax=ax, projection='cyl',
            llcrnrlat=-44, urcrnrlat=-10,
            llcrnrlon=112, urcrnrlon=156)
        x, y = map_axes(*meshgrid(lon, lat))
        cs = map_axes.contourf(x, y, eofs[count,:,:].squeeze(),
            levels=levs, cmap=plt.cm.RdBu_r)
        map_axes.drawparallels(parallels, labels=[True,False,False,False])
        map_axes.drawmeridians(meridians, labels=[False,False,False,True])
        map_axes.drawcoastlines()
        ax.set_title(string+str(count+1))
        count = count+1

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cb = plt.colorbar(cs, cax=cbar_ax, orientation="vertical")

    plt.savefig('eofs_0.5d.v0.3.eps', format='eps')

def plot_eigenvalues(eigens):
    '''
    plot_eigenvalues(eigens) makes a scree plot of the first 
    20 eigenvalues stored in the variable eigen.
    '''
    import matplotlib.pyplot as plt
    plt.figure()
    neig = range(len(eigens))
    plt.plot(range(1,20), eigens[0:19], "bs")
    plt.xlabel("PC")
    plt.ylabel('Explained variance')
    plt.title("Scree Plot")
    plt.savefig("scree_0.5d.v0.3.eps")

if __name__ == "__main__":
    from eofs.standard import Eof
    import sys
    from numpy import dot, load
    import rotate

    if sys.argv[1] == '--reload':
        # Load the t_max data.
        print 'Loading data'
        filename = 'AWAP_TX_1911-2011_0.5deg.nc'
        maskname = 'AWAP_Land-Sea-Mask_0.5deg.nc'
        t_max, lon, lat, time = load_data(filename, maskname)
        # We are not interested in the seasonal cycle so this will be
        # removed with a +-7 (15 day) moving window average.
        print 'Deseasoning'
        t_max = deseason(t_max,time,n=7)
        # Save the masked and deseasoned data to file.
        t_max.dump("masked_deseasoned_tmax_0.5d")
        lon.dump("lons_0.5d")
        lat.dump("lats_0.5d")
        time.dump("times_0.5d")
    elif sys.argv[1] == '--continue':
        print 'Loading masked and deseasoned data.'
        t_max = load('masked_deseasoned_tmax_0.5d')
        lon = load("lons_0.5d")
        lat = load("lats_0.5d")
        time = load("times_0.5d")
    else:
        print "Specify --load or --continue"
        sys.exit()

    # Set up the EOF solver.
    print 'Solving eofs, pcs and eigenvalues'
    solver = Eof(t_max)
    # Find the first four pc series scaled to unit variance.
    pcs = solver.pcs(npcs=4, pcscaling=1)
    explained_variance = solver.varianceFraction()
    eigenvalues = solver.eigenvalues()
    # Find the first four eofs.
    eofs = solver.eofs(neofs=4, eofscaling=1)
    eofs_covariance = solver.eofsAsCovariance(neofs=4)
    eofs_correlation = colver.eofsAsCorrelation(neofs=4)

    # Apply rotation
    if sys.argv[2] == '--rotate':
        pcs, R = rotate.varimax(pcs)
        nmaps, ny, nx = eofs.shape
        channels = nx*ny
        eofs2d = eofs.reshape([nmaps, channels])
        rot_eofs = dot(R, eofs2d)
        eofs = rot_eofs.reshape([nmaps,ny,nx])

    # Plotting.
    print 'Plotting'
    plot_eigenvalues(explained_variance)
    plot_eofs(eofs,lon,lat)
    plot_pcs(pcs,time)
    print 'Done!'
