"""
Calculate the eofs of tmax from AWAP data.
This is a crappy test script.
"""

def deseason(data,n=15):
    """
    deseason(data,n) deseasonalizes a numpy matrix of data with a
    moving window that is +-n timesteps wide. At the moment it 
    converts slices of a numpy matrix into a pandas series and then
    performs the moving average. A numpy matrix of deseasonalized
    data is returned. This is probably not the fastest way to do 
    this but, whatever man...
    """
    from pandas import Series, rolling_mean
    from numpy import array
    tdim, ydim, xdim = data.shape
    for x in range(0,xdim,1):
        for y in range(0,ydim,1):
            a_series = Series(data[:,y,x])
            roll_mean = rolling_mean(a_series, window=2*n+1).shift(-n)
            deseasoned = a_series - roll_mean
            np_deseasoned = array(deseasoned)
            data[:,y,x] = np_deseasoned
    return data

def load_data(filename):
    """
    load_data(filename) loads t_max and dimension data from the 
    file in filename.
    """
    from netCDF4 import Dataset
    from numpy import arange
    # Load the data from file.
    ncin = Dataset(filename, 'r')
    temp = ncin.variables['tmax'][:]
    # No need to remove the mean from the data. 
    # Apparently the EOF solver does this automatically.
    #temp_mean = temp.mean(axis=0)
    #temp_anom = temp - temp_mean
    lons = ncin.variables['lon'][:]
    lats = ncin.variables['lat'][:]
    nctime = ncin.variables['time'][:]
    # The time axis in the nc file is all screwy so 
    # I'll just use a fake time axis for now.
    nctime = arange(nctime.size)
    ncin.close()
    return temp, lons, lats, nctime

def plot_pcs(pcs,time):
    """
    plot_pcs(pcs,time) plots the principle component time series.
    """
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(time, pcs)
    plt.xlim([1,20200])
    plt.savefig('pc1_v0.2.eps', format='eps')

def plot_eofs(eofs,lon,lat):
    """
    plot_eofs(eofs,lon,lat) plots the eof patterns.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    # Use an equidistant cylyndrical map projection.
    m = Basemap(projection='cyl', 
                llcrnrlat=-44, urcrnrlat=-10, 
                llcrnrlon=112, urcrnrlon=156)
    x, y = m(*np.meshgrid(lon, lat))
    #levs = np.arange(-8.,8.,1.)
    m.contourf(x, y, eof1.squeeze(), cmap=plt.cm.RdBu_r)
    parallels = np.arange(-45., -10., 10.)
    m.drawparallels(parallels, labels=[True,False,False,False])
    meridians = np.arange(115., 156., 10.,)
    m.drawmeridians(meridians, labels=[False,False,False,True])
    m.drawcoastlines()
    plt.title("EOF 1")
    cb = plt.colorbar(orientation='horizontal')
    #plt.show()
    plt.savefig('eof1_v0.2.eps', format='eps')


if __name__ == "__main__":
    from eofs.standard import Eof

    # Load the t_max data.
    filename = 'AWAP_tmax_1951-2009_1deg.nc'
    t_max, lon, lat, time = load_data(filename)

    # We are not interested in the seasonal cycle so this will be
    # removed with a +-15 moving window average.
    t_max2 = deseason(t_max)
    t_max3 = t_max2[16:len(time)-16,:,:]

    # Set up the EOF solver.
    solver = Eof(t_max3)
    # Find the first eof.
    eof1 = solver.eofsAsCovariance(neofs=1)
    # Find the first pc series scaled to unit variance.
    pc1 = solver.pcs(npcs=1, pcscaling=1)

    # Plotting.
    plot_eofs(eof1,lon,lat)
    plot_pcs(pc1,time[16:len(time)-16])

