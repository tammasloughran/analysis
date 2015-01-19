'''Plot EOFs, PCs and eigenvalue scree plots.

These functions are intended to plot EOFs, PCs and eigenvalues caluclated
from the eofs python package.
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
    if frequency == 'M':
        string = str(int(time[0]))
        start_date = string[6:] + '/' + string[4:6] + '/' + string[0:4]
        string = str(time[-1])
        end_date = string[6:] + '/' + string[4:6] + '/' + string[0:4]
    elif frequency == 'A'
        start_date = str(time[0])
        end_date = str(time[-1]+1)
    dates = date_range(start_date, end_date, freq=frequency)
    return dates

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


def plot_pcs(pcs, time, yearmean=False):
    '''Plot the 1st and 2nd principle component time series.

    Arguments
    pcs -- matrix of PCs.
    time -- array of time axis.
    '''
    import matplotlib.pyplot as plt
    from pandas import date_range, Series
    dates = get_dates(time, frequency='A')
    pc1 = Series(pcs[:,0], index=dates)
    pc2 = Series(pcs[:,1], index=dates)
    if yearmean:
        pc1 = pc1.resample('A', how='mean')
        pc2 = pc2.resample('A', how='mean')
    plt.figure()
    fig, axes = plt.subplots(nrows=2, sharex=True, squeeze=True)
    pc1.plot(ax=axes[0], style='k', title='PC 1', lw=0.5)
    pc2.plot(ax=axes[1], style='k', title='PC 2', lw=0.5)
    fig.text(0.06, 0.5, 'Normalized Score', ha='center', 
            va='center', rotation='vertical')
    plt.xlabel('Date')
    plt.savefig('pcs_0.6.eps',format='eps')

