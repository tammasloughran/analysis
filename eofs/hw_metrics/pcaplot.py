"""Plot EOFs, PCs and eigenvalue scree plots.

These functions are intended to plot EOFs, PCs and eigenvalues caluclated
from the eofs python package.
"""

def get_dates(time, frequency='M'):
    """Create a date_range object from time axis array.

    Arguments
    time -- the array containing the time axis.
    frequency -- the required frequency of the date_range object.

    Returns
    dates -- a date_range object.
    """
    from pandas import date_range
    if frequency == 'M':
        string = str(int(time[0]))
        start_date = string[6:] + '/' + string[4:6] + '/' + string[0:4]
        string = str(time[-1])
        end_date = string[6:] + '/' + string[4:6] + '/' + string[0:4]
    elif frequency == 'A':
        start_date = str(time[0])
        end_date = str(time[-1]+1)
    dates = date_range(start_date, end_date, freq=frequency)
    return dates

def plot_eigenvalues(eigens, errors, name):
    """Make a scree plot of the first 20 eigenvalues.
    
    It also plots the error bars from the Nort test. 
    PC/EOF pairs whose error pars do not overlap are significant.

    Arguments
    eigens -- array of eigenvalues.
    errors -- the corresponding array of errors from a North test.
    """
    import matplotlib.pyplot as plt
    plt.figure()
    neig = range(len(eigens))
    plt.errorbar(range(1,20), eigens[0:19], yerr=errors[0:19], fmt='.')
    plt.xlabel("PC")
    plt.ylabel("Eigenvalue")
    plt.title("Scree Plot")
    rev = svnversion()
    plotfilename = name+"_scree_r"+rev+".eps"
    plt.savefig(plotfilename, format='eps')


def plot_eofs(eofs, lon, lat, name):
    """Plots the eof patterns and saves them in .eps format.

    Arguments
    eofs -- matrix of EOFs.
    lon -- array of longitudes from the netcdf file.
    lat -- array of latitudes from the netcdf rile.
    name -- output plot file name.
    """
    from numpy import arange, meshgrid
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    # Define the scale of the plot.
    maxval = eofs[0,:,:].max()
    absminval = abs(eofs[0,:,:].min())
    if absminval>maxval:
        maxval = absminval
    step = (maxval*2.)/20.
    levs = arange(-maxval, maxval+step, step)
    # Define the parallels and meridians
    parallels = arange(-40., -9., 10.)
    meridians = arange(120., 160., 10.,)
    # make the figure
    plt.figure()
    fig, axes = plt.subplots(nrows=2, ncols=2)
    string = "EOF "
    count = 0
    for ax in axes.flat:
        # Use an equidistant cylyndrical map projection.
        map_axes = Basemap(ax=ax, projection='cyl',
            llcrnrlat=-44, urcrnrlat=-10,
            llcrnrlon=112, urcrnrlon=156)
        # Create the basemap grid 
        x, y = map_axes(*meshgrid(lon, lat))
        # Plot
        cs = map_axes.contourf(x, y, eofs[count, :, :].squeeze(),
            levels=levs, cmap=plt.cm.RdBu_r)
        map_axes.drawparallels(parallels, labels=[True,False,False,False])
        map_axes.drawmeridians(meridians, labels=[False,False,False,True])
        map_axes.drawcoastlines()
        # Labels
        ax.set_title(string + str(count+1))
        count = count + 1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.12, 0.05, 0.76])
    cb = plt.colorbar(cs, cax=cbar_ax, orientation='vertical')
    # Save
    rev = svnversion()
    plotfilename = name+"_r"+rev+".eps"
    plt.savefig(plotfilename, format='eps')

def plot_lags(rhos, ps, name):
    """Plot the laged correlations betweena a pc and index.

    Arguments
    rhos -- correlations for each PC against each month. PCs on second axis.
    ps -- same as rhos but p-values.
    """
    import matplotlib.pyplot as plt
    from numpy import nan, arange
    plt.figure()
    mrange = arange(1,13,1)
    months = ["Jan","Feb","Mar","Apr","May","Jun",
              "Jul","Aug","Sep","Oct","Nov","Dec"]
    for pc in range(0,4,1):
        plt.plot(mrange, rhos[:,pc], label="PC%s"%(pc+1))
        #Plot dots for significance. p<=0.05
        plt.plot(mrange[ps[:,pc]<=0.05], rhos[ps[:,pc]<=0.05,pc], 'ko')
    plt.xticks(mrange, months, size='small')
    plt.title(name+" Lag Correlations")
    plt.xlabel("Month")
    plt.ylabel(r'$\rho$')
    plt.ylim((-0.8,0.8))
    plt.plot([0,12],[0,0],'k')
    plt.legend()
    rev = svnversion()
    plotfilename = name+"_lagrho_r"+rev+".eps"
    plt.savefig(plotfilename,format='eps')

def plot_pcs(pcs, mode, time, name, yearmean=False):
    """Plot the 1st and 2nd principle component time series.

    Arguments
    pcs -- matrix of PCs.
    time -- array of time axis.
    """
    import matplotlib.pyplot as plt
    from pandas import date_range, Series
    # Create pandas timestamps for PCs.
    dates = get_dates(time, frequency='A')
    pc1 = Series(pcs[:,0], index=dates)
    pc2 = Series(pcs[:,1], index=dates)
    # If PCs are not yearly data resample to yearly.
    if yearmean:
        pc1 = pc1.resample('A', how='mean')
        pc2 = pc2.resample('A', how='mean')
    # Make the figure.
    plt.figure()
    fig, axes = plt.subplots(nrows=2, sharex=True, squeeze=True)
    pc1.plot(ax=axes[0], style='k', title="PC 1", lw=0.5)
    axr = axes[0].twinx()
    mode.plot(ax=axr, style='b--', lw=0.6)
    # The right axis needs to be realigned to match the left.
    miny, maxy = axes[0].get_ylim() 
    axr.set_ylim(miny, maxy)
    # Repeat for the second PC.
    pc2.plot(ax=axes[1], style='k', title="PC 2", lw=0.5)
    axr = axes[1].twinx()
    mode.plot(ax=axr, style='b--', lw=0.6)
    miny, maxy = axes[1].get_ylim()
    axr.set_ylim(miny, maxy)
    # Add labels.
    fig.text(0.06, 0.5, "Normalized Score", ha='center', 
            va='center', rotation='vertical')
    plt.xlabel("Date")
    # Get the repository version and save in eps format.
    rev = svnversion()
    plotfilename = name+"_pcs_r"+rev+".eps"
    plt.savefig(plotfilename,format='eps')

def svnversion():
    """Return the current svn revision.
    """
    import subprocess
    p = subprocess.Popen('svnversion', shell=True, \
       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()
    return stdout[:-1]
