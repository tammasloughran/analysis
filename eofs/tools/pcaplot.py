"""Plot EOFs, PCs and eigenvalue scree plots.

These functions are intended to plot EOFs, PCs and eigenvalues caluclated
from the eofs python package.
"""
import matplotlib.pyplot as plt


def eofscatter(eofs):
    """Creates an interractive 3D plot of the first 3 EOFs in state space.
    
    Arguments:
    eofs -- the eofs to be plotted.
    """
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(eofs[0,:,:], eofs[1,:,:], eofs[2,:,:], marker='.')
    ax.plot([-1,1],[0,0],[0,0])
    ax.plot([0,0],[-1,1],[0,0])
    ax.plot([0,0],[0,0],[-1,1])
    ax.set_xlabel('EOF1')
    ax.set_ylabel('EOF2')
    ax.set_zlabel('EOF3')
    plt.show()
    plt.close()

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
    plt.figure()
    neig = range(len(eigens))
    plt.errorbar(range(1,20), eigens[0:19], yerr=errors[0:19], fmt='.')
    plt.xlabel("PC")
    plt.ylabel("Eigenvalue")
    plt.title(name+" Scree Plot")
    rev = svnversion()
    plotfilename = name+"_scree_r"+rev+".eps"
    plt.savefig(plotfilename, format='eps')
    plt.close()


def plot_eofs(eofs, lon, lat, name, single_eof=None, head=''):
    """Plots the eof patterns and saves them in .eps format.

    Arguments
    eofs -- matrix of EOFs.
    lon -- array of longitudes from the netcdf file.
    lat -- array of latitudes from the netcdf rile.
    name -- output plot file name.
    single_eof -- plot the single_eofth EOF in eofs on its own. (Default None)
    rot -- sting indicating rotation.
    """
    from numpy import arange, meshgrid
    import matplotlib
    from mpl_toolkits.basemap import Basemap
    import math
    from netCDF4 import Dataset
    # Define the scale of the plot based on minima and maxima.
    maxval = eofs[0,:,:].max()
    absminval = abs(eofs[0,:,:].min())
    if absminval>maxval:
        maxval = absminval
    maxval = math.ceil(maxval)
    if maxval<5:
        #If the maximum is very small round steps to nearest decimal point.
        step = round((maxval)/4.,1)
    else:
        step = math.ceil((maxval)/4.)
    levs = arange(-(5.*step), (5.*step)+step, step)
    # Define the parallels and meridians
    parallels = arange(-40., -9., 10.)
    meridians = arange(120., 160., 10.,)
    # make the figure
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    plt.figure()
    if single_eof!=None:
        fig, axes = plt.subplots(nrows=1, ncols=1)
        count = single_eof
        subfig = ''
    else:
        fig, axes = plt.subplots(nrows=2, ncols=2)
        count = 0
        subfig = 'a) '
    string = "EOF "
    replace = ['b) ','c) ','d) ', '']
    for ax in axes.flat:
        try:
            dummy = eofs[count, :, :]
        except:
            break
        # Use an equidistant cylyndrical map projection.
        map_axes = Basemap(ax=ax, projection='cyl',
            llcrnrlat=-44, urcrnrlat=-10,
            llcrnrlon=112, urcrnrlon=156)
        # Create the basemap grid 
        x, y = map_axes(*meshgrid(lon, lat))
        # Plot
        if eofs.shape[1:] == (68,88):
            msknc = Dataset('/srv/ccrc/data35/z5032520/AWAP/mask/varmask.nc')
            grey = msknc.variables['mask'][:]
            map_axes.contourf(x, y, grey, (0,0.5,1), colors=('1','0.5'))
        contours = map_axes.contour(x, y, eofs[count, :, :].squeeze(),
                linewidths=0.4, colors='k',levels=levs)
        cs = map_axes.contourf(x, y, eofs[count, :, :].squeeze(),
                levels=levs, cmap=plt.cm.RdBu_r)
        map_axes.drawparallels(parallels, labels=[True,False,False,False])
        map_axes.drawmeridians(meridians, labels=[False,False,False,True])
        map_axes.drawcoastlines()
        # Labels
        ax.set_title(subfig + string + str(count+1))
        subfig = replace[count]
        count = count + 1
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.12, 0.05, 0.76])
    cb = plt.colorbar(cs, cax=cbar_ax, orientation='vertical')
    fig.text(0.5,0.975, head+' of '+name, horizontalalignment='center',
                   verticalalignment='top')
    # Save
    rev = svnversion()
    plotfilename = name+"_r"+rev+".eps"
    plt.savefig(plotfilename, format='eps')
    plt.close()

def plot_three_eofs(eofs, lon, lat, name, single_eof=None, head=''):
    from numpy import arange, meshgrid
    import matplotlib
    from mpl_toolkits.basemap import Basemap
    import math
    from netCDF4 import Dataset
    import matplotlib as mpl
    # Define the scale of the plot based on minima and maxima.
    maxval = eofs[0,:,:].max()
    absminval = abs(eofs[0,:,:].min())
    if absminval>maxval:
        maxval = absminval
    maxval = math.ceil(maxval)
    if maxval<5:
        #If the maximum is very small round steps to nearest decimal point.
        step = round((maxval)/4.,1)
    else:
        step = math.ceil((maxval)/4.)
    levs = arange(-(5.*step), (5.*step)+step, step)
    # Define the parallels and meridians
    parallels = arange(-40., -9., 10.)
    meridians = arange(120., 160., 10.)
    # make the figure
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    plt.figure()
    if single_eof!=None:
        fig, axes = plt.subplots(nrows=1, ncols=1)
        count = single_eof
        subfig = ''
    else:
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,4.5))
        count = 0
        subfig = 'a) '
    string = "EOF "
    replace = ['b) ','c) ','d) ', '']
    for ax in axes.flat:
        try:
            dummy = eofs[count, :, :]
        except:
            break
        # Use an equidistant cylyndrical map projection.
        map_axes = Basemap(ax=ax, projection='cyl',
            llcrnrlat=-44, urcrnrlat=-10,
            llcrnrlon=112, urcrnrlon=156)
        # Create the basemap grid 
        x, y = map_axes(*meshgrid(lon, lat))
        # Plot
        if eofs.shape[1:] == (68,88):
            msknc = Dataset('/srv/ccrc/data35/z5032520/AWAP/mask/varmask.nc')
            grey = msknc.variables['mask'][:]
            map_axes.contourf(x, y, grey, (0,0.5,1), colors=('1','0.5'))
        contours = map_axes.contour(x, y, eofs[count, :, :].squeeze(),
                linewidths=0.4, colors='k',levels=levs)
        cs = map_axes.contourf(x, y, eofs[count, :, :].squeeze(),
                levels=levs, cmap=plt.cm.RdBu_r)
        map_axes.drawparallels(parallels, labels=[True,False,False,False],fontsize=8)
        map_axes.drawmeridians(meridians, labels=[False,False,False,True],fontsize=8)
        map_axes.drawcoastlines()
        # Labels
        ax.set_title(string + str(count+1))
        subfig = replace[count]
        count = count + 1
    #fig.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.12, 0.05, 0.76])
    cbar_ax, kw = mpl.colorbar.make_axes([ax for ax in axes.flat], orientation='horizontal')
    cb = plt.colorbar(cs, cax=cbar_ax, orientation='horizontal')
    if name=='HWF': units='days'
    if name=='HWN': units='events'
    if name=='HWD': units='days'
    if name=='HWA': units='$^\circ C^2$'
    if name=='HWM': units='$^\circ C^2$'
    if name=='HWT': units=''
    cb.set_label(units)
    #fig.text(0.5,0.975, name+' '+head, horizontalalignment='center',
    #               verticalalignment='top')
    # Save
    #plt.tight_layout()
    rev = svnversion()
    plotfilename = name+'block'+"_r"+rev+".eps"
    plt.savefig(plotfilename, format='eps')
    plt.close()

def plot_lags(rhos, ps, name):
    """Plot the laged correlations betweena a pc and index.

    Arguments
    rhos -- correlations for each PC against each month. PCs on second axis.
    ps -- same as rhos but p-values.
    """
    from numpy import nan, arange
    plt.figure()
    mrange = arange(0,-25,-1)
    months = ["Jan","Dec","Nov","Oct","Sep","Aug","Jul","Jun","May","Apr","Mar","Feb"]
    months += months+["Jan"]
    for pc in range(0,4,1):
        plt.plot(mrange, rhos[:,pc], label="PC%s"%(pc+1))
        #Plot dots for significance. p<=0.05
        plt.plot(mrange[ps[:,pc]<=0.05], rhos[ps[:,pc]<=0.05,pc], 'ko')
    plt.xticks(mrange, months, size='small', rotation=25)
    plt.title(name+" Lag Correlations")
    plt.xlabel("Lag Month")
    plt.ylabel(r'$\rho$')
    plt.ylim((-0.8,0.8))
    plt.plot([-25,0],[0,0],'k')
    plt.legend(loc=2)
    rev = svnversion()
    plotfilename = name+"_lagrho_r"+rev+".eps"
    plt.savefig(plotfilename,format='eps')
    plt.close()

def plot_pcs(pcs, mode, time, name, yearmean=False, head=''):
    """Plot the 1st and 2nd principle component time series.

    Arguments
    pcs -- matrix of PCs.
    time -- array of time axis.
    """
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
    pc1.plot(ax=axes[0], style='k', title=head+"PC 1", lw=0.5)
    axr = axes[0].twinx()
    mode.plot(ax=axr, style='b--', lw=0.6)
    # The right axis needs to be realigned to match the left.
    miny, maxy = axes[0].get_ylim() 
    axr.set_ylim(miny, maxy)
    # Repeat for the second PC.
    pc2.plot(ax=axes[1], style='k', title=head+"PC 2", lw=0.5)
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
    plt.close()

def station_scatter(stations,lats,lons,string,reverse=True):
    """Plot a field of station EOF.

    Arguments
    stations -- the station data eof to plot.
    lats -- the latitudes of each station.
    lons -- the longitudes of each station.
    string -- the title of the plot.
    reverse -- reverse the signs.
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import pdb
    plt.figure()
    # Create map
    mapping = Basemap(projection='cyl',llcrnrlat=-44,
            urcrnrlat=-10,llcrnrlon=112, urcrnrlon=156)
    mapping.drawcoastlines()
    # Reverse sign if desired.
    if reverse:
        stations = stations*-1
    # Separate positive and negative values since 
    # I want different markers for each.
    plus, circ = stations.copy(), stations.copy()
    plus[plus<0] = 0
    circ[circ>0] = 0
    circ = abs(circ)
    # Plot data
    mapping.scatter(lons,lats,s=plus*(300/0.9),marker='+',c='r')
    mapping.scatter(lons,lats,s=circ*(300/0.9),marker='o',
            edgecolors='b',facecolors='none')
    # Plot key
    mapping.scatter([126,126,126],[-34,-35.75,-37.5],[300,200,100],
            marker='+',c='r')
    mapping.scatter([126,126,126],[-39.25,-41,-42.75],[100,200,300],
            marker='o',edgecolors='b',facecolors='none')
    plt.text(128,-34-0.38,'= 0.9')
    plt.text(128,-35.75-0.38,'= 0.6')
    plt.text(128,-37.5-0.38,'= 0.3')
    plt.text(128,-39.25-0.38,'= -0.3')
    plt.text(128,-41-0.38,'= -0.6')
    plt.text(128,-42.75-0.38,'= -0.9')
    plt.title(string)
    plt.savefig(string+'.eps',format='eps')
    plt.close()

def svnversion():
    """Return the current svn revision.
    """
    import subprocess
    p = subprocess.Popen('svnversion', shell=True, \
       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()
    return stdout[:-1]
