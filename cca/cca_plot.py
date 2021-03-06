import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from mpl_toolkits.basemap import Basemap

def get_levels(data):
    """Creates a sequence of centred levels for plotting anomaly contours.
    
    Inputs
    data -- data
    
    The range is calculated as the largest value from data symmetric 
    about zero. The steps depend on how large the maximum is.
    """
    maxval = np.nanmax(data)
    absminval = abs(np.nanmin(data))
    if absminval>maxval:
        maxval = absminval
    maxval = math.ceil(maxval)
    if maxval<5:
        #If the maximum is very small round steps to nearest decimal point.
        step = round((maxval)/5.,1)
    else:
        step = math.ceil((maxval)/5.)
    return np.arange(-(5.*step), (5.*step)+step, step)

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

def plot_cp(cp, lons, lats, name, mask=None):
    """Plot a single cannonical pattern.

    Inputs
    cp -- the cannonical pattern
    lons -- the longitudes of the field
    lats -- the latitudees of the field
    """
    fig = plt.figure()
    levels = get_levels(cp)
    map_axes = Basemap(projection='cyl',
            llcrnrlat=lats[-1], urcrnrlat=lats[0],
            llcrnrlon=lons[0], urcrnrlon=lons[-1])
    x, y = map_axes(*np.meshgrid(lons, lats))
    if cp.shape == (68,88):
        from netCDF4 import Dataset
        msknc = Dataset('/srv/ccrc/data35/z5032520/AWAP/mask/varmask.nc')
        grey = msknc.variables['mask'][:]
        #grey = np.ma.array(grey,mask=np.logical_not(grey))
        map_axes.contourf(x, y, 0.5*grey, vmin=0, vmax=1, cmap='Greys')
    map_axes.contour(x, y, cp, levels, colors='k', linewidths=0.4)
    cs = map_axes.contourf(x, y, cp, levels, cmap=plt.cm.RdBu_r)
    if mask!=None:
        plt.contourf(x, y, mask, 1, colors='none', 
            hatches=[None,'.'], extend='lower')
    map_axes.drawcoastlines()
    parallels = np.arange(lats[-1], lats[0], 10.)
    meridians = np.arange(lons[0], lons[-1], 20.)
    map_axes.drawparallels(parallels, 
            labels=[True,False,False,False], linewidth=0.0)
    map_axes.drawmeridians(meridians, 
            labels=[False,False,False,True], linewidth=0.0)
    plt.title(name+' Canonical Pattern')
    cb = plt.colorbar(cs, orientation='horizontal')
    filename = 'cp_'+name+'_'+svnversion()+'.eps'
    plt.savefig(filename ,format='eps')
    plt.close()

def plot_coefs(left, right, years, llabel='HWF', rlabel='SST'):
    """Plot the left and right canonical coefficients.

    Inputs
    left -- the left coefficients
    right -- the right coefficients
    years -- list of years for x axis
    """
    plt.figure()
    dates = pd.date_range(str(years[0]), str(years[-1]), freq='A')
    left = pd.Series(left, index=dates)
    right = pd.Series(right, index=dates)
    fig = left.plot()
    fig = right.plot()
    lines, labels = fig.get_legend_handles_labels()
    fig.legend(lines, [llabel,rlabel], loc='best')
    plt.xlabel('Year')
    plt.ylabel('Score')
    plt.title('Temporal Expansion Coefficients')
    filename = 'expcoef_'+llabel+'-'+rlabel+'_'+svnversion()+'.eps'
    plt.savefig(filename ,format='eps')
    plt.close()

def svnversion():
    """Return the current svn revision.
    """
    import subprocess
    p = subprocess.Popen('svnversion', shell=True, \
       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = p.communicate()
    return stdout[:-1]
