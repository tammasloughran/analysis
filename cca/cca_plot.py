import numpy as np
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

def plot_cp(cp, lons, lats):
    """Plot a single cannonical pattern.

    Inputs
    cp -- the cannonical pattern
    lons -- the longitudes of the field
    lats -- the latitudees of the field
    """
    levels = get_levels(cp)
    map_axes = Basemap(projection='cyl',llcrnrlat=lats[-1],
            urcrnrlat=lats[0],llcrnrlon=lons[0], urcrnrlon=lons[-1])
    x, y = map_axes(*np.meshgrid(lons, lats))
    map_axes.contourf(x,y,cp,levels,cmap=plt.cm.RdBu_r)
    map_axes.drawcoastlines()
    cb = plt.colorbar(orientation='horizontal')
    plt.show()

def plot_coefs(left, right):
    """Plot the left and right canonical coefficients.

    Inputs
    left -- the left coeficients.
    right -- the right coefficients.
    """
    plt.plot(left)
    plt.plot(right)
    plt.xlim([0,len(left)])
    plt.show()
