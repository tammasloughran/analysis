"""
Created on Wed Mar 15 14:12:39 2017

@author: Tammas Loughran
"""
import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from subprocess import call
import glob
import scipy.stats as stats
import matplotlib.collections as mcoll

def colorline(x, y, 
        z=None, cmap=plt.get_cmap('viridis'), norm=plt.Normalize(0.0, 1.0),
        linewidth=2, alpha=1.0):
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
    z = np.asarray(z)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)
    return lc

def plot_tracks(ax, lons, lats, levs):
    mp = Basemap(ax=ax,projection='mill',
                llcrnrlon=40.,llcrnrlat=-70.,
                urcrnrlon=200.,urcrnrlat=0.,
                resolution='l')
    lons[lons<40] = -9999
    for traj in xrange(lats.shape[0]):
        x,y = mp(lons[traj],lats[traj])
        lc = colorline(x, y, levs[traj]/100000., linewidth=1)
        ax.add_collection(lc)
    cbar = mp.colorbar(lc, location='right')
    cbar.set_label('hPa')
    cbar.ax.set_yticklabels([0,100,200,300,400,500,600,700,800,900,1000])
    cbar.ax.invert_yaxis()
    mp.drawcoastlines()


def make_ibox(ur_box,size, lons, lats):
    cnrs = (ur_box[0],ur_box[0]+size,ur_box[1],ur_box[1]-size)
    ibox = (lons[:,0]>=cnrs[0])&(lons[:,0]<=cnrs[1])&(lats[:,0]<=cnrs[2])&(lats[:,0]>=cnrs[3])
    # Remove stationary trajectories that get stuck
    delta = np.sqrt((lons[:,0]-lons[:,-1])**2 + (lats[:,0]-lats[:,-1])**2)
    ibox[delta<10] = False
    return ibox


cwd = os.getcwd()
phase = 'lanina'
ntimes = 241
# Load the trajectories
def get_traj(cphase):
    trajdir = '/srv/ccrc/data35/z5032520/traj3d/work/'
    os.chdir(trajdir+cphase)
    trajfiles = glob.glob('?????_heatwave_trajectories.nc')
    lons = np.ones((1,ntimes))*np.nan
    lats = np.ones((1,ntimes))*np.nan
    levs = np.ones((1,ntimes))*np.nan
    temp = np.ones((1,ntimes))*np.nan
    q = np.ones((1,ntimes))*np.nan
    for ifile in trajfiles:
        ncfile = nc.Dataset(ifile,'r')
        lons = np.concatenate((lons, np.squeeze(ncfile.variables['lon'][...])),axis=0)
        lats = np.concatenate((lats, np.squeeze(ncfile.variables['lat'][...])),axis=0)
        levs = np.concatenate((levs, np.squeeze(ncfile.variables['lev'][...])),axis=0)
        temp = np.concatenate((temp, np.squeeze(ncfile.variables['TEMP'][...])),axis=0)
        q = np.concatenate((q, np.squeeze(ncfile.variables['Q'][...])),axis=0)
    return lons, lats, levs, temp, q

lons, lats, levs, temp, q = get_traj('elnino')
alons, alats, alevs, atemp, aq = get_traj('lanina')

os.chdir(cwd)

size = 7.5
rnames = ['North','Northeast','Central-East','Southeast']
regions = [(129.,-12),(139.,-18),(145.5,-24),(141,-31)]
f, axes = plt.subplots(nrows=4, ncols=2,figsize=(5.5,7.75),sharex=True,sharey=True)
for i,ur_box in enumerate(regions):
    ibox = make_ibox(ur_box,size, lons, lats)
    aibox = make_ibox(ur_box,size, alons, alats)
    num = ibox.sum()
    anum = aibox.sum()
    lab = [['a','b'],['c','d'],['e','f'],['g','h']]
    plot_tracks(axes[i][0],lons[ibox],lats[ibox],levs[ibox])
    axes[i][0].set_title(lab[i][0]+')', loc='left')
    plot_tracks(axes[i][1],alons[aibox],alats[aibox],alevs[aibox])
    axes[i][1].set_title(lab[i][1]+')', loc='left')
#plt.savefig('tracks.eps',format='eps')
plt.show()
