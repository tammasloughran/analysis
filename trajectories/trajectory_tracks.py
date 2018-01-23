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

def plot_tracks(ax, lons, lats, levs, lab=[1,0,0,1]):
    mp = Basemap(ax=ax,projection='mill',
                llcrnrlon=40.,llcrnrlat=-70.,
                urcrnrlon=200.,urcrnrlat=0.,
                resolution='l')
    mp.drawmeridians([60,90,120,150,180],labels=lab,fontsize=8,linewidth=0.03,dashes=[5,10000])
    mp.drawparallels([-10,-30,-50,-70],labels=lab,fontsize=8,linewidth=0.03,dashes=[5,10000])
    lons[lons<40] = -9999
    for traj in xrange(lats.shape[0]):
        x,y = mp(lons[traj],lats[traj])
        lc = colorline(x, y, levs[traj]/100., norm=None, linewidth=1, cmap='nipy_spectral')
        ax.add_collection(lc)
    #cbar = mp.colorbar(lc, location='right')
    #cbar.set_label('hPa')
    #cbar.ax.set_yticklabels([0,100,200,300,400,500,600,700,800,900,1000])
    #cbar.ax.set_yticklabels([1000,900,800,700,600,500])
    #cbar.ax.invert_yaxis()
    mp.drawcoastlines()
    return lc


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
levs[:,0] = 45000.
levs[:,1] = 100000.
alevs[:,0] = 45000.
alevs[:,1] = 100000.


os.chdir(cwd)

size = 7.5
rnames = ['North','Northeast','Central-East','Southeast']
regions = [(129.,-12),(139.,-18),(145.5,-24),(141,-31)]
f, axes = plt.subplots(nrows=4, ncols=2,figsize=(5.5,7.75))
#for i,ur_box in enumerate(regions):
#    ibox = make_ibox(ur_box,size, lons, lats)
#    aibox = make_ibox(ur_box,size, alons, alats)
#    num = ibox.sum()
#    anum = aibox.sum()
#    lab = [['a','b'],['c','d'],['e','f'],['g','h']]
#    plot_tracks(axes[i][0],lons[ibox],lats[ibox],levs[ibox])
#    axes[i][0].set_title(lab[i][0]+')', loc='left')
#    plot_tracks(axes[i][1],alons[aibox],alats[aibox],alevs[aibox])
#    axes[i][1].set_title(lab[i][1]+')', loc='left')
# a)
ibox = make_ibox(regions[0],size,lons,lats)
num = ibox.sum()
plot_tracks(axes[0][0],lons[ibox],lats[ibox],levs[ibox],lab=[1,0,0,0])
axes[0][0].set_title('a)', loc='left')
# b)
aibox = make_ibox(regions[0],size,alons,alats)
anum = aibox.sum()
plot_tracks(axes[0][1],alons[aibox],alats[aibox],alevs[aibox],lab=[0,0,0,0])
axes[0][1].set_title('b)', loc='left')
# c)
ibox = make_ibox(regions[1],size,lons,lats)
num = ibox.sum()
plot_tracks(axes[1][0],lons[ibox],lats[ibox],levs[ibox],lab=[1,0,0,0])
axes[1][0].set_title('c)', loc='left')
# d)
aibox = make_ibox(regions[1],size,alons,alats)
anum = aibox.sum()
plot_tracks(axes[1][1],alons[aibox],alats[aibox],alevs[aibox],lab=[0,0,0,0])
axes[1][1].set_title('d)', loc='left')
# e)
ibox = make_ibox(regions[2],size,lons,lats)
num = ibox.sum()
plot_tracks(axes[2][0],lons[ibox],lats[ibox],levs[ibox],lab=[1,0,0,0])
axes[2][0].set_title('e)', loc='left')
# f)
aibox = make_ibox(regions[2],size,alons,alats)
anum = aibox.sum()
plot_tracks(axes[2][1],alons[aibox],alats[aibox],alevs[aibox],lab=[0,0,0,0])
axes[2][1].set_title('f)', loc='left')
# g)
ibox = make_ibox(regions[3],size,lons,lats)
num = ibox.sum()
plot_tracks(axes[3][0],lons[ibox],lats[ibox],levs[ibox],lab=[1,0,0,1])
axes[3][0].set_title('g)', loc='left')
# h)
aibox = make_ibox(regions[3],size,alons,alats)
anum = aibox.sum()
lc = plot_tracks(axes[3][1],alons[aibox],alats[aibox],alevs[aibox],lab=[0,0,0,1])
axes[3][1].set_title('h)', loc='left')

axes[0][0].annotate('El Nino               La Nina',(0.19,0.92),textcoords='figure fraction',fontsize=20)

f.subplots_adjust(bottom=0.1)
cbar_ax = f.add_axes([0.12, 0.06, 0.78, 0.02])
cbar = plt.colorbar(lc,cax=cbar_ax, orientation='horizontal')
cbar.set_label('hPa')
cbar.ax.invert_xaxis()
plt.savefig('tracks.eps',format='eps')
#plt.show()

