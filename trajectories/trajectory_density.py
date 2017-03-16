# -*- coding: utf-8 -*-
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

# Load the trajectories 
trajdir = '/home/nfs/z5032520/traj3d/work/'
os.chdir(trajdir+'lanina')
trajfiles = glob.glob('?????_heatwave_trajectories.nc')
ntimes = 241
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

# Calculate potential temperature
Rcp = 0.286
theta = temp*(100000./levs)**Rcp
theta = theta - 273.15
theta[theta>60] = np.nan
theta[theta<0] = np.nan

# Identify the trajectories that start within a given region.
ibox = (lons[:,0]>=140.)&(lons[:,0]<=150.)&(lats[:,0]>=-30.)&(lats[:,0]<=-20.)

# Remove stationary trajectories that get stuck
delta = np.sqrt((lons[:,0]-lons[:,-1])**2 + (lats[:,0]-lats[:,-1])**2)
ibox[delta<10] = False

# Plot theta
fig = plt.figure()
axes = fig.gca()
plt.boxplot(theta[ibox,2::6])
axes.set_xticks(np.arange(4,41,4))
axes.set_xticklabels(np.arange(1,11,1))
plt.xlabel('Days before heatwave')
plt.ylabel(r'$\theta$ ($^{\circ}$C)')
axes.spines['right'].set_visible(False)
axes.yaxis.set_ticks_position('left')
axes.spines['top'].set_visible(False)
axes.xaxis.set_ticks_position('bottom')
plt.savefig('theta.eps', format='eps')

# Create an array for the parcel density
dlats = np.arange(-70.,1.,2.)
dlons = np.arange(90.,201.,2.)
density = np.ones((ntimes,len(dlats),len(dlons)))*np.nan
for t in xrange(ntimes):
    tlons = lons[ibox,t]
    tlats = lats[ibox,t]
    for y,dlat in enumerate(dlats):
        for x,dlon in enumerate(dlons):
            density[t,y,x] = ((tlons>=dlon)&(tlons<=dlon+5)&(tlats>=dlat)&(tlats<=dlat+5)).sum()
dens_levs = np.arange(0,70,5)

#xx,yy = np.meshgrid(dlons+2,dlats+2)
# Plot all the trajectory points
#for i in xrange(lons.shape[1]):
#    fig = plt.figure()
#    m = Basemap(projection='mill',
#                llcrnrlon=90.,llcrnrlat=-70.,
#                urcrnrlon=200.,urcrnrlat=0.,
#                resolution='l')
#    xxx,yyy=m(xx,yy)
#    cont = m.contourf(xxx,yyy,density[i],levels=dens_levs,cmap='gist_heat_r',extend='max')
#    levels = [-1, 1]
#    plt.contourf(xxx, yyy, density[i], levels=levels, colors='w')
#    m.colorbar(cont)
#    x, y = m(lons[ibox,i],lats[ibox,i])
#    m.scatter(x,y,marker='.')
#    m.drawcoastlines()
#    plt.savefig(str(i).zfill(3), dpi=900)
#    plt.close()

#call(['ffmpeg', '-i %03d.png output.mov'])