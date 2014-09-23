"""
Calculate the eofs of tmax from AWAP data.
This is a crappy test script.
"""
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from eofs.standard import Eof

# Load the data from file.
filename = 'AWAP_tmax_1951-2009_1deg.nc'
ncin = Dataset(filename, 'r')
temp = ncin.variables['tmax'][:]
# No need to remove the mean from the data. 
# Apparently the EOF solver does this automatically.
#temp_mean = temp.mean(axis=0)
#temp_anom = temp - temp_mean
lon = ncin.variables['lon'][:]
lat = ncin.variables['lat'][:]
nctime = ncin.variables['time'][:]
# The time axis in the nc file is all screwy so 
# I'll just use a fake time axis for now.
time = np.arange(nctime.size)

# Set up the EOF solver.
solver = Eof(temp)
# Find the first eof.
eof1 = solver.eofsAsCovariance(neofs=1)
# Find the first pc series scaled to unit variance.
pc1 = solver.pcs(npcs=1, pcscaling=1)

# Plotting.
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
plt.savefig('eof1.eps', format='eps')

plt.figure()
plt.plot(time, pc1)
plt.savefig('pc1.eps', format='eps')
