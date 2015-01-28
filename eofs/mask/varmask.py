'''Create a mask from the variance of the first 40 years of AWAP data.
'''
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from numpy import arange, meshgrid
import matplotlib.pyplot as plt
from numpy import zeros
from scipy import stats
import numpy as np

# Load data.
hwdset = Dataset('CCRC_NARCliM_1911-2014_EHFheatwaves_summer_AWAP0.5deg.nc',
                 'r')
hwn = hwdset.variables['HWN_EHF'][:]
lon = hwdset.variables['lon'][:]
lat = hwdset.variables['lat'][:]

# Calculate variance of first 40 years.
hwn_var = hwn[:40,:,:].var(axis=0)

# Calculate 5 yearly variance trend.
var_series = zeros([hwn.shape[0]-4,hwn.shape[1],hwn.shape[2]])
for yearc in range(var_series.shape[0]):
    if yearc>2 &  yearc<hwn.shape[0]-2:
        var_chunk = hwn[yearc-2:yearc+2,:,:].var(axis=0)
        var_series[yearc,:,:] = var_chunk
trend = zeros(hwn_var.shape)
for ilon in range(lon.size):
    for ilat in range(lat.size):
        # over the entire period
        #reg = stats.linregress(arange(var_series.shape[0]), 
        #    var_series[:40,ilat,ilon])
        #or over the firest 40 years.
        reg = stats.linregress(arange(40), var_series[:40,ilat,ilon])
        trend[ilat,ilon] = reg[0]

# Make the mask
mask = (hwn_var<0.8) & (trend<0.015)
in1 = np.where(lat<-20.)
in2 = np.where(lat<-31.)
in3 = np.where(lon>120.)
in4 = np.where(lon>132.)
mask[:in1[0][0],:] = False
mask[in2[0][0]:,:] = False
mask[:,:in3[0][0]] = False
mask[:,in4[0][0]:] = False

# Plot
plt.figure()
map_axes = Basemap(projection='cyl',llcrnrlat=-44, urcrnrlat=-10,
        llcrnrlon=112, urcrnrlon=156)
x, y = map_axes(*meshgrid(lon, lat))
vplt = plt.contourf(x, y, hwn_var)
plt.colorbar(vplt, orientation='vertical')
map_axes.drawcoastlines()
plt.savefig('variance40.eps', format='eps')

# plot trend
plt.figure()
cs = plt.contourf(x, y, trend, levels=arange(-0.13,0.14,0.01), 
        cmap=plt.cm.RdBu_r)
map_axes.drawcoastlines()
plt.colorbar(cs, orientation='vertical')
plt.savefig('trend.eps', format='eps')

# Plot mask
plt.figure()
cs = plt.pcolor(x, y, mask)
map_axes.drawcoastlines()
plt.savefig('mask.eps', format='eps')

# Save mask to a netCDF file.
masknc = Dataset('varmask.nc', 'w')
masknc.createDimension('lat', lat.shape[0])
masknc.createDimension('lon', lon.shape[0])
latout = masknc.createVariable('lat', 'f8', ('lat'), fill_value=-99.99)
setattr(latout, 'long_name', 'Latitude')
setattr(latout, 'units', 'degrees_north')
setattr(latout, 'axis', 'Y')
setattr(latout, 'standard_name', 'latitude')
lonout = masknc.createVariable('lon', 'f8', ('lon'), fill_value=-99.99)
setattr(lonout, 'long_name', 'Longitude')
setattr(lonout, 'units', 'degrees_east')
setattr(latout, 'axis', 'X')
setattr(latout, 'standard_name', 'longitude')
maskout = masknc.createVariable('mask', 'i', ('lat','lon'), 
        fill_value=-99.99)
latout[:] = lat[:]
lonout[:] = lon[:]
maskout[:] = mask[:]
setattr(masknc, "author", "Tammas Loughran, CCRC, UNSW, Australia")
setattr(masknc, "contact", "t.loughran@student.unsw.edu.au")
setattr(masknc, "script", "varmask.py")
