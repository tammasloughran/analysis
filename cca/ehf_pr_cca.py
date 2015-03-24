from pyclimate.bpcca import BPCCA
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from numpy import arange, meshgrid
from mpl_toolkits.basemap import Basemap
import load_data
import numpy

def timevsspace(data):
    ntimes = data.shape[0]
    npoints = data.shape[1]*data.shape[2]
    reshaped_data = data.reshape([ntimes,npoints])
    return reshaped_data

def removenans(data):
    index = numpy.where(numpy.isnan(data) == False)[0]
    no_nan_data = data[index]
    return no_nan_data, index

def restorenans(data, shape, index):
    with_nans = numpy.ones(shape)
    with_nans[index] = data
    return with_nans

def mask_to_zero(data):
    fill = numpy.zeros(data.shape)
    index = numpy.where(data.mask == False)
    fill[index] = data[index]
    return fill

prfile = ('/home/nfs/z5032520/heatwaves-svn/dataproc/summer_mslp_1911-2013.nc')
hwfile = ('/media/Jupiter/reanalysis/AWAP/yearly/ehfhw/CCRC_NARCliM_1911-2014_EHFheatwaves_summer_AWAP0.5deg_detrended.nc')
maskfile = ('/media/Jupiter/reanalysis/AWAP/mask/varmask.nc')
#masknc = Dataset(maskfile, 'r')
#hwnc = Dataset(hwfile,'r')
hwf, hwn, hwd, hwa, hwm, hwt, hw_lats, hw_lons, times\
        = load_data.load_heat_waves(hwfile)
hwf = mask_to_zero(hwf)
prnc = Dataset(prfile,'r')
#mask = masknc.variables['mask'][:]
pr = prnc.variables['mslp'][:]
pr_lats = prnc.variables['lat'][:]
pr_lons = prnc.variables['lon'][:]
#hw = hwnc.variables['HWF_EHF'][:]
#hw_lats = hwnc.variables['lat'][:]
#hw_lons = hwnc.variables['lon'][:]
CCA = BPCCA(pr, hwf[:-3,:,:], (4,4))
L = CCA.leftPatterns()
R = CCA.rightPatterns()
r = CCA.correlation()
a = CCA.rightExpCoeffs()
b = CCA.leftExpCoeffs()
print (r)
map_axes = Basemap(projection='cyl',llcrnrlat=pr_lats[-1], urcrnrlat=pr_lats[0],llcrnrlon=pr_lons[0], urcrnrlon=pr_lons[-1])
x, y = map_axes(*meshgrid(pr_lons, pr_lats))
map_axes.contourf(x,y,L[:,:,0],cmap=plt.cm.RdBu_r)
map_axes.drawcoastlines()
cb = plt.colorbar()
plt.show()
map_axes = Basemap(projection='cyl',llcrnrlat=hw_lats[-1], urcrnrlat=hw_lats[0],llcrnrlon=hw_lons[0], urcrnrlon=hw_lons[-1])
x, y = map_axes(*meshgrid(hw_lons, hw_lats))
map_axes.contourf(x,y,R[:,:,0],cmap=plt.cm.RdBu_r)
map_axes.drawcoastlines()
cb = plt.colorbar()
plt.show()
plt.plot(a[:,0])
plt.plot(b[:,0])
plt.show()
