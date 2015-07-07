import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import load_data
from netCDF4 import Dataset
from pyclimate.bpcca import BPCCA
from mpl_toolkits.basemap import Basemap
from scipy.signal import detrend

# Load variables
sstfile = ('/media/Jupiter/observations/HadISST/pacific_SST.nc')
hwfile = ('/media/Jupiter/reanalysis/AWAP/yearly/ehfhw/'
        'CCRC_NARCliM_1911-2014_EHFheatwaves_summer_AWAP0.5deg_detrended.nc')
hwf, hwn, hwd, hwa, hwm, hwt, hw_lats, hw_lons, times\
        = load_data.load_heat_waves(hwfile)
hadisst = Dataset(sstfile,'r')
sst = hadisst.variables['sst'][:]
sst_time = hadisst.variables['time'][:]
sst_lats = hadisst.variables['lat'][:]
sst_lons = hadisst.variables['lon'][:]

# Define dates
sst_dates = np.array([dt.datetime(1,1,1) + dt.timedelta(hours=hrs) 
        for hrs in sst_time])
sst_years = np.array([sst_dates[i].year 
        for i in np.arange(sst_dates.size)])
sst_months = np.array([sst_dates[i].month 
        for i in np.arange(sst_dates.size)])

# Slice the desired period from 1910-2013
period = np.where(sst_years>=1911)[0]
sst_dates = sst_dates[period]
sst_months = sst_months[period]
sst_period = sst[period,...]

# Select the desired months to average over
sel = (sst_months==1)|(sst_months==2)|\
        (sst_months==3)|(sst_months==11)|(sst_months==12)
sst_sel = sst_period[sel,...]

# Calculate annual summer means
sst_annual = np.ma.zeros((102,)+sst_sel.shape[-2:])
for year in range(102):
    sst_annual[year,...] = sst_sel[3+5*year:3+5+5*year,...].mean(axis=0)

# Calculate anomalies with respect to 1961-1990 base period
years = np.arange(1911,2013)
base = np.where((years>=1961)&(years<=1990))[0]
ave = sst_annual[base,...].mean(axis=0)
ssta = sst_annual - ave

# Detrend
ssta_detrended = np.ma.array(detrend(ssta, axis=0), mask=ssta.mask)
hwf_detrended = np.ma.array(detrend(hwf, axis=0), mask=hwf.mask)

# Perform CCA

CCA = BPCCA(ssta_detrended, hwf_detrended[:-2,:,:])
L = CCA.leftPatterns()
R = CCA.rightPatterns()
r = CCA.correlation()
a = CCA.rightExpCoeffs()
b = CCA.leftExpCoeffs()
print (r)

# Plot
# Left paterns
map_axes = Basemap(projection='cyl',
        llcrnrlat=sst_lats[-1], urcrnrlat=sst_lats[0],
        llcrnrlon=sst_lons[0], urcrnrlon=sst_lons[-1])
x, y = map_axes(*np.meshgrid(sst_lons, sst_lats))
map_axes.contourf(x,y,L[:,:,0],cmap=plt.cm.RdBu_r)
map_axes.drawcoastlines()
cb = plt.colorbar()
plt.show()
# Right paterns
map_axes = Basemap(projection='cyl',
        llcrnrlat=hw_lats[-1], urcrnrlat=hw_lats[0],
        llcrnrlon=hw_lons[0], urcrnrlon=hw_lons[-1])
x, y = map_axes(*np.meshgrid(hw_lons, hw_lats))
map_axes.contourf(x,y,R[:,:,0],cmap=plt.cm.RdBu_r)
map_axes.drawcoastlines()
cb = plt.colorbar()
plt.show()
# Coeffs
plt.plot(a[:,0])
plt.plot(b[:,0])
plt.show()
