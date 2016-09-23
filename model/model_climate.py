import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import pandas as pd
import os
import glob


# Define directory and ensembles
data_directory = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'
elnino_ens = os.listdir(data_directory+'elnino')
lanina_ens = os.listdir(data_directory+'lanina')
control = []
for year in xrange(1970,2000):
    for month in xrange(1,13):
        control += [data_directory+'vamrb/vamrba.pa'+str(year)+'-'+"%02d"%(month,)+'.nc']

# load climatology from control
controlnc = nc.MFDataset(control)
times = controlnc.variables['t'][:]
dates = nc.num2date(times,units=controlnc.variables['t'].units)
date_range = pd.period_range(dates[0],dates[-1], freq='M')
summer = (date_range.month==12)|(date_range.month==1)|(date_range.month==2)
# Pressure
pr = controlnc.variables['p'][:]
pr_summer = pr[summer]
pr_summer_clim = np.mean(pr_summer, axis=0)
# Temperature
temp = controlnc.variables['temp'][:]
temp_summer = temp[summer]
temp_summer_clim = np.mean(temp_summer, axis=0)
# Rain
rain = controlnc.variables['precip'][:]
rain_summer = rain[summer]
rain_summer_clim = np.mean(rain_summer, axis=0)
# wind
u = controlnc.variables['u'][:]
u_summer = u[summer]
u_summer_clim = np.mean(u_summer, axis=0)
v = controlnc.variables['u'][:]
v_summer = u[summer]
v_summer_clim = np.mean(v_summer, axis=0)

# loop over all ensembles
mslp_ensmean = np.zeros(pr_summer_clim.shape)
temp_ensmean = np.zeros(temp_summer_clim.shape)
rain_ensmean = np.zeros(rain_summer_clim.shape)
u_ensmean = np.zeros(u_summer_clim.shape)
v_ensmean = np.zeros(v_summer_clim.shape)
for ens in elnino_ens:
    ens_dir = data_directory+'elnino/'+ens+'/'
    summer_files = glob.glob(ens_dir+'*.pe2000-12.nc') # December
    summer_files += glob.glob(ens_dir+'*.pe2001-01.nc') # Janurary
    summer_files += glob.glob(ens_dir+'*.pe2001-02.nc') # February
    summernc = nc.MFDataset(summer_files)
    mslp = summernc.variables['p_1'][:]
    mslp_ensmean += np.mean(mslp, axis=0)-pr_summer_clim
    tave = summernc.variables['temp_2'][:]
    temp_ensmean += np.mean(tave, axis=0)-temp_summer_clim
    precip = summernc.variables['precip'][:]
    rain_ensmean += np.mean(precip, axis=0)-rain_summer_clim
    wu = summernc.variables['u'][:]
    u_ensmean += np.mean(wu, axis=0)#-u_summer_clim
    wv = summernc.variables['v'][:]
    v_ensmean += np.mean(wv, axis=0)#-v_summer_clim
mslp_ensmean = mslp_ensmean/30.
temp_ensmean = temp_ensmean/30.
rain_ensmean = rain_ensmean/30.
u_ensmean = u_ensmean.squeeze()/30.
v_ensmean = v_ensmean.squeeze()/30.

lats = summernc.variables['latitude'][:]
lats1 = summernc.variables['latitude_1'][:]
lons = summernc.variables['longitude'][:]


def plot_map(data,lats,lons,units,rng,fname):
    mp = Basemap(projection='robin', lon_0=180)
    x,y = np.meshgrid(lons,lats)
    xx,yy = mp(x,y)
    mp.pcolormesh(xx,yy,data.squeeze(),
                  vmin=-rng,vmax=rng, cmap='bwr_r')
    cb = mp.colorbar(location='bottom')
    cb.set_label(units)
    mp.drawcoastlines()
    plt.savefig(fname,format='eps')


def plot_stream(u,v,lats,lons,units,rng,fname):
    mp = Basemap(projection='cyl', lon_0=180.)
    velo = np.sqrt(np.square(u)+np.square(v))
    x,y = np.meshgrid(lons,lats)
    xx,yy = mp(x,y)
    mp.pcolormesh(xx,yy,velo,cmap='viridis',vmin=0,vmax=rng)
    mp.streamplot(xx,yy,u,v, color='k')
    cb = mp.colorbar(location='bottom',pad=0.3)
    mp.drawcoastlines()
    mp.drawparallels([-80,-40,-60,-20,0,20,40,60,80],labels=[1, 0, 0, 0],linewidth=0)
    mp.drawmeridians([0,40,80,120,160,200,240,280,320],labels=[0, 0, 0, 1],linewidth=0)
    cb.set_label(units)
    plt.show(fname,format='eps')

plot_map(mslp_ensmean,lats,lons,'Pa',500,'mslp_nino_ensmean.eps')
plot_map(temp_ensmean,lats,lons,'${}^{\circ}C$',5,'temp_nino_ensmean.eps')
plot_map(rain_ensmean*86400,lats,lons,'$mm/day$', 7,'rain_nino_ensmean.eps')
plot_stream(u_ensmean,v_ensmean,lats1,lons,'$ms^{-1}$',12,'winds_nino_ensmean.eps')