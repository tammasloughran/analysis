import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4 as nc
import pandas as pd
import glob


print "Loading control data"
# Define directory and ensembles
data_directory = '/srv/ccrc/data46/z5032520/modelout/ACCESS/'
elnino_ens = glob.glob('/srv/ccrc/data46/z5032520/modelout/ACCESS/elnino/*')
lanina_ens = glob.glob('/srv/ccrc/data46/z5032520/modelout/ACCESS/lanina/*')
modoki_ens = glob.glob('/srv/ccrc/data48/z5032520/modelout/ACCESS/modokielnino/va*')
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
pr_summer_clim = np.squeeze(np.mean(pr_summer, axis=0))
# Temperature
temp = controlnc.variables['temp'][:]
temp_summer = temp[summer]
temp_summer_clim = np.squeeze(np.mean(temp_summer, axis=0))
# Rain
rain = controlnc.variables['precip'][:]
rain_summer = rain[summer]
rain_summer_clim = np.squeeze(np.mean(rain_summer, axis=0))
# wind
u = controlnc.variables['u'][:]
u_summer = u[summer]
u_summer_clim = np.squeeze(np.mean(u_summer, axis=0))
v = controlnc.variables['u'][:]
v_summer = u[summer]
v_summer_clim = np.squeeze(np.mean(v_summer, axis=0))
# 500hPa geopotential
hgt500 = controlnc.variables['hgt_1'][:,5,...]

# Load latitude and longitude
lats = controlnc.variables['latitude'][:]
lats1 = controlnc.variables['latitude_1'][:]
lons = controlnc.variables['longitude'][:]


# Function to loop and mean over all ensembles
def mean_ensembles(experiment_ens):
    mslp_ensmean = np.zeros(pr_summer_clim.shape)
    temp_ensmean = np.zeros(temp_summer_clim.shape)
    rain_ensmean = np.zeros(rain_summer_clim.shape)
    u_ensmean = np.zeros(u_summer_clim.shape)
    v_ensmean = np.zeros(v_summer_clim.shape)
    for ens in experiment_ens:
        summer_files = glob.glob(ens+'/*.pe2000-12.nc') # December
        summer_files += glob.glob(ens+'/*.pe2001-01.nc') # Janurary
        summer_files += glob.glob(ens+'/*.pe2001-02.nc') # February
        summernc = nc.MFDataset(summer_files)
        mslp = summernc.variables['p_1'][:]
        mslp_ensmean += np.squeeze(np.mean(mslp, axis=0))-pr_summer_clim
        tave = summernc.variables['temp_2'][:]
        temp_ensmean += np.squeeze(np.mean(tave, axis=0))-temp_summer_clim
        precip = summernc.variables['precip'][:]
        rain_ensmean += np.squeeze(np.mean(precip, axis=0))-rain_summer_clim
        wu = summernc.variables['u'][:]
        u_ensmean += np.squeeze(np.mean(wu, axis=0))#-u_summer_clim
        wv = summernc.variables['v'][:]
        v_ensmean += np.squeeze(np.mean(wv, axis=0))#-v_summer_clim
        summernc.close()
    n = float(len(experiment_ens))
    mslp_ensmean = mslp_ensmean/n
    temp_ensmean = temp_ensmean/n
    rain_ensmean = rain_ensmean/n
    u_ensmean = u_ensmean.squeeze()/n
    v_ensmean = v_ensmean.squeeze()/n
    return mslp_ensmean, temp_ensmean, rain_ensmean, u_ensmean, v_ensmean


# Define plotting functions
def plot_map(data,lats,lons,units,rng,fname):
    plt.clf()
    mp = Basemap(projection='robin', lon_0=180)
    x,y = np.meshgrid(lons,lats)
    xx,yy = mp(x,y)
    if units!='$mm/day$':
        colours = 'bwr'
        plot_cont = True
    else: 
        colours = 'bwr_r'
        plot_cont = False
    clrmsh = mp.pcolormesh(xx,yy,data.squeeze(),vmin=-rng,vmax=rng, cmap=colours)
    if plot_cont:
        step = 2.*rng/10.
        levels = np.arange(-rng,rng+step,step)
        mp.contour(xx,yy,data.squeeze(),colors='k',levels=levels,linewidths=0.5)
    cb = mp.colorbar(clrmsh, location='bottom', pad=0.3)
    cb.set_label(units)
    mp.drawcoastlines()
    mp.drawparallels([-80,-40,-60,-20,0,20,40,60,80],labels=[1, 0, 0, 0],linewidth=0)
    mp.drawmeridians([0,80,160,240,320],labels=[0, 0, 0, 1],linewidth=0)
    plt.savefig(fname,format='eps')
    

def plot_stream(u,v,lats,lons,units,rng,fname):
    plt.clf()
    mp = Basemap(projection='cyl', lon_0=180.)
    velo = np.sqrt(np.square(u)+np.square(v))
    x,y = np.meshgrid(lons,lats)
    xx,yy = mp(x,y)
    mp.pcolormesh(xx,yy,velo,cmap='viridis',vmin=0,vmax=rng)
    mp.streamplot(xx,yy,u,v, color='k')
    cb = mp.colorbar(location='bottom',pad=0.3)
    mp.drawcoastlines()
    mp.drawparallels([-80,-40,-60,-20,0,20,40,60,80],labels=[1, 0, 0, 0],linewidth=0)
    mp.drawmeridians([0,80,160,240,320],labels=[0, 0, 0, 1],linewidth=0)
    cb.set_label(units)
    plt.savefig(fname,format='eps')


# Plot the meaning function and plot the results
print "Meaning El Nino"
mslp_ninomean, temp_ninomean, rain_ninomean, u_ninomean, v_ninomean = mean_ensembles(elnino_ens)
plot_map(mslp_ninomean,lats,lons,'Pa',500,'mslp_nino_ensmean.eps')
plot_map(temp_ninomean,lats,lons,'${}^{\circ}C$',3,'temp_nino_ensmean.eps')
plot_map(rain_ninomean*86400,lats,lons,'$mm/day$', 7,'rain_nino_ensmean.eps')
plot_stream(u_ninomean,v_ninomean,lats1,lons,'$ms^{-1}$',12,'winds_nino_ensmean.eps')

print "Meaning La Nina"
mslp_ninamean, temp_ninamean, rain_ninamean, u_ninamean, v_ninamean = mean_ensembles(lanina_ens)
plot_map(mslp_ninamean,lats,lons,'Pa',500,'mslp_nina_ensmean.eps')
plot_map(temp_ninamean,lats,lons,'${}^{\circ}C$',3,'temp_nina_ensmean.eps')
plot_map(rain_ninamean*86400,lats,lons,'$mm/day$', 7,'rain_nina_ensmean.eps')
plot_stream(u_ninamean,v_ninamean,lats1,lons,'$ms^{-1}$',12,'winds_nina_ensmean.eps')

print "Meaning Modoki"
mslp_modokimean, temp_modokimean, rain_modokimean, u_modokimean, v_modokimean = mean_ensembles(modoki_ens)
plot_map(mslp_modokimean,lats,lons,'Pa',500,'mslp_modoki_ensmean.eps')
plot_map(temp_modokimean,lats,lons,'${}^{\circ}C$',3,'temp_modoki_ensmean.eps')
plot_map(rain_modokimean*86400,lats,lons,'$mm/day$', 7,'rain_modoki_ensmean.eps')
plot_stream(u_modokimean,v_modokimean,lats1,lons,'$ms^{-1}$',12,'winds_modoki_ensmean.eps')