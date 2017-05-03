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
pacnino_ens = glob.glob('/srv/ccrc/data48/z5032520/modelout/ACCESS/pac_nino/va*')
pacnina_ens = glob.glob('/srv/ccrc/data48/z5032520/modelout/ACCESS/pac_nina/va*')
indpiod_ens = glob.glob('/srv/ccrc/data48/z5032520/modelout/ACCESS/ind_piod/va*')
indniod_ens = glob.glob('/srv/ccrc/data48/z5032520/modelout/ACCESS/ind_niod/va*')
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
pr = pr[summer]
pr_summer_clim = np.squeeze(np.mean(pr, axis=0))
# Temperature
temp = controlnc.variables['temp'][:]
temp = temp[summer]
temp_summer_clim = np.squeeze(np.mean(temp, axis=0))
# Rain
rain = controlnc.variables['precip'][:]
rain = rain[summer]
rain_summer_clim = np.squeeze(np.mean(rain, axis=0))
# wind
u = controlnc.variables['u'][:]
u = u[summer]
u_summer_clim = np.squeeze(np.mean(u, axis=0))
v = controlnc.variables['u'][:]
v = v[summer]
v_summer_clim = np.squeeze(np.mean(v, axis=0))
# 500hPa geopotential
hgt500 = controlnc.variables['ht_1'][:,5,...]
hgt500 = hgt500[summer]
hgt500_summer_clim = np.squeeze(np.mean(hgt500, axis=0))

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
    hgt500_ensmean = np.zeros(hgt500_summer_clim.shape)
    for ens in experiment_ens:
        summer_files = glob.glob(ens+'/*.pe2000-12.nc') # December
        summer_files += glob.glob(ens+'/*.pe2001-01.nc') # Janurary
        summer_files += glob.glob(ens+'/*.pe2001-02.nc') # February
        summernc = nc.MFDataset(summer_files)
        mslp = summernc.variables['p_1'][:]
        mslp_ensmean += np.squeeze(np.mean(mslp, axis=0))-pr_summer_clim
#        tave = summernc.variables['temp_2'][:]
#        temp_ensmean += np.squeeze(np.mean(tave, axis=0))-temp_summer_clim
#        precip = summernc.variables['precip'][:]
#        rain_ensmean += np.squeeze(np.mean(precip, axis=0))-rain_summer_clim
#        wu = summernc.variables['u'][:]
#        u_ensmean += np.squeeze(np.mean(wu, axis=0))#-u_summer_clim
#        wv = summernc.variables['v'][:]
#        v_ensmean += np.squeeze(np.mean(wv, axis=0))#-v_summer_clim
#        hgt = summernc.variables['ht_1'][:,3,...]
#        hgt500_ensmean += np.squeeze(np.mean(hgt, axis=0))-hgt500_summer_clim
#        summernc.close()
    n = float(len(experiment_ens))
    mslp_ensmean = mslp_ensmean/n
#    temp_ensmean = temp_ensmean/n
#    rain_ensmean = rain_ensmean/n
#    u_ensmean = u_ensmean.squeeze()/n
#    v_ensmean = v_ensmean.squeeze()/n
#    hgt500_ensmean = hgt500_ensmean/n
    return mslp_ensmean, temp_ensmean, rain_ensmean, u_ensmean, v_ensmean, hgt500_ensmean


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
mslp_ninomean, temp_ninomean, rain_ninomean, u_ninomean, v_ninomean, hgt_ninomean = mean_ensembles(elnino_ens)
#plot_map(mslp_ninomean,lats,lons,'Pa',500,'mslp_nino_ensmean.eps')
#plot_map(temp_ninomean,lats,lons,'${}^{\circ}C$',3,'temp_nino_ensmean.eps')
#plot_map(rain_ninomean*86400,lats,lons,'$mm/day$', 7,'rain_nino_ensmean.eps')
#plot_stream(u_ninomean,v_ninomean,lats1,lons,'$ms^{-1}$',12,'winds_nino_ensmean.eps')
#plot_map(hgt_ninomean,lats1,lons,'m', 80,'hgt_nino_ensmean.eps')

#print "Meaning La Nina"
#mslp_ninamean, temp_ninamean, rain_ninamean, u_ninamean, v_ninamean, hgt_ninamean = mean_ensembles(lanina_ens)
#plot_map(mslp_ninamean,lats,lons,'Pa',500,'mslp_nina_ensmean.eps')
#plot_map(temp_ninamean,lats,lons,'${}^{\circ}C$',3,'temp_nina_ensmean.eps')
#plot_map(rain_ninamean*86400,lats,lons,'$mm/day$', 7,'rain_nina_ensmean.eps')
#plot_stream(u_ninamean,v_ninamean,lats1,lons,'$ms^{-1}$',12,'winds_nina_ensmean.eps')
#plot_map(hgt_ninamean,lats1,lons,'m', 80,'hgt_nina_ensmean.eps')

print "Meaning Modoki"
mslp_modokimean, temp_modokimean, rain_modokimean, u_modokimean, v_modokimean, hgt_modokimean = mean_ensembles(modoki_ens)
#plot_map(mslp_modokimean,lats,lons,'Pa',500,'mslp_modoki_ensmean.eps')
#plot_map(temp_modokimean,lats,lons,'${}^{\circ}C$',3,'temp_modoki_ensmean.eps')
#plot_map(rain_modokimean*86400,lats,lons,'$mm/day$', 7,'rain_modoki_ensmean.eps')
#plot_stream(u_modokimean,v_modokimean,lats1,lons,'$ms^{-1}$',12,'winds_modoki_ensmean.eps')
#plot_map(hgt_modokimean,lats1,lons,'m', 80,'hgt_modoki_ensmean.eps')

fig,(ax1,ax2) = plt.subplots(nrows=2, ncols=1,figsize=(6,8))
units='Pa'

mp = Basemap(ax=ax1,projection='robin', lon_0=180)
x,y = np.meshgrid(lons,lats)
xx,yy = mp(x,y)
colours = 'bwr'
#clrmsh = mp.pcolormesh(xx,yy,data.squeeze(),vmin=-rng,vmax=rng, cmap=colours)
clrmsh = mp.contourf(xx,yy,mslp_ninomean.squeeze(),levels=range(-300,300+1,60), cmap=colours,extend='both')
step = 2.*300/10.
levels = np.arange(-300,300+step,step)
cont = mp.contour(xx,yy,mslp_ninomean.squeeze(),colors='k',levels=levels,linewidths=0.5)
for c in cont.collections:
    if c.get_linestyle() == [(None, None)]:
       continue
    else:
        c.set_dashes([(0, (2.0, 2.0))])
eqtr = mslp_ninomean[(lats>=-5)&(lats<=5)]
lt = lats[(lats>=-5)&(lats<=5)]
eqmax = np.where(eqtr==np.max(eqtr))
eqmax = (lons[eqmax[1][0]],lt[eqmax[0][0]])
eqmin = np.where(eqtr==np.min(eqtr))
eqmin = (lons[eqmin[1][0]],lt[eqmin[0][0]])
mp.scatter(eqmax[0],eqmax[1],s=35,c='grey',marker='D',latlon=True)
mp.scatter(eqmin[0],eqmin[1],s=35,c='grey',marker='o',latlon=True)
#mp.scatter(eqmin,'k',latlon=True)
cb = mp.colorbar(clrmsh, location='bottom', pad=0.3)
cb.set_label(units)
mp.drawcoastlines()
mp.drawparallels([-80,-40,-60,-20,0,20,40,60,80],labels=[1, 0, 0, 0],linewidth=0)
mp.drawmeridians([0,80,160,240,320],labels=[0, 0, 0, 1],linewidth=0)
ax1.set_title('a) El Nino', loc='center')

mp = Basemap(ax=ax2,projection='robin', lon_0=180)
x,y = np.meshgrid(lons,lats)
xx,yy = mp(x,y)
colours = 'bwr'
#clrmsh = mp.pcolormesh(xx,yy,data.squeeze(),vmin=-rng,vmax=rng, cmap=colours)
clrmsh = mp.contourf(xx,yy,mslp_modokimean.squeeze(),levels=range(-300,300+1,60), cmap=colours,extend='both')
step = 2.*300/10.
levels = np.arange(-300,300+step,step)
cont = mp.contour(xx,yy,mslp_modokimean.squeeze(),colors='k',levels=levels,linewidths=0.5)
for c in cont.collections:
    if c.get_linestyle() == [(None, None)]:
       continue
    else:
        c.set_dashes([(0, (2.0, 2.0))])
eqtr = mslp_modokimean[(lats>=-5)&(lats<=5)]
lt = lats[(lats>=-5)&(lats<=5)]
eqmax = np.where(eqtr==np.max(eqtr))
eqmax = (lons[eqmax[1][0]],lt[eqmax[0][0]])
eqmin = np.where(eqtr==np.min(eqtr))
eqmin = (lons[eqmin[1][0]],lt[eqmin[0][0]])
mp.scatter(eqmax[0],eqmax[1],s=35,c='grey',marker='D',latlon=True)
mp.scatter(eqmin[0],eqmin[1],s=35,c='grey',marker='o',latlon=True)
#mp.scatter(eqmin,'k',latlon=True)
cb = mp.colorbar(clrmsh, location='bottom', pad=0.3)
cb.set_label(units)
mp.drawcoastlines()
mp.drawparallels([-80,-40,-60,-20,0,20,40,60,80],labels=[1, 0, 0, 0],linewidth=0)
mp.drawmeridians([0,80,160,240,320],labels=[0, 0, 0, 1],linewidth=0)
ax2.set_title('b) Modoki', loc='center')

plt.savefig('ACCESS_mslp_shift.eps',format='eps')

#print "Meaning Pacific El nino"
#mslp_ninomean, temp_ninomean, rain_ninomean, u_ninomean, v_ninomean, hgt_ninomean = mean_ensembles(pacnino_ens)
#plot_map(mslp_ninomean,lats,lons,'Pa',500,'mslp_pacnino_ensmean.eps')
#plot_map(temp_ninomean,lats,lons,'${}^{\circ}C$',3,'temp_pacnino_ensmean.eps')
#plot_map(rain_ninomean*86400,lats,lons,'$mm/day$', 7,'rain_pacnino_ensmean.eps')
#plot_stream(u_ninomean,v_ninomean,lats1,lons,'$ms^{-1}$',12,'winds_pacnino_ensmean.eps')
#plot_map(hgt_ninomean,lats1,lons,'m', 80,'hgt_pacnino_ensmean.eps')
#
#print "Meaning Pacific La nina"
#mslp_ninamean, temp_ninamean, rain_ninamean, u_ninamean, v_ninamean, hgt_ninamean = mean_ensembles(pacnina_ens)
#plot_map(mslp_ninamean,lats,lons,'Pa',500,'mslp_pacnina_ensmean.eps')
#plot_map(temp_ninamean,lats,lons,'${}^{\circ}C$',3,'temp_pacnina_ensmean.eps')
#plot_map(rain_ninamean*86400,lats,lons,'$mm/day$', 7,'rain_pacnina_ensmean.eps')
#plot_stream(u_ninamean,v_ninamean,lats1,lons,'$ms^{-1}$',12,'winds_pacnina_ensmean.eps')
#plot_map(hgt_ninamean,lats1,lons,'m', 80,'hgt_pacnina_ensmean.eps')
#
#print "Meaning Indian PIOD"
#mslp_piodmean, temp_piodmean, rain_piodmean, u_piodmean, v_piodmean, hgt_piodmean = mean_ensembles(indpiod_ens)
#plot_map(mslp_piodmean,lats,lons,'Pa',500,'mslp_indpiod_ensmean.eps')
#plot_map(temp_piodmean,lats,lons,'${}^{\circ}C$',3,'temp_indpiod_ensmean.eps')
#plot_map(rain_piodmean*86400,lats,lons,'$mm/day$', 7,'rain_indpiod_ensmean.eps')
#plot_stream(u_piodmean,v_piodmean,lats1,lons,'$ms^{-1}$',12,'winds_indpiod_ensmean.eps')
#plot_map(hgt_piodmean,lats1,lons,'m', 80,'hgt_indpiod_ensmean.eps')
#
#print "Meaning Indian NIOD"
#mslp_niodmean, temp_niodmean, rain_niodmean, u_niodmean, v_niodmean, hgt_niodmean = mean_ensembles(indniod_ens)
#plot_map(mslp_niodmean,lats,lons,'Pa',500,'mslp_indniod_ensmean.eps')
#plot_map(temp_niodmean,lats,lons,'${}^{\circ}C$',3,'temp_indniod_ensmean.eps')
#plot_map(rain_niodmean*86400,lats,lons,'$mm/day$', 7,'rain_indniod_ensmean.eps')
#plot_stream(u_niodmean,v_niodmean,lats1,lons,'$ms^{-1}$',12,'winds_indniod_ensmean.eps')
#plot_map(hgt_niodmean,lats1,lons,'m', 80,'hgt_indniod_ensmean.eps')