# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 16:28:44 2017

@author: Tammas Loughran
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.basemap import Basemap
import scipy.stats as stats
import matplotlib

def plot_map(data,lats,lons,units,rng,ax):
    plt.clf()
    mp = Basemap(ax=ax,projection='robin', lon_0=180)
    x,y = np.meshgrid(lons,lats)
    xx,yy = mp(x,y)
    if units!='$mm/day$':
        colours = 'bwr'
        plot_cont = True
    else: 
        colours = 'bwr_r'
        plot_cont = False
    #clrmsh = mp.pcolormesh(xx,yy,data.squeeze(),vmin=-rng,vmax=rng, cmap=colours)
    clrmsh = mp.contourf(xx,yy,data.squeeze(),levels=range(-rng,rng+1,60), cmap=colours,extend='both')
    if plot_cont:
        step = 2.*rng/10.
        levels = np.arange(-rng,rng+step,step)
        cont = mp.contour(xx,yy,data.squeeze(),colors='k',levels=levels,linewidths=0.5)
    for c in cont.collections:
        if c.get_linestyle() == [(None, None)]:
           continue
        else:
            c.set_dashes([(0, (2.0, 2.0))])
    eqtr = data[(lats>=-5)&(lats<=5)]
    lt = [4,2,0,-2,-4]
    eqmax = np.where(eqtr==np.max(eqtr))
    eqmax = (lons[eqmax[1][0]],lt[eqmax[0][0]])
    eqmin = np.where(eqtr==np.min(eqtr))
    eqmin = (lons[eqmin[1][0]],lt[eqmin[0][0]])
    mp.scatter(eqmax[0],eqmax[1],s=35,c='grey',marker='D',latlon=True)
    mp.scatter(eqmin[0],eqmin[1],s=35,c='grey',marker='*',latlon=True)
    #mp.scatter(eqmin,'k',latlon=True)
    cb = mp.colorbar(clrmsh, location='bottom', pad=0.3)
    cb.set_label(units)
    mp.drawcoastlines()
    mp.drawparallels([-80,-40,-60,-20,0,20,40,60,80],labels=[1, 0, 0, 0],linewidth=0)
    mp.drawmeridians([0,80,160,240,320],labels=[0, 0, 0, 1],linewidth=0)
    #plt.savefig(fname,format='eps')


# Define modoki years
modoki_years = [ 1986,
                1990, 1991, 1992, 1994, 2002, 2004, 2009]
                #1923, 1929, 1940, 1946, 1958, 1963, 1977,
elnino_years = [1911, 1913, 1914, 1918, 1925, 1930, 1941, 1951,
                1957, 1965, 1969, 1972, 1976, 1982, 1987, 1997, 2006]
                
# Load the pressure data
twntycrdir = '/media/Jupiter/reanalysis/20crv2/prmsl'
ncs = nc.MFDataset(twntycrdir+'/monthly_prmsl.????.nc','r')
lons = ncs.variables['lon'][:]
lats = ncs.variables['lat'][:]
dates = pd.date_range(start='1871-01-01',end='2012-12-31',freq='M')
prmsl = ncs.variables['prmsl'][...]
base_summer = ((dates.month==12)|(dates.month==1)|(dates.month==2))&(dates.year>=1961)
prmsl_clim = prmsl[base_summer].mean(axis=0)
for t in xrange(0,prmsl.shape[0]):
    prmsl[t,...] = prmsl[t,...] - prmsl_clim



shift_year = np.roll(dates.year,2)
shift_year[0:2] = 0
modoki_pr = np.zeros((len(modoki_years),len(lats),len(lons)))
for i,year in enumerate(modoki_years):
    summer = ((dates.month==12)|(dates.month==1)|(dates.month==2))&(shift_year==year)
    modoki_pr[i,...] = prmsl[summer].mean(axis=0)
modoki_pr = modoki_pr.mean(axis=0)
#plot_map(modoki_pr,lats,lons,'Pa',300,ax2)

elnino_pr = np.zeros((len(elnino_years),len(lats),len(lons)))
for i,year in enumerate(elnino_years):
    summer = ((dates.month==12)|(dates.month==1)|(dates.month==2))&(shift_year==year)
    elnino_pr[i,...] = prmsl[summer].mean(axis=0)
elnino_pr = elnino_pr.mean(axis=0)
#plot_map(elnino_pr,lats,lons,'Pa',300,ax1)

fig,(ax1,ax2) = plt.subplots(nrows=2, ncols=1,figsize=(6,8))
units='Pa'

mp = Basemap(ax=ax1,projection='robin', lon_0=180)
x,y = np.meshgrid(lons,lats)
xx,yy = mp(x,y)
colours = 'bwr'
#clrmsh = mp.pcolormesh(xx,yy,data.squeeze(),vmin=-rng,vmax=rng, cmap=colours)
clrmsh = mp.contourf(xx,yy,elnino_pr.squeeze(),levels=range(-300,300+1,60), cmap=colours,extend='both')
step = 2.*300/10.
levels = np.arange(-300,300+step,step)
cont = mp.contour(xx,yy,elnino_pr.squeeze(),colors='k',levels=levels,linewidths=0.5)
for c in cont.collections:
    if c.get_linestyle() == [(None, None)]:
       continue
    else:
        c.set_dashes([(0, (2.0, 2.0))])
eqtr = elnino_pr[(lats>=-5)&(lats<=5)]
lt = [4,2,0,-2,-4]
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
clrmsh = mp.contourf(xx,yy,modoki_pr.squeeze(),levels=range(-300,300+1,60), cmap=colours,extend='both')
step = 2.*300/10.
levels = np.arange(-300,300+step,step)
cont = mp.contour(xx,yy,modoki_pr.squeeze(),colors='k',levels=levels,linewidths=0.5)
for c in cont.collections:
    if c.get_linestyle() == [(None, None)]:
       continue
    else:
        c.set_dashes([(0, (2.0, 2.0))])
eqtr = modoki_pr[(lats>=-5)&(lats<=5)]
lt = [4,2,0,-2,-4]
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

plt.savefig('observed_mslp_shift.eps',format='eps')
