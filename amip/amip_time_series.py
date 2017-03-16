# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 15:55:21 2017

@author: Tammas Loughran
"""
import os
import glob
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MultipleLocator
# User specified parameters
aspect = 'HWF'


# Time series of model spread and uncertainty against ENSO background
# Define some variables
amipdir = '/srv/ccrc/data48/z5032520/amip/'
models = os.listdir(amipdir)
elnino_years = [1982, 1987, 1991, 1997, 2002]
lanina_years = [1984, 1988, 1998, 1999, 2007]


# Prepare figures
fig = plt.figure()
axes = fig.gca()
minorLocator = MultipleLocator(1)
for oyear in elnino_years:
    rect = patches.Rectangle((oyear-.5,0),1,45,color='lightcoral')
    axes.add_patch(rect)
for ayear in lanina_years:
    rect = patches.Rectangle((ayear-.5,0),1,45,color='lightblue')
    axes.add_patch(rect)


colors = ['darkred','red','yellow','olive','yellowgreen','darkblue','cyan','peru','green','orange','indigo','gray','lightgreen','khaki','teal']
colors = ['gray','red','gray','gray','gray','gray','gray','gray','gray','gray','gray','gray','gray','gray','gray']
nyears = 29
mmm = np.zeros(nyears)
all_series = np.ones((len(models),nyears))*np.nan
for i,model in enumerate(models):
    print model
    times = 0
    hwfiles = glob.glob(amipdir+model+'/EHF*yearly*')
    ifile = hwfiles[0]
    model_color = colors[i]
    hwnc = nc.Dataset(ifile)
    times = hwnc.variables['time'][:]
    period = (times>=1979)&(times<=2007)
    times = times[period]
    var = hwnc.variables[aspect+'_EHF'][period,...]
    lons = hwnc.variables['lon'][:]
    lats = hwnc.variables['lat'][:]
    # Select a region to average over.
    ibox = (lons>=140.)&(lons<=150.)
    jbox = (lats>=-30.)&(lats<=-20.)
    var_series = np.mean(np.mean(var[:,jbox,:][:,:,ibox],axis=1),axis=1)
    all_series[i,:] = var_series
    # Plot
    plt.plot(times.astype(int), var_series, color=model_color, label=model)
    #aggregate mmm
    mmm += var_series
mmm /= float(len(models))

awapnc = nc.Dataset('/srv/ccrc/data35/z5032520/AWAP/yearly/ehfhw/EHF_heatwaves_AWAP_bp1979-2008_yearly_summer.nc')
times = awapnc.variables['time'][:]
period = (times>=1979)&(times<=2007)
var = awapnc.variables[aspect+'_EHF'][period,...]
lons = awapnc.variables['lon'][:]
lats = awapnc.variables['lat'][:]
times = times[period]
ibox = (lons>=140.)&(lons<=150.)
jbox = (lats>=-30.)&(lats<=-20.)
awap = np.mean(np.mean(var[:,jbox,:][:,:,ibox],axis=1),axis=1)
plt.plot(times.astype(int), awap, color='black', linestyle='dotted', label='AWAP')
plt.plot(times.astype(int), mmm, color='k', linewidth=2, label='MMM')
plt.legend(fontsize=6, loc=2, ncol=3)
plt.xlim((times[0]-1,times[-1]+1))
plt.xlabel('Time')
axes.xaxis.set_minor_locator(minorLocator)
axes.spines['right'].set_visible(False)
axes.yaxis.set_ticks_position('left')
axes.spines['top'].set_visible(False)
axes.xaxis.set_ticks_position('bottom')
plt.ylabel('Days',rotation=0.)
plt.title('AMIP Australian Heatwave Frequency')
plt.savefig('amip_hwf.eps',format='eps')

fig = plt.figure()
axes = fig.gca()
nino_box = [4,9,13,19,24]
nina_box = [6,10,20,21,29]
for oyear in nino_box:
    rect = patches.Rectangle((oyear-.5,0),1,45,color='lightcoral')
    axes.add_patch(rect)
for ayear in nina_box:
    rect = patches.Rectangle((ayear-.5,0),1,45,color='lightblue')
    axes.add_patch(rect)
bp = plt.boxplot(all_series,whis=99)
plt.plot(np.arange(1,30,1),awap, color='black', linestyle='-', label='AWAP')
axes.set_xticks([2,7,12,17,22,27])
axes.set_xticklabels(['1980','1985','1990','1995','2000','2005'])
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black',linestyle='-')
axes.xaxis.set_minor_locator(minorLocator)
axes.spines['right'].set_visible(False)
axes.yaxis.set_ticks_position('left')
axes.spines['top'].set_visible(False)
axes.xaxis.set_ticks_position('bottom')
plt.xlabel('Time')
plt.ylabel('Days',rotation=0.)
plt.title('AMIP Spread of HWF in Northeastern Australia')
plt.savefig('AMIP_HWF_boxplot.eps',format='eps')