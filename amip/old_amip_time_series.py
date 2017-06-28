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

# Time series of model spread and uncertainty against ENSO background
# Define some variables
amipdir = '/srv/ccrc/data48/z5032520/amip/'
models = os.listdir(amipdir)
elnino_years = [1982, 1987, 1991, 1997, 2002]
lanina_years = [1984, 1988, 1998, 1999, 2007]


# Prepare figures
fig_hwf = plt.figure()
ax_hwf = fig_hwf.gca()
minorLocator = MultipleLocator(1)
for oyear in elnino_years:
    rect = patches.Rectangle((oyear-.5,0),1,45,color='lightcoral')
    ax_hwf.add_patch(rect)
for ayear in lanina_years:
    rect = patches.Rectangle((ayear-.5,0),1,45,color='lightblue')
    ax_hwf.add_patch(rect)

fig_hwa = plt.figure()
ax_hwa = fig_hwa.gca()
minorLocator = MultipleLocator(1)
for oyear in elnino_years:
    rect = patches.Rectangle((oyear-.5,0),1,45,color='lightcoral')
    ax_hwa.add_patch(rect)
for ayear in lanina_years:
    rect = patches.Rectangle((ayear-.5,0),1,45,color='lightblue')
    ax_hwa.add_patch(rect)

colors = ['darkred','red','yellow','olive','yellowgreen','darkblue','cyan','peru','green','orange','indigo','gray','lightgreen','khaki','teal']
nyears = 30
hwf_mmm = np.zeros(nyears)
hwd_mmm = np.zeros(nyears)
hwn_mmm = np.zeros(nyears)
hwa_mmm = np.zeros(nyears)
hwm_mmm = np.zeros(nyears)
for i,model in enumerate(models):
    print model
    
    hwfiles = glob.glob(amipdir+model+'/EHF*yearly*')
    model_color = colors[i]
    hwf_ensm = np.zeros(nyears)
    hwd_ensm = np.zeros(nyears)
    hwn_ensm = np.zeros(nyears)
    hwa_ensm = np.zeros(nyears)
    hwm_ensm = np.zeros(nyears)
    for ens,ifile in enumerate(hwfiles):
        hwnc = nc.Dataset(ifile)
        time = hwnc.variables['time'][:]
        period = (time>=1979)&(time<=2008)
        hwf = hwnc.variables['HWF_EHF'][period,...]
        hwd = hwnc.variables['HWD_EHF'][period,...]
        hwn = hwnc.variables['HWN_EHF'][period,...]
        hwa = hwnc.variables['HWA_EHF'][period,...]
        hwm = hwnc.variables['HWM_EHF'][period,...]
        lons = hwnc.variables['lon'][:]
        lats = hwnc.variables['lat'][:]
        time = hwnc.variables['time'][period]

        # Select a region to average over.
        ibox = (lons>=140.)&(lons<=150.)
        jbox = (lats>=-30.)&(lats<=-20.)
        hwf_series = np.mean(np.mean(hwf[:,jbox,:][:,:,ibox],axis=1),axis=1)
        hwd_series = np.mean(np.mean(hwd[:,jbox,:][:,:,ibox],axis=1),axis=1)
        hwn_series = np.mean(np.mean(hwn[:,jbox,:][:,:,ibox],axis=1),axis=1)
        hwa_series = np.mean(np.mean(hwa[:,jbox,:][:,:,ibox],axis=1),axis=1)
        hwm_series = np.mean(np.mean(hwm[:,jbox,:][:,:,ibox],axis=1),axis=1)
        
        # Agregate ensemble mean
        hwf_ensm += hwf_series
        hwd_ensm += hwd_series
        hwn_ensm += hwn_series
        hwa_ensm += hwa_series
        hwm_ensm += hwm_series
    hwf_ensm /= float(ens+1)
    hwd_ensm /= float(ens+1)
    hwn_ensm /= float(ens+1)
    hwa_ensm /= float(ens+1)
    hwm_ensm /= float(ens+1)

    # Plot
    ax_hwf.plot(time.astype(int), hwf_ensm, color=model_color, label=model)
    ax_hwa.plot(time.astype(int), hwa_ensm, color=model_color, label=model)
    
    #aggregate mmm
    hwf_mmm += hwf_ensm
    hwd_mmm += hwd_ensm
    hwn_mmm += hwn_ensm
    hwa_mmm += hwa_ensm
    hwm_mmm += hwm_ensm
hwf_mmm /= 15.
hwd_mmm /= 15.
hwn_mmm /= 15.
hwa_mmm /= 15.
hwm_mmm /= 15.

plt.figure(1)
plt.plot(time.astype(int), hwf_mmm, color='k', linewidth=2, label='MMM')
plt.legend(fontsize=6, loc=9)
plt.xlim((time[0]-1,time[-1]+1))
plt.xlabel('Time')
ax_hwf.xaxis.set_minor_locator(minorLocator)
ax_hwf.spines['right'].set_color('none')
ax_hwf.spines['top'].set_color('none')
plt.ylabel('Days',rotation=0.)
plt.title('AMIP Australian Heatwave Frequency')
plt.savefig('amip_hwf.eps',format='eps')

plt.figure(2)
plt.plot(time.astype(int), hwa_mmm, color='k', linewidth=2, label='MMM')
plt.legend(fontsize=6, loc=9)
plt.xlim((time[0]-1,time[-1]+1))
plt.xlabel('Time')
ax_hwa.xaxis.set_minor_locator(minorLocator)
ax_hwa.spines['right'].set_color('none')
ax_hwa.spines['top'].set_color('none')
plt.ylabel('$^{\circ}C^{2}$',rotation=0.)
plt.title('AMIP Australian Heatwave Amplitude')
plt.savefig('amip_hwa.eps',format='eps')