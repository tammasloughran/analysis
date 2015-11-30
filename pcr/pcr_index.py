import numpy as np
from netCDF4 import Dataset
from scipy.stats import linregress
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import datetime as dt
import load_data
from pandas import concat

# Principle components from PCA of heatwaves
hwn_pcs = np.load('HWN_rotated_pcs.npy')[:-1]
hwf_pcs = np.load('HWF_rotated_pcs.npy')[:-1]
hwd_pcs = np.load('HWD_rotated_pcs.npy')[:-1]
hwa_pcs = np.load('HWA_rotated_pcs.npy')[:-1]
hwm_pcs = np.load('HWM_rotated_pcs.npy')[:-1]
hwt_pcs = np.load('HWT_rotated_pcs.npy')[:-1]

# Load indices
directory = ('/srv/ccrc/data35/z5032520/')
fname = 'indices/Nino3.4/Nino3.4_monthly_1870_2014_hadisst.dat'
nino34 = load_data.load_index(directory+fname)
fname = 'indices/DMI/DMI_monthly_1870_2014_hadisst.dat'
dmi = load_data.load_index(directory+fname)
fname = 'indices/SOI/SOI_monthly_1866_2014_UCAR.txt'
soi = load_data.load_index2(directory+fname)
fname = 'indices/SAM/SAM_index_monthly_visbeck.txt'
sam1 = load_data.load_index2(directory+fname)
fname = 'indices/SAM/SAM_monthly_1957_2014_BAS.txt'
sam2 = load_data.load_index2(directory+fname)
sam = concat([sam1, sam2])
fname = 'indices/STRH/STRH_monthly_1911_2012_20CRv2.txt'
strh = load_data.load_index2(directory+fname)

start_year = 1911
end_year = 2012
year_range = np.arange(1911, 2014+1, 1)

# Take spring (9,10,11)/summer (11,12,1,2,3) (A)nnual (mean) (S)tarting in
# (JUL)y. i.e. 'AS-JUL'i
# SAM
samslice = concat([sam1,sam2])
samslice = samslice['%s-07'%(start_year):'%s-06'%(end_year)]
axis = samslice.index.month
months = (axis==11)|(axis==12)|(axis==1)|(axis==2)|(axis==3)
samslice = samslice[months]
samslice = samslice.resample('AS-JUL', how='mean')
# Nino3.4
ninoslice = nino34['%s-07'%(start_year):'%s-06'%(end_year)]
axis = ninoslice.index.month
months = (axis==11)|(axis==12)|(axis==1)|(axis==2)|(axis==3)
ninoslice = ninoslice[months]
ninoslice = ninoslice.resample('AS-JUL', how='mean')
# DMI
dmislice = dmi['%s-07'%(start_year):'%s-06'%(end_year)]
dmislice = dmislice[months]
dmislice = dmislice.resample('AS-JUL', how='mean')
# SOI
soislice = soi['%s-07'%(start_year):'%s-06'%(end_year)]
soislice = soislice[months]
soislice = soislice.resample('AS-JUL', how='mean')
# STRH
strhslice = strh['%s-07'%(start_year):'%s-06'%(end_year)]
strhslice = strhslice[months]
strhslice = strhslice.resample('AS-JUL', how='mean')

def plot_pcscatter(index, pc, slope, intercept, correlation, p, index_name, pc_name, fname):
    fig, ax = plt.subplots()
    ax.axhline(y=0, color='k')
    ax.axvline(x=0, color='k')
    plt.scatter(index,pc,color='k')
    plt.plot(index,slope*index+intercept,color='k')
    ax.set_xlabel(index_name)
    ax.set_ylabel(pc_name)
    if p<0.01:
        p = '< 0.01'
    else:
        p = '= %s'%(round(p,3))
    ax.text(-3, 2, 'slope = %s\nr = %s\np '%(round(slope,2),round(correlation,2))+p, fontsize=15)
    plt.savefig(fname, format='eps')
    plt.close()

slope, intercept, correlation, p, error = linregress(ninoslice,hwn_pcs[:,0])
plot_pcscatter(ninoslice, hwn_pcs[:,0], slope, intercept, correlation, p, 'Nino3.4 index ($^\circ$C)', 'HWN rPC1', 'nino3.4_HWNrpc1.eps')
slope, intercept, correlation, p, error = linregress(ninoslice,hwn_pcs[:,1])
plot_pcscatter(ninoslice, hwn_pcs[:,1], slope, intercept, correlation, p, 'Nino3.4 index ($^\circ$C)', 'HWN rPC2', 'nino3.4_HWNrpc2.eps')
slope, intercept, correlation, p, error = linregress(ninoslice,hwn_pcs[:,2])
plot_pcscatter(ninoslice, hwn_pcs[:,2], slope, intercept, correlation, p, 'Nino3.4 index ($^\circ$C)', 'HWN rPC3', 'nino3.4_HWNrpc3.eps')
slope, intercept, correlation, p, error = linregress(ninoslice,hwn_pcs[:,3])
plot_pcscatter(ninoslice, hwn_pcs[:,3], slope, intercept, correlation, p, 'Nino3.4 index ($^\circ$C)', 'HWN rPC4', 'nino3.4_HWNrpc4.eps')

slope, intercept, correlation, p, error = linregress(ninoslice,hwd_pcs[:,0])
plot_pcscatter(ninoslice, hwd_pcs[:,0], slope, intercept, correlation, p, 'Nino3.4 index ($^\circ$C)', 'HWD rPC1', 'nino3.4_HWDrpc1.eps')
slope, intercept, correlation, p, error = linregress(ninoslice,hwd_pcs[:,1])
plot_pcscatter(ninoslice, hwd_pcs[:,1], slope, intercept, correlation, p, 'Nino3.4 index ($^\circ$C)', 'HWD rPC2', 'nino3.4_HWDrpc2.eps')
slope, intercept, correlation, p, error = linregress(ninoslice,hwd_pcs[:,2])
plot_pcscatter(ninoslice, hwd_pcs[:,2], slope, intercept, correlation, p, 'Nino3.4 index ($^\circ$C)', 'HWD rPC3', 'nino3.4_HWDrpc3.eps')
slope, intercept, correlation, p, error = linregress(ninoslice,hwd_pcs[:,3])
plot_pcscatter(ninoslice, hwd_pcs[:,3], slope, intercept, correlation, p, 'Nino3.4 index ($^\circ$C)', 'HWD rPC4', 'nino3.4_HWDrpc4.eps')

slope, intercept, correlation, p, error = linregress(samslice,hwn_pcs[:,0])
plot_pcscatter(samslice, hwa_pcs[:,0], slope, intercept, correlation, p, 'SAM', 'HWN rPC1', 'SAM_HWNrpc1.eps')
slope, intercept, correlation, p, error = linregress(samslice,hwn_pcs[:,1])
plot_pcscatter(samslice, hwa_pcs[:,1], slope, intercept, correlation, p, 'SAM', 'HWN rPC2', 'SAM_HWNrpc2.eps')
slope, intercept, correlation, p, error = linregress(samslice,hwn_pcs[:,2])
plot_pcscatter(samslice, hwa_pcs[:,2], slope, intercept, correlation, p, 'SAM', 'HWN rPC3', 'SAM_HWNrpc3.eps')
slope, intercept, correlation, p, error = linregress(samslice,hwn_pcs[:,3])
plot_pcscatter(samslice, hwa_pcs[:,3], slope, intercept, correlation, p, 'SAM', 'HWN rPC4', 'SAM_HWNrpc4.eps')

slope, intercept, correlation, p, error = linregress(strhslice,hwa_pcs[:,0])
plot_pcscatter(strhslice, hwa_pcs[:,0], slope, intercept, correlation, p, 'STRH', 'HWA rPC1', 'STRH_HWArpc1.eps')
slope, intercept, correlation, p, error = linregress(strhslice,hwa_pcs[:,1])
plot_pcscatter(strhslice, hwa_pcs[:,1], slope, intercept, correlation, p, 'STRH', 'HWA rPC2', 'STRH_HWArpc2.eps')
slope, intercept, correlation, p, error = linregress(strhslice,hwa_pcs[:,2])
plot_pcscatter(strhslice, hwa_pcs[:,2], slope, intercept, correlation, p, 'STRH', 'HWA rPC3', 'STRH_HWArpc3.eps')
slope, intercept, correlation, p, error = linregress(strhslice,hwa_pcs[:,3])
plot_pcscatter(strhslice, hwa_pcs[:,3], slope, intercept, correlation, p, 'STRH', 'HWA rPC4', 'STRH_HWArpc4.eps')
