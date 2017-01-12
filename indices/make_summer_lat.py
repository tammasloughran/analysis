"""Make latitude of summer Tasman high index from 20C data.
"""
from netCDF4 import Dataset, MFDataset
import glob
import datetime as dt
import numpy as np
import pdb
import matplotlib.pyplot as plt

# Data directory
twentyc_dir = '/media/Jupiter/reanalysis/20crv2/prmsl/'

# Load time and grid data
nc = MFDataset(twentyc_dir+"prmsl.*.nc", 'r')
time = nc.variables['time'][:]
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]

# Define the Tasman region
y_tasman = np.where((lat<=0)&(lat>=-60))[0]
x_tasman = np.where((lon>=150)&(lon<=166))[0]
lat_tasman = lat[y_tasman]
lon_tasman = lon[x_tasman]

# Make time arrays
deltas = np.array([dt.timedelta(hours=time[i]) for i in range(time.size)])
reference = dt.datetime(1800,1,1,0,0)
dates = reference + deltas
months = np.array([dates[i].month for i in range(dates.size)])
years = np.array([dates[i].year for i in range(dates.size)])

# Select summer months
summer_months = (months==11)|(months==12)|(months==1)|(months==2)|(months==3)
pr_months = months[summer_months]
pr_years = years[summer_months]

# Shift years array by a few months so that summer is contained within the same year
shift_months = (pr_months==1)|(pr_months==2)|(pr_months==3)
pr_years[shift_months] = pr_years[shift_months]-1

# Load summer pressure
pr = nc.variables['prmsl'][pr_months]

# Cut out Tasman region
pr = pr[:,y_tasman,:]
pr = pr[:,:,x_tasman]

# Average over summer months
years = range(years[0],years[-1],1)
summer_pr = np.zeros((len(years),pr.shape[1],pr.shape[2]))
for iyear, year in enumerate(years):
    summer_pr[iyear,...] = pr[pr_years==year,...].mean(axis=0)

# Average over longitudes
summer_pr = summer_pr.mean(axis=2)

# Identify the maximum
maximum = np.amax(summer_pr, axis=1)

# Identify the latitude of the max
max_lat = np.zeros(len(years))
for iyear in xrange(0,len(years)):
    ilat = np.where(summer_pr[iyear,:]==maximum[iyear])[0]
    max_lat[iyear] = lat_tasman[ilat]

# Save to file
np.savetxt('tasman_max_lats.txt',max_lat)
fout = open('fileout.txt','w')
for iyear, year in enumerate(years):
    fout.write(year+' '+max_lat[iyear])
