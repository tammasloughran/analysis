"""Make sub tropical ridge Tasman high (STRH) index from 20C data.
"""
from netCDF4 import Dataset
import glob
import datetime as dt
import numpy as np
import pdb

# Load grid data
twentyc_dir = '/media/Jupiter/reanalysis/20crv2/prmsl/'
file_list = glob.glob(twentyc_dir+'prmsl.????.nc')
ncget = Dataset(file_list[0])
lat = ncget.variables['lat'][:]
lon = ncget.variables['lon'][:]

# Define base period
base_start = dt.datetime(1961,1,1,0,0)
base_end = dt.datetime(1991,1,1,0,0)
base_length = base_end - base_start

# Define the Tasman region
y_tasman = np.where((lat<=0)&(lat>=-60))[0]
x_tasman = np.where((lon>=150)&(lon<=166))[0]
lat_tasman = lat[y_tasman]
lon_tasman = lon[x_tasman]

# Define time stamps
deltas = np.array([dt.timedelta(days=i) for i in range(base_length.days)])
dates = base_start + deltas
months = np.array([dates[i].month for i in range(dates.size)])
years = np.array([dates[i].year for i in range(dates.size)])

# Load and cut out Tasman pressure data over the base period
pr_base = np.zeros([base_length.days,lat_tasman.size,lon_tasman.size])
for cyear in range(1961,1991,1):
    ncdata = Dataset(twentyc_dir+'prmsl.%s.nc'%(cyear))
    pr_cyear = ncdata.variables['prmsl'][:]
    pr_cyear = pr_cyear[:,y_tasman,:]
    pr_cyear = pr_cyear[:,:,x_tasman]
    pr_base[years==cyear,:,:] = pr_cyear

# Calculate the climatology over the base period
pr_clim = np.zeros([12,30,lat_tasman.size,lon_tasman.size])
for cyear in range(1961,1991,1):
    iyear = cyear-1961
    for cmonth in range(1,13,1):
        imonth = cmonth-1
        pr_clim[imonth,iyear,:,:] = \
                pr_base[(years==cyear)|(months==cmonth),:,:].mean(axis=0)
pr_clim = pr_clim.mean(axis=1) # Average across years axis
pr_clim = pr_clim.mean(axis=2) # Average across tasman longitudes

# Find the latitude of maximum as the subtropical ridge.
maximum = np.amax(pr_clim, axis=1)
ridge_lat = np.zeros((12))
for m in range(0,12,1):
    iridge = np.where(pr_clim[m,:]==maximum[m])[0]
    ridge_lat[m] = lat_tasman[iridge]
box_lats_upper = ridge_lat + 8.
box_lats_lower = ridge_lat - 8.
box_lons = [150., 166.]

# Load and calculate summer anomalies.
STRH = np.zeros([102,12])
reference = dt.datetime(1800,1,1,0,0)
for cyear in range(1911,2013,1):
    iyear = cyear-1911
    # Load
    ncdata = Dataset(twentyc_dir+'prmsl.%s.nc'%(cyear))
    time = ncdata.variables['time'][:]
    deltas = np.array([dt.timedelta(hours=time[i]) for i in range(time.size)])
    dates = reference + deltas
    months = np.array([dates[i].month for i in range(dates.size)])
    pr_cyear = ncdata.variables['prmsl'][:]
    pr_monthly = np.zeros([12,lat.size,lon.size])
    # Monthly means
    for imonth in range(0,12,1):
        cmonth = imonth+1
        pr_monthly[imonth,:,:] = pr_cyear[months==cmonth].mean(axis=0)
        y1 = box_lats_lower[imonth]
        y2 = box_lats_upper[imonth]
        y_box = np.where((lat>=y1)&(lat<=y2))[0]
        x_box = np.where((lon>=150)&(lon<=165))[0]
        box = pr_monthly[imonth,y_box,:]
        box = box[:,x_box]
        STRH[iyear, imonth] = box.mean()

# Calculate the STRH climatology
years = np.arange(1911,2013)
iyears = np.where((years>=1961)&(years<=1990))[0]
STRH_base = STRH[iyears]
STRH_clim = STRH_base.mean(axis=0)

# Express STRH as anomaly and save to a file
outfile = open('STRH_monthly_1911_2012_20CRv2.txt','w')
for year in range(years.size):
    outfile.write(str(years[year])+" ")
    for month in range(0,12,1):
        STRH[year,month] = STRH[year,month]-STRH_clim[month]
        outfile.write(str(round(STRH[year,month],3))+" ")
    outfile.write("\n")
outfile.close()
