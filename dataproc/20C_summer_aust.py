"""Turns daily 20C reanalysis mean sea level pressure data into 
summer anomalies.
"""
# Contact: Tammas Loughran, t.loughran@student.unsw.edu.au
#
# Modules:
from netCDF4 import Dataset
import glob
import datetime as dt
import numpy as np

# Start code
# Define data and dates
twentyc_dir = '/media/Jupiter/reanalysis/20crv2/prmsl/'
file_list = glob.glob(twentyc_dir+'aus_prmsl.????.nc')
ncget = Dataset(file_list[0])
lat = ncget.variables['lat'][:]
lon = ncget.variables['lon'][:]
base_start = dt.datetime(1961,1,1,0,0)
base_end = dt.datetime(1991,1,1,0,0)
base_length = base_end - base_start
pr_base = np.zeros([base_length.days,lat.size,lon.size])
deltas = np.array([dt.timedelta(days=i) for i in range(base_length.days)])
dates = base_start + deltas
months = np.array([dates[i].month for i in range(dates.size)])
years = np.array([dates[i].year for i in range(dates.size)])

# Get base period data (1961-1990)
for cyear in range(1961,1991,1):
    ncdata = Dataset(twentyc_dir+'aus_prmsl.%s.nc'%(cyear))
    pr_cyear = ncdata.variables['prmsl'][:]
    pr_base[years==cyear,:,:] = pr_cyear

# Calculate climatology over base period 
pr_clim = np.zeros((12,30,lat.size,lon.size))
for cyear in range(1961,1991,1):
    iyear = cyear-1961
    for cmonth in range(1,13,1):
        imonth = cmonth-1
        pr_clim[imonth,iyear,:,:] = \
                pr_base[(years==cyear)|(months==cmonth),:,:].mean(axis=0)
pr_clim = pr_clim.mean(axis=1)

# Load and calculate summer anomalies.
pr_summer = np.zeros([101,lat.size,lon.size])
pr_monthly = np.zeros([12,lat.size,lon.size])
reference = dt.datetime(1800,1,1,0,0)
for cyear in range(1911,2013,1):
    iyear = cyear-1911
    # Load
    ncdata = Dataset(twentyc_dir+'aus_prmsl.%s.nc'%(cyear))
    time = ncdata.variables['time'][:]
    deltas = np.array([dt.timedelta(hours=time[i]) for i in range(time.size)])
    dates = reference + deltas
    months = np.array([dates[i].month for i in range(dates.size)])
    pr_cyear = ncdata.variables['prmsl'][:]
    # Monthly means
    for imonth in range(0,12,1):
        cmonth = imonth+1
        pr_monthly[imonth,:,:] = pr_cyear[months==cmonth].mean(axis=0)
    # Anomalies
    pr_anom = pr_monthly - pr_clim
    # Summer mean (NOV of previous year to MAR of the current year)
    if cyear!=1911:
        pr2 = pr_anom[0:3,:,:]
        pr3 = np.append(pr1,pr2,axis=0)
        pr_summer[iyear-1,:,:] = pr3.mean(axis=0)
    # NOV-DEC Pressure data for next year
    pr1 = pr_anom[10:12,:,:]

# Output file
outfile = Dataset("summer_mslp_1911-2013.nc",'w')
# Dims
outfile.createDimension("time", size=pr_summer.shape[0])
outfile.createDimension("lat", size=pr_summer.shape[1])
outfile.createDimension("lon", size=pr_summer.shape[2])
# Vars
otime = outfile.createVariable("time",'i',dimensions=('time'))
setattr(otime,"long_name","Time")
setattr(otime,"units","year")
setattr(otime,"standard_name","time")
olat = outfile.createVariable("lat",'f',dimensions=('lat'))
setattr(olat,"long_name","Latitude")
setattr(olat,"units","degrees_north")
setattr(olat,"standard_name","latitude")
olon = outfile.createVariable("lon",'f',dimensions=('lon'))
setattr(olon,"long_name","Longitude")
setattr(olon,"units","degrees_east")
setattr(olon,"standard_name","longitude")
omslp = outfile.createVariable("mslp",'f',dimensions=('time','lat','lon'))
setattr(omslp,"long_name","Summer Mean Sea Level Pressure Anomaly")
setattr(omslp,"units","Pa")
setattr(omslp,"standard_name","air_pressure_anomaly")
setattr(omslp,"missing_value",-9.96921e+36)
otime[:] = range(1911,2012,1)
olat[:] = lat
olon[:] = lon
omslp[:] = pr_summer
outfile.close()
