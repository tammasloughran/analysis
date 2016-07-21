"""nonlinear_lanina.py uses HadISST SSTs to construct a composite of sst anomalies for
weak Eastern Pacific La nina years where the Nino1+2 is greater than 1 standard deviation.
"""
import numpy as np
import netCDF4 as nc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap as bm
from scipy import signal

# Load data
hadisst_file = '/media/Jupiter/observations/HadISST/sst/HadISST_sst.nc'
hadisstnc = nc.Dataset(hadisst_file)
sst_ncvar = hadisstnc.variables['sst']
sst = sst_ncvar[:]
time_ncvar = hadisstnc.variables['time']
time = time_ncvar[:]
lats = hadisstnc.variables['lat'][:]
lons = hadisstnc.variables['lon'][:]

# Dates
time[-1] = time[-2] + (time[-2]-time[-3])
dates = nc.netcdftime.num2date(time, time_ncvar.units)
dates2 = pd.date_range(dates[0],dates[-1]+(dates[1]-dates[0]),freq='M')

# Select base period
base_period_sst = sst[(dates2.year>=1971)&(dates2.year<=2000),...]
base_period_dates = dates2[(dates2.year>=1971)&(dates2.year<=2000)]

# Calculate climatology
sst_clim = np.ones((12,)+base_period_sst.shape[-2:])*np.nan
for imonth in xrange(0,12):
    sst_clim[imonth,...] = base_period_sst[base_period_dates.month==(imonth+1),...].mean(axis=0)

# Subtract climatology
month = 0
sst_absolute = sst.copy()
for imonth in xrange(0,dates2.shape[0]):
    sst[imonth,...] = sst[imonth,...] - sst_clim[month,...]
    if month<11: month += 1
    else: month = 0

# Detrend
satera = dates2.year>=1970
sst = sst[satera] # Select 1970 and later
dates2 = dates2[satera]
mask = np.sum(sst.mask,axis=0)>0 # Create a mask for where data extsts for all time.
mask = np.broadcast_to(mask,sst.shape) # Broacast to shape of sst
sst = signal.detrend(sst,axis=0) # Detrend
sst = np.ma.array(sst,mask=mask) # Reapply mask

# Calculate Nino12 index
nino12_sst = sst[:,(lats<=0)&(lats>=-10),:]
nino12_sst = nino12_sst[:,:,(lons>=-90)&(lons<=-80)]
nino12 = np.zeros(nino12_sst.shape[0])
for i in xrange(nino12_sst.shape[0]): nino12[i] = np.mean(nino12_sst[i])

# Plot Nino3
fig, ax = plt.subplots()
rnge = 528
ax.plot_date(dates[-rnge:], nino12[-rnge:], '-')
ax.plot_date(dates[-rnge:], np.zeros(rnge), 'k')
ax.plot_date(dates[-rnge:], np.ones(rnge)*nino12.std(), 'r--')
ax.plot_date(dates[-rnge:], -np.ones(rnge)*nino12.std(), 'g--')
ax.xaxis.set_major_locator(YearLocator())
ax.xaxis.set_major_formatter(DateFormatter('%Y'))
ax.autoscale_view()
ax.grid(True)
ax.set_xlabel('Date')
ax.set_ylabel('Nino 1+2 ($^\circ$C)')
fig.autofmt_xdate()
plt.show()

# Define the years that have nonlinear lanina summers. The given year is the year containing January.
nlln_years = [2014, 2011, 2008, 2004, 1989, 1986, 1976, 1974]
nnllnyears = len(nlln_years)
nlln_composite = np.zeros((24,)+sst_absolute.shape[1:])
for year in nlln_years:
    sst_iyear = sst[(dates2.year==year)|(dates2.year==(year-1))]
    nlln_composite += sst_iyear
nlln_composite /= nnllnyears

# Plot a JJA composite
nllnpic = nlln_composite[5:8,...].mean(axis=0)
nllnpic[nllnpic < -5.] = np.nan
m = Basemap(projection='robin', lon_0=180.)
m.drawcoastlines()
parallels = np.arange(-90., 120., 30.)
meridians = np.arange(0., 360., 30.)
m.drawparallels(parallels,labels=[True,False,False,False], linewidth=0.3)
meridians = m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=0.3)
for mr in meridians:
    try: meridians[mr][1][0].set_rotation(50)
    except: pass
nllnpic, lons = bm.shiftgrid(1, nllnpic, lons)
xx,yy = np.meshgrid(lons, lats)
x,y = m(xx,yy)
fcont = m.contourf(x,y,nllnpic,cmap='seismic',levels=np.arange(-1.5,1.7,0.2))
m.contour(x,y,nllnpic,colors='k',levels=[0])
m.fillcontinents(color='gray',lake_color='gray')
cbar = plt.colorbar(fcont, orientation='horizontal')
cbar.ax.get_xaxis().set_ticks([])
for j, lab in enumerate([str(i) for i in np.arange(-1.5,1.7,0.2)]):
    cbar.ax.text((1/15.)*j, -0.5, lab, ha='center', va='center')
plt.title('JJA SSTA Eastern Pacific La Nina Composite')
plt.show()

# Load the low resolution anomalies and create the model compatible composite
low_nc = nc.Dataset('/media/Jupiter/observations/HadISST/HadISST_sst_data_N96_filled_anomaly.nc','r')
sst_low = low_nc.variables['sst'][:]
sst_low = sst_low[satera]
sst_low = signal.detrend(sst_low,axis=0) # Detrend
lats = low_nc.variables['lat'][:]
lons = low_nc.variables['lon'][:]
sst_iyear = 0
nlln_composite_low = np.zeros((24,)+sst_low.shape[1:])
for year in nlln_years:
    sst_iyear = sst_low[(dates2.year==year)|(dates2.year==(year-1))]
    nlln_composite_low += sst_iyear
nlln_composite_low /= nnllnyears

# Save to file
nllnnc = nc.Dataset('nllnsst.nc','w')
nllnnc.createDimension('time',24)
nllnnc.createDimension('lat',lats.shape[0])
nllnnc.createDimension('lon',lons.shape[0])
tout = nllnnc.createVariable('time','f8',('time'))
setattr(tout, 'units', 'month')
tout[:] = range(1,25)
latout = nllnnc.createVariable('lat','f8',('lat'))
setattr(latout, 'units', 'degrees_north')
latout[:] = lats
lonout = nllnnc.createVariable('lon','f8',('lon'))
setattr(lonout, 'units', 'degrees_north')
lonout[:] = lons
sstout = nllnnc.createVariable('sst','f8',('time','lat','lon'))
setattr(sstout, 'units', 'degC')
sstout[:] = nlln_composite_low
nllnnc.close()
