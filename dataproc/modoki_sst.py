"""modoki_sst.py uses HadISST SSTs to construct a composite of sst anomalies for
modoki years where the EMI is greater than 1 standard deviation.
"""
import numpy as np
import netCDF4 as nc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap as bm
from scipy import signal
import sys

# Load data
hadisst_file = '/media/Jupiter/observations/HadISST/sst/HadISST_sst.nc'
hadisstnc = nc.Dataset(hadisst_file)
sst_ncvar = hadisstnc.variables['sst']
sst = sst_ncvar[:]
time_ncvar = hadisstnc.variables['time']
time = time_ncvar[:]
lats = hadisstnc.variables['lat'][:]
lons = hadisstnc.variables['lon'][:]


ninofile = '/srv/ccrc/data35/z5032520/ancilforge/modal_anomalies.nc'
ninonc = nc.Dataset(ninofile,'r')
nino = ninonc.variables['elnino_sst'][:]
ninopic = nino[11:14,...].mean(axis=0)
lons2 = ninonc.variables['lon'][:]
lats2 = ninonc.variables['lat'][:]

# Dates
time[-1] = time[-2] + (time[-2]-time[-3])
dates = nc.netcdftime.num2date(time, time_ncvar.units)
dates2 = pd.date_range(dates[0],dates[-1]+(dates[1]-dates[0]),freq='M')
hadisstnc.close()
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
for imonth in xrange(0,sst.shape[0]):
    sst[imonth,...] = sst[imonth,...] - sst_clim[month,...]
    if month<11: month += 1
    else: month = 0

# Detrend
sel_from = dates2.year>=1911
sst = sst[sel_from] # Select 1911 and later
dates2 = dates2[sel_from]
mask = np.sum(sst.mask,axis=0)>0 # Create a mask for where data extsts for all time.
mask = np.broadcast_to(mask,sst.shape) # Broacast to shape of sst
sst = signal.detrend(sst,axis=0) # Detrend
sst = np.ma.array(sst,mask=mask) # Reapply mask

def plt_dec(data,year):
    m = Basemap(projection='robin', lon_0=180.)
    m.drawcoastlines()
    parallels = np.arange(-90., 120., 30.)
    meridians = np.arange(0., 360., 30.)
    m.drawparallels(parallels,labels=[True,False,False,False], linewidth=0.3)
    meridians = m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=0.3)
    for mr in meridians:
        try: meridians[mr][1][0].set_rotation(50)
        except: pass
    fafa, ulons = bm.shiftgrid(1, data, lons)
    xx,yy = np.meshgrid(ulons, lats)
    x,y = m(xx,yy)
    fcont = m.contourf(x,y,fafa,cmap='seismic',levels=np.arange(-2,2.2,0.2),extend='both')
    m.contour(x,y,fafa,colors='k',levels=[0])
    m.fillcontinents(color='gray',lake_color='gray')
    cbar = plt.colorbar(fcont, orientation='horizontal')
    #cbar.ax.get_xaxis().set_ticks([])
    #for j, lab in enumerate([str(i) for i in np.arange(-2,2.2,0.2)]):
    #    cbar.ax.text((1/15.)*j, -0.5, lab, ha='center', va='center')
    plt.title(str(year))
    plt.savefig(str(year)+'_dec.eps',format='eps')
    plt.close()

#i = -1
#for iyear in xrange(1911,2015,1):
#    i += 12
#    plt_dec(sst[i,...],iyear)
#sys.exit()
# Calculate EMI
# Regions
sst_a = sst[:,(lats<=10)&(lats>=-10),:]
sst_a = sst_a[:,:,(lons<=-140)|(lons>=165)]
emia = np.zeros(sst_a.shape[0])
for i in xrange(emia.shape[0]): emia[i] = np.mean(sst_a[i])
sst_b = sst[:,(lats<=5)&(lats>=-15),:]
sst_b = sst_b[:,:,(lons>=-110)&(lons<=-70)]
emib = np.zeros(sst_b.shape[0])
for i in xrange(emib.shape[0]): emib[i] = np.mean(sst_b[i])
sst_c = sst[:,(lats<=20)&(lats>=-10),:]
sst_c = sst_c[:,:,(lons>=125)&(lons<=145)]
emic = np.zeros(sst_c.shape[0])
for i in xrange(emib.shape[0]): emic[i] = np.mean(sst_c[i])
# EMI
emi = emia - 0.5*emib - 0.5*emic
# Normaisation from Ashok. Do not use.
#emi2 = emi[(dates2.year>=1980)&(dates2.year==2005)]
#emi = emi/0.54

exceed = emi>1.5
for xx in dates[exceed]:
    print xx

# Plot EMI
#fig, ax = plt.subplots()
##rnge = 528
#rnge = 1248
#ax.plot_date(dates[-rnge:], emi[-rnge:], '-')
#ax.plot_date(dates[-rnge:], np.zeros(rnge), 'k')
#ax.plot_date(dates[-rnge:], np.ones(rnge)*emi.std(), 'r--')
#ax.plot_date(dates[-rnge:], -np.ones(rnge)*emi.std(), 'g--')
#ax.xaxis.set_major_locator(YearLocator(base=5))
#ax.xaxis.set_minor_locator(YearLocator(base=1))
##ax.xaxis.set_major_formatter(DateFormatter('%Y'))
#ax.autoscale_view()
#ax.grid(True)
#ax.xaxis.grid(True, which="minor")
#ax.set_xlabel('Date')
#ax.set_ylabel('EMI')
##fig.autofmt_xdate()
#plt.show()

# This just constructs the anomalies for visual inspection of each year. 
#climnc = nc.Dataset('/media/Jupiter/observations/HadISST/HadISST_sst_1971-2000_clim_N96_filled.nc')
#datanc = nc.Dataset('/media/Jupiter/observations/HadISST/HadISST_sst_data_N96_filled.nc')
#clim_sst_low = climnc.variables['sst'][:]
#sst_low = datanc.variables['sst'][:]
#month = 0
#for imonth in xrange(sst_low.shape[0]):
#    sst_low[imonth,...] = sst_low[imonth,...] - clim_sst_low[month,...]
#    if month<11: month += 1
#    else: month = 0
#outnc = nc.Dataset('/media/Jupiter/observations/HadISST/HadISST_sst_data_N96_filled_anomaly.nc','r+')
#sstout = outnc.variables['sst']
#sstout[:] = sst_low[:]
#outnc.close()

# Define the years that have modoki summers. The given year is the year containing January.
# 1988 is questionable and not in Ashok et al. (2007).
#modoki_years = [1973, 1978, 1980, 1987, 1988, 1991, 1992, 1993, 1995, 2003, 2005, 2010]
modoki_years = [1980, 1987, 1988, 1991, 1992, 1993, 1995, 2003, 2005, 2010]
nmyears = len(modoki_years)
modoki_composite = np.zeros((24,)+sst_absolute.shape[1:])
for year in modoki_years:
    sst_iyear = sst[(dates2.year==year)|(dates2.year==(year-1))]
    modoki_composite += sst_iyear
modoki_composite /= nmyears
modoki_composite[modoki_composite<-9] = np.nan
modoki_composite[modoki_composite>9] = np.nan

del sst
del sst_absolute
# Plot a DJF composite

i, j = np.where(lons==-165.5)[0][0], np.where(lats==0.5)[0][0]
mdk_series = modoki_composite[:,j,i]
i, j = 189, 96
nno_series = nino[:,j,i]

fig, axes = plt.subplots(nrows=3,figsize=(7,10))

modokipic = modoki_composite[11:14,...].mean(axis=0)
modokipic[modokipic < -5.] = np.nan
m = Basemap(ax=axes[0],projection='robin', lon_0=180.)
m.drawcoastlines()
parallels = np.arange(-90., 120., 30.)
meridians = np.arange(0., 360., 30.)
par = m.drawparallels(parallels,labels=[True,False,False,False], linewidth=0.3)
for pr in par:
    try:
        par[pr][1][0].set_fontsize(9)
    except: pass
meridians = m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=0.3)
for mr in meridians:
    try: 
        meridians[mr][1][0].set_rotation(50)
        meridians[mr][1][0].set_fontsize(9)
    except: pass
modokipic, lons = bm.shiftgrid(1, modokipic, lons)
xx,yy = np.meshgrid(lons, lats)
x,y = m(xx,yy)
fcont = m.contourf(x,y,modokipic,cmap='seismic',levels=np.arange(-1.6,1.8,0.2))
#m.contour(x,y,modokipic,colors='k',levels=[0])
m.fillcontinents(color='gray',lake_color='gray')
#cbar = plt.colorbar(fcont, orientation='horizontal')
cbar = m.colorbar(fcont, location='right', pad=0.3)
#cbar.ax.get_xaxis().set_ticks([])
#for j, lab in enumerate([str(i) for i in np.arange(-1.5,1.7,0.2)]):
#    cbar.ax.text((1/15.)*j, -0.6, lab, ha='center', va='center')
plt.title('a)', loc='left')

m = Basemap(ax=axes[1], projection='robin', lon_0=180.)
m.drawcoastlines()
par = m.drawparallels(parallels,labels=[True,False,False,False], linewidth=0.3)
for pr in par:
    try:
        par[pr][1][0].set_fontsize(9)
    except: pass
meridians = m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=0.3)
for mr in meridians:
    try: 
        meridians[mr][1][0].set_rotation(50)
        meridians[mr][1][0].set_fontsize(9)
    except: pass
ninopic, lons2 = bm.shiftgrid(1, ninopic, lons2)
xx,yy = np.meshgrid(lons2, lats2)
x,y = m(xx,yy)
fcont = m.contourf(x,y,ninopic,cmap='seismic',levels=np.arange(-2.2,2.3,0.2))
m.fillcontinents(color='gray',lake_color='gray')
cbar = m.colorbar(fcont, location='right', pad=0.3)
plt.title('b)', loc='left')

plt.sca(axes[2])
plt.plot(nno_series,'k-',label='El Nino')
plt.scatter(np.arange(0,24), nno_series, marker='x',color='k')
plt.plot(mdk_series,'k',label='Modoki', linestyle='dotted')
plt.scatter(np.arange(0,24), mdk_series, marker='x',color='k')
labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
labels.extend(labels)
plt.xticks(np.arange(0,24), labels, rotation=45)
plt.axhline(y=0,color='k')
plt.xlim(0, 24)
plt.grid(True)
plt.xlabel('Month')
plt.ylabel('SST $^{\circ}$C')
plt.legend(loc='upper right')
plt.title('c)', loc='left')


plt.savefig('grl_modoki_forcing.eps',format='eps')



# Define years that have modoki la nina summers. make composite
#modoki_nina_years = [2012, 2011, 2008, 2006, 2001, 2000, 1999, 1989, 1976, 1974]
#nmayears = len(modoki_nina_years)
#modoki_nina_composite = np.zeros((24,)+sst_absolute.shape[1:])
#for year in modoki_nina_years:
#    sst_iyear = sst[(dates2.year==year)|(dates2.year==(year-1))]
#    modoki_nina_composite += sst_iyear
#modoki_nina_composite /= nmayears


# Plot modoki nina composite
#f = plt.figure()
#modokipic2 = modoki_nina_composite[11:14,...].mean(axis=0)
#modokipic2[modokipic2 < -5.] = np.nan
#modokipic2[modokipic2 > 5.] = np.nan
#m = Basemap(projection='robin', lon_0=180.)
#m.drawcoastlines()
#parallels = np.arange(-90., 120., 30.)
#meridians = np.arange(0., 360., 30.)
#m.drawparallels(parallels,labels=[True,False,False,False], linewidth=0.3)
#meridians = m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=0.3)
#for mr in meridians:
#    try: meridians[mr][1][0].set_rotation(50)
#    except: pass
#modokipic2, lonsblank = bm.shiftgrid(180, modokipic2, lons)
#xx,yy = np.meshgrid(lons, lats)
#x,y = m(xx,yy)
#fcont = m.contourf(x,y,modokipic2,cmap='seismic',levels=np.arange(-1.5,1.7,0.2),extend='both')
#m.contour(x,y,modokipic2,colors='k',levels=[0])
#m.fillcontinents(color='gray',lake_color='gray')
#cbar = plt.colorbar(fcont, orientation='horizontal')
#cbar.ax.get_xaxis().set_ticks([])
#for j, lab in enumerate([str(i) for i in np.arange(-1.5,1.7,0.2)]):
#    cbar.ax.text((1/15.)*j, -0.5, lab, ha='center', va='center')
#plt.title('DJF SSTA La Nina Modoki Composite')
#plt.show()

# Create composites of the ACCESS compatible low resolution anomalyfile
#low_nc = nc.Dataset('/media/Jupiter/observations/HadISST/HadISST_sst_data_N96_filled_anomaly.nc','r')
#sst_low = low_nc.variables['sst'][:]
#sst_low = sst_low[satera]
#sst_low = signal.detrend(sst_low,axis=0)
#lats = low_nc.variables['lat'][:]
#lons = low_nc.variables['lon'][:]
#sst_iyear = 0
#modoki_composite_low = np.zeros((24,)+sst_low.shape[1:])
#for year in modoki_years:
#    sst_iyear = sst_low[(dates2.year==year)|(dates2.year==(year-1))]
#    modoki_composite_low += sst_iyear
#modoki_composite_low /= nmyears
#modoki_nina_composite_low = np.zeros((24,)+sst_low.shape[1:])
#for year in modoki_nina_years:
#    sst_iyear = sst_low[(dates2.year==year)|(dates2.year==(year-1))]
#    modoki_nina_composite_low += sst_iyear
#modoki_nina_composite_low /= nmayears

## Save to file
#modokinc = nc.Dataset('modokisst_dt.nc','w')
#modokinc.createDimension('time',24)
#modokinc.createDimension('lat',lats.shape[0])
#modokinc.createDimension('lon',lons.shape[0])
#tout = modokinc.createVariable('time','f8',('time'))
#setattr(tout, 'units', 'month')
#tout[:] = range(1,25)
#latout = modokinc.createVariable('lat','f8',('lat'))
#setattr(latout, 'units', 'degrees_north')
#latout[:] = lats
#lonout = modokinc.createVariable('lon','f8',('lon'))
#setattr(lonout, 'units', 'degrees_north')
#lonout[:] = lons
#sstout = modokinc.createVariable('sst','f8',('time','lat','lon'))
#setattr(sstout, 'units', 'degC')
#sstout[:] = modoki_composite_low
##sstout = modokinc.createVariable('sst_modoki_nina','f8',('time','lat','lon'))
##setattr(sstout, 'units', 'degC')
##sstout[:] = modoki_nina_composite_low
#modokinc.close()
