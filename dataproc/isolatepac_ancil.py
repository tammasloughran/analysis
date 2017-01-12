"""isolatepac_ancil.py uses HadISST SSTs to construct a composite of sst anomalies for
ENSO and IOD excluding each basin.
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

# Calculate Nino34 index
nino34_sst = sst[:,(lats<=5)&(lats>=-5),:]
nino34_sst = nino34_sst[:,:,(lons>=-170)&(lons<=-120)]
nino34 = np.zeros(nino34_sst.shape[0])
wght = lats[(lats<=5)&(lats>=-5)]
wght = np.broadcast_to(wght, (50,10)).T
wght = np.broadcast_to(wght, (540,10,50))
wght = np.cos(np.deg2rad(wght))
nino34_sst = nino34_sst*wght
for i in xrange(nino34_sst.shape[0]): 
    nino34[i] = np.sum(nino34_sst[i])/np.sum(wght[0])

# Plot Nino34
fig, ax = plt.subplots()
rnge = 528
ax.plot_date(dates[-rnge:], nino34[-rnge:], '-')
ax.plot_date(dates[-rnge:], np.zeros(rnge), 'k')
ax.plot_date(dates[-rnge:], np.ones(rnge)*nino34.std(), 'r--')
ax.plot_date(dates[-rnge:], -np.ones(rnge)*nino34.std(), 'g--')
ax.xaxis.set_major_locator(YearLocator())
ax.xaxis.set_major_formatter(DateFormatter('%Y'))
ax.autoscale_view()
ax.grid(True)
ax.set_xlabel('Date')
ax.set_ylabel('Nino 3.4 ($^\circ$C)')
fig.autofmt_xdate()
#plt.show()

# Define the years that have nonlinear lanina summers. The given year is the year containing January.
nlln_years = [2014, 2011, 2008, 2004, 1989, 1986, 1976, 1974]
nnllnyears = len(nlln_years)
#nlln_composite = np.zeros((24,)+sst_absolute.shape[1:])
#for year in nlln_years:
#    sst_iyear = sst[(dates2.year==year)|(dates2.year==(year-1))]
#    nlln_composite += sst_iyear
#nlln_composite /= nnllnyears

# Plot a JJA composite
#nllnpic = nlln_composite[5:8,...].mean(axis=0)
#nllnpic[nllnpic < -5.] = np.nan
#m = Basemap(projection='robin', lon_0=180.)
#m.drawcoastlines()
#parallels = np.arange(-90., 120., 30.)
#meridians = np.arange(0., 360., 30.)
#m.drawparallels(parallels,labels=[True,False,False,False], linewidth=0.3)
#meridians = m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=0.3)
#for mr in meridians:
#    try: meridians[mr][1][0].set_rotation(50)
#    except: pass
#nllnpic, lons = bm.shiftgrid(1, nllnpic, lons)
#xx,yy = np.meshgrid(lons, lats)
#x,y = m(xx,yy)
#fcont = m.contourf(x,y,nllnpic,cmap='seismic',levels=np.arange(-1.5,1.7,0.2))
#m.contour(x,y,nllnpic,colors='k',levels=[0])
#m.fillcontinents(color='gray',lake_color='gray')
#cbar = plt.colorbar(fcont, orientation='horizontal')
#cbar.ax.get_xaxis().set_ticks([])
#for j, lab in enumerate([str(i) for i in np.arange(-1.5,1.7,0.2)]):
#    cbar.ax.text((1/15.)*j, -0.5, lab, ha='center', va='center')
#plt.title('JJA SSTA Eastern Pacific La Nina Composite')
#plt.show()

# Load the anomalies
ensonc = nc.Dataset('/srv/ccrc/data35/z5032520/ancilforge/modal_anomalies.nc','r')
elnino = np.array(ensonc.variables['elnino_sst'][:])
elnino[elnino>100] = 0
lanina = np.array(ensonc.variables['lanina_sst'][:])
lanina[lanina>100] = 0
lats = ensonc.variables['lat'][:]
lons = ensonc.variables['lon'][:]



pac_mask = np.ones(elnino.shape)
# zero anomalies outside pacific +-30N of equator
pac_mask[:,lats>=30,:] = 0
pac_mask[:,lats<=-30,:] = 0
pac_mask[:,:,lons<=120] = 0
pac_mask[:,:,lons>=290] = 0
# Linear damping north south and westward
for i,j in enumerate(np.where((lats<30)&(lats>=20))[0]):
    n = float(np.sum((lats<30)&(lats>=20)))
    factor = (n-i)/n
    pac_mask[:,j,(lons>120)&(lons<290)] = pac_mask[:,j,(lons>120)&(lons<290)]*factor
for i,j in enumerate(np.where((lats>-30)&(lats<=-20))[0]):
    n = float(np.sum((lats<30)&(lats>=20)))
    factor = 1+(1./n) - (n-i)/n
    pac_mask[:,j,(lons>120)&(lons<290)] = pac_mask[:,j,(lons>120)&(lons<290)]*factor
for ii,i in enumerate(np.where((lons>120)&(lons<=140))[0]):
    n = float(np.sum((lons>120)&(lons<=140)))
    factor = 1+(1./n) - (n-ii)/n
    pac_mask[:,(lats>-30)&(lats<30),i] = pac_mask[:,(lats>-30)&(lats<30),i]*factor
pac_mask[:,114:,208:] = 0 # remove caribean 
pac_mask[:,105:,221:] = 0 # remove caribean

# Indo-pacific region
indopac_mask = np.ones(elnino.shape)
# zero anomalies outside pacific +-30N of equator
indopac_mask[:,lats>=30,:] = 0
indopac_mask[:,lats<=-30,:] = 0
indopac_mask[:,:,lons<=30] = 0
indopac_mask[:,:,lons>=290] = 0
# Linear damping north south and westward
for i,j in enumerate(np.where((lats<30)&(lats>=20))[0]):
    n = float(np.sum((lats<30)&(lats>=20)))
    factor = (n-i)/n
    indopac_mask[:,j,(lons>30)&(lons<290)] = indopac_mask[:,j,(lons>30)&(lons<290)]*factor
for i,j in enumerate(np.where((lats>-30)&(lats<=-20))[0]):
    n = float(np.sum((lats<30)&(lats>=20)))
    factor = 1+(1./n) - (n-i)/n
    indopac_mask[:,j,(lons>30)&(lons<290)] = indopac_mask[:,j,(lons>30)&(lons<290)]*factor
indopac_mask[:,114:,208:] = 0 # remove caribean 
indopac_mask[:,105:,221:] = 0 # remove caribean

indopac_nino = elnino*indopac_mask
indopac_nina = lanina*indopac_mask

elnino = elnino*pac_mask
# There are still positive anomalies in the far western pacific 
# that need to be removed.
for t in xrange(0,17):
    for y in xrange(0,150):
        for x in xrange(0,120):
            if elnino[t,y,x]>0.: elnino[t,y,x] = 0
lanina = lanina*pac_mask
# And remove negative anomalies in far western Pacific.
for t in xrange(0,17):
    for y in xrange(0,150):
        for x in xrange(0,120):
            if lanina[t,y,x]<0.: lanina[t,y,x] = 0

# Isolated indian ocean ssts
piod = np.array(ensonc.variables['piod_sst'][:])
piod[piod>100] = 0
niod = np.array(ensonc.variables['niod_sst'][:])
niod[niod>100] = 0
ind_mask = np.ones(elnino.shape)
# zero anomalies outside indian ocean +-30N of equator
ind_mask[:,lats>=30,:] = 0
ind_mask[:,lats<=-30,:] = 0
ind_mask[:,:,lons<=30] = 0
ind_mask[:,:,lons>=130] = 0
# Linear damping north south and eastward
for i,j in enumerate(np.where((lats<30)&(lats>=20))[0]):
    n = float(np.sum((lats<30)&(lats>=20)))
    factor = (n-i)/n
    ind_mask[:,j,(lons>30)&(lons<130)] = ind_mask[:,j,(lons>30)&(lons<130)]*factor
for i,j in enumerate(np.where((lats>-30)&(lats<=-20))[0]):
    n = float(np.sum((lats<30)&(lats>=20)))
    factor = 1+(1./n) - (n-i)/n
    ind_mask[:,j,(lons>30)&(lons<130)] = ind_mask[:,j,(lons>30)&(lons<130)]*factor
for ii,i in enumerate(np.where((lons>110)&(lons<=130))[0]):
    n = float(np.sum((lons>110)&(lons<=130)))
    factor = (n-ii)/n
    ind_mask[:,(lats>-30)&(lats<30),i] = ind_mask[:,(lats>-30)&(lats<30),i]*factor
ind_mask[:,91:,83:] = 0 # Remove south China sea
ind_mask[:,93:,81:] = 0 # South China sea
ind_mask[:,89:,85:] = 0 # South China
piod = piod*ind_mask
niod = niod*ind_mask

# Save to file
isonc = nc.Dataset('isolated_anomalies.nc','w')
isonc.createDimension('time',24)
isonc.createDimension('lat',lats.shape[0])
isonc.createDimension('lon',lons.shape[0])
tout = isonc.createVariable('time','f8',('time'))
setattr(tout, 'units', 'month')
tout[:] = range(1,25)
latout = isonc.createVariable('lat','f8',('lat'))
setattr(latout, 'units', 'degrees_north')
latout[:] = lats
lonout = isonc.createVariable('lon','f8',('lon'))
setattr(lonout, 'units', 'degrees_east')
lonout[:] = lons
elninoout = isonc.createVariable('elnino_sst','f8',('time','lat','lon'))
setattr(elninoout, 'units', 'degC')
elninoout[:] = elnino
laninaout = isonc.createVariable('lanina_sst','f8',('time','lat','lon'))
setattr(laninaout, 'units', 'degC')
laninaout[:] = lanina
piodout = isonc.createVariable('piod_sst','f8',('time','lat','lon'))
setattr(piodout, 'units', 'degC')
piodout[:] = piod
niodout = isonc.createVariable('niod_sst','f8',('time','lat','lon'))
setattr(niodout, 'units', 'degC')
niodout[:] = niod
ipelout = isonc.createVariable('indopac_nino','f8',('time','lat','lon'))
setattr(ipelout, 'units', 'degC')
ipelout[:] = indopac_nino
iplaout = isonc.createVariable('indopac_nina','f8',('time','lat','lon'))
setattr(iplaout, 'units', 'degC')
iplaout[:] = indopac_nina
isonc.close()
