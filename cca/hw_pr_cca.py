import numpy as np
import load_data
import cca_plot
from netCDF4 import Dataset
from pyclimate.bpcca import BPCCA
from scipy.signal import detrend

# Load variables
prfile = ('/home/nfs/z5032520/heatwaves-svn/dataproc/summer_mslp_1911-2013.nc')
hwfile = ('/media/Jupiter/reanalysis/AWAP/yearly/ehfhw/'
        'CCRC_NARCliM_1911-2014_EHFheatwaves_summer_AWAP0.5deg_detrended.nc')
hwf, hwn, hwd, hwa, hwm, hwt, hw_lats, hw_lons, times\
        = load_data.load_heat_waves(hwfile)
prnc = Dataset(prfile,'r')
pr = prnc.variables['mslp'][:]
pr_lats = prnc.variables['lat'][:]
pr_lons = prnc.variables['lon'][:]
years = range(1911,2013)

# Detrend
pr_detrended = detrend(pr, axis=0)
pr_detrended = (pr_detrended - pr_detrended.mean(axis=0))/pr_detrended.std()
hwa_detrended = np.ma.array(detrend(hwa, axis=0), mask=hwa.mask)
hwa_detrended = np.ma.array((hwa_detrended - hwa_detrended.mean(axis=0))/hwa_detrended.std(), mask=hwa.mask)

# Perform CCA
CCA = BPCCA(hwa_detrended[:-3,...], pr_detrended, (5,5))
L = CCA.leftPatterns()
R = CCA.rightPatterns()
r = CCA.correlation()
a = CCA.rightExpCoeffs()
b = CCA.leftExpCoeffs()

if __name__=='__main__':
    lmask, rmask = CCA.MCTestMask(500,100,0.05)
    # Print corelations
    print (r)
    # Plot
    pattern = 0
    # Left paterns
    cca_plot.plot_cp(L[...,pattern], hw_lons, hw_lats, 'HWA', lmask[...,pattern])
    # Right paterns
    cca_plot.plot_cp(R[...,pattern], pr_lons, pr_lats, 'Pressure', rmask[...,pattern])
    # Plot coeficients
    cca_plot.plot_coefs(a[:,pattern], b[:,pattern], years,
            llabel='HWA', rlabel='Pressure')
