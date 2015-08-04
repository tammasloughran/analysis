import numpy as np
import load_data
import cca_plot
from netCDF4 import Dataset
from pyclimate.bpcca import BPCCA
from scipy.signal import detrend
from eofs.standard import Eof

# Load variables
prfile = ('/home/nfs/z5032520/heatwaves-svn/dataproc/summer_mslp_1911-2011.nc')
hwfile = ('/media/Jupiter/reanalysis/AWAP/yearly/ehfhw/'
        'CCRC_NARCliM_1911-2014_EHFheatwaves_summer_AWAP0.5deg_detrended.nc')
hwf, hwn, hwd, hwa, hwm, hwt, hw_lats, hw_lons, times\
        = load_data.load_heat_waves(hwfile)
prnc = Dataset(prfile,'r')
pr = prnc.variables['mslp'][:]
pr_lats = prnc.variables['lat'][:]
pr_lons = prnc.variables['lon'][:]
years = range(1911,2013)

# Detrend and standardize.
pr_detrended = detrend(pr, axis=0)
pr_detrended = pr_detrended - pr_detrended.mean(axis=0)
pr_detrended = pr_detrended/pr_detrended.std(axis=0)
hwa_detrended = np.ma.array(detrend(hwa, axis=0), mask=hwa.mask)
hwa_detrended = hwa_detrended - hwa_detrended.mean(axis=0)
hwa_detrended = hwa_detrended/hwa_detrended.std(axis=0)

# Calculate significant number of eofs using North test
def nsigpcs(varfrac, north):
    upper = varfrac + north
    lower = varfrac - north
    upper = upper[1:]
    lower = lower[:-1]
    clear = upper<lower
    return np.argmax(clear==False)

pca = Eof(pr_detrended)
pr_sig = nsigpcs(pca.varianceFraction(), pca.northTest(vfscaled=True))
#import pdb
#pdb.set_trace()
pca = Eof(hwa_detrended)
hw_sig = nsigpcs(pca.varianceFraction(), pca.northTest(vfscaled=True))
print hw_sig, pr_sig
# Perform CCA
CCA = BPCCA(hwa_detrended[:-3,...], pr_detrended, (hw_sig,pr_sig))
L = CCA.leftPatterns()
R = CCA.rightPatterns()
a = CCA.rightExpCoeffs()
b = CCA.leftExpCoeffs()

if __name__=='__main__':
    lmask, rmask = CCA.MCTestMask(250,100,0.05)
    # Print corelations
    print CCA.correlation()
    # Variance fraction
    #print CCA.varianceFractions()
    # Plot
    pattern = 0
    # Left paterns
    cca_plot.plot_cp(L[...,pattern], hw_lons, hw_lats, 'HWA', 
            lmask[...,pattern])
    # Right paterns
    cca_plot.plot_cp(R[...,pattern], pr_lons, pr_lats, 'Pressure', 
            rmask[...,pattern])
    # Plot coeficients
    cca_plot.plot_coefs(a[:,pattern], b[:,pattern], years,
            llabel='HWA', rlabel='Pressure')
