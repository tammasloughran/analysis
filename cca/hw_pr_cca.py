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
def standardize(hwmet):
    """Standardize a heatwave metric by subtracting the long term mean and 
    dividing by the long term standard deviation.

    Arguments
    hwmet -- heatwave metric.

    Returns
    hwmet -- heatwave metric standardized along the first dimension.
    """
    hwmet = hwmet - hwmet.mean(axis=0)
    hwmet = hwmet/hwmet.std(axis=0)
    return hwmet
pr_detrended = detrend(pr, axis=0)
pr_detrended = standardize(pr_detrended)
hwa_detrended = np.ma.array(detrend(hwa, axis=0), mask=hwa.mask)
hwa_detrended = standardize(hwa_detrended)
hwm_detrended = np.ma.array(detrend(hwm, axis=0), mask=hwa.mask)
hwm_detrended = standardize(hwm_detrended)
hwf_detrended = np.ma.array(detrend(hwf, axis=0), mask=hwa.mask)
hwf_detrended = standardize(hwf_detrended)
hwd_detrended = np.ma.array(detrend(hwd, axis=0), mask=hwa.mask)
hwd_detrended = standardize(hwd_detrended)
hwn_detrended = np.ma.array(detrend(hwd, axis=0), mask=hwa.mask)
hwn_detrended = standardize(hwd_detrended)
hwt_detrended = np.ma.array(detrend(hwt, axis=0), mask=hwa.mask)
hwt_detrended = standardize(hwt_detrended)


# Calculate significant number of eofs using North test
def nsigpcs(varfrac, north):
    upper = varfrac + north
    lower = varfrac - north
    upper = upper[1:]
    lower = lower[:-1]
    clear = upper<lower
    return np.argmax(clear==False)

if __name__=='__main__':
    pca = Eof(pr_detrended)
    pr_sig = nsigpcs(pca.varianceFraction(), pca.northTest(vfscaled=True))
    metric_dict = {"HWN":hwn,"HWF":hwf,"HWD":hwd,"HWA":hwa,"HWM":hwm,"HWT":hwt}
    for metric_name in metric_dict.keys():
        pca = Eof(metric_dict[metric_name])
        hw_sig = nsigpcs(pca.varianceFraction(), pca.northTest(vfscaled=True))
        if hw_sig<2: hw_sig = 2
        # Perform CCA
        CCA = BPCCA(metric_dict[metric_name][:-3,...], pr_detrended, (hw_sig,4))
        L = CCA.leftPatterns()
        R = CCA.rightPatterns()
        a = CCA.rightExpCoeffs()
        b = CCA.leftExpCoeffs()
    
        #lmask, rmask = CCA.MCTestMask(250,100,0.05)
        # Print corelations
        print metric_name
        print CCA.correlation()
        # Variance fraction
        print np.array(CCA.varianceFractions())
        # Plot
        for pattern in range(2):
            # Left paterns
            cca_plot.plot_cp(L[...,pattern], hw_lons, hw_lats, 
                    '%d %s'%(pattern+1,metric_name))
            #        lmask[...,pattern])
            # Right paterns
            cca_plot.plot_cp(R[...,pattern], pr_lons, pr_lats,
                    '%d Pressure %s'%(pattern+1,metric_name))
            #        rmask[...,pattern])
            # Plot coeficients
            cca_plot.plot_coefs(a[:,pattern], b[:,pattern], years,
                    llabel='%d %s'%(pattern+1,metric_name), rlabel='Pressure')
