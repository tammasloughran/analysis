"""Calculate the EOFs ond PCs of EHF heatwave metrics.

EHF heat wave metrics have been calculated from AWAP Tmax and Tmin data using a
90th percentile with respect to the base period 1961-1990. Metrics are annual
heat wave characteristics. They are HWF (frequency), HWD (duration), HWA 
(amplitude), HWM (magnitude), HWN (number) and HWT (timing).
"""
import __init__
import numpy as np
import scipy.stats as stats
from scipy.signal import detrend
from pandas import concat
from netCDF4 import Dataset
from eofs.standard import Eof
from tools import rotate
from tools import pcaplot
from tools import load_data

def detrend_kendal(data):
    from scipy.stats.mstats import theilslopes
    import numpy as np
    slope = np.ones(data.shape[1:])*np.nan
    intercept = np.ones(data.shape[1:])*np.nan
    import pdb
    for y in xrange(0,data.shape[1],1):
        for x in xrange(0,data.shape[2],1):
            if not data.mask[:,y,x].any():
                slope[y,x], intercept[y,x], _, _ = theilslopes(data.data[:,y,x])
    trend = np.ones(data.shape)*np.nan
    for t in xrange(0,data.shape[0],1):
        trend[t,...] = slope*t # + intercept
    detrended = data - trend
    return detrended

if __name__ == '__main__':
    # Load the heat wave metrics.
    directory = ('/srv/ccrc/data35/z5032520/')
    fname = (directory+'AWAP/yearly/ehfhw/EHF_heatwaves____yearly_summer.nc')
    hwf, hwn, hwd, hwa, hwm, hwt, lat, lon, times\
            = load_data.load_heat_waves(fname)
    start_year = 1911
    end_year = 2012
    year_range = np.arange(1911, 2014+1, 1)
    period = np.where((year_range>=start_year)&(year_range<=end_year))[0]
    hwf = hwf[period,...]
    hwn = hwn[period,...]
    hwd = hwd[period,...]
    hwa = hwa[period,...]
    hwm = hwm[period,...]
    hwt = hwt[period,...]

    # Load variability indices.
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
    strh = load_data.load_index2(directory+fname, standardize=True)

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

    # Perform PCA on all metrics.
    def nsigpcs(varfrac, north):
        upper = varfrac + north
        lower = varfrac - north
        upper = upper[1:]
        lower = lower[:-1]
        clear = upper<lower
        return np.argmax(clear==False)

    # Calculate weightings.
    coslat = np.cos(np.deg2rad(lat)).clip(0.,1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]
    hwn = detrend_kendal(hwn)
    hwf = detrend_kendal(hwf)
    hwd = detrend_kendal(hwd)
    hwa = detrend_kendal(hwa)
    hwm = detrend_kendal(hwm)
    hwt = detrend_kendal(hwt)
    metric_dict = {"HWN":hwn,"HWF":hwf,"HWD":hwd,"HWA":hwa,"HWM":hwm,"HWT":hwt}
    print "HW?    pc1        pc12        pc3        pc4"
    for metric_name in metric_dict.keys():
        # Set up solver
        solver = Eof(metric_dict[metric_name], weights=wgts)
        explained_variance = solver.varianceFraction()
        total_variance = solver.totalAnomalyVariance()
        errors = solver.northTest(vfscaled=True)
        retain = nsigpcs(explained_variance, errors)
        if retain==1: retain = 4
        pcs = solver.pcs(pcscaling=0, npcs=retain)
        eigens = solver.eigenvalues()
        eofs = solver.eofs(eofscaling=2, neofs=retain)
        eofs_covariance = solver.eofsAsCovariance(pcscaling=1, neofs=retain)
        eofs_correlation = solver.eofsAsCorrelation(neofs=retain)

        np.save(metric_name+'_simple_pcs', pcs)

        # Apply rotation to EOFs.
        eofs2 = eofs
        pcs, eofs, R = rotate.do_rotation(pcs, eofs)

        # Get explained variances of rotated EOFs
        expvar = rotate.expvar_from_eofs(eofs,total_variance)
        print metric_name+' ', expvar

        # Scale the rotated PCs by dividing by the square root of the variances.
        pcs = pcs/np.sqrt(expvar*total_variance)

        np.save(metric_name+'_rotated_pcs', pcs)

        #pcaplot.eofscatter(eofs2)
        #pcaplot.eofscatter(eofs)

        # Plotting.
        years = np.arange(start_year,end_year+1)
        pcaplot.plot_eigenvalues(explained_variance, errors, metric_name)
        pcaplot.plot_eofs(eofs, lon, lat, "%s_Rotated_EOFs"%(metric_name), head='VARIMAX EOFs')
        #pcaplot.plot_eofs(eofs2, lon, lat, "%s_EOFs"%(metric_name), head='Simple EOFs')
        pcaplot.plot_three_eofs(eofs, expvar, lon, lat, metric_name, head='Rotated EOFs')
        #pcaplot.plot_eofs(eofs_covariance, lon, lat, 
        #        "%s_EOFs_Covariance"%(metric_name))
        #pcaplot.plot_eofs(eofs_correlation, lon, lat, 
        #        "%s_EOFs_Correlation"%(metric_name))
        pcaplot.plot_pcs(pcs, ninoslice, year_range[period], metric_name, head='VARIMAX ')

        # Correlations for Summer/Winter
        outfile = open("%s_correlations"%(metric_name),'w')
        outfile.write("      Nino3.4       SOI          DMI           SAM           STRH\n")
        for pc in range(pcs.shape[1]):
            outfile.write("PC%0.f: "%(pc+1))
            for mode in [ninoslice, soislice, dmislice, samslice, strhslice]:
                rho, p = stats.spearmanr(mode, pcs[:-1,pc])
                outfile.write("%+.2f (%.3f) "%(rho, p))
            outfile.write("\n")
        outfile.close()

        # Lag Correlations
        mds = ["Nino 3.4","SOI","DMI","SAM","STRH"]
        mdsn = 0
        for mode in [nino34, soi, dmi, sam, strh]:
            mode = mode['%s-01'%(start_year-2):'%s-12'%(end_year)]
            rho_lag = np.zeros((25, pcs.shape[1]))
            p_lag = np.zeros((25, pcs.shape[1]))
            for lag in range(0,25,1):
                mode_lag = mode.shift(lag)
                axis = mode_lag.index.month
                mode_lag = mode_lag[(axis==12)|(axis==1)|(axis==2)]
                mode_lag = mode_lag.resample('AS-JUL', how='mean')
                mode_lag = mode_lag['%s-01'%(start_year):'%s-12'%(end_year-1)]
                for pc in range(pcs.shape[1]):
                    rho_lag[lag,pc], p_lag[lag,pc] = \
                        stats.mstats.spearmanr(mode_lag, pcs[:-1,pc])
            pcaplot.plot_lags(rho_lag, p_lag, metric_name+" PCs & %s"%(mds[mdsn]))
            mdsn+=1
