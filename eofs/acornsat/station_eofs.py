import __init__
import numpy as np
from pandas import concat
from netCDF4 import Dataset
from eofs.standard import Eof
from tools import load_data
from tools import rotate
from tools import pcaplot
from scipy.signal import detrend
import scipy.stats as stats

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

def nsigpcs(varfrac, north):
    """Calculate the number of significant PCs using the north test.

    Arguments
    varfrac -- eigenvalues expressed as a fraction of the total variance
    north -- errors from a north test scaled by the variance fraction.

    Returns
    retain -- the number of significant pcs to retain. 
    """
    upper = varfrac + north
    lower = varfrac - north
    upper = upper[1:]
    lower = lower[:-1]
    clear = upper<lower
    return np.argmax(clear==False)

if __name__ == '__main__':
    # Load the data
    directory = '/srv/ccrc/data35/z5032520/ACORNSAT/yearly/EHF/'
    filename = 'EHF_heatwaves_ACORNSAT_1911_2014_summer.nc'
    ncin = Dataset(directory+filename, 'r')
    hwf = ncin.variables['HWF'][:]
    hwn = ncin.variables['HWN'][:]
    hwd = ncin.variables['HWD'][:]
    hwm = ncin.variables['HWM'][:]
    hwa = ncin.variables['HWA'][:]
    hwt = ncin.variables['HWT'][:]
    lat = ncin.variables['lat'][:]
    lon = ncin.variables['lon'][:]
    stationids = ncin.variables['stationid'][:]

    # Select years for which mode indices exist
    start_year = 1911
    end_year = 2012
    year_range = np.arange(start_year, end_year+1, 1)
    period = np.where((year_range>=start_year)&(year_range<=end_year))[0]
    hwf = hwf[period,...]
    hwn = hwn[period,...]
    hwd = hwd[period,...]
    hwa = hwa[period,...]
    hwm = hwm[period,...]
    hwt = hwt[period,...]

    # Quality control, detrend, standardize
    good_stations = np.load('../../dataproc/qqtest.npy')
    hwf = detrend(hwf[:,good_stations], axis=0)
    hwf = standardize(hwf)
    hwn = detrend(hwn[:,good_stations], axis=0)
    hwn = standardize(hwn)
    hwd = detrend(hwd[:,good_stations], axis=0)
    hwd = standardize(hwd)
    #hwm = detrend(hwm[:,good_stations], axis=0)
    #hwm = standardize(hwm)
    #hwa = detrend(hwa[:,good_stations], axis=0)
    #hwa = standardize(hwa)
    hwt = detrend(hwt[:,good_stations], axis=0)
    hwt = standardize(hwt)
    lat = lat[good_stations]
    lon = lon[good_stations]
    stationids = stationids[good_stations]

    # Load variability indices.
    directory = '/srv/ccrc/data35/z5032520/'
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
    strh = load_data.load_index2(directory+fname)

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

    # Perform PCA
    metric_dict = {"HWN":hwn,"HWF":hwf,"HWD":hwd,"HWT":hwt}
    for metric_name in metric_dict.keys():
        solver = Eof(metric_dict[metric_name])
        explained_variance = solver.varianceFraction()
        errors = solver.northTest(vfscaled=True)
        eigens = solver.eigenvalues()
        retain = nsigpcs(explained_variance, errors)
        if retain<2:
            retain = 2
        pcs = solver.pcs(pcscaling=1, npcs=retain)
        eofs = solver.eofs(eofscaling=2, neofs=retain)
        #eofs_covariance = solver.eofsAsCovariance(pcscaling=1, neofs=retain)
        #eofs_correlation = solver.eofsAsCorrelation(neofs=retain)

        # Rotate. No need to reshape. 
        eofs_r, R = rotate.varimax(eofs.T, normalize=False)
        pcs_r = np.dot(pcs, R)

        # Corelations for summer
        outfile = open("%s_correlations"%(metric_name),'w')
        outfile.write("      Nino3.4       SOI          DMI           SAM           STRH\n")
        for pc in range(pcs.shape[1]):
            outfile.write("PC%0.f: "%(pc+1))
            for mode in [ninoslice, soislice, dmislice, samslice, strhslice]:
                rho, p = stats.spearmanr(mode, pcs[:-1,pc])
                outfile.write("%+.2f (%.3f) "%(rho, p))
            outfile.write("\n")
        outfile.close()

        # Make picures
        pcaplot.plot_eigenvalues(explained_variance, errors, metric_name)
        pcaplot.station_scatter(eofs[:,0],lat,lon,
            metric_name+'_Simple_EOF1',reverse=False)
        pcaplot.station_scatter(eofs[:,1],lat,lon,
            metric_name+'_Simple_EOF2',reverse=False)
        pcaplot.station_scatter(eofs_r[:,0],lat,lon,
            metric_name+'_Rotated_EOF1',reverse=False)
        pcaplot.station_scatter(eofs_r[:,1],lat,lon,
            metric_name+'_Rotated_EOF2',reverse=False)
        pcaplot.plot_pcs(pcs_r, ninoslice, year_range, 
            metric_name, yearmean=False, head='')
