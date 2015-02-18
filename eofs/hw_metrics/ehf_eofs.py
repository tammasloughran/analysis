"""Calculate the EOFs ond PCs of EHF heatwave metrics.

EHF heat wave metrics have been calculated from AWAP Tmax and Tmin data using a
90th percentile with respect to the base period 1961-1990. Metrics are annual
heat wave characteristics. They are HWF (frequency), HWD (duration), HWA 
(amplitude), HWM (magnitude), HWN (number) and HWT (timing).
"""

if __name__ == '__main__':
    from numpy import arange, cos, sqrt, deg2rad, newaxis, dot
    from pandas import concat
    from eofs.standard import Eof
    import rotate
    import pcaplot
    from netCDF4 import Dataset
    import load_data
    import scipy.stats as stats

    # Load the heat wave metrics.
    fname = ('/srv/ccrc/data35/z5032520/AWAP/yearly/ehfhw/CCRC_NARCliM_1911-'
             '2014_EHFheatwaves_summer_AWAP0.5deg_detrended.nc')
    hwf, hwn, hwd, hwa, hwm, hwt, lat, lon, times\
            = load_data.load_heat_waves(fname)
    start_year = 1911
    end_year = 2014

    # Load variability indices.
    nino34 = load_data.load_index('../../indices/NINO_3.4_monthly_index.csv',
                                    standardize=False)
    dmi = load_data.load_index('../../indices/DMI_monthly_index.csv', 
                                    standardize=False)
    soi = load_data.load_index2('../../indices/SOI_UCAR_stdzd_1886_2013.txt', 
                                    standardize=False)
    sam1 = load_data.load_index2('../../indices/SAM_index_monthly_visbeck.txt', 
                                    standardize=False)
    sam2 = load_data.load_index2('../../indices/SAM_1957_2013.txt', 
                                    standardize=False)
    sam = concat([sam1['1911-07':'1956-12'],sam2['1957-01':'2013-06']])

    # Take spring (9,10,11)/summer (11,12,1,2,3) (A)nnual (mean) (S)tarting in (JUL)y. i.e. 'AS-JUL'i
    mnth = sam.index.month
    # SAM
    sam = sam[(mnth==11)|(mnth==12)|(mnth==1)|(mnth==2)|(mnth==3)]
    sam = sam.resample('AS-JUL', how='mean')
    # Nino3.4
    ninoslice = nino34['1911-07':'2013-06']
    ninoslice = ninoslice[(mnth==11)|(mnth==12)|(mnth==1)|(mnth==2)|(mnth==3)]
    ninoslice = ninoslice.resample('AS-JUL', how='mean')
    # DMI
    dmislice = dmi['1911-07':'2013-06']
    dmislice = dmislice[(mnth==11)|(mnth==12)|(mnth==1)|(mnth==2)|(mnth==3)]
    dmislice = dmislice.resample('AS-JUL', how='mean')
    # SOI
    soislice = soi['1911-07':'2013-06']
    soislice = soislice[(mnth==11)|(mnth==12)|(mnth==1)|(mnth==2)|(mnth==3)]
    soislice = soislice.resample('AS-JUL', how='mean')

    # Perform PCA on all metrics.
    # Calculate weightings.
    coslat = cos(deg2rad(lat)).clip(0.,1.)
    wgts = sqrt(coslat)[..., newaxis]
    metric_dict = {"HWF":hwf,"HWN":hwn,"HWD":hwd,"HWA":hwa,"HWM":hwm,"HWT":hwt}
    for metric_name in metric_dict.keys():
        # Set up solver
        retain = 4
        solver = Eof(metric_dict[metric_name], weights=wgts)
        pcs = solver.pcs(pcscaling=1, npcs=retain)
        explained_variance = solver.varianceFraction()
        errors = solver.northTest(vfscaled=True)
        eigens = solver.eigenvalues()
        eofs = solver.eofs(eofscaling=2, neofs=retain)
        eofs_covariance = solver.eofsAsCovariance(pcscaling=1, neofs=retain)
        eofs_correlation = solver.eofsAsCorrelation(neofs=retain)

        # Apply rotation to PCs and EOFs.
        eofs2 = eofs
        pcs, eofs = rotate.do_rotation(pcs, eofs, space='state')

        # Plotting.
        years = arange(start_year,end_year+1)
        pcaplot.plot_eigenvalues(explained_variance, errors, metric_name)
        pcaplot.plot_eofs(eofs, lon, lat, "%s_Rotated_EOFs"%(metric_name))
        pcaplot.plot_eofs(eofs2, lon, lat, "%s_EOFs"%(metric_name))
        pcaplot.plot_eofs(eofs_covariance, lon, lat, 
                "%s_EOFs_Covariance"%(metric_name))
        pcaplot.plot_eofs(eofs_correlation, lon, lat, 
                "%s_EOFs_Correlation"%(metric_name))
        pcaplot.plot_pcs(pcs, ninoslice, years, metric_name)

        # Correlations
        outfile = open("%s_correlations"%(metric_name),'w')
        outfile.write("      Nino3.4       SOI          DMI           SAM\n")
        for pc in [0,1,2,3]:
            outfile.write("PC%0.f: "%(pc+1))
            for mode in [ninoslice, soislice, dmislice, sam]:
                rho, p = stats.spearmanr(mode, pcs[:-2,pc])
                outfile.write("%+.2f (%.3f) "%(rho, p))
            outfile.write("\n")
        outfile.close()
