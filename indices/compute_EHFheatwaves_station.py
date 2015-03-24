#!/usr/bin/env python
"""compute_EHFheatwaves.py 
Author: Daniel Argueso (d.argueso@unsw.edu.au) at CoECSS, UNSW, and CCRC, 
UNSW. Australia. Contains the function to calculate Excess Heat Factor 
indices from tmax, tmin, dates. It also takes nwindow: number of days to 
calculate the percentiles, and thres_file if a previously calculated 
threshold file is available. It gives back 5 characteristics of the heatwaves

Yearly maximum heatwave intensity 
Yearly average heatwave intensity 
Yearly number of heatwave days
Yearly number of heatwaves
Duration of yearly longest heatwave

It also provides: 

90th percentile of mean temperature for each calendar day using a n-day window

Script based on Sarah's matlab script EHF_index.m
To be used with EHF_WRF.py or EHF_AWAP.pu

Created: 09 September 2013
Modified: 15 October 2013
"""
import numpy as np
import datetime as dt
import netCDF4 as nc
from constants import const
import glob as glob
import sys


def calc_percentile(tave, nyears, nwindow=15, modified=True, thres_file=None):
    if modified==False:
        # Original method - 95th percentile over the entire period
        pct_calc = np.percentile(tave, 95, axis=0)
    else:
        # Modified method - 90th percentile for each calendar day using a 
        # n-day window (default 15-day window)
        if thres_file==None:
            # No percentile file is provided and thus they are calculated from
            # the given data
            windowrange=np.zeros((365,), dtype=np.bool)
            windowrange[:np.ceil(nwindow/2)+1] = True
            windowrange[-np.floor(nwindow/2):] = True
            windowrange = np.tile(windowrange, nyears)
            pct_calc = np.ones((365,)+tave.shape[1:], 
                    np.float64)*const.missingval
            for d in xrange(365):
                tave_window = tave[windowrange==True,...]
                for station in range(tave.shape[1]):
                    tave_no_nan = tave_window[:,station]
                    tave_no_nan = tave_no_nan[np.isnan(tave_no_nan)==False]
                    pct_calc[d,station] = np.percentile(tave_no_nan,
                            90,
                            axis=0)
        	windowrange = np.roll(windowrange, 1)
        else:
            # A percentile file is provided and it contains a 
            # PRCTILE90 variable
            pct_file = nc.Dataset(thres_file, 'r')
            pct_calc = pct_file.variables['PRCTILE90'][:].astype('float64')
    return pct_calc


def compute_EHF(
        tmax, tmin, dates=None, modified=True, nwindow=15,
        thres_file=None, season='yearly'):
    """Function to calculate Excess Heat Factor (EHF) heatwaves from tmax
    and tmin. It removes the leapdays and makes all calculations as if lead
    days didn't exist. It requires that the dates include complete years 
    (not portions of them). It gives an error otherwise.

    tmax: daily maximum temperature [degC] - dimensions=[time,lat,lon]
    tmin: daily minimum temperature [degC]
    dates: a list with the dates for which we have data (tmax and tmin) 
        - basically the range
    modified: if true, it uses the window to calculate the 90th percentile, 
        otherwise uses the entire period to calculate the 95th percentile, 
        which was the original definition (no winter warm spells then)
    nwindow: number of days of the window to calculate the 90th 
        percentile (default=15 days)
    thres_file: Provides a file where to retrieve the percentiles from. In 
        case different datasets are to be compared, the same thresholds 
        should be used for comparability. If thres_file=None, then it 
        calculates the percentile from the input data.
    season: takes yearly (leave as is, heat waves are calculated over the 
        entire year) or summer (heat waves are calcualted over the period 
        NOV-MARCH - austral summer)
    ---
    PRCTILE90:90th percentile for each calendar day using n-day 
        window (nwindow)
    HWA_EHF: Peak of the hottest heatwave per year - yearly maximum of each 
        heatwave peak [degC2] 
    HWM_EHF: Average magnitude of the yearly heatwave - yearly average of 
        heatwave magnitude [degC2]
    HWN_EHF: Number of heatwaves per year
    HWF_EHF: Number of heatwave days - expressed as the percentage relative 
        to the total number of days [%]
    HWD_EHF: Duration of the longest heatwave per year [days]
    HWT_EHF: First heat wave day [days from 1st nov]

    Script based on Sarah's matlab script EHF_index.m
    (Perkins and Alexander, 2012 JCLIM; Perkins et al.,2012 GRL)
    """


    months_all = np.asarray([dates[i].month for i in xrange(len(dates))])    
    days_all = np.asarray([dates[i].day for i in xrange(len(dates))])

    if (months_all[0]!=1) or (days_all[0]!=1):
        sys.exit(
                '''The period must start on 1st of January and the dates 
                provided start on %s/%s''' %(days_all[0],months_all[0]))
    elif (months_all[-1]!=12) or (days_all[-1]!=31):
        sys.exit(
                ('''The period must finish on 31st of December and the dates 
                provided end on %s/%s''') %(days_all[-1],months_all[-1]))

    # Removing leap days
    dates = dates[((months_all==2) & (days_all==29))==False]
    tmax = tmax[((months_all==2) & (days_all==29))==False,...]
    tmin = tmin[((months_all==2) & (days_all==29))==False,...]

    years = np.asarray([dates[i].year for i in xrange(len(dates))])
    months = np.asarray([dates[i].month for i in xrange(len(dates))])  
    days = np.asarray([dates[i].day for i in xrange(len(dates))])

    # Taking starting and ending years
    syear = np.min(years)
    eyear = np.max(years)
    nyears = eyear-syear+1

    # Start year: financial year so that summers belong to the same year 
    # (not split in DEC/JAN)
    years[months<7] -= 1
    shift_pct = np.argmax(years==syear) # Because the percentile is 
    # calculated using natural years (1st day corresponds to 1st Jan), 
    # so it must be shifted when comparing with temp.

    # Calculate average temperature 
    if tmax.shape!=tmin.shape:
        sys.exit(
                '''Maximum and minimum temperature arrays do not have the 
                same dimensions. They must have the same dimensions to 
                calculate daily mean temperature''')
    tave = (tmax+tmin)/2.

    # Calculate the percentiles
    pct = calc_percentile(tave, nyears, nwindow, modified, thres_file)

    ## DEFINING ARRAYS HOLDING METRICS
    HWA_EHF_yearly = np.zeros((nyears,)+tave.shape[1:],np.float64)
    HWM_EHF_yearly = np.zeros((nyears,)+tave.shape[1:],np.float64)
    HWN_EHF_yearly = np.zeros((nyears,)+tave.shape[1:],np.float64)
    HWF_EHF_yearly = np.zeros((nyears,)+tave.shape[1:],np.float64)
    HWD_EHF_yearly = np.zeros((nyears,)+tave.shape[1:],np.float64)
    HWT_EHF_yearly = np.zeros((nyears,)+tave.shape[1:],np.float64)

    ###### CALCULATING EHF ##########
    ## Calculating EHI acclimatization (EHI_accl)
    ndays = tave.shape[0]
    nstat = tave.shape[1]
    #nlon = tave.shape[2]

    for year in range(nyears):
        this_year = year+syear
        print "Processing year: %s" %(this_year)

        # Selecting variables for the running year
        tave_y = tave[years==this_year,...].copy()
        ndays_y = tave_y.shape[0]
        months_y = months[years==this_year]
        years_y = years[years==this_year]
        EHIaccl = np.zeros(tave_y.shape,dtype=np.float64)

        for t in xrange(0, ndays_y):
            EHIaccl[t,...] = np.mean(tave_y[t+30:t+33,:,:], axis=0) \
                    -np.mean(tave_y[t:t+30,:,:], axis=0) 

        ###### CALCULATING Significance (EHIsig) ########
        EHIsig = np.zeros(tave_y.shape, dtype=np.float64)
        if modified==False:
            for t in xrange(3, ndays_y):
                EHIsig[t,...] = np.mean(tave_y[t-2:t+1,...], axis=0)-pct
        else:
            for t in xrange(3, ndays_y):
                EHIsig[t,...] = np.mean(tave_y[t-2:t+1,...],axis=0) \
                        -pct[(t+shift_pct)%365,...]

        ##### CALCULATING EHF AND EHF_EXCEED ########
        EHF = np.maximum(1, EHIaccl)*EHIsig
        EHF[EHF<0] = 0
        EHF_exceed = np.zeros(tave_y.shape, dtype=np.bool)
        EHF_exceed[EHF>0] = True

        ###### ZEROING DAYS NOT BELONGING TO SUMMER considered as 
        ###### NOV,DEC,JAN,FEB,MAR
        if season=='summer':
            EHF_exceed[(months_y>=4) & (months_y<=10),...] = False
            years_y[(months_y>=4) & (months_y<=10)] =- 99
            shift_start_year = (dt.datetime(syear, 11, 01) \
                    -dt.datetime(syear, 07, 01)).days
        elif season=='yearly':
            shift_start_year = 0

        ########### IDENTIFYING HEAT WAVES, 
        # AVERAGE MAGNITUDE AND PEAK ##############
        spell = np.zeros(tave_y.shape, dtype=np.int)

        # If the period starts with a heatwave, set that in the 
        # spell counter variable
        spell[0,EHF_exceed[0,...]==True] = 1

        # Calculate the spells for the rest of the period
        for t in xrange(1,ndays_y):
            # Count the number of consecutive heatwave days 
            spell[t,EHF_exceed[t,...]==True] = \
                    spell[t-1,EHF_exceed[t,...]==True]+1

        # Count the length of each of the spells
        spell_final = np.zeros(tave_y.shape, dtype=np.int)
        # Identify negative changes in spell to locate the end of the heatwave
        spell_final[np.diff(spell,axis=0)<0] = spell[np.diff(spell, axis=0)<0]
        # The period might end in a heatwave day
        spell_final[-1,...] = spell[-1,...]

        # Define HEATWAVE variables
        heatwave_EHF_avg = np.zeros(tave_y.shape, dtype=np.float64)
        heatwave_EHF_peak = np.zeros(tave_y.shape, dtype=np.float64)
        heatwave_EHF = np.zeros(tave_y.shape, dtype=np.float64)

        # Calculate the peak and the average of the heatwaves
        #for i in xrange(nlon):
        for j in xrange(nstat):
            for t in xrange(ndays_y):
                # For each time that there is a heatwave, locate the 
                # beginning and compute the mean and the peak
                spell_l = spell_final[t,j] 
                if spell_l>=3:
                    heatwave_EHF_avg[t,j] = \
                            np.mean(EHF[t-spell_l:t+1,j])
                    heatwave_EHF_peak[t,j] = \
                            np.max(EHF[t-spell_l:t+1,j])
                    heatwave_EHF[t-spell_l:t+1,j] = \
                            EHF_exceed[t-spell_l:t+1,j]
 
        ########## PULLING OUT THE 5 CHARACTERISTICS BASED ON THE 3 ARRAYS
        # CALCUALTED ABOVE #########
        hwavg_masked = np.ma.masked_equal(heatwave_EHF_avg, 0) 
        hwpeak_masked = np.ma.masked_equal(heatwave_EHF_peak, 0)
        HWA_EHF_yearly[year,...] = np.max(heatwave_EHF_peak, axis=0)
        HWM_EHF_yearly[year,...] = np.ma.mean(hwavg_masked, axis=0)
        HWF_EHF_yearly[year,...] = np.sum(heatwave_EHF, axis=0)*100./ \
                                    np.float(np.sum(years_y==syear+year))
        HWN_EHF_yearly[year,...] = np.sum(spell_final>=3,axis=0)
        HWD_EHF_yearly[year,...] = np.max(spell_final,axis=0)
        HWD_EHF_yearly[HWD_EHF_yearly<3] = 0
        for j in xrange(nstat):
            #for i in xrange(nlon):
            # Find the days heat wave events occur.
            event_index = np.where(spell_final[:,j]>=3)
            if event_index[0].size:
                # Find the duration of the first event -1.
                duration = spell_final[event_index[0][0],j]-1
                # Find the starting day of this event.
                HWT_EHF_yearly[year,j] = event_index[0][0] \
                                    -duration-shift_start_year
            else:
                # If the array is empty there were no heat waves 
                # for this location.
                HWT_EHF_yearly[year,j] = 0
        HWM_EHF_yearly[HWM_EHF_yearly==0]=const.missingval
        HWA_EHF_yearly[HWA_EHF_yearly==0]=const.missingval
        HWT_EHF_yearly[HWT_EHF_yearly<0]=0
    return HWA_EHF_yearly, HWM_EHF_yearly, HWF_EHF_yearly, \
            HWN_EHF_yearly, HWD_EHF_yearly, HWT_EHF_yearly, pct
