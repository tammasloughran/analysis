"""EHF_ACORNSAT.py calculates heatwave stats from ACORNSAT station data 
"""
import glob
import numpy as np
import pandas as pd
import datetime as dt
import os
import compute_EHFheatwaves_station as comp
import netCDF4


def dstr_to_dt(string):
    """Convert a string do a datetime object with the format %Y%m%d
    """
    return dt.datetime.strptime(string, '%Y%m%d')


def load_acornsat(fnamelist):
    """Load the acornsat data as a pandas dataframe object.
    """
    station_no = 1
    dataframe = pd.read_table('acorn.sat.maxT.002012.daily.txt',
            sep=' ',
            skipinitialspace=True,
            skiprows=1,
            index_col=0,
            header=None,
            names=['Dates',station_no],
            na_values=99999.9,
            parse_dates=True,
            date_parser=dstr_to_dt)
    dates = dataframe.index
    data = pd.DataFrame(index=dates)
    for file_name in fnamelist:
        station_no = file_name[15:21]
        dataframe = pd.read_table(file_name,
                sep=' ', 
                skipinitialspace=True,
                skiprows=1,
                index_col=0,
                header=None,
                names=['Dates',station_no],
                na_values=99999.9,
                parse_dates=True,
                date_parser=dstr_to_dt)
        data[station_no] = dataframe
    return data


def do_detrend(data):
    '''do_detrend detrends station data using an OLS linear 
    regression while ignoring missing values.

    Arguments
    data -- station data with time on the 0 axis and station on the 1 axis

    Returs
    data -- detrended data
    '''
    from numpy import array, arange, nan, ones, empty, where, isnan
    from scipy.stats import linregress
    for station in range(data.shape[1]):
        # Select the station to detrend
        a_series = data[:,station]
        # Create an empty array and independant variable
        result = ones(len(a_series))*nan
        x = arange(len(a_series))
        # Locate nans and remove from the series and independant variable
        nan_index = isnan(a_series)
        a_series = a_series[where(nan_index==False)]
        x = x[where(nan_index==False)]
        # Detrend
        slope, intcpt, r, p, err = linregress(x, a_series)
        trend = x*slope + intcpt
        result[nan_index==False] = a_series - trend
        # Return series to data matrix
        data[:,station] = result
    return data


if __name__ == '__main__':
    syear = 1911
    eyear = 2014

    # Define the file names.
    os.chdir('/srv/ccrc/data35/z5032520/ACORNSAT/')
    file_list = sorted(glob.glob('acorn.sat.maxT.??????.daily.txt'))
    thres_file = 'EHF_heatwaves_ACORNSAT_1961_1990_summer.nc'
    #thres_file=None

    # Load the files.
    tmax = load_acornsat(file_list)

    # Select the analysis period
    dates = tmax.index
    tmax = tmax[(dates.year>=syear)&(dates.year<=eyear)]
    dates = tmax.index

    # Find the latitude and longitudes for each station
    latslons = np.genfromtxt('BoM_STN_latlon.txt')
    station_ids = tmax.columns.values
    lats = np.zeros(len(station_ids))
    lons = np.zeros(len(station_ids))
    for i, station in enumerate(station_ids):
        station = float(station)
        lats[i] = latslons[np.where(latslons[:,0]==station)[0],1]
        lons[i] = latslons[np.where(latslons[:,0]==station)[0],2]

    # Convert to a regular numpy array for compute_EHF
    tmax = np.array(tmax)
    # Detrend
    tmax = do_detrend(tmax)

    # Repeat for tmin
    file_list = sorted(glob.glob('acorn.sat.minT.??????.daily.txt'))
    tmin = load_acornsat(file_list)
    dates = tmin.index
    tmin = tmin[(dates.year>=syear)&(dates.year<=eyear)]
    dates = tmin.index
    tmin = np.array(tmin)
    tmin = do_detrend(tmin)
#    import numpy as np
#    np.save('tminfile', tmin)
#    np.save('tmaxfile', tmax)
    # Caclulate heatwave stats
    HWA, HWM, HWF, HWN, HWD, HWT, pct = comp.compute_EHF(tmax,
            tmin,
            dates=dates,
            modified=True,
            thres_file=None,
            season='summer',
            bsyear=1961,
            beyear=1990)

    # Save to a netCDF file.
    years = range(syear,eyear+1,1)
    season = 'summer'
    ofile_name = 'EHF_heatwaves_ACORNSAT_%s_%s_%s.nc' %(syear, eyear, season)
    outfile = netCDF4.Dataset(ofile_name, mode='w')

    # Create dimentions
    outfile.createDimension('time',len(HWA[:,0]))
    outfile.createDimension('day',365)
    outfile.createDimension('station', len(HWA[0,:]))

    # Create variables
    otime = outfile.createVariable("time", 'i', dimensions=('time'))
    setattr(otime, "long_name", "Time")
    setattr(otime, "units", "year")
    setattr(otime, "standard_name", "time")

    oday = outfile.createVariable("day", 'i', dimensions=('day'))
    setattr(oday, "long_name", "Day of Year")
    setattr(oday, "units", "days since 0000-01-01")

    olat = outfile.createVariable("lat",'f', dimensions=('station'))
    setattr(olat, "long_name", "Latitude")
    setattr(olat, "units", "degrees_north")
    setattr(olat, "standard_name", "latitude")

    olon = outfile.createVariable("lon",'f',dimensions=('station'))
    setattr(olon, "long_name", "Longitude")
    setattr(olon, "units", "degrees_east")
    setattr(olon, "standard_name", "longitude")

    ostation = outfile.createVariable("station", 'f', dimensions=('station'))
    setattr(ostation, "long_name", "Station ID")

    oHWA = outfile.createVariable("HWA", 'f', dimensions=(['time','station']))
    setattr(oHWA, "long_name"," Heat Wave Amplitude")
    setattr(oHWA, "units", "C2")
    setattr(oHWA, "description", "Peak of the hottest heatwave")

    oHWM = outfile.createVariable("HWM", 'f', dimensions=(['time','station']))
    setattr(oHWM, "long_name", "Heat Wave Magnitude")
    setattr(oHWM, "units", "C2")
    setattr(oHWM, "description", "Average of heatwave magnitude")

    oHWF = outfile.createVariable("HWF", 'f', dimensions=(['time','station']))
    setattr(oHWF, "long_name", "Heat Wave Frequency")
    setattr(oHWF, "units", "days")
    setattr(oHWF, "description", "Number of heat wave days")

    oHWN = outfile.createVariable("HWN", 'f', dimensions=(['time','station']))
    setattr(oHWN, "long_name", "Heat Wave Number")
    setattr(oHWN, "units", "heatwaves")
    setattr(oHWN, "description", "Number of heat waves")

    oHWD = outfile.createVariable("HWD", 'f', dimensions=(['time','station']))
    setattr(oHWD, "long_name", "Heat Wave Duration")
    setattr(oHWD, "units", "days")
    setattr(oHWD, "description", "Duration of the longest heat wave")

    oHWT = outfile.createVariable("HWT", 'f', dimensions=(['time','station']))
    setattr(oHWT, "long_name","Heat Wave Timing")
    if season == 'summer':
        setattr(oHWT, "description", 
                "First heat wave day of the year from 1st of Nov")
    elif season == 'yearly':
        setattr(oHWT, "description", 
                "First heat wave day of the year from 1st of Jul")
    setattr(oHWT, "units", "day")

    opct = outfile.createVariable("PRCTILE90", 'f', dimensions=(['day','station']))
    setattr(opct, "long_name", "90th Percentile")
    setattr(opct, "units", "C")
    setattr(opct, "description", "90th percentile used as threshold")

    ostids = outfile.createVariable("stationid", 'i', 'station')
    setattr(opct, "long_name", "Station ID")
    setattr(opct, "cf_role", "timeseries_id")

    # Transfer data
    otime[:] = years
    oday[:] = range(1,366,1)
    olat[:] = lats 
    olon[:] = lons
    ostation[:] = station_ids
    oHWA[:] = HWA
    oHWM[:] = HWM
    oHWF[:] = HWF
    oHWN[:] = HWN
    oHWD[:] = HWD
    oHWT[:] = HWT
    opct[:] = pct
    ostids[:] = station_ids

    # Global attributes
    setattr(outfile, "author", "Tammas Loughran")
    setattr(outfile, "contact", "t.loughran@student.unsw.edu.au")
    setattr(outfile, "featureType", "timeSeries")
    setattr(outfile, "date", dt.datetime.today().strftime('%Y-%m-%d'))
    setattr(outfile, "script", "EHF_ACORNSAT.py")

    outfile.close()
