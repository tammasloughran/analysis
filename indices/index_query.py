# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 14:55:39 2017

@author: Tammas Loughran
"""
from pandas import concat
def load_index(fname, standardize=False):
    """Load an index from a file in a single column. 
    
    Arguments
    fname -- name of csv file containing monthly nino3.4 data.
    standardize -- if True the series is standardized.

    Returns
    series -- pandas time series of index.
    """
    from numpy import genfromtxt, nan
    from pandas import date_range, Series
    data = genfromtxt(fname)
    years = data[:,0].astype(int)
    dates = date_range(str(years[0]), str(years[-1]+1), freq='M')
    series = Series(data[:,-1], index=dates)
    series = series.replace(-99.9, nan)
    series = series.replace(-999, nan)
    if standardize == True:
        mu = series.mean()
        sigma = series.std()
        series = (series - mu)/sigma
    return series

def load_index2(fname, standardize=False):
    """Load an index with timestamps from a text file with multiple
    columns by month.

    Arguments
    fname -- name of whitespace delimited file.
    standardize -- if Ture the series is standardized.

    Returns
    series -- pandas time series of index.
    """
    from numpy import genfromtxt, nan
    from pandas import date_range, Series
    data = genfromtxt(fname)
    years = data[:,0].astype(int)
    index = data[:,1:]
    dates = date_range(str(years[0]), str(years[-1]+1), freq='M')
    series = index.reshape(index.shape[0]*index.shape[1])
    series = Series(series, index=dates)
    series = series.replace(-99.9, nan)
    series = series.replace(-999, nan)
    if standardize == True:
        mu = series.mean()
        sigma = series.std()
        series = (series - mu)/sigma
    return series
    
import numpy as np
start_year = 1911
end_year = 2012
year_range = np.arange(1911, 2014+1, 1)
period = np.where((year_range>=start_year)&(year_range<=end_year))[0]
directory = ('/srv/ccrc/data35/z5032520/')
fname = 'indices/SAM/SAM_index_monthly_visbeck.txt'
sam1 = load_index2(directory+fname)
fname = 'indices/SAM/SAM_monthly_1957_2014_BAS.txt'
sam2 = load_index2(directory+fname)
sam = concat([sam1, sam2])
samslice = concat([sam1,sam2])
samslice = samslice['%s-07'%(start_year):'%s-06'%(end_year)]
axis = samslice.index.month
months = (axis==11)|(axis==12)|(axis==1)|(axis==2)|(axis==3)
samslice = samslice[months]
samslice = samslice.resample('AS-JUL', how='mean')

modoki_years = [1986, 1990, 1991, 1992, 1994, 2002, 2004, 2009]
elnino_years = [1911, 1913, 1914, 1918, 1925, 1930, 1941, 1951,1957, 1965, 1969, 1972, 1976, 1982, 1987, 1997, 2006]

print 'Modoki SAM:'
mean = 0
for year in modoki_years:
    print year, ': ', float(samslice[str(year)])
    mean += float(samslice[str(year)])
mean /= len(modoki_years)
print 'Modoki Mean: ', mean

print 'EP SAM:'
mean = 0
for year in elnino_years:
    print year, ': ', float(samslice[str(year)])
    mean += float(samslice[str(year)])
mean /= len(modoki_years)
print 'EP Mean: ', mean

fname = 'indices/DMI/DMI_monthly_1870_2014_hadisst.dat'
dmi = load_index(directory+fname)
# DMI
axis = dmi.index.month
months = (axis==11)|(axis==12)|(axis==1)|(axis==2)|(axis==3)
dmislice = dmi[months]
dmislice = dmislice['%s-07'%(start_year):'%s-06'%(end_year)]
dmislice = dmislice.resample('AS-JUL', how='mean')
print 'Modoki IOD:'
mean = 0
for year in modoki_years:
    print year, ': ', float(dmislice[str(year)])
    mean += float(dmislice[str(year)])
mean /= len(modoki_years)
print 'Modoki Mean: ', mean

print 'EP IOD:'
mean = 0
for year in elnino_years:
    print year, ': ', float(dmislice[str(year)])
    mean += float(dmislice[str(year)])
mean /= len(modoki_years)
print 'EP Mean: ', mean