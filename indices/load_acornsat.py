"""load_acornsat.py loads ACORN-SAT data from the BoM in a pandas dataframe
"""
import glob
import numpy as np
import pandas as pd
import datetime as dt


def dstr_to_dt(string):
    """Convert a string do a datetime object witht he format %Y%m%d
    """
    return dt.datetime.strptime(string, '%Y%m%d')


def load_acornsat(fnamelist):
    """Load the acornsat data as a pandas dataframe object.
    """
    first = True
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
        if first:
            data = dataframe
            first = False
        else:
            data[station_no] = dataframe
    return data


file_list = glob.glob('acorn.sat.maxT.??????.daily.txt')
tmax = load_acornsat(file_list)
file_list = glob.glob('acorn.sat.minT.??????.daily.txt')
tmin = load_acornsat(file_list)

import pdb
pdb.set_trace()

tave = (tmax+tmin)/2.



#for file_name in file_listi:
#    data = np.genfromtxt(file_list, skip_header=1, missing_values=99999.9)
#    start_date = str(data[0,0])
#    end_date = str(data[-1,0])
#
#    tmax[no,:] = data[:,1]
