#!/bin/python
import numpy as np
import load_data
from pandas import concat
import scipy.stats as stats

directory = ('/srv/ccrc/data35/z5032520/')
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

start_year = 1911
end_year = 2012
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

a, p = stats.pearsonr(ninoslice, dmislice)
b, p = stats.pearsonr(ninoslice, samslice)
c, p = stats.pearsonr(ninoslice, strhslice)
d, p = stats.pearsonr(dmislice, samslice)
e, p = stats.pearsonr(dmislice, strhslice)
f, p = stats.pearsonr(strhslice, samslice)

print a, b, c
print '\t\t', d, e
print '\t\t', '\t\t', f
