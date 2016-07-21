#!/bin/bash
cdo add HadISST_sst_1970-2014_clim_N96_filled_24.nc modokisst_dt.nc temp.nc
ncrcat clim_december.nc temp.nc -o temp1.nc
ncrcat temp1.nc clim_january.nc -o modoki_model_in.nc

cdo add HadISST_sst_1970-2014_clim_N96_filled_24.nc nllnsst.nc temp.nc
ncrcat clim_december.nc temp.nc -o temp1.nc
ncrcat temp1.nc clim_january.nc -o nlln_model_in.nc

rm temp.nc
rm temp1.nc

