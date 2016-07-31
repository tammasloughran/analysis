cdo remapbil,grid1.txt isolated_anomalies.nc isolated_anomalies_N96.nc

ncks -O -v elnino_sst isolated_anomalies_N96.nc -o pac_nino.nc
ncks -O -v lanina_sst isolated_anomalies_N96.nc -o pac_nina.nc
ncks -O -v piod_sst isolated_anomalies_N96.nc -o ind_piod.nc
ncks -O -v niod_sst isolated_anomalies_N96.nc -o ind_niod.nc

ncrename -v elnino_sst,sst pac_nino.nc
ncrename -v lanina_sst,sst pac_nina.nc
ncrename -v piod_sst,sst ind_piod.nc
ncrename -v niod_sst,sst ind_niod.nc

cdo -O add HadISST_sst_1970-2014_clim_N96_filled_24.nc pac_nino.nc pac_nino_ancil.nc
cdo -O add HadISST_sst_1970-2014_clim_N96_filled_24.nc pac_nina.nc pac_nina_ancil.nc
cdo -O add HadISST_sst_1970-2014_clim_N96_filled_24.nc ind_piod.nc ind_piod_ancil.nc
cdo -O add HadISST_sst_1970-2014_clim_N96_filled_24.nc ind_niod.nc ind_niod_ancil.nc

ncrcat -O clim_december.nc pac_nino_ancil.nc -o pac_nino_ancil.nc
ncrcat -O clim_december.nc pac_nina_ancil.nc -o pac_nina_ancil.nc
ncrcat -O clim_december.nc ind_piod_ancil.nc -o ind_piod_ancil.nc
ncrcat -O clim_december.nc ind_niod_ancil.nc -o ind_niod_ancil.nc

ncrcat -O pac_nino_ancil.nc clim_january.nc -o pac_nino_ancil.nc
ncrcat -O pac_nina_ancil.nc clim_january.nc -o pac_nina_ancil.nc
ncrcat -O ind_piod_ancil.nc clim_january.nc -o ind_piod_ancil.nc
ncrcat -O ind_niod_ancil.nc clim_january.nc -o ind_niod_ancil.nc
