#!/usr/bin/env bash

#Review these variables before running: 
description='test'
namelist_file='options.nml'
spinup_years=20
only_spinup='False'

mkdir ./results/${description}
cp ./src/paramMod.f90 ./results/${description}/parameters_${description}_${site}.txt
cp ./src/mycmimMod.f90  ./results/${description}/main_${description}_${site}.txt

for site in BOR2
do
  echo '/home/ecaas/${site}_hist_BGC/lnd/hist/'${site}'_hist_BGC.clm2.'
  cp options.nml ./results/${description}/options_${description}_${site}.txt
  clm_data_file='/home/ecaas/'${site}'_hist_BGC/lnd/hist/'${site}'_hist_BGC.clm2.'
  clm_surface_file='/home/ecaas/surfdata_0.9x1.25_hist_16pfts_Irrig_CMIP6_simyr2000_'${site}'_c221027.nc'
  ./run_script_dev $site $description $clm_data_file $clm_surface_file $namelist_file "./results/"${description}"/" $spinup_years $only_spinup
done
