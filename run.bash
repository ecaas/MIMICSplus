#!/usr/bin/env bash

#Review these variables before running: 
description='test_fixes2'
namelist_file='options.nml'
spinup_years=10
only_spinup='False'

mkdir ./results/${description}
cp ./src/paramMod.f90 ./results/${description}/parameters_${description}_${site}.txt
cp ./src/mycmimMod.f90  ./results/${description}/main_${description}_${site}.txt

for site in 31464_Hurdal
do
  cp options.nml ./results/${description}/options_${description}_${site}.txt
  clm_data_file='./test/'${site}'_historical/lnd/hist/'${site}'_historical.clm2.'
  clm_surface_file='./test/surfdata_'${site}'_simyr2000.nc'
  ./run_script_dev $site $description $clm_data_file $clm_surface_file $namelist_file "./results/"${description}"/" $spinup_years $only_spinup
done
