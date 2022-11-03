#!/usr/bin/env bash

description="test_fixes2"
namelist_file='options.nml'
mkdir ./results/${description}
cp ./src/paramMod.f90  ./results/${description}/parameters_${description}_${site}.txt
cp ./src/fluxMod.f90  ./results/${description}/flux_${description}_${site}.txt

for site in 31539_Modum

do
  cp options.nml ./results/${description}/options_${description}_${site}.txt
  clm_data_file='./test/'${site}'_historical/lnd/hist/'${site}
  clm_surface_file='./test/surfdata_'${site}'_simyr2000.nc'
  ./run_script_dev $site $description $clm_data_file $clm_surface_file $namelist_file "./results/"${description}"/"
done