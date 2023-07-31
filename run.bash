#!/usr/bin/env bash

#Review these variables before running: 
description='equal_input'
namelist_file='options.nml'
spinup_years=1000
only_spinup='True'
CLM_version='new' #"old" or "new", "old" uses old variable names.
mkdir Spinup_values
mkdir ./results/${description}
mkdir ./results/${description}/nml_options
mkdir ./results/${description}/Spinup_values
cp ./src/paramMod.f90 ./results/${description}/parameters_${description}_${site}.txt
cp ./src/mycmimMod.f90  ./results/${description}/main_${description}_${site}.txt
make

for site in langtjern #31464_Hurdal #31650_Valle 32103_Halden 32087_Dovre NR31927_Os 31464_Hurdal NR31585_Flekkefjord NR32361_Lyngdal NR31682_Tysvar 32288_Sortland NR32249_Vik NR32182_Stryn NR31881_Sande NR31578_Kvinesdal 32379_Hemne 32441_Sel 32258_Maaselv 32032_VestreToten 31984_Namdalseid 31976_Namdalseid 31461_Nittedal 31513_Nes 31539_Modum 31652_Bygland 31767_Kongsvinger 31780_Vaaler 31941_Roeyrvik 31997_Verdal 32088_Lesja  32139_Rennebu 32374_Saltdal 32404_Vinje 32409_Vang 32438_Porsanger 31519_Nissedal NR32170_Selbu NR31801_Trysil NR31883_Trysil NR32174_Selbu NR32173_Selbu NR31579_Kvinesdal NR31577_Kvinesdal NR31799_Trysil NR31677_Suldal 31714_Flaa 32085_Skjaak NR31906_Voss 32124_Engerdal NR32485_Stord 32246_SoerVaranger
do
  echo ${site}
  cp options.nml ./results/${description}/nml_options/options_${description}_${site}.txt
  clm_data_file='./test/'${site}'_historical/'${site}'_historical.clm2.' 
  clm_mortality_file='./test/'${site}'_mortality/mort_'${site}'_'
  clm_surface_file='./test/surfdata_'${site}'_simyr2000.nc'
  #clm_surface_file=/home/elisacw/soil_decomp/test/surfdata_langtjern_simyr2000.nc

  ./run_script_dev $site $description $clm_data_file $clm_mortality_file $clm_surface_file $namelist_file "./results/"${description}"/" $spinup_years $only_spinup $CLM_version #> log.txt
done
mv Spinup_values/* ./results/${description}/Spinup_values/
rm -r Spinup_values

