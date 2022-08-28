#!/usr/bin/env bash

#site_names=31463_Hurdal 31464_Hurdal NR31585_Flekkefjord NR31581_Lyngdal NR32361_Lyngdal 
#NR31682_Tysvar 32288_Sortland NR32249_Vik NR32182_Stryn NR31881_Sande NR31578_Kvinesdal 
#32087_Dovre 32379_Hemne 32441_Sel 32258_Maaselv 32032_VestreToten 31984_Namdalseid 31976_Namdalseid
# 31461_Nittedal 31513_Nes 31539_Modum 31652_Bygland 31767_Kongsvinger 31780_Vaaler 31941_Roeyrvik 
#31997_Verdal 32088_Lesja 32103_Halden 32139_Rennebu 32374_Saltdal 32404_Vinje 32409_Vang  
#32438_Porsanger 31519_Nissedal NR32170_Selbu NR31801_Trysil NR31883_Trysil NR32174_Selbu 
#NR32173_Selbu NR31579_Kvinesdal NR31577_Kvinesdal NR31799_Trysil NR31677_Suldal 31650_Valle 31714_Flaa 32085_Skjaak NR31906_Voss 32124_Engerdal
description="test_normscale"
namelist_file='options.nml'
mkdir /home/ecaas/soil_decomp/results/${description}
cp /home/ecaas/soil_decomp/src/paramMod.f90  /home/ecaas/soil_decomp//results/${description}/parameters_${description}_${site}.txt
for site in 31464_Hurdal 31714_Flaa 31650_Valle 31984_Namdalseid 32087_Dovre 31463_Hurdal 32438_Porsanger 32103_Halden NR32182_Stryn #NR31906_Voss
do
  cp options.nml /home/ecaas/soil_decomp/results/${description}/options_${description}_${site}.txt
  clm_data_file='/home/ecaas/Documents/sites/'${site}'_historical/lnd/hist/'${site}
  clm_surface_file='/home/ecaas/Documents/surface_files/'${site}'/surfdata_'${site}'_simyr2000.nc'
  ./run_script_dev $site $description $clm_data_file $clm_surface_file $namelist_file
done