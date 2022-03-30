#!/usr/bin/env bash

#site_names=31463_Hurdal 31464_Hurdal NR31585_Flekkefjord
description="test"
for site in NR32485_Stord #NR31908_Ullensvang 31463_Hurdal 31464_Hurdal NR31585_Flekkefjord
do

  clm_data_file='/home/ecaas/nird/'${site}'_historical/lnd/hist/'${site}
  clm_surface_file='/home/ecaas/nird/surface_files/'${site}'/surfdata_'${site}'_simyr2000.nc'
  echo $clm_data_file
  echo $clm_surface_file
  ./run_script_dev $site $description $clm_data_file $clm_surface_file
done
