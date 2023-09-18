#!/usr/bin/env bash

#
#requires package nco to run. Install to conda environment with "conda install -c conda-forge nco"
mkdir ../test/tmp_folder
mv ../test/langtjern_historical/lnd/hist/*1901* ../test/tmp_folder
for f in ../test/langtjern_historical/lnd/hist/*.nc
do
        ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,NPP_NACTIVE,LEAFC_TO_LITTER,QDRAI,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,nbedrock,T_SCALAR,NPP_NNONMYC,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN $f -O $f
done
mv ../test/tmp_folder/*1901* ../test/langtjern_historical/lnd/hist/
rm -r ../test/tmp_folder