CDF      
      time       levdcmp       lndgrid       levsoi        natpft        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Fri Nov  4 12:00:27 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,HR,HR_vr,NPP,NPP_GROWTH,NPP_NACTIVE,NPP_NACTIVE_NO3,NPP_NACTIVE_NH4,NPP_NAM,NPP_NECM,NPP_NFIX,NECM,NAM,NEE,NEP,NACTIVE,LITFALL,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,T_SCALAR,NPP_NNONMYC,LEAFC_TO_LITTER_FUN,GPP,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN,SMIN_NO3_LEACHED 31539_Modum_historical.clm2.all.1968.nc -O 31539_Modum_historical.clm2.all.1968.nc
Sun Jan  9 16:25:51 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1968.nc /nird/home/ecaas/31539_Modum_historical/lnd/hist/31539_Modum_historical.clm2.all.1968.nc
created on 12/14/21 11:28:24       source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31539_Modum_hist_for_decomp    Surface_dataset       !surfdata_31539_Modum_simyr2000.nc      Initial_conditions_dataset        -31539_Modum_Spinup.clm2.r.1201-01-01-00000.nc      #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         9./31539_Modum_hist_for_decomp.clm2.h0.1901-02-01-00000.nc      Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:25:51 2022: Appended file /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1968.nc had following "history" attribute:
created on 12/14/21 11:28:24
   NCO       `netCDF Operators version 4.8.1 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)       (   CWDC_TO_LITR2C_vr                         	long_name         .decomp. of coarse woody debris C to litter 2 C     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d      5�   CWDC_TO_LITR3C_vr                         	long_name         .decomp. of coarse woody debris C to litter 3 C     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d      6X   CWDN_TO_LITR2N_vr                         	long_name         .decomp. of coarse woody debris N to litter 2 N     units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d      6�   CWDN_TO_LITR3N_vr                         	long_name         .decomp. of coarse woody debris N to litter 3 N     units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d      7    FROOTC_TO_LITTER                   	long_name         fine root C litterfall     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             7�   FROOTN_TO_LITTER                   	long_name         fine root N litterfall     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             7�   GPP                    	long_name         gross primary production   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             7�   HR                     	long_name         total heterotrophic respiration    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             7�   HR_vr                         	long_name         3total vertically resolved heterotrophic respiration    units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P      7�   LEAFC_TO_LITTER_FUN                    	long_name         leaf C litterfall used by FUN      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             7�   LEAFN_TO_LITTER                    	long_name         leaf N litterfall      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             7�   LITFALL                    	long_name         "litterfall (leaves and fine roots)     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             7�   NACTIVE                    	long_name         Mycorrhizal N uptake flux      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             7�   NAM                    	long_name         AM-associated N uptake flux    units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             7�   NDEP_TO_SMINN                      	long_name         *atmospheric N deposition to soil mineral N     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             7�   NECM                   	long_name         ECM-associated N uptake flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             7�   NEE                    	long_name         �net ecosystem exchange of carbon, includes fire and hrv_xsmrpool (latter smoothed over the year), excludes landuse and harvest flux, positive for source   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8    NEP                    	long_name         Unet ecosystem production, excludes fire, landuse, and harvest flux, positive for sink      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8   NPP                    	long_name         net primary production     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8   
NPP_GROWTH                     	long_name         Total C used for growth in FUN     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8   NPP_NACTIVE                    	long_name         Mycorrhizal N uptake used C    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8   NPP_NACTIVE_NH4                    	long_name         Mycorrhizal N uptake use C     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8   NPP_NACTIVE_NO3                    	long_name         Mycorrhizal N uptake used C    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8   NPP_NAM                    	long_name         AM-associated N uptake used C      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8   NPP_NECM                   	long_name         ECM-associated N uptake used C     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8    NPP_NFIX                   	long_name         Symbiotic BNF uptake used C    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8$   NPP_NNONMYC                    	long_name         Non-mycorrhizal N uptake used C    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8(   PCT_NAT_PFT                       	long_name         =% of each PFT on the natural vegetation (i.e., soil) landunit      units         %      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       <      8,   QDRAI                      	long_name         sub-surface drainage   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8h   SMIN_NO3_LEACHED                   	long_name         soil NO3 pool loss to leaching     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             8l   SOILICE                       	long_name         4soil ice (natural vegetated and crop landunits only)   units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       P      8p   SOILLIQ                       	long_name         =soil liquid water (natural vegetated and crop landunits only)      units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       P      8�   TSOI                      	long_name         <soil temperature (natural vegetated and crop landunits only)   units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       d      9   T_SCALAR                      	long_name         'temperature inhibition of decomposition    units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P      9t   W_SCALAR                      	long_name         .Moisture (dryness) inhibition of decomposition     units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P      9�   levdcmp                	long_name         2coordinate levels for soil decomposition variables     units         m      axis      Y         d      5@   levsoi                 	long_name         Dcoordinate soil levels (equivalent to top nlevsoi levels of levgrnd)   units         m      axis      Y         P      5�   mcdate                  	long_name         current date (YYYYMMDD)             :   time                standard_name         time   	long_name         time   bounds        time_bounds    units         days since 1850-01-01 00:00:00     calendar      noleap     axis      T               :   time_bounds                    	long_name         history time interval endpoints             :<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�747	=�7��,76Z�<5��5J�4��$                                                                    68�\6-[u6��6)a�5�:A4۟34@^�3�,�                                                                    2�82���3$t#2�$11=��0���08�                                                                    1�o1���2O�1��L1��0o_�/ѭ�/iqR                                                                    6��3��6]��7"z�7n��7��B8Z|-7�)s7R�/6�Φ60�5�b�                                                7b�3-7���        2���    7�������H                                    B��@Q�5                @Qе?	�                ?�G�    5�*O�:A�1�A�y                                                                        ?���Av�>AГ�B�8B:��Bq%�B���B� <�<�<�<�<�<�<�<�<�<�<�<�C��|C��$C��LC��+C��SC�+�C�^�C��YC�̯C��PC�3�C�k�C���C�ߺC�pC�T�C���C��=C���C��MC��hC�v6C�hoC�r�C�v�>��3>�v`>��}>�Ȗ>��r>��>��v>���>�cW>��>���>��>�O4>���>�J�>���>��
>Ä!>Č>Ģ>���?D��?F[Z?I�?M;�?T-5?_�>?u��                                                ,K�G(] @��    @��    8.�k8n7�?�736W{5�
�5CQ�5
A                                                                    7\Ľ7$]�6�?6(@v5��4�y\4v��4.�                                                                    3�y�3�Ϗ3$��2�*m1�h�1H��0���0��2                                                                    2��2�!!2P6�1�^1h0}a�0u�/�Y                                                                    6�-)3��6ٌ�7a��8dql8��t8[P�7��y7O��6���6f�G6JK                                                7�+38+7�ɓ        2���    7���������i�                                    B��@Q�5                @Qе?	�                ?�G�    6>l�-�bAe��?='                                                                        @H�AxL�AѣxB��B;��BrxdB��2Bŧ]<�<�<�<�<�<�<�<�<�<�<�<�C���C��C���C���C���C�UC�?�C�j�C��+C���C��+C��C�G#C�|�C���C��C�(�C�WC�yxC���C���C�|�C�jJC�r�C�v�>�f
>�/z>��>��6>�d�>�=>��>�<">�eC>�z>��v>�WI>��>�	�>�Zq>��Y>�O>��z>��>� �?!��?Ew?G[�?JQ?Nc?U�?a��?y1                                                ,L-G(y @��    @�     8��98=��7���73�6i��6 C�5��p54�K                                                                    8�7o��6�^/67i�5��$5"�4��=4dw                                                                    4�ȡ3Τ�34�$2�?
1��:1��1A��0�#�                                                                    3��m3��2d,c1���1 ۹0���0u/��                                                                    6��B3�a�7��7��T9V�9�8pX�7� B7a��7��6ϓ�6M��                                                7�43��7���3/&�-OT�2�d|1*�r7	��	� 7-��6hJ�1!L�/��-[��        *�#c1M�    B��@Q�5                @Qе?	�                ?�G�    6�>w0��@+�q                                                                            @���A��A�B BJ-	B��qB���B�g�<�<�<�<�<�<�<�<�<�<�<�<�C��\C�ԥC��C���C��C�&UC�A�C�aUC��C��#C��?C���C��C�BC�v�C��C�� C��C�B>C�bC�C�}�C�l.C�rtC�v�>���>��.>��6>��
>���>��>�!>���>�e>��>��Y>��:>�4�>��>���>��'>�]>���>��c>�t5?V�b?S��?U��?X&?\{?ba?m`?~�q                                                ,L�G(� @�     @�     9#2�8}'o7�r-7"��6\��5�$5,�4��D                                                                    8N$�7��7 ~6Mn65�VP4�cr46_3���                                                                    4��H4	��3]��2�=�1�p�17B�0�%B0)�V                                                                    3�33.>�2�?1��1ۏ0g|�/��/Vs�                                                                    6�=3�\d8w+7��9K�91^C8���8 ��7TZ>6���6'��5��                                                7�3`Y7�mP4�VZ/�z12�ͯ2��۷�b�7�b�8&��7��q2���1x�N(>�_        ,�Y(2��    B��@Q�5                @Qе?	�                ?�G�    7��f)鑡9V,                                                                            A�zvB�:B X�B7ǠBv�B�
�B�CB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C�%�C��C���C��1C���C��C��UC��C��@C���C���C��0C��!C�"�C�N�C��C��C�� C��C�8�C�bpC�yC�m�C�r\C�vk>�#�>�QG>�f�>��'>�>�\a>���>�US>��<>�0�>�">�m�>�
�>��>�wf>�7�>�"?>��>���>��Y?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�                                                  ,L�G(� @�     @��    9C�>8���7���7,T6[*�5��5a�4��                                                                    8w�7��I7�=6Yi�5�k�4��4+F3�l�                                                                    4Ռ�4!�3xx 2���1���1-]�0���0#s^                                                                    4�~3L}�2��1��Y1۷0Z�$/�d�/Nv�                                                                    6��33�?�8���7���9tG�9P�l8�~�8>�7S$�6���6 >5�C^                                                6�N43p7�j�5�0B�2�,2�q�*�8*�8l�~8>�F3��2:>�,z�        .�N�3RlH    B��@Q�5                @Qе?	�                ?�G�    7t�u,�O                                                                                A?4:A⅀B%B3��Bt��B�
�B�CB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C���C�u�C�=^C�
�C��2C��*C�k�C�5uC�	>C��C�ôC��RC��PC��+C���C���C���C��QC��AC�_C�F�C�pMC�n�C�rPC�vW>��->�cn>�^>�~
>�
>Öh>�)>��0>�2d>�/>��>���>�fs>��`>���>�F�>���>���>���>��s?w�<?yd�?{i�?}�?~�f?��?�  ?�                                                  ,MYG(� @��    @��    9W�^8���8	%b7l��6��A6'. 5�)D4�@�                                                                    8�H7�~J7-<�6��P5�R5S,�4���3�6�                                                                    4�(-40mP3�u�3�2Mƫ1�5�1$�V01��                                                                    4�*3^��2���2"��1���0�(�0P�/`�                                                                    6�	43��-8�
�7��9��9bp�8���87Le7�=677�6�`5���                                                6�2���7��4�R /�9�2�Ӭ2����Կ�7��F8a��8Gʛ6h=�541�0>        5�6B�    B��@Q�5                @Qе?	�                ?�G�    7�/,2i                                                                                @��A�E�A�.�B��BK8-B��uB�z�B�;)<�<�<�<�<�<�<�<�<�<�<�<�C��C���C�,<C��C���C�$�C��DC�O0C��C���C�?MC���C��sC�b`C�*�C�mC���C��JC��'C��C�3nC�e�C�nlC�rGC�vB?m�?��?w>�W>�:�>�U>� �>�6>�1�>�B7>͆�>�?�>à�>��o>�+�>��Z>���>�rE>�#}>�ef?R.o?Sal?V6�?Zb?aa�?m�?zv?�                                                  ,M�G(� @��    @�`    9N`j8���8��7g�6��S6J_5���4��<                                                                    8�W�7Ž�7(�y6�v5��5?�5Ψ4�                                                                    4��4*�3���2���2Jg�1�<�1���0���                                                                    4p3WL2�ۊ2�%1�510���/�:                                                                    6��:3��'8�M67��h9� 9Z��8� U80�E7�P�7U�97�26k6                                                6��2�Q37��x4�T1/M}h2�1�2M6|����7��_8:�c82<*6��25Q�2��        5VU�6�H    B��@Q�5                @Qе?	�                ?�G�    6���+1�                                                                                @oA��0A�#B
1GBB�?B~7�B�DB��$<�<�<�<�<�<�<�<�<�<�<�<�C��C��9C�^�C�%@C��C���C�_�C�oC�̣C��C�B�C���C��*C�Z�C�OC��#C�w�C�IC�.�C�'TC�1�C�\C�m�C�r:C�v.?�A?}!?s�?��>���>��G>�U>�1�>��>�.>��>ٸb>ԃ2>�)�>���>���>���>��>�_
>��V?J�?K��?N"?Q�[?Wf2?a@i?q�q?��                                                ,N!G) @�`    @�"@    97�58���7��!7Q��6���64�5��5���                                                                    8hY7��u7�6���5�$5d��5{4�r                                                                    4�vA4#�3�:�2�27�1�6S1�o�1P
*                                                                    3�73Ap�2��b2g1g330�30�I�0�d�                                                                    6�+3��~8���7�'|9en49Er>8�{�8��7��[7:�7�6ݘ�                                                6���3�!7��04��R/U;Y2�k�2LR�����6���8�8��6h}5�2
n        4��6B��    B��@Q�5                @Qе?	�                ?�G�    2g��%�w�                                                                                @RğAa�dA�r�A���B/�Bb-tB���B���<�<�<�<�<�<�<�<�<�<�<�<�C��EC���C�|	C�YC�3�C�C���C��C�P�C��C�ԞC��NC�D�C��)C��9C�H�C���C���C��|C�^�C�F�C�W]C�l"C�r!C�v?"?ȉ?�?>?�@? ��>���>��`>���>�`w>�~�>�Z�>��>ٌ,>Ӵ�>���>���>Ĉ>�Pe>�A5?6>�?7�?8l!?:�$?>�?D�?Lc�?Y,                                                ,N�G)1 @�"@    @�&     9#mB8}�7���7@3/6�/�6&C�5�
\5�}�                                                                    8No7�щ7
�&6r�d5���5R�5�4�h�                                                                    4�4	�3o}�2�vu2(
$1�71s�1>0I                                                                    3���3.+�2�B2J�1TB�0��*0���0p<�                                                                    6�ϸ4�)8@a 7�n�9L�Q92ow8���8�R7��=7)��6��t6��3                                                6�X�3)?�7��4�|/�R2��2%�-6��\���;7� �7��P5�$�4i:�1��        4F95��c    B��@Q�5                @Qе?	�                ?�G�                                                                                            @XeAc�fA�{�A�F�B-��B]��B�SWB���<�<�<�<�<�<�<�<�<�<�<�<�C�;C�FJC�L�C�PC�P�C�M�C�DC�1�C�&C���C�٭C���C�z�C�>yC���C���C�Z�C�.C��tC���C�jC�Y�C�j�C�q�C�v>���>��>��H>��!>��>�	>�/�>�>�W>씍>��>�N>���>ޤ�>ٲ�>�i�>�(�>�P�>�5�>��?8�?8%5?9_�?:YE?<-R??w�?D�o?MlV                                                ,N�G)O @�&     @�)�    9	d88Z��7��c7-F�6���6VK5�R�5tT                                                                    8-�7�a6��G6Z�35�F�5ID�5��4�P                                                                    4���3�<g3S4�2�׬2�W1���1t�v1'�                                                                    3�#/3v�2�d�11Dz0�^90�mX0(2=                                                                    6��4��7��h7�-9-&�9��8��8�O7�@�7"�c6��6�                                                6�¹3��7�D�3ͭ.Pk\2��y1�*�7������{6\|�7#�3)��1�δ-t��        /;ݹ3sp    B��@Q�5                @Qе?	�                ?�G�    6� ',��0;؄                                                                            @~K�A��A��B	
�B@C�By��B�܄B��:<�<�<�<�<�<�<�<�<�<�<�<�C�C�W�C��BC��,C��3C�&C�F�C�s)C��@C��7C���C���C��C��C�ϺC��`C�~C�EZC�
C���C���C�cqC�j>C�q�C�u�>���>�m>>�~>�;�>� >��}>��>��M>�6>��>֫->��9>�XN>�(�>��p>Զ>>с�>Ϳ_>��>�|*?JU�?J�?MP?O��?T?"?\K?i�8?|��                                                ,OMG)n @�)�    @�-�    8�]�8I�p7��w7"�?6���6U�5>`x4�,                                                                    8�7~�=6�0�6MGx5�$5=�J4py�3ވ�                                                                    4���3��3D�.2�P2�W1��=0ς�0@k                                                                    3��3
�2x��1߹=19o+0���0(/r��                                                                    6���4!6�(7��9f09��8{\<7�#"7y|�7��6M =5��                                                6��@3287���        2�H�    7���������                                    B��@Q�5                @Qе?	�                ?�G�    7/�-��#=�N                                                                            @��A�T�A��B�BQ��B�JB��[B�Ph<�<�<�<�<�<�<�<�<�<�<�<�C�M]C���C��RC��C�S�C���C���C�!�C�]#C���C��qC��HC��C�0�C�GXC�P�C�KC�6nC��C��C���C�q`C�j�C�q�C�u�>��?>�Ѓ>�K!>�w>��>���>��V>��T>�C>�Z�>�X>��>�s>�s>��}>΁�>��>��;>ʭp>�,�?W�?X}�?[�K?_�?e��?o�?}�(?�                                                  ,O�G)� @�-�    @�1�    7��7��6�$�6��M6Z�-5�' 51�4�T�                                                                    6Ho�6%o�5�I5���5�4��4F��3�'�                                                                    2���2��2T�f2!E1�E�1B.0�XI09�                                                                    1�rg1�M1�R�1K��1|�0uG�/�o�/j��                                                                    6��4��6cl�6�9x7n+�7���7�]7c�67La�6��%6.z5�K�                                                7 "�3S17��*        2�j�    7�PH��PH���                                    B��@Q�5                @Qе?	�                ?�G�    6���,��A���A���A��A�%                                                                ?Z�@��A�fnA���BBI�B}/�B��.B�K�<�<�<�<�<�<�<�<�<�<�<�<�C�)�C�O�C���C���C���C�?�C���C��C�&�C�g9C��qC���C�)C�e�C���C��QC���C��C���C��jC���C��C�lC�q�C�u�>���>�K�>���>��E>�o�>��>�C>�_�>�M�>��>��u>�h�>�->���>�&v>��>� �>��>ȿ@>Ƿ2=�)2>U>�>��?<��?S�O?\M�?g%�?u[�                                                ,ruG)� @�1�    @�5`    