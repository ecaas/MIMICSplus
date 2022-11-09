CDF      
      time       levdcmp       lndgrid       levsoi        natpft        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Fri Nov  4 12:00:26 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,HR,HR_vr,NPP,NPP_GROWTH,NPP_NACTIVE,NPP_NACTIVE_NO3,NPP_NACTIVE_NH4,NPP_NAM,NPP_NECM,NPP_NFIX,NECM,NAM,NEE,NEP,NACTIVE,LITFALL,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,T_SCALAR,NPP_NNONMYC,LEAFC_TO_LITTER_FUN,GPP,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN,SMIN_NO3_LEACHED 31539_Modum_historical.clm2.all.1955.nc -O 31539_Modum_historical.clm2.all.1955.nc
Sun Jan  9 16:25:50 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1955.nc /nird/home/ecaas/31539_Modum_historical/lnd/hist/31539_Modum_historical.clm2.all.1955.nc
created on 12/14/21 10:58:55       source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31539_Modum_hist_for_decomp    Surface_dataset       !surfdata_31539_Modum_simyr2000.nc      Initial_conditions_dataset        -31539_Modum_Spinup.clm2.r.1201-01-01-00000.nc      #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         9./31539_Modum_hist_for_decomp.clm2.h0.1901-02-01-00000.nc      Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:25:50 2022: Appended file /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1955.nc had following "history" attribute:
created on 12/14/21 10:58:55
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�4o�5?�6S�`6���6U&�5��5xL4�N�                                                                    3�Tz4qWY5��/5��E5��>4��4Db3�H�                                                                    0~}0�91��2)91�.R1A�"0�[�01��                                                                    /$ի0p$1�1Awq1��0u�/��/`��                                                                    6�-83���6��?6�}�5d�6�B7O<7\��7D��6�t�6-��5�<�                                                7 D�3*7��        2$-h    7_�O�_�O��V                                    B��@Q�5                @Qе?	�                ?�G�    6C=f+��uA��EA���A���@�E                                                                ?��@�IA�x�A�/�B;]�BtR"B��B�j�<�<�<�<�<�<�<�<�<�<�<�<�C�?eC�Y�C��WC���C���C��C�=�C�}�C���C���C�(�C�c�C��5C��yC�C�U�C��PC���C��_C��oC�ɜC���C��C�x&C�o�>��[>��>��+>�Q�>�b�>��>��1>�<8>�U�>�<�>�_�>��|>��>���>�5�>���>���>�'�>ţv>�4�<X�=�5G>��C?=�n?M�"?V�?bX�?sF
                                                *O�G� @ⶠ    @⺀        4|�5/�=5��p6J{�5��5|4��                                                                        3�;�4^�4�$�5��4��w4>�3�B,                                                                        0	P0��81Pϻ1ܐ�1=�0��10+׏                                                                        /-r�/��0��1M�0o��/��/YI                                                                    6��4��6Ĭ$6s��    5�76��6���7=;M6�e�6+�*5�Ȇ                                                6� 32S
7�m^        29��    73w�3w���                                    B��@Q�5                @Qе?	�                ?�G�            A��VA�1�A� +A�<                                                                ?qf@��A9�A� bB3rBf�B���B���<�<�<�<�<�<�<�<�<�<�<�<�C�fLC�oC�}C��C��C�صC��C�8�C�h�C���C�ƝC��(C�2C�m�C��UC��C�,�C�a�C���C��C���C��QC���C�x\C�p>�Q]>��4>�e�>�uN>��>��>�b>��7>�%1>�xW>��>��>���>�2e>��@>���>�C�>�m>��/>À    <�$�>W�?=��?C��?IT�?Q�?_��                                                *P]G� @⺀    @�         4��H5H=:5���6E�5�H�5r@4�X�                                                                        3���4|�4��5x��4���47��3�pd                                                                        0 ��0�%1N�k1ֳu18h�0�q�0%�v                                                                        /K5�0	�20���1��0h��/�#�/Q�                                                                    6��83�C�7��,6s@    5Д6 �g6�Xw7:<�6�ۀ6)[f5�2                                                6��3�7�֨2�?]    2B��0�^���o�6�o�7�&5��0pa�/;>�,�E�        )��0@sw    B��@Q�5                @Qе?	�                ?�G�            A��A���A�҆ALS                                                                ?	�@��jA>:tA�ߥB/�3Ba^�B�N@B��X<�<�<�<�<�<�<�<�<�<�<�<�C���C���C��DC��#C��yC��hC���C��C�F�C�lxC��JC�ôC�� C�)�C�e�C���C��LC�JC�N�C�t�C���C��C��lC�x�C�p,>��?>�ԛ>�S>�vT>���>�Mi>�(?>�P�>�c�>�T!>���>���>�~�>�i�>���>�OH>�� >�T�>�Jd>��k    <�2&>-�S?:��?>��?CR�?JK�?Uy                                                *P�G @�     @���    7�H6�Cr6Q!6�b�6>��5�05	�4��,                                                                    6���5��u5F�R5�z5qA4�oe4-�3��n                                                                    3"��28`�1�]�2Q�C1��	1+!�0��M0"!                                                                    2Mg�1h��0�v"1�[�1Kd0X*�/�u�/L��                                                                    6�7�3�3�8P�u6�~�86�7��p7�c7�f75��6�ۖ6#6�5�F�                                                6���3Vg7��~4io    2k=�2Y¨���7��8�Y7���2W�.1&��.N�        +���2-8    B��@Q�5                @Qе?	�                ?�G�    2�5k))�AgqA���A��!?�                                                                @Z�A�hA���A�t#B.��B_ahB��bB��<�<�<�<�<�<�<�<�<�<�<�<�C��LC��TC���C���C��aC���C��C�]C�6�C�VJC�y�C��C��mC��!C�1�C�l�C��PC��#C��C�F�C�y)C���C��vC�x�C�pW>�E�>�	�>�|n>� �>�9>��+>��>��+>���>�/=>��>�F>�P�>��>��>�#\>���>��E>��>���>�&�>���>�މ?:j?<��?@��?G2N?P�p                                                *Q%G- @���    @�Š    9Q�8h�7�~�7-�6���6�45>J4�,                                                                    8B�7��D6��&6Z��5��5?�A4p�3���                                                                    4�3�Y3X92���2�1�.[0�
T0T�7                                                                    3�V�3 [2���1�1_1A`0ЦW0�'/�H�                                                                    6�g/3�jt8�T�7���9?��9%o�8��8��7��97$36d��5��                                                6��3�7�
5�    2�62�_I�j�8j�8`�l81Cl3�`2>c�-��        .���3Y1    B��@Q�5                @Qе?	�                ?�G�    7b_-�-;��                                                                            @���A��A�xB��BIp�B�]KB���B�IH<�<�<�<�<�<�<�<�<�<�<�<�C�[&C���C���C�Q�C�
�C��C��[C�Q�C�.�C�|C�kC��C��C�C�/�C�VfC��C��ZC��C�3C�YC���C���C�y?C�p�>��>�%>���>�!z>��,>���>���>�Ҫ>�٫>��>��>�V�>���>�M�>�>��3>��F>���>���>�b=?P�U?R�?TƠ?X�V?_~�?k�?{C?�                                                  *Q�GL @�Š    @�ɀ    9/QM8�ُ7�s7B�j6��6%n5�l�5��                                                                     8]t+7��7��6v=(5��5P��5�4��y                                                                    4��?4�,3ta2�Zb2*��1�7�1j�-1-5�                                                                    3�<35l2�kf2#1W��0�q0�8�0Z��                                                                    6�m�3�E<8�S�7�P9Y��9;b8��d8 +7��-7/Q96�7d6�l&                                                6��H2�o,7�vE4���-�ح2_�2m����#�7�$�8U:O8:r*6j�	5�r1�L$        5�M6C��    B��@Q�5                @Qе?	�                ?�G�    5ߋ�+�                                                                                @g))Aqu@A�;B\B8n"Bn�?B�B��s<�<�<�<�<�<�<�<�<�<�<�<�C�aC���C�A8C���C���C�MC��C��uC�2C��C���C�\�C�"nC��C�ɄC��XC���C��NC��*C��C�=�C�~�C���C�y�C�p�?t�>���>�(�>�2�>�C>���>�j>�C>̘<>�ۉ>�T�>�:>��/>��T>�u�>�9>�A>��>��>��l?@0?A'0?B�,?E�?J?�?Qv6?]	�?q��                                                *Q�Gj @�ɀ    @��@    9-��8���7���7<mA6�G&6v5� 95��                                                                     8[l#7�vD7��6nE5��5H)J4�(�4�BC                                                                    4�:�4�3l8)2�C%2$�81���1\1&�                                                                    3��30�~2�0�2��1P*b0�R0���0R��                                                                    6���3��8�N�7�:M9Wۼ97,8���8 E7�HI7% �6��f6���                                                6��#2��H7�p�4�m�.���2F��2)��^�67^��8 �8v�6��5.,}2��        5/Y|6c��    B��@Q�5                @Qе?	�                ?�G�                                                                                            @*1AB�zA��A�:}B ��BOF�Bz��B�n�<�<�<�<�<�<�<�<�<�<�<�<�C�jIC���C�C/C��LC���C�)C���C�&�C��/C�^�C��C���C�N#C���C��hC�b�C�0 C��C�
:C�-C�2RC�qC���C�y�C�p�?��?�o?r?
0a?�?�>��>�֔>�\>�P>ڇ�>�<�>�h�>��J>��>��j>�w�>���>�4I>�� ?��?�_? �3?$w?(qk?.�?5+?=2O                                                *RQG� @��@    @��     8��
8/6�7�*�7�&6_�l5�q�5��5r��                                                                    87]R�6��=6*w�5�p5t4��S4�ne                                                                    4t�3�ފ3#��2��1��B1��1-781T!                                                                    3���2�2N�g1���1�0���0Z�|0'&�                                                                    6�X3�Q�8��@7�o#9�>8���8gB7�^�7O�96���6�#�6�                                                6� 2�ŕ7���4\��.o��2MNE2�r�,E~7,O�7��7�6\U=5�3F�        5�683l    B��@Q�5                @Qе?	�                ?�G�                                                                                            ?�7�A��A��0A�kVB�mB,��BQ"�B�v�<�<�<�<�<�<�<�<�<�<�<�<�C��rC�[6C��C��C��%C�xC�,�C�ӡC���C�:�C��C���C�BC��C���C�"�C���C��vC�\�C�EC�?C�g�C��JC�y�C�q?u�?��?�e?
t-?o?Y�?44>�"�>��>>�,z>�G
>��>��x>�s�>��&>ˊ�>� �>���>�"�>���>�G�>� >�hh>�,>��4? �?w�?L                                                *R�G� @��     @��     8�v�8h�7���6�65�|5�;F5rv�5:tz                                                                    8/�7F�&6���6�5e�u4��4�"�4k��                                                                    4gtX3�zE3�2z�1��<1N(e1�0�!�                                                                    3�.�2ؚs2;�1��0��0�4v0&�t0 KB                                                                    6�33���8E�57��797�8�G08RJR7���7+�6���6�0p6R�                                                6�c92�^�7���4,a�.��2I��1�$���f�5�b7��M7���6��4���2�        4�.�5��g    B��@Q�5                @Qе?	�                ?�G�                                                                                            @O�A��A���A��vA��B ��B=Be��<�<�<�<�<�<�<�<�<�<�<�<�C�{C�q{C�oOC�l[C�gXC�^ZC�OpC�7�C�HC� �C��"C���C�|�C�=eC���C��FC�IoC���C��nC�� C�`�C�e�C��C�y�C�q8>��Z>���>�{�>�,�>��>�>��s>��{>���>��v>�"7>��T>��>ސ�>�M>ӡ?>�;>��>��A>��]>�zF>���>�6S>�R>��>ܸ
>�J�>儼                                                *SG� @��     @���    8��8/�7��p7
�6ab85�J�5���5Zk                                                                    8	��7^?6���6/o$5�X�5<�4��^4���                                                                    4m��3��T3*�F2�M�1���1��\1 �"0��Q                                                                    3�\2��B2Ww�1�1i0���0KTa0J                                                                    6���4k27��07��9��8���8r��7��7V��6���6�c�6u�w                                                6Į�3y�)7���3�-���2DҮ1gZ�7i��i6{��6�&�4y�Y37�1*\        2m�4PO�    B��@Q�5                @Qе?	�                ?�G�            >�u�:4R6                                                                        @2ajAG�1A�A�i_B#�BCk�Ba��B�k<�<�<�<�<�<�<�<�<�<�<�<�C�1�C���C��lC��C�O�C���C���C��C��C��C�$�C�%(C��C�3C���C��C�{C�;�C��1C��hC���C�lC�}�C�y�C�qd>�B�>�\�>�`�>��_>�6 >Ҷ�>��>��>�>�B>��>��|>�;�>��>�h�>�(k>�FF>�!�>�)�>���?"�-?$�?$�"?#�?!�<?�?�?b�                                                *S}G� @���    @�ܠ    8�x�8(�H7��A7�6_�5�5�,�5���                                                                    8��7Ux[6�J6*�5�<5;c4˗ 4��_                                                                    4eKL3�3$C2�)�1��1�*\1/��19                                                                    3��=2�/2ON�1��+1��0��#0]��02c                                                                    6���3��7��7���9�38�ݙ8j[7��G7WR�72�6�|`6�Y`                                                6��3��7�J�10=    26)M.�L�7�_۷�_�����4kv0A�.�5-7�        +,�0��    B��@Q�5                @Qе?	�                ?�G�            ?7Qo<;�^                                                                        @X�Ae�-A��]A��)B.��B_ �B��B�� <�<�<�<�<�<�<�<�<�<�<�<�C�g�C���C���C�C�H�C���C��CC��C�U�C���C��EC��dC�	C�?�C�X�C�a�C�W�C�=jC�[C��C���C�xbC�|gC�y�C�q�>�r >���>���>�<�>�+}>���>��d>��>���>��>�N�>�g>��>�o#>��>ϝ:>���>�;R>��q>��?8?9h#?:��?;�p?=QL?@l�?E�_?N��                                                *S�G @�ܠ    @��`    7Ʉ�7G$�6�yE6�?6@p 5�ͭ5�O�5h-                                                                    6���6{��6��5���5s65f�4�\4���                                                                    3[��2���2c2�1Ѩ�1gڋ1$0�ϟ                                                                    2���2	E1�q�14��1j�0�o0E:�0��                                                                    6�AY4`k6v7��7��8�47���7V�=7:`6�
6�\6~��                                                6��3$ 87���        2#��    7��#�����G                                    B��@Q�5                @Qе?	�                ?�G�            @��'A��lAD�A>�F@.l(                                                            ?��@���Af�aA�8�B(�BX�B�!;B���<�<�<�<�<�<�<�<�<�<�<�<�C��C�OC��,C�� C��C�P:C���C��]C�G�C��9C��EC��C�C�C�|�C���C���C���C��C���C��C���C��vC�|�C�y�C�q�>�X>���>�h�>�j�>�C�>��>�I�>� 2>�-0>��(>�z9>�(�>��>��>�P5>���>ȯ>�d{>��>��H>3�7>c>��>�x�?3��?9��?@_?Iӿ                                                *v�G" @��`    @��@    