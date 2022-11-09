CDF      
      time       levdcmp       lndgrid       levsoi        natpft        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Fri Nov  4 12:00:20 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,HR,HR_vr,NPP,NPP_GROWTH,NPP_NACTIVE,NPP_NACTIVE_NO3,NPP_NACTIVE_NH4,NPP_NAM,NPP_NECM,NPP_NFIX,NECM,NAM,NEE,NEP,NACTIVE,LITFALL,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,T_SCALAR,NPP_NNONMYC,LEAFC_TO_LITTER_FUN,GPP,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN,SMIN_NO3_LEACHED 31539_Modum_historical.clm2.all.1876.nc -O 31539_Modum_historical.clm2.all.1876.nc
Sun Jan  9 16:25:45 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1876.nc /nird/home/ecaas/31539_Modum_historical/lnd/hist/31539_Modum_historical.clm2.all.1876.nc
created on 12/13/21 20:40:18       source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31539_Modum_hist_for_decomp    Surface_dataset       !surfdata_31539_Modum_simyr2000.nc      Initial_conditions_dataset        -31539_Modum_Spinup.clm2.r.1201-01-01-00000.nc      #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         9./31539_Modum_hist_for_decomp.clm2.h0.1850-02-01-00000.nc      Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:25:45 2022: Appended file /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1876.nc had following "history" attribute:
created on 12/13/21 20:40:18
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�78z�6��7LR�6�J�6>�15�f�5�4�,q                                                                    6i�5סh6��6��5p��4��433˖t                                                                    2�RH2:J-2���2z<1��14�'0��0/�;                                                                    1�L�1kP82��1� 1O)0d�/�fS/^/                                                                    6�H�3�a	6�B�6� �7Lz�7W��8��7���7+x�6��p6Kk5�k                                                6��3
�n7��2        1�x    7�������                                    B��@Q�5                @Qе?	�                ?�G�            A�ԙA��?w�"                                                                    ?�dvA=A�F�A�ɱB#n�BP��By��B�rR<�<�<�<�<�<�<�<�<�<�<�<�C���C���C��C��5C��C�dC�3�C�e�C���C��.C��C�+C�P�C��;C��uC� �C�8�C�g^C��C���C���C��C�gsC�\�C�V�>�c�>�oq>�H^>���>�W�>�;>>�kC>���>�j>���>��>���>���>���>�#">���>��>��>�զ>��I>ʅ�>�*�?)�*?+�?-05?0K�?4�r?:v|                                                B	F� @     @�    8K(�7�_7�\N6��?6HB�5�q52�4苚                                                                    7�O�6�96��p6�y5|��4仁4a�t4��                                                                    3ݵ�3,�z3��2��1ڏ�1E��0�)�0}�=                                                                    3 2ZO�2>��1�̹1
	�0y��/��^/�I�                                                                    6��k3߮�6�~a77�!8W��8C�$8W��7��!75��6�)�6B�5�zL                                                6��3�7��        1O �    7R~�R}��|�                                    B��@Q�5                @Qе?	�                ?�G�            AL��@`�$                                                                        @��%AZ��A�HCA�Z�B)V�BW�eB�'B�G<�<�<�<�<�<�<�<�<�<�<�<�C��3C��9C���C��=C��=C��C�'�C�P�C�v�C��qC���C��C�!C�GPC�~1C���C��_C�$C�OC�nbC���C��"C�iC�\�C�W>�w>��5>�f>�>�>��s>�e>��&>��>�ܽ>��>��d>��>�X�>��>��>�f�>��M>���>�N>�4Y?Pk�?2_�?31�?4�?5��?8s�?<�0?CD�                                                BmF4 @�    @¦�    8��8/��7���7ڍ6\5�5�S�5d�                                                                    8)7^B�6��(6'в5��5_�4�Ռ4�                                                                    4x�+3��3%��2���1��1|�1�70��l                                                                    3�4�2� 2Q��1�%�1�d0��M0I�|06_                                                                    6�~�3�a7�pk7�6`8�
�8�h�8mf�7�h�7Hb�6��6�ı6k�                                                6�7S3��7���1�f�    1H<�/�s'7K���K��6��o5N�0!*�.���-�        )�Mt0�    B��@Q�5                @Qе?	�                ?�G�    5J�0��@SM�                                                                            @sӲAp�A�w�B��B5{<Bh�B�C�B��<�<�<�<�<�<�<�<�<�<�<�<�C���C���C��DC���C���C��C�)�C�IHC�hC��,C���C��C���C��C�N_C���C���C���C��C�DC�l�C�~�C�jFC�]JC�W6>��G>��M>���>�v�>�i�>���>��>���>��>��~>�T;>�=�>�S9>��m>�q�>�s>���>���>�\�>��?E�=?@��?B=!?C�z?FF�?J�A?R��?`�D                                                B�F� @¦�    @¶     9N8k̩7��f7!�?6rUH5�a�5Lj5 8�                                                                    8=��7��6���6L~�5�q5�4��4!��                                                                    4�b4 �o3O;+2��2=�1`)�0��0���                                                                    3�33"�q2�� 1�.1'
�0���0�$/��)                                                                    6��3ٕv8N�7�)79Z 9�8���7�J>7\'6��#6\��6��                                                6���3��7���4t�T/X9�1�na2J���^��7^��8��7���2��1���/a�        ,��2��4    B��@Q�5                @Qе?	�                ?�G�    7c�q1(�:K�                                                                            ABP�A�
�B��B-�[Bm^�B���B���B�Ph<�<�<�<�<�<�<�<�<�<�<�<�C�n�C�f�C�[�C�SbC�NoC�N6C�THC�a�C�t�C��<C��FC��C��$C�C�/.C�_C���C��9C���C��C�OkC�u�C�j�C�]�C�WU>��>��U>��>�R�>��>��C>��>��>���>���>�"�>���>�^�>�\�>��%>�_�>�8�>�O>��>�Z~?s
�?s��?u�S?wu�?zk?}>�?�d?�                                                  C5F( @¶     @��     9PON8�Y7�/�7��6]L�5��Z5	��4�P                                                                    8��g7��i6���6I�k5�Ė4�#d4-��3���                                                                    4�V;4+��3Ya~2�xr1�p14��0��0(i�                                                                    4��3X�2�K1�b@1�o0dC�/��0/T��                                                                    6���3�_*8�^�7��/9\b9BF^8���7���7I�6��f6yZ5��                                                6��l3=�7�GV4���/�m1���2�h%��/�7�/�8Q8(?3�zT2�� ,R�        /�3̓
    B��@Q�5                @Qе?	�                ?�G�    7}�)�ڙ                                                                                AZIIA�9XB�B7ǠBv�B�
�B�CB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C�	�C��1C���C�['C��C�ٙC���C�H�C��C�ڪC���C��cC�o:C�_�C�]�C�kTC���C��C�գC���C�43C�i�C�j�C�]�C�Wv>ܼ�>��'>�c�>�iq>�^�>��~>�>��>��>��'>���>�ۉ>�Z>�s�>�Ry>�8>��>��>��>��?|�!?}�3?��?��?�  ?�  ?�  ?�                                                  C�F� @��     @�Ԁ    9]�h8�y�8�7B�Y6r/@5�Ǩ5U34���                                                                    8��O7�V/7%��6u�5��l4��@4*�\3��D                                                                    4���4:o3��2�U2)�120��40&�5                                                                    4�3k �2���2�1&�b0a�E/��O/R~�                                                                    6�U3٘�8�l27��'9k&9S.&8�;�8��7\)�6��16�^5���                                                6���2�f�7��<4{��.(m�1���2]շ�%=7�%S89�=8!�6��f5)�N0�d�        5E�(6V�+    B��@Q�5                @Qе?	�                ?�G�    7e��*;�^                                                                                @���A���Bo{B/]Bs�sB�
�B�CB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C��1C�GrC���C���C�z�C�1aC���C���C�K�C��C���C��C�`C�.uC��C��C��AC�� C�ݐC��pC� qC�\�C�i^C�^7C�W�>��>�1%>��>�>�.0>��>؁>��(>�G>�9�>�N�>§F>�`j>�a�>���>���>��>���>���>���?n?p�.?vu?y�?}��?�t?�  ?�                                                  C�F @�Ԁ    @��    9hn�8�3�8F�7��L6Ӎ�6$u�5]�k4�E�                                                                    8�̎7�Հ7<�{6��6��5O�*4�5�3Ʃ                                                                    4���4Aco3��W3pB2f��1��'0�N�0+��                                                                    4 6	3tG�2�ɲ21ea1��}0�0	?/X�i                                                                    6��3��o8�0k7�(9w��9\�?8�8?��7��s7*��6})L5��                                                6��Y2�-�7��y4�J�.\�1�(#2>���S{N7S{�8,�8%0u6��!59��1�ag        5B�J6l��    B��@Q�5                @Qе?	�                ?�G�    7�a+�s                                                                                @�_ A�r�A�9B̯BOj�B���B���B�K<�<�<�<�<�<�<�<�<�<�<�<�C��OC�CC��RC���C�F�C��PC���C�C���C�f"C��C���C�cbC��C��HC�{MC�DtC��C�C��C�OC�Q�C�gVC�^jC�W�?�;?�X?
��?�??h�>��>�$�>��2>���>Ꮗ>�V>�Z�>��#>�i�>�Z�>�#>���>�}�>�nm>�]B?T9?U��?XѾ?]��?ehf?q��?|F�?�                                                  DaF� @��    @��     9?�L8�$�7�G7^J6��x6ARL5�h5�C$                                                                    8r'7�dO7!:�6�d�5���5t25�54�ۋ                                                                    4�8�4"��3�M�2�\2Hn�1�L1��1=�                                                                    4#�3M�w2��?2:U1}-�1D�0��P0n�G                                                                    6�8�3��_8}Z�7�a9MM9:�08��8#�B7�~7B��7�6���                                                6�o2��7��4��.pY 1��2M�-����6���8�8+�67�4�*1�u&        4���6�r    B��@Q�5                @Qе?	�                ?�G�    5�U*��C                                                                                @k��At�jA�/ZB3�B9ƜBp$�B�	�B�cH<�<�<�<�<�<�<�<�<�<�<�<�C�r�C�SC�9C� �C�
C��#C���C��rC�GC��C�˚C��cC�4XC��C���C�'3C���C��=C�]�C�@cC�0eC�K1C�d�C�^�C�W�?Y�?!�?�w?�? � >�T3>�ؙ>���>��>���>��->�|t>��>��>�ϖ>��g>�~�>�9�>�2G>�m?BX�?C+�?D��?G��?L�?S?^��?s��                                                D�F @��     @��    9u 8}(:7��A7A�96��A6,{5چ�5�~                                                                    8Ik_7��7
�6t��5�X�5Yc�5
(4�[�                                                                    4��4
%%3o�02�b�20��1��<1n�,1:�                                                                    3��3.�2��G2��1_	�0�HL0���0k�                                                                    6���3�=�82��7�Rm9+^�9��8��8�7�=e7)�-6錳6�&�                                                6�0�3	n�7�_�4�d�.G��1�Z�2+b�6>�A�>��7�57�ҟ5W5�4�00x�        3`��53��    B��@Q�5                @Qе?	�                ?�G�                                                                                            @e)�AoǄA�۶B2�B6͋Bk��B�DB���<�<�<�<�<�<�<�<�<�<�<�<�C�`C�~�C��xC���C���C���C��fC���C���C���C���C�pC�F�C�'C�қC���C�;�C���C��1C���C�S5C�LC�b|C�^�C�W�>�F>��>壖>��>��>�S�>�>�c�>�n>�S>��>�.�>�?�>ۛ�>�)�>�2�>��>�b�>�au>�`&??ZY?@*�?AҀ?DH'?H*�?Niy?XVL?i��                                                E)F� @��    @��    9�38Q�7���7$=�6��)6*5�\�5�q�                                                                    8%@7���6�4�6Ov(5�V�5=��4��Q4���                                                                    4�� 3���3I�2�C"2��1�̍1T'1'y�                                                                    3�0 3�g2~�P1�o�1?a�0��o0���0S��                                                                    6�%4/�7���7��9��9�8�[7��7r
�7�6��56���                                                6�m�3E��7��|3���-�d�1�V�1Is�7������h6,�6��M3R�2�-�I        /�^�3/w    B��@Q�5                @Qе?	�                ?�G�    68�G,=��>&�                                                                            @z�{A~~�A��]B��B>i]Bv%�B��B���<�<�<�<�<�<�<�<�<�<�<�<�C�h3C���C��eC��C�P�C��dC��C��C�6�C�[�C�{�C��(C���C���C��8C�}C�S�C� {C��>C���C�{�C�TKC�`�C�^{C�X>�. >���>��e>��>��p>�G)>��>���>��.>�[�>�k$>���>��W>�A>�;�>�j�>ήV>�[�>��(>ļ�?H�U?H�w?J��?M@�?Q8Y?X�?c �?v�                                                E�F @��    @�!     8�E8n��7���7/�X6~�J5��85:�R4���                                                                    8��7��U7�F6^5���5 �14k��3�_                                                                    4xZ�46;3a�x2�߅2
��1^k20�ۿ0A�                                                                    3��3$zf2���1�]�1/xa0�y�0 ��/s�                                                                    6��3��7��7�o�8�h�9$T8�~8�u7_��6�"�6E��5���                                                6�A�3�/7�m�27I,tq]1j?G0XR7Ћ��Ћ���o�5t%0��R/��<)�        +�4�0�l=    B��@Q�5                @Qе?	�                ?�G�    7{�v*���@%��                                                                            AV��A��B�=B6��Bv�B�
�B�CB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C��6C�
�C�7C�`�C��IC�� C���C� �C�MC�s�C���C��AC��#C��}C��C�C��C��C��)C��/C���C�aC�`3C�^bC�X?>��x>���>���>�u>�[�>���>���>��#>�E�>��3>��1>��>��>���>�MN>�'>��>�
�>�U<>�0,?|�?|� ?~I�?C?�+?�  ?�  ?�                                                  E�F� @�!     @�0     8:�A7���7�|i7k�6Y�=5��5
?4���                                                                    7k�z7��6���66m`5���4֫�4APO3�us                                                                    3˴3�;�3�p2��Q1��H19�50�N0>J                                                                    3 ��2��u29j=1�#1Iq0jT�/��/p��                                                                    6��3��c6�v�7F�u8I�|8�{�8F�7��p7@��6�=z6%�5��W                                                6���3�R7���        1*�o    7������`�                                    B��@Q�5                @Qе?	�                ?�G�    7C��*��>A�%A+b]@|E�                                                                    A	.yA��iB/(B%MvBb3�B�"B���B���<�<�<�<�<�<�<�<�<�<�<�<�C���C��yC��kC�$C�:cC�sC��iC��C�0�C�i[C��C�ޙC��C�MdC���C���C��wC��C��C���C��C�nC�`�C�^LC�X]>��>>��}>�3>�F�>���>���>���>�z�>��a>��9>�N2>���>��P>�;,>�Vl>��>Ŷs>Ɛ>�~'>ũn?SB?H?d�&?mm�?pu�?t�?{)�?�                                                  h�F� @�0     @�?�    