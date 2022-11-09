CDF      
      time       levdcmp       lndgrid       levsoi        natpft        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Fri Nov  4 12:00:23 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,HR,HR_vr,NPP,NPP_GROWTH,NPP_NACTIVE,NPP_NACTIVE_NO3,NPP_NACTIVE_NH4,NPP_NAM,NPP_NECM,NPP_NFIX,NECM,NAM,NEE,NEP,NACTIVE,LITFALL,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,T_SCALAR,NPP_NNONMYC,LEAFC_TO_LITTER_FUN,GPP,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN,SMIN_NO3_LEACHED 31539_Modum_historical.clm2.all.1920.nc -O 31539_Modum_historical.clm2.all.1920.nc
Sun Jan  9 16:25:48 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1920.nc /nird/home/ecaas/31539_Modum_historical/lnd/hist/31539_Modum_historical.clm2.all.1920.nc
created on 12/14/21 09:36:13       source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31539_Modum_hist_for_decomp    Surface_dataset       !surfdata_31539_Modum_simyr2000.nc      Initial_conditions_dataset        -31539_Modum_Spinup.clm2.r.1201-01-01-00000.nc      #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         9./31539_Modum_hist_for_decomp.clm2.h0.1901-02-01-00000.nc      Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:25:48 2022: Appended file /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1920.nc had following "history" attribute:
created on 12/14/21 09:36:13
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�5���5��63��6��6)�25��u4���4�\�                                                                    4��4�z5b��6 F5Vi�4�t�4�i3�Z                                                                    1\�1�o1�Ҷ2]u�1��1 ��0�F30{i                                                                    0���0I��0�[1��|0��"0KU9/�f&/F�                                                                    6�Kj3�4f6[}�6�F�6<�M6[ac7	K7��47!�>6��~6�5��0                                                6�\�3�Z7�Mu        1�]�    7^�˷^�˷7g                                    B��@Q�5                @Qе?	�                ?�G�            A��A�ȦA��                                                                    ?.?4@���A���A��B_�B>�IBb��B�(�<�<�<�<�<�<�<�<�<�<�<�<�C�|�C��qC��xC���C�ԊC���C�0uC�h�C���C�ǧC��tC�*'C�^�C���C��JC�C�A�C�laC���C��3C���C�l�C�XBC�[�C�\{>�a>��	>�n>��7>��j>��>�@�>�#�>���>�,>��>�r�>�Y�>���>��@>�`2>��;>�%>��>¤�=]r�=�zQ?R�?�8?��?��?��? zy                                                $��F�� @��    @��@    7&:�6��7.�?6ѱ6-�5���4��m4�j?                                                                    6Q�05��/6\�]6o�5[p�4�M4^`3�]�                                                                    2�CA2n�2��2d��1�m11ɗ0�0,�                                                                    1���1C�1�ݮ1�h	0�F�0I�8/���/A{�                                                                    6�)-3�n7��6��n7F�~79��8_7�8j7&�c6���6��5��[                                                6��3(�7�p�        1���    75zŷ5zŶ~P7                                    B��@Q�5                @Qе?	�                ?�G�            A���A���?���                                                                    ?̂pA�A�V[AװlBZ<B@s�Bc�YB�2R<�<�<�<�<�<�<�<�<�<�<�<�C���C���C���C���C��*C��YC��C�K(C�t0C��	C���C���C��C�J�C���C���C��`C�&+C�N(C�i�C�~�C�pfC�Z-C�[�C�\w>�o�>�r#>�B>�e�>�Ȕ>�m�>�Xz>��>��i>��$>���>�8>�{�>�A�>�a:>���>��n>���>�@�>���>�=>�zl?��?�a?�G?�S?ޗ? �                                                $�-F� @��@    @�@    8�%7U$B7��g6�i.6=��5�bm5�N4���                                                                    7:�!6���6���6'|5o{L4��+43M3ֲ                                                                    3�F�2�i�3�2|R1ι1+��0��@09h�                                                                    2˷�2�b26׮1�\E1��0Xɯ/Í�/j3k                                                                    6�A�3�1v7��k7'��8#?�8��8^3�7�H�76�6�Al6&ޑ5���                                                6��3	�t7�/3)��*��01ݨ�1I��\l2�^�7'�/6b+�1M�0c�+ 0<        *�01%�    B��@Q�5                @Qе?	�                ?�G�    5LV,GnvAI�'@�'9                                                                        @��AO>[A��zA�VkB!!'BLBr��B��P<�<�<�<�<�<�<�<�<�<�<�<�C���C��YC���C��.C��LC���C�C�=RC�`'C�C��mC���C���C��C�NC�� C��hC���C�C�>�C�cvC�nC�[�C�[�C�\s>�LI>�c�>�3->�(>�#W>�q�>�V>��*>��?>�J�>�s>��>�8	>��^>�m�>�s�>���>��,>�2}>�Ux?\�?)F?)?(Վ?(�?)�O?+�?->P                                                $��F�P @�@    @�
     9�8w�t7��71R6��5���5Fſ56&                                                                    8H��7���7��6_��5�Z�4��4{�4f	�                                                                    4�~�4,�3hU�2�V�2:$1D׮0��l0ƨ                                                                    3�'#3*��2���1�71@Iq0x��0�/��?                                                                    6�ݿ3։�80Z�7���9*�b9��8�%�8�H7�6�6��B6i��6L�                                                6�M39*7�{42e/]
1���2_ȶ���6��>7�Ӄ7m�*2tBt1<g+�	        ,�]K2E$    B��@Q�5                @Qе?	�                ?�G�    7G�w-���<��I                                                                            @���A���A��B �qBd��B���B�CB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C��C���C��TC�2C�a�C�D�C�*C��C��C��	C���C� C��C�6C�=C�dRC���C���C��C��C�FOC�f�C�]	C�[�C�\o>�[�>�q�>�(�>�Z�>��`>�i>���>�J�>��	>��>��i>�,T>��>��:>�}N>��>�=>���>���>���?aɺ?c�?h4�?m�?uNs?}�u?�  ?�                                                  $��FȌ @�
     @��    9F_^8��R8 4�7Y�F6�et5�AX4� �4���                                                                    8z�i7�ph7!��6��:5�!�4�R�4�a3�K                                                                    4�L�4%@p3�ɖ2�u52X�1:��0��0ԧ                                                                    4�U3P�/2���2�/1E}�0k�/��/C�Z                                                                    6�+�3��8�N7��9U�T9=�n8�28*��7�?X6�6��5���                                                6�+�3�d7���5 ^�/��l2�~2����K-7�KO8Q��8*��3�d2�5�-�R|        /�L'3���    B��@Q�5                @Qе?	�                ?�G�    7KP-�/:��                                                                            @���A�sbA��B ��Bi͔B���B�CB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C���C�8C��(C��C�V�C��C���C�w�C�=�C�GC���C�øC���C��bC���C���C���C��@C��oC���C�,QC�\bC�]VC�[�C�\l>�a,>��D>ٗc>�rH>υ�>�b>Ť�>��>�vv>���>�)z>�%�>��>���>��>�:l>�>��
>�g=>�h�?b�?dH?h��?n�?w%?~,�?�  ?�                                                  $�YF�� @��    @�@    9Q�18�7 8�K7fG�6��6Bk�5�^F52z�                                                                    8���7�57(�6�p�5��C5u�k4�4ar�                                                                    4�֝4-��3��,2�72N(�1��1C�p0¯                                                                    4��3[H�2��R2�>1�4�1�l0w�/��                                                                    6�2�3��8�^7�k9b�9G��8�]}81
�7�%�7Q��6�W6L�                                                6�0�2�67��4�+�.nj)2 u2^����vj7�v�8H��80q
6��s57�2��        546i^�    B��@Q�5                @Qе?	�                ?�G�    6�|�+c�J                                                                                @}@A���A�GOB	�cBBuPB%�B��;B�˗<�<�<�<�<�<�<�<�<�<�<�<�C��+C�Z�C��C���C�j-C��C���C�T_C��C���C�^|C�5C���C��qC�B�C��C���C��C��C��0C��C�QC�\�C�[�C�\i?
#?��? ��>�2J>��w>�>=>�V~>�J�>�,�>�γ>�y�>�c\>ż�>�kI>�� >��_>�׃>��>�^�>�P�?I�
?K�?M`�?Q X?V�?`��?ol�?~h�                                                $��F� @�@    @� �    9V5G8�f8
�7pW�6Š�6N$|6 ܌5�8t                                                                    8�J7�>f7/U$6�ˁ5��l5�1�5"ž4���                                                                    4��42�A3�V�3�2Ww�1��1��1g�                                                                    4��3a��2�)�2%1��1�c0��l08��                                                                    6�3�L�8�67�� 9h�9N�W8�%�85��7���7U�M7,�6�.�                                                6�1�2�I�7��X4�T.; �1�/2Y
L�o9�7o:�86��8-��6�7k5N��2*�        5Bv66�+    B��@Q�5                @Qе?	�                ?�G�    6{�^*��0                                                                                @x\�A}U�A�.Bx<B>6bBv_�B�ZbBΐb<�<�<�<�<�<�<�<�<�<�<�<�C�/�C��C��lC��?C�I�C�sC���C�h�C��C���C�xC�!�C��/C�qC�pC��C�{�C�G�C�(C��C� GC�G�C�[^C�[�C�\g?�8?
��?#�?�?w�? ��>�h�>���>�]>���>�ʙ>ܯ�>ִ�>а�>ʿ�>�eb>��>��i>��.>�=?GϞ?Hf|?Jq"?M,�?Q�4?Yz ?gn?|��                                                $�!F�D @� �    @�(�    9D28�"x8b7c/26�w�6Fq�5�x�5��                                                                    8w�27�+�7$DF6�|5��5z��51f4��?                                                                    4�ǡ4%��3���2���2Mz�1�`�1��o1��                                                                    4�3Q�{2��2o1�ƨ1��0��V0FK�                                                                    6�]�3� ]8��a7�39U<9@K�8�e�8)��7�y7G��74�6���                                                6�`�2��7�V�4�5�.A"1�2Y���Ϡ6��}8d�8��6V�5�}1ŭL        4ݭ63(�    B��@Q�5                @Qе?	�                ?�G�    6��y+ �                                                                                @x)�A}4dA���B`�B>Bu��B���B�6f<�<�<�<�<�<�<�<�<�<�<�<�C�3EC�fC��5C��C��C���C�|�C�M�C��C���C���C�tC�?pC���C���C�RC� C���C�~?C�WZC�9?C�DdC�Y�C�[�C�\d?�$?��? 9q>�Y�>�L>�o=>�J	>�2>�J>낐>�">�F^>޶�>ٳ9>�2�>Ε�>�J|>ĸH>�-o>��>?G�L?H\�?Jc4?M\?Q|?Y8?fs?}��                                                $��Fɂ @�(�    @�0@    94�8��17��7X@�6��06B�5�[@4��                                                                    8cd	7��7�
6���5��5v0�4���3�˄                                                                    4�FS4o�3�~�2��2F)u1Ԃ�1\�o0B�f                                                                    3��3C�2��2�1zOP17�0�`/vI�                                                                    6�5�3��8:3�7�Ǒ9D��93�8��s8 Z7��7@&66�s�5��                                                6�F"3��7���4��.EX2��270�6~E:�~Bb7��D7��}5I�_4Uc0(��        3?�5(a�    B��@Q�5                @Qе?	�                ?�G�    6ҍ*+�XQ                                                                                @�0�A� mA�LB�IBD��B��zB���B���<�<�<�<�<�<�<�<�<�<�<�<�C��C�ݒC��IC��C���C�ݣC��C��-C��C���C��	C�o�C�HVC�C�ܫC���C�PpC�	xC�ɗC���C�_8C�HgC�X5C�[�C�\a>�Pu>�5>ꋹ>ꂿ>�q�>�>%>��>���>��(>�Z~>�zL>�">�V�>��>��}>�@{>�y3>��m>��'>	?Mͫ?N��?Q?To�?ZJ?d ?t�6?�                                                  $��Fɾ @�0@    @�7�    9|�8VH�7���7+��6���6!��5ɸ�4�o�                                                                    8&�7�VZ6��q6X��5�Z�5LA�4�Ϋ3�/	                                                                    4�\.3�3O!�2�*2 �1�P1[�+01-�                                                                    3�3��2���1�k1K+�0޶0��/_�%                                                                    6� n4��7��7���9[�9i�8���7�";7�χ75�6�{"5��6                                                6�#�3H��7�z�3�=	-�<1��N1s�\7r���r��6��`6���3N�2׊-�X        /�xX3+��    B��@Q�5                @Qе?	�                ?�G�    6�Z/,�I�?.,�                                                                            @|�5A��lA�B
3|BB�_B}ïB��uB��v<�<�<�<�<�<�<�<�<�<�<�<�C�.8C��7C�ݤC�&^C�m�C��&C��cC�@PC�o�C���C��C��C��HC��C��oC���C�h�C�3C��=C��C��pC�R�C�W�C�[�C�\^>���>��>���>���>��L>�Yv>ɶc>���>�ʜ>��>���>�1�>ַ->�x>�B�>�c>��>̎�>��\>Ůc?I_?K�??N8�?Q�U?WH�?`�^?q>�?�                                                  $�MF�� @�7�    @�?�    8���8B��7��7ƴ6�86I�5�N,4��5                                                                    85�7u��6׍�6@��5�XG52x74��3��C                                                                    4�at3��3:�2���2��1��1J�0//�                                                                    3��]3�2kc1�d�13,�00B�/]I�                                                                    6��'3��%7�H7�YR9Q�8���8��P7� 7e117
�k6��5���                                                6�g3��7��1j�},�1��E/%�]7ǭ#�ǭ��N�4��(0,e.� g*K��        *��/��=    B��@Q�5                @Qе?	�                ?�G�    6��K-�,�?���                                                                            @x�A���A׻�B	InBA#�B{)�B��gB�|<�<�<�<�<�<�<�<�<�<�<�<�C�P�C��C��wC��1C�WC�CC�}�C���C��&C�'gC�[�C��C���C��aC��C�%�C�+C�C�gC���C���C�a2C�W�C�[�C�\Z>�b�>���>�y>�p}>��(>���>�_>���>���>��>��>�;`>�-�>��Z>�B�>˹�>�>�5~>�n(>��?G�?J��?L��?P6/?Ui�?^K~?m��?�                                                  $��F�8 @�?�    @�G     8C��7a�6���6���6Y�%5֧�5zo�4���                                                                    7wk�6�+A6{�5��V5��5�f4�+�3ĒX                                                                    3Տ�2�m*2�l2�;1���1j�1�0)�                                                                    3�2�1��|17�1�0��0,}�/ViI                                                                    6��`3�&�6o�77	�g8Y�A8 g7��R7IZ�7AW�6���6���5�B#                                                6��3�7�*I        1��i    7��ַ��ַ�D                                    B��@Q�5                @Qе?	�                ?�G�    5�-�-6��AE�<A�|�A^�9A9q:�Q                                                            @A�dA�A���A�sMB5Q�Bj��B��AB�oR<�<�<�<�<�<�<�<�<�<�<�<�C�0+C�W�C���C���C��C�1-C�x�C��NC�2C�H"C���C��C� �C�;C�r�C��"C���C��C��?C���C���C�o4C�YnC�[nC�\U>��8>���>�o>�w�>��>�U>�z>��>���>�$�>���>�5>���>�"r>�}�>�e>Ō�>��6>���>�<�>���>���>ⴥ?��?F�?M/�?X�?g�                                                %uF�v @�G     @�N�    