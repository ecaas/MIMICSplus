CDF      
      time       levdcmp       lndgrid       levsoi        natpft        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Fri Nov  4 12:00:21 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,HR,HR_vr,NPP,NPP_GROWTH,NPP_NACTIVE,NPP_NACTIVE_NO3,NPP_NACTIVE_NH4,NPP_NAM,NPP_NECM,NPP_NFIX,NECM,NAM,NEE,NEP,NACTIVE,LITFALL,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,T_SCALAR,NPP_NNONMYC,LEAFC_TO_LITTER_FUN,GPP,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN,SMIN_NO3_LEACHED 31539_Modum_historical.clm2.all.1898.nc -O 31539_Modum_historical.clm2.all.1898.nc
Sun Jan  9 16:25:46 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1898.nc /nird/home/ecaas/31539_Modum_historical/lnd/hist/31539_Modum_historical.clm2.all.1898.nc
created on 12/13/21 20:53:26       source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31539_Modum_hist_for_decomp    Surface_dataset       !surfdata_31539_Modum_simyr2000.nc      Initial_conditions_dataset        -31539_Modum_Spinup.clm2.r.1201-01-01-00000.nc      #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         9./31539_Modum_hist_for_decomp.clm2.h0.1850-02-01-00000.nc      Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:25:46 2022: Appended file /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1898.nc had following "history" attribute:
created on 12/13/21 20:53:26
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�                2j�|4<Q�4�:�4���                                                                                    1�D�3m�3�!�3�uL                                                                                    . 3�/ͱ&0-�U0	�                                                                                    -!�\/�%/[�C/E{                                                                    6��*3��y6mG~55�b                3X�5D��5���5�_                                                6�.d3Ⱦ7�ȳ        1��E    7�H��H� k"                                    B��@Q�5                @Qе?	�                ?�G�            A�<	A֝YA�MA��BރA�G�@��6                                                    =tIl@M�+A��AGpDA��>B�Bw��B�8�<�<�<�<�<�<�<�<�<�<�<�<�C���C��XC��KC��HC�"�C�nyC��xC���C�>HC�|�C���C��C�J�C���C��#C� �C�Z�C��C���C���C���C�l�C�\TC�^�C�[]>�k�>�]�>��D>��>�-�>���>��Z>�Y|>��4>�1�>��.>�o>�@�>�L{>���>��3>�l>���>� m>Ù�                <���>�]?2??9�v                                                !�iF� @�     @�#�                        3��!4���4���                                                                                        3�4��3�                                                                                        /�u!0�	0g1                                                                                        .�^/���/?>�                                                                    6���3�$6��M5Q�&                    5
�j6+35�LV                                                6��F3
	B7��7        1�qh    6��ᶝ�ᶃj                                    B��@Q�5                @Qе?	�                ?�G�            A��A֌GA�A�gBL�B�=�@N                                                    >��@Nd~A�6AFk�A���B	�Bu��B�|�<�<�<�<�<�<�<�<�<�<�<�<�C�b�C�g C�nC�vmC���C���C���C���C�C�=�C�n�C���C���C�KC�c\C��^C��C�*�C�W�C�t�C���C�rcC�]�C�^�C�[m>�%�>�[�>��F>��>���>�B!>��*>�Op>���>��>�t3>�B�>�N_>���>��t>���>��>�">��V>���                    >�.?0e�?5��                                                !��F�V @�#�    @�*�                        4k^N4��4�R�                                                                                        3��Z4�3�h\                                                                                        0 ��0�-�0�B                                                                                        /"^&/��R/7�S                                                                    6���3��77���5hX�                    5�n6�"5�ѧ                                                6�� 3��7���2�D    1��A0�ɵ���,6��,7W"5�q0�0�/� �,$_�        *O�0��Z    B��@Q�5                @Qе?	�                ?�G�            A���A�t�AꤩA�rB��A�G                                                        >��d@O5�A��AF��A��jB/�nBko�B���<�<�<�<�<�<�<�<�<�<�<�<�C�C���C���C���C��=C��bC���C�ݖC�4C�*TC�R�C�4C���C��C�!�C�b�C��AC���C�4C�@�C�iC�qAC�_�C�^�C�[}>��>��9>���>��>�!B>�l�>��o>��>�<>��>�Y>�L�>�Β>���>���>���>�I>��H>���>�o                    ?	 ?%��?*�                                                !�1F�� @�*�    @�2�    8�*7��a6�w�6��3���4�<�4�KM4|�                                                                    7�k
7{f6��5'� 3!��3�$,4/�3�3                                                                    4Lh�3��32�9�1���/��m0-�@0y\�0	��                                                                    3��2�,v1�'0�$t.�j8/[�l/�}�/-�                                                                    6�B3�M8+�7�8��8��7�6�[�5	85��L6�	5��@                                                6��"3w^7���4,�    1�j,2��pE7�pE7�-7eh�2�!g1Q��,��        ,z��2g��    B��@Q�5                @Qе?	�                ?�G�            @��IAB�VA��A���B	�A!�                                                        A,�FA�0�A�� A���A�׋BR�BtuB��K<�<�<�<�<�<�<�<�<�<�<�<�C�xC�*BC��kC���C��/C��3C��kC��C��cC��C�=dC�dJC��FC���C��%C�/�C�mC��nC��	C�iC�F�C�jC�`�C�^�C�[�>��n>�H>�3}>��>�v�>�w>��y>���>���>��>�� >���>�9>��$>��k>���>�(�>��y>���>��F?-��?i�>�|�>��>t?1ĵ?/2 ?2o�                                                !��F�� @�2�    @�:     92�B8���7�e27;8�6�j-6g�5���5 �                                                                    8a�7��d72p6l}�5�5F5??w4�s[4=y�                                                                    4�2�4M3p�R2�y�2!b1�_1:N0���                                                                    3��37Mo2�w2$�1K`0���0kU/��                                                                    6�Z3�J�8�΁7���9>x�9$�n8�.�8 I�7�71��6��63ְ                                                6��2�Ŧ7���4�q�)TYO2�t2����N7�Nd8J��8!۲5�54�P�16cr        4�<5�&�    B��@Q�5                @Qе?	�                ?�G�    6��,�x:<��8� �        >���                                                            @z�A���A�̿B	C�B?�Bx��B��ABά�<�<�<�<�<�<�<�<�<�<�<�<�C���C�3�C��_C�?C�ʹC�c�C���C��WC�_0C�.eC�=C��!C��"C��&C��sC�RC�K�C���C��vC���C�%�C�^[C�`�C�^�C�[�>�_�>�9�>֌�>���>Ǔ#>��$>�u:>���>��Q>�	>�Ն>�`#>��0>�Ւ>��2>���>�OR>�F1>�S�>�/�?Hs?I�_?K�?OD?S�|?[�<?i��?{��                                                !��F� @�:     @�A�    91�p8��K7�O87DZ�6�`@6"�K5�J<5��                                                                    8`՛7�by72	6x�5�Q(5M��4�5Y4�w�                                                                    4�d}4P3yW�2�s2.�1�=1W~�1 ^]                                                                    3��i3:U2�z�2q1[�P0��0�,0J�%                                                                    6�E.3�`�8���7�U9=�9'c�8�г8"�7��w77}�6�+�6�b                                                6��2�a7�E4T��-PS�1ޥZ2HO��Y7�Y�8C��8*3�6�}�5Wz2?��        5�'6�.p    B��@Q�5                @Qе?	�                ?�G�                                                                                            @Z�Af�A��"A��cB1E+Bc��B�B�)<�<�<�<�<�<�<�<�<�<�<�<�C��C�kC�#�C��
C���C�YzC��C��{C�^RC�C��uC���C�D2C��C��"C��AC��C��#C��-C��!C�
�C�P=C�_�C�^�C�[�>�~<>��Z>���>�ߧ>���>��O>�޷>��(>�l�>��s>�B�>��>���>��>>��6>�l*>�n�>���>��>���?:,�?:f�?;ǐ?=�?@�-?E��?Mq�?Y��                                                !�]F�J @�A�    @�I@    9N��8���8<(7`�Z6�J�68w�5ޓ�5�9�                                                                    8��7�� 7%�?6��+5���5i>5�M4���                                                                    4���4*�3�S)2�?�2F1�{�1s�15��                                                                    4��3Vߔ2�
�2� 1zG0��p0���0eT�                                                                    6���3�V�8�k�7�6	9\��9AȐ8�׶84�o7�l�7H�6��76��                                                6��92��7�*4�}0.C�1��2,�U�Q�7Q��8+�|8%,/6��)5V.2QK(        5]�]6�[>    B��@Q�5                @Qе?	�                ?�G�    3��1(R�>                                                                                @c��AmK+A��B ��B4M�Bg�}B���B�]�<�<�<�<�<�<�<�<�<�<�<�<�C�k C���C���C�OC���C���C�A3C��{C�|=C�+�C�ؾC��C�5�C��C���C�RaC�|C��dC��C��pC��C�B�C�^C�^�C�[�?��?P?�?��? W�>�a#>�ı>��>�#�>�h�>ר�>�T>�͂>Ǭ�>���>��[>�<�>�>�,e>�V.?>_�?>l�??��?A��?D��?I�?R5�?`��                                                !��F�� @�I@    @�Q     9H8�b8c�7h��6�Tj6F�X5�� 5� 5                                                                    8|¡7��|7(~	6��?5��5{�5&4�6(                                                                    4ډD4)h3��2��M2P��1��1�(�1L>�                                                                    4
�3U��2��2 f�1���1	�0�v�0��2                                                                    6�ܐ3�\�8��	7�c9V�j9A�8���86Hz7��y7O��7)26��$                                                6�2�D�7�i�4��X.u~�1�[j2L���U)Z6U-�828��6L�m5��1�/        4�c@6+"�    B��@Q�5                @Qе?	�                ?�G�    5���*9>I                                                                                @n�9AvF9A�: B��B:)BpF�B��pB��Q<�<�<�<�<�<�<�<�<�<�<�<�C���C��[C�{oC�]�C�;*C��C�ձC��-C�K�C�HC���C�k{C��C��C�\C���C���C�h�C�9�C��C�)C�:�C�[�C�^]C�[�?	��?X�?��?l)?�T? ��>�\I>���>�z+>�Y�>��>���>��>ը?>�B�>�8>���>���>��>�w�?C�?D$?E��?Hh�?Lt&?S<�?^��?s�l                                                !�%F�� @�Q     @�X�    9-�8�$7�u�7U}F6���6-�n5��5ML                                                                    8Z��7�7ފ6���5��,5[6�4Ֆ>4��T                                                                    4��/4�%3�N%2�)r2D8e1��$18�0�:�                                                                    3�q3=�2���2B�1w�s0�n�0iKc0�8                                                                    6�e�3��8/A7�w59:r�9+�w8���8$��7�0�70��6���6T�                                                6���3}7���4�-�.;��1�}2'H76��=���7��7�S5s�4 	0m�        3�.�5J�    B��@Q�5                @Qе?	�                ?�G�    6洭+~��                                                                                @�n+A�o�A���B�BH��B���B��!B�h)<�<�<�<�<�<�<�<�<�<�<�<�C�5C�LhC�_'C�n-C�{zC���C��C��4C��vC�x�C�d}C�G�C�"RC��C���C�i�C�2C���C��C�b�C�8#C�:>C�Y	C�^&C�[�>�T3>�)�>�9>�3�>��>��>�V�>�Tk>���>���>�]�>�Q&>ܬ@>�H>�C>�$>�>�[>�e;>��?Q�N?Rcu?Ut?X��?^y�?h�_?t�?~�k                                                !��F� @�X�    @�`@    9(Ͽ8��7��7Vɭ6���5�!c5�v4��H                                                                    8U<O7�t�7��6���5�4��42�z3�L&                                                                    4�\�4_�3�y�2�{2'�16��0��0+vB                                                                    3��&3:'�2�Ϫ2'�1S�e0f�u/�q/X�`                                                                    6��:4{�7�Z�7��96fj9)މ8��8$C(7���6���6}W5�e@                                                6�@N3DA7�%�3���-�KY1�_1h4�7��k���M6 ��6�x�3"�1ܺ�,�3�        /�$�3<�    B��@Q�5                @Qе?	�                ?�G�    7Is,*/9ZXu                                                                            @���A���A��B��BiE�B��!B�CB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C�oC�%qC�@�C�X�C�pPC���C���C���C���C��gC��QC�ĘC���C��UC���C�g�C�8kC��C��_C��%C�a�C�A�C�W!C�]�C�[�>�=u>̘&>�K>�g�>��K>�\^>��f>�^>��>�Z�>�z'>�8X>Ռ>�_k>҉9>�a>���>�t�>�B>��?a�,?c�7?h�?n,r?vȝ?~��?�  ?�                                                  !��F�@ @�`@    @�h     9�28d�<7ͫ+7.�6h��5��	5�%4�2                                                                    8*�7���7�_6\��5��4��47�J3���                                                                    4��3���3`�n2��1�!�1;��0��0-�<                                                                    3��K3�2�ܱ1�[.1 ��0mw/���/[q�                                                                    6���3�V�7�7�P|9��9c�8���8T�7T��6���6!\5��Q                                                6�#C3��7��w1u�Z    1�9/8��7ͼ}�ͼu�ɯ�4��=/�N�.�	(Jܜ        *���/�m�    B��@Q�5                @Qе?	�                ?�G�    7V��-�?��Q=,�                                                                        @�ŦA�|NA���B'q�Bpy�B���B�CB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C��4C���C�&�C�^�C��uC���C��C�L~C��|C���C�՜C��TC��C�*�C�4�C�2"C�!{C�C��C��xC���C�N"C�VBC�]�C�[�>� �>�.�>���>���>���>��>��>�A�>�e�>��>Ƨ�>��>ʳ�>�j>̳>́r>�l�>ɘ�>�I�>��l?Y��?h9:?mL�?sr�?z��?)x?�  ?�                                                  !�QF�| @�h     @�o�    85��7ƣ�7���7M'6]�E5�� 5
�i4�L�                                                                    7e�k6��6���6*�25�o4ɾ�4/0�3�nD                                                                    3ƺL3X�3��2��1�7�1.rP0�}�0*��                                                                    2�E2�2FK�1��1��0\Z�/�[Y/W�J                                                                    6�5�3�D6d�_7AFl8G �8{�H8\ٱ7β�7Kn�6�3�6��5���                                                6�x�3�|7���        1iM�    7�!Z��!Z��H                                    B��@Q�5                @Qе?	�                ?�G�    6勤,�T�Ag\rA<y=>9��                                                                    @d5Ah��A�l�B��BL[zB��fB���B�k<�<�<�<�<�<�<�<�<�<�<�<�C�q�C���C�� C��C�2VC�o�C��TC��_C�<NC�xQC���C��zC�+�C�bC���C���C��5C���C�ωC���C��!C�[�C�V�C�]=C�[�>��>��k>��>��>�]>���>��>��+>�z�>���>�U�>��;>�7�>�}�>2>��>�S>ƾ
>�6>��I>��.?1�?RX�?W�*?_j?i?s�?}�                                                !�F�� @�o�    @�w@    