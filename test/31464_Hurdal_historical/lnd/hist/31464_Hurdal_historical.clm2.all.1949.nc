CDF      
      time       levdcmp       lndgrid       natpft        levsoi        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Sun Nov 13 13:48:41 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,NPP_NACTIVE,LEAFC_TO_LITTER,QDRAI,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,nbedrock,T_SCALAR,NPP_NNONMYC,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1949.nc -O ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1949.nc
Sun Jan  9 16:23:30 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1949.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1949.nc
created on 12/10/21 17:01:00    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:30 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1949.nc had following "history" attribute:
created on 12/10/21 17:01:00
     NCO       `netCDF Operators version 5.1.1 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)          CWDC_TO_LITR2C_vr                         	long_name         .decomp. of coarse woody debris C to litter 2 C     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d      $<   CWDC_TO_LITR3C_vr                         	long_name         .decomp. of coarse woody debris C to litter 3 C     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d      $�   CWDN_TO_LITR2N_vr                         	long_name         .decomp. of coarse woody debris N to litter 2 N     units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d      %   CWDN_TO_LITR3N_vr                         	long_name         .decomp. of coarse woody debris N to litter 3 N     units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d      %h   FROOTC_TO_LITTER                   	long_name         fine root C litterfall     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             %�   FROOTN_TO_LITTER                   	long_name         fine root N litterfall     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown             %�   LEAFC_TO_LITTER                    	long_name         leaf C litterfall      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             %�   LEAFN_TO_LITTER                    	long_name         leaf N litterfall      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             %�   NDEP_TO_SMINN                      	long_name         *atmospheric N deposition to soil mineral N     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             %�   NPP_NACTIVE                    	long_name         Mycorrhizal N uptake used C    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             %�   NPP_NNONMYC                    	long_name         Non-mycorrhizal N uptake used C    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             %�   PCT_NAT_PFT                       	long_name         =% of each PFT on the natural vegetation (i.e., soil) landunit      units         %      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       <      %�   QDRAI                      	long_name         sub-surface drainage   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown             &$   SOILICE                       	long_name         4soil ice (natural vegetated and crop landunits only)   units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       P      &(   SOILLIQ                       	long_name         =soil liquid water (natural vegetated and crop landunits only)      units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       P      &x   TSOI                      	long_name         <soil temperature (natural vegetated and crop landunits only)   units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       d      &�   T_SCALAR                      	long_name         'temperature inhibition of decomposition    units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P      ',   W_SCALAR                      	long_name         .Moisture (dryness) inhibition of decomposition     units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P      '|   levdcmp                	long_name         2coordinate levels for soil decomposition variables     units         m      axis      Y         d      #�   levsoi                 	long_name         Dcoordinate soil levels (equivalent to top nlevsoi levels of levgrnd)   units         m      axis      Y         P      #�   mcdate                  	long_name         current date (YYYYMMDD)             '�   nbedrock               	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����            $8   time                standard_name         time   	long_name         time   bounds        time_bounds    units         days since 1850-01-01 00:00:00     calendar      noleap     axis      T               '�   time_bounds                    	long_name         history time interval endpoints             '�<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�   7S�~7��7�t�7*U?6e.�5x|4�%O3���                                                                    6��B6��#6�]�6W(P5��X4���3��;2�P�                                                                    2��3&3K�2�P�1�D�1��0#N^/MU�                                                                    2��2<f2���1�"1am0$�/NH
.���                                                                    6��3�A6A-31��                B��                    A�                @�p�    6��A��@3�                                                                        ?�]AAweA��*B$�B:��Bq�B�YB��<�<�<�<�<�<�<�<�<�<�<�<�C��;C���C��JC��KC��<C�^C�=1C�kFC��aC��C���C� C�S9C���C��*C�hC�7C�b-C�~�C���C��fC�O�C�&zC�"�C�!`>�>��>�t>�~u>�~>��?>��>�DW>���>�Į>�9>��>���>���>�<r>��O>���>�w�>�2�>��>�}h?D��?F��?IS�?M�?T��?`p�?w.�                                                )e�GF @��    @��    7a_n8~�7��'7'��6e�5dz_4>�v3�pK                                                                    6�W871w�6���6T 5��4�MJ3qAl2�=                                                                     2�=�3���3H�!2���1�h0���/�
�/AK                                                                    2�
2��n2}�1�;U1�O0�f.��.D�                                                                    6�353��m6>��3^r1׆�                B��                    A�                @�p�        A�Q#?+XE                                                                        ?�xgAr��A�� B�B8�=Bn�B��KB�J�<�<�<�<�<�<�<�<�<�<�<�<�C��QC���C���C���C��C�WC�*�C�QiC�vKC��C���C��:C�UC�E�C�{�C��C��uC��C�A�C�[{C�mC�ScC�)�C�"�C�!g>�ه>���>�$�>�a�>��T>�:�>��r>��=>��t>��n>���>���>�M�>��>���>�0�>�V�>�1E>�z>��>�f�?A�l?C��?F�?J��?Qְ?\�[?q�                                                )e�Gb @��    @�@    5�P6��6��Z7G�6[�L5]��46�%3�lB                                                                    5�@50M�5�M6?�5�ҁ4���3g!2�ٖ                                                                    1x�1��2J��2�B�1�el0�.�/���/�                                                                    0��p0���1��1�,21%0E�.�J#.8�                                                                    6���48n6u'�3azM2 /�0��0���        B��                    A�                @�p�        A���A���Ab                                                                    ?3ޱ@�3�A��iA� B+�~B\�YB�ɮB�kL<�<�<�<�<�<�<�<�<�<�<�<�C�h�C�|\C���C��lC��$C��C��C�+C�O�C�p�C�� C��?C��C�C�E�C�{XC���C���C�C�/�C�PPC�P�C�,eC�#C�!o>�t>�_�>�}�>�r>���>�E�>���>��q>��J>���>�rQ>��y>���>�-�>��>���>� w>��>�ne>�l�=Z�>��?(a�?5J�?9?>Tl?E�L?O��                                                )faG� @�@    @�     9�v8F��7�%�76��6���5��4�43Ƽ5                                                                    8<)�7{N�6�n6f��5��4��.3պ�2��                                                                    4���3�jM3'�2�L�2pT1(PR02'/QE�                                                                    3��3C#2S��1���1+�0T�Z/`��.�+�                                                                    6�;�3�;�6Q�Q3-��2#�2i!1��a        B��                    A�                @�p�    6�?@2�-@��@�Ð                                                                    @��A���A�R�B��B6��Bj�YB�2�B�w<�<�<�<�<�<�<�<�<�<�<�<�C��(C�?ZC�|C���C���C��:C��7C��uC��VC��C���C��C���C��uC��C�LnC��C��mC���C�/C�2PC�I�C�.�C�#rC�!x>���>�7�>��o>���>��q>�k>�:P>���>�c#>��b>�2=>�0w>��C>�[=>���>�V<>�8>�Q>��>��?K�?B?Aã?D��?G��?L�D?U��?gs2                                                )f�G� @�     @��    9M�~8���8
H7i��6�n�5�A4�@�4e�                                                                    8�'7�۠7.u�6���5�:�5#�4��3(�i                                                                    4��D4+��3�d�2�R�2/6$1l��0�<�/�x.                                                                    4��3X�)2���2��1]Q�0��</�	.�ob                                                                    6�1�3���66�
3|42B��3V��3,q�        B��                    A�                @�p�    6�x;� �                                                                            @m�
Auh-AϣBOJB9� Bo� B�kB�Y�<�<�<�<�<�<�<�<�<�<�<�<�C���C�jyC�?C��bC��uC�:�C��}C���C�R�C�@C��C��	C���C�y�C�lC�n�C��wC��0C��2C��C�(C�>�C�/�C�#�C�!�>���>�e�>�ý>מR>Ҩj>�l�>�%x>���>��c>�PL>�.�>��.>�ko>��>��>�@
>�@�>��>�  >��?C5H?C�?Eh;?G�D?K��?R��?]��?re                                                )g)G� @��    @��    9^�8���8�7Nq6�t�5��R5
�4h                                                                    8��`7�j67=��6�?5�x�5~�4/S�37>h                                                                    4ꉃ4:,Q3��Z3cw2@1�kz0��/��}                                                                    4 �3k*�2ǅ2)�1r�M0���/���.���                                                                    6���3ԻJ67}�2��2(�[6Z��66�m        B��                    A�                @�p�    6���                                                                                @ch�An�Aʾ�B�RB7)�Bm]B���B�X<�<�<�<�<�<�<�<�<�<�<�<�C�GC��AC���C�BMC���C��2C�SeC��5C���C�k�C�'C��C��C�dUC�)�C��mC��=C��HC��C��C��C�2�C�/�C�$4C�!�?u�>�!�>�ɞ>��V>�Z>�W?>�6�>�>��t>�N�>��>ǃ�>�z�>��b>��>�J>�u�>��>��>��?=�9?>�I?@]?C|?HYH?O�?[̞?l�D                                                )g�G� @��    @Ề    9C�8��8 �7[�6��B5��4��3�b�                                                                    8vp7�$x7"~c6�VJ5���5�y4��3r                                                                    4�`B4 �3�lV2��2$�-1_G�0u�2/v��                                                                    4�3J@O2��2�M1Pa�0�/�HV.��>                                                                    6�3ش06:��2��w236��G6�ľ        B��                    A�                @�p�                                                                                        @��A3}jA�B+A�ƔBD�BC��Bl,�B���<�<�<�<�<�<�<�<�<�<�<�<�C���C�)EC���C��bC�V�C��C���C�=�C���C���C�;�C��C��C�<�C��C���C�X�C�)C�4C� �C�`C�(�C�/C�$�C�!�?� ?��?	�m?��?��? ��>��f>�t>�]�>�$>�>،�>��>�:q>��J>���>��>�	d>�RY>��B?'�?��?9�?q�?;?��?&e-?-"                                                )g�G� @Ề    @�`    9+��8���7��7J�6�nq5�̈́4�|3�?                                                                    8Y<o7���7ff6�?5���4���4t�3	y                                                                    4�	4�c3y�2�r32�1O.
0cn�/d�<                                                                    3䯠37#�2�G�2��1Af�0�ٵ/��_.�\\                                                                    6�H3�s`6>�2�2�%6���6c��        B��                    A�                @�p�                                                                                        @fKA3HA���Aӛ_BngB?E�Bc�B���<�<�<�<�<�<�<�<�<�<�<�<�C�F�C�C���C��MC��C��bC�[�C�)KC���C��oC���C�^|C�NC�ՕC���C�.:C��wC��6C�]�C�8�C��C�$iC�-�C�$�C�!�?>�? �>���>�Jq>��u>���>�>��>�0F>��L>�~>��D>�b�>�^�>��g>�A�>��>1>�2]>��o?G�?%m?�?h?G�?E�?��?!*�                                                )hUG @�`    @��@    9'�.8���7��	7B'�6��}5���4��3� w                                                                    8S��7��L786u?�5���4���3��L2��9                                                                    4�c4��3r�2�iC2�1?�&0O@�/N�{                                                                    3���31�L2��215�0q��/��.�q�                                                                    6�ގ3�,�6E�3�2�a5���5��n        B��                    A�                @�p�                                                                                        @�3A2>A��A�|BF,B7�]BXv�B���<�<�<�<�<�<�<�<�<�<�<�<�C��	C��C��1C���C���C�f�C�CYC��C���C��QC��wC�q�C�>�C��C��>C�y�C�/�C��rC���C�u�C�?�C�'^C�,�C�%C�!�?�>��}>��w>�Ο>��>���>��>�p>�qe>��>婬>�I]>ީ3>ڜ�>��>�59>�Y�>�ȫ>���>���?&?�N?�'?yN?��?u]?m?��                                                )h�G8 @��@    @��     9�8p��7��79;�6�A�5��B4��3���                                                                    85\�7��7�6i��5���4��S4 ��3��                                                                    4�'V3�-3`Q�2��2Z1A;W0V��/\�z                                                                    3��R3�2���1�X312��0t/���.��M                                                                    6�{�4,��6�Fo3�}=2Kc39��3Vm        B��                    A�                @�p�        ;�D                                                                            @=m7AM�A�KA���B ��BKQ�Bp6B�©<�<�<�<�<�<�<�<�<�<�<�<�C�xC�̊C��C�8�C�h�C���C��4C��;C��gC�C�gC�dC���C���C�ėC��C�[�C��C��tC��?C�f�C�0�C�,/C�%=C�!�>��>�'�>ʁC>�zE>�vy>ӈ�>�d�>��P>�X >� a>�gQ>�=>�5�>ء�>�6�>��>�8�>�/�>�F�>�҇?)K�?'��?(��?(�?(�?)�?*`?,�F                                                )iGW @��     @���    9��8p!�7֠@7=ϛ6���5�ل4�8�4��                                                                    8$�K7���7��6o��5� :4�q4��3#�$                                                                    4�x3��)3a��2��2^N1Rg�0w��/�~�                                                                    3��43�2���1�p�1;i�0��/�e_.�jd                                                                    6� D3祖6N�I3$�2�/�K�/��        B��                    A�                @�p�    7P?���?e�                                                                        @e*A~Z A���B�RB>+hBu��B��xB�qd<�<�<�<�<�<�<�<�<�<�<�<�C�M�C��&C�ԭC��C�N�C���C��DC�&�C�fC���C�� C���C�";C�?IC�P2C�P�C�@fC� �C��VC��sC���C�>�C�,�C�%iC�!�>�&f>��>�i>�5�>���>�qt>��>���>��>�A>�D>�%,>ˈ�>�e�>�{:>΃
>�m�>�d4>ȹ�>��'?6�?H�_?K}?M�]?Q��?X��?e�?v;�                                                )i�Gu @���    @�Π    7���7���7���7&9�6mAH5��l4��&3�V�                                                                    6��96�l%6�'6Q�%5��c4�E3�n36�                                                                    30��3�3E�M2�R1�ή1��0/a�/p{                                                                    2_ 2GiU2y��1�11��08G/]��.���                                                                    6��34�M6q"3Rm1�4�                B��                    A�                @�p�    6.}A�4�AU�O                                                                        ?��8Ac�=A̩�B��B8�Bm��B� �B�*�<�<�<�<�<�<�<�<�<�<�<�<�C��SC��!C��@C��dC�'uC�f�C��XC��C�< C�x�C���C���C�/DC�g�C��SC��C��8C��jC��LC��fC���C�M�C�.dC�%�C�!�>��c>�v2>�v>�	>��6>�
�>��M>��l>�u�>��I>�e>��v>�p
>���>�6>Ţ=>�E|>�ȟ>�45>��P>��D?6�G?B�\?E}�?I�?P�]?[s�?n��                                                )�EG� @�Π    @�Ҁ    