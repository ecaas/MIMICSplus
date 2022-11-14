CDF      
      time       levdcmp       lndgrid       natpft        levsoi        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Sun Nov 13 13:48:36 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,NPP_NACTIVE,LEAFC_TO_LITTER,QDRAI,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,nbedrock,T_SCALAR,NPP_NNONMYC,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1902.nc -O ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1902.nc
Sun Jan  9 16:23:26 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1902.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1902.nc
created on 12/10/21 16:37:19    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:26 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1902.nc had following "history" attribute:
created on 12/10/21 16:37:19
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�                   1��@5P��4q��3tf�                                                                                    0�|�4��N3���2�[�                                                                                    -7KH0�<�/���/ F                                                                                    ,g��0
wT/ #�."�                                                                    6�3�h�63��3+�1hQ�                B��                    A�                @�p�        A*7�A���AP7hA���A�5A*��                                                        ?B@@@O0@��.A@��A���B)��BN�By�<�<�<�<�<�<�<�<�<�<�<�<�C�U�C�\RC�gyC�t�C��C���C�ÏC�{C�99C�hSC��HC��C�zC�M�C��IC���C� C�=C�Y�C�dCC�\C�*MC�9C��C�d>��m>��a>�_r>� i>���>�s)>�љ>���>���>��>���>�Ȇ>��I>�h�>�6C>��>��>�;z>��I>��                :�P;>��?�*?D^                                                "9�F�� @҉     @Ґ�                    4i%L5]��4wM�3{�g                                                                                    3�?�4��z3�12�i                                                                                    /��)0��0��/#�                                                                                    /��0�/#��.&��                                                                    6�J�3Ǔ�61\�3x�1n�;                B��                    A�                @�p�        A*7�A���AP7hA���A�5@�Y�                                                        ?<�@@_V@�ȌA@��A��B.aaBQ�\B~Z�<�<�<�<�<�<�<�<�<�<�<�<�C�C�(KC�<%C�S�C�lC��fC��'C���C��C�@ C�j�C���C�ɭC���C�;�C�zC��[C���C�5C�2_C�E�C�/�C�C��C�o>�޾>�d�>�Rb>�q%>���>�Ts>�h3>���>��>��>�;�>���>�9^>�!8>�i�>��u>�R�>�k�>��_>���                =5�?a?	\�?�                                                ":F�� @Ґ�    @җ�    8��            5�15jK%4y7�3y��                                                                    76��            44��4���3�f�2��i                                                                    3���            0�'�0���0Ш/)0                                                                    2��K            /��|0]�/%=u.%�I                                                                    6�)�3�T16.�}3 d!1���0��>0�t�        B��                    A�                @�p�        A*'�A���AP'�A���A��@EP�                                                        ?�њ@Ad�@���AA9�A�׈B3IwBR�B}u3<�<�<�<�<�<�<�<�<�<�<�<�C�wrC�zC�~�C���C���C���C��YC���C�C�(C�L C�sqC���C��mC�UC�<�C�wC���C�ݮC�C�(C�.JC��C��C�y>�"�>�B�>�{�>��I>��>���>�^/>�]M>�7Z>�� >��;>���>��p>�li>�S>�{�>���>���>��>�ǘ>Gc            =���?
RA?
F�?�                                                ":qF�� @җ�    @ҟ�    8*�6�g         5�x5_��4tO�3y�                                                                    7V�86�        4���4�j3�MW2�O�                                                                    3���2\��        1	g0��0 <�/��                                                                    2�1���        0-�0y|/!��.%!                                                                    6��/3�z�6+W*2��1��3U�Q30/        B��                    A�                @�p�        A&�	A�ˊANܧA���Aa��                                                            @���@�dgA _�AJ�LA�U�B.T�BO�IB|q5<�<�<�<�<�<�<�<�<�<�<�<�C��3C��%C��C��,C��3C���C���C���C��C�"0C�@C�aC���C���C��C��C�G$C�|�C���C��7C��C�'VC�&C��C��>�w>�vb>�u�>�v�>�w>��>��>��f>�1�>���>�>���>��>���>�Bp>��>�5>�B>��4>�F�?�9>��l        >�~L?YH?ږ?
1                                                ":�F�8 @ҟ�    @ҧ     93ɖ8�7���7 �6c��5��4�W3�t�                                                                    8c�7��Q7 $�6J>N5��54�o3��2���                                                                    4���4Z3U�2�#1��1->�06��/5�                                                                    3�3)�%2��j1�b-1�0Z��/f��.d�                                                                    6~�3��.6)�2���26XJ7642�        B��                    A�                @�p�        >��
?���?�C�@]�h?��                                                            @��A/�A�
jA�?�B$eFBO�LBv@B���<�<�<�<�<�<�<�<�<�<�<�<�C���C�0C��6C�TC�FC��%C��+C���C��MC�~�C�}C���C���C���C��TC��iC�&C�W�C��IC��HC���C��C��C��C��>��>��>���>���>�}[>��.>�q>�
>��6>�]�>�=3>��,>�V>���>�]@>���>�7x>��>>��k>�3n?Emq?A��?3CY?��?-x�?.�?/�U?2�3                                                ";9F�v @ҧ     @Ү�    91&�8�4�7�[7O#6�U�5���4��?3�D\                                                                    8_�/7�]�7NT6���5�_4��[4�k3 a                                                                    4�*4N�3{��2ل�2��1Ke�0X�/UY@                                                                    3��3;U�2��2	ad1D��0�v2/�t�.��                                                                    6}<�3�
?6(�2���1Ƙ�6�'�6�r�        B��                    A�                @�p�        :��                                                                            @��A-W�A�2�A�|�B�BDJUBm_B�o<�<�<�<�<�<�<�<�<�<�<�<�C���C���C�8C��SC�m�C���C��CC��C���C�`C�C��C��C�k�C�JkC�;�C�@pC�U�C�v%C��WC��fC��C��C��C��?�v>��>�K�>�I)>�>ڳ%>ҫ�>ʠ�>�Q�>���>�(>��>�ό>�"�>�A�>�l�>��?>���>��>��/?mc?��?��?��?�? �K?'�-?.�D                                                ";�F�� @Ү�    @Ҷ@    8�$�8+��7�l�6�ng6.�P5yw�4�4�3�;Q                                                                    7��7Y a6�8�6�5\��4���3�B�2˩>                                                                    4;�3�m�3��2�c�1�~q1�0�/)9�                                                                    3lUH2���2C��1�p60��"0%j-/IѸ.U��                                                                    6~ơ3�>�6)T2��i1�P�6��6��P        B��                    A�                @�p�                                                                                        ?��(@�XAmܰA��A��B��B?A�Bp�y<�<�<�<�<�<�<�<�<�<�<�<�C�O�C��9C���C�N�C�C��C��C�4�C��}C��C�]�C�|C�ȶC�|2C�0�C��C���C��eC��OC���C���C��C�gC�C��?L>���>��e>��>��>� >�0�>���>؅�>�^>�^>>ʋU>�ϼ>�w>��y>���>���>� �>���>�MO>�_�>�˶>�Y�>�¢>�t�>��>�&}>��+                                                "<F�� @Ҷ@    @Ҿ     98M�7�w7�|6EXB5x��4��]3��                                                                    8&ֆ7�z�6��A66��5yG4�3�\�2��#                                                                    4��J3�O�39QH2��s1�4?1~'0	j�/
l                                                                    3�8^3�|2j�1���1ݡ0$�>/-�T..�E                                                                    6��3�-�6+ā2�)J1���6é�6�z'        B��                    A�                @�p�                                                                                        ?�`\A��A���A� 3A�sB�B4]B[�~<�<�<�<�<�<�<�<�<�<�<�<�C�]�C�C��}C�ƦC��>C�~VC�Q�C�cC��C���C���C�kC�6fC���C���C�x�C�:�C��C��[C��IC���C��UC�C�C��>��>�Oe>��>苩>��>�86>��>�S�>�u>�o�>�vL>�7�>���>��>���>���>��>��>��>���>�->�@�>��>��>�i�>���>�� >�9�                                                "<eF�. @Ҿ     @���    9-�8��(8 �b7`\6��5�N�4�N�3�?�                                                                    8Z��7��7"��6���5�N�5u4L�3^                                                                    4��D4�)3�>�2�D�2(�V1]�}0lr/iA                                                                    3垧3A"�2���2�1Uu�0�/�]�.�Q�                                                                    6�q3�_�6.Ni2�w�1ڲ6j��6E��        B��                    A�                @�p�        8�#                                                                            @8�nAM�A�tA�fdB#
�BO��Bv��B�z<�<�<�<�<�<�<�<�<�<�<�<�C���C��{C��wC��>C���C���C���C�޹C�ͣC��nC���C��"C�[�C�/�C��>C��*C���C�SBC�&�C�gC��DC���C��C��C��>���>�]�>؟k>���>�D6>�T�>��>�!d>���>Ր8>���>ѰA>�9j>�UL>���>�`>���>���>���>��.?'�(?(F?))�?*�'?,{m?/�?1�?1�3                                                "<�F�j @���    @��@    9�=8r�7ڂ�7?�z6�D25��<4ʢ�3�V�                                                                    82Ξ7�k�7
�6rk�5��Z4�Rg3���3ؕ                                                                    4��`3��3ev{2Ʉ�27�1B��0T��/\�                                                                    3��z3!�2��1���18�?0u�H/�Uk.�k�                                                                    6�г3�j6V
3D(1� �4˸�4�U�        B��                    A�                @�p�        >|�                                                                            @GC�AW{A��AB'B�BT�!B}�@B�3+<�<�<�<�<�<�<�<�<�<�<�<�C��dC�1�C�l�C��%C���C��C�6DC�iRC��aC��:C�ŮC���C��=C��#C���C��C���C�}�C�U~C�2C�
uC��C��C��C��>�*>�e�>�}f>�7>��{>���>�>�j>�h*>�*>Ū+>��k>�k�>ǆ�>��>Ŏ�>È�>� �>���>���?.��?/L�?0Pt?1	�?2V@?4�1?8:�?>9�                                                "=-F�� @��@    @��     8�0�8!F�7��[7(�z6��5�V4���3�h�                                                                    7���7K��6Ҙ6U1�5�d:4⇘3���31                                                                    4	��3�^�3/�2�8j2
L1<B�0R�1/b�                                                                    3.v2��2]%01��y1.�x0m��/�?.��                                                                    6�1^3�=
62��3/J1qm�0��0�'�        B��                    A�                @�p�        A�A��@yGm>�\                                                                 ?��MA1�	A�3WA�6B.�Ba[�B��B�x�<�<�<�<�<�<�<�<�<�<�<�<�C�mtC���C��C�#�C�\dC��gC�ߜC�%rC�`EC��fC���C��C��C�=�C�\%C�o�C�u�C�n�C�]hC�F\C�#RC��C�qC��C��>��>�r�>���>���>��<>��>�}`>�C�>���>�D�>��>���>��>�IQ>��>�I�>��2>�:>�,'>��&>��?��? @J?7�?=i�?C6�?KG�?WB�                                                "=�F�� @��     @�܀                2]�3���4���4xw|3��"                                                                                1���3M^3�ى3��&2���                                                                                -�!b/v�50��0e(/��                                                                                -��.��/-</$��.F8�                                                                    6�jS41�6�o]3���1a'#                B��                    A�                @�p�        Ak�GA���Aʑ8A�m�A�ŹA��,A-��                                                    <g:e@�I@�}�A9 A�A�ipBT|�B�fC<�<�<�<�<�<�<�<�<�<�<�<�C���C�խC��C�h�C��_C�FC���C���C�HmC��UC�ƚC�C�E�C���C���C��C�UC�4+C�>HC�<C�,QC��C�C�DC��>�-�>���>�?�>��V>���>�ҍ>���>��">���>��>�R>���>��U>�t�>��E>���>�8�>��>�Lj>�*
            ;�ړ=[N�>�B?
A/?#�0                                                "`UF�" @�܀    @��@    