CDF      
      time       levdcmp       lndgrid       natpft        levsoi        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Sun Nov 13 13:48:33 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,NPP_NACTIVE,LEAFC_TO_LITTER,QDRAI,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,nbedrock,T_SCALAR,NPP_NNONMYC,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1872.nc -O ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1872.nc
Sun Jan  9 16:23:24 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1872.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1872.nc
created on 12/10/21 16:12:21    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1850-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:24 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1872.nc had following "history" attribute:
created on 12/10/21 16:12:21
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�   4��0        4]*�4�P4uȳ4z5<3���                                                                    3��        3��U4��3�;d3��2�(�                                                                    0@��        /�b0r�20Ij0��/-R�                                                                    /s4�        /��/�6q/#OB/&:
.Z�y                                                                    6�?�3�K�60��3".0���                B��                    A�                @�p�        A�'�A��A�D�A饛A��YA��+AC�>�C�                                                =ڋ�@%�`@�5�A[sA�B�BY]zB���<�<�<�<�<�<�<�<�<�<�<�<�C�vhC���C��C�I]C��PC�4kC���C��/C�;C�t�C�� C��SC�*�C�h�C��>C��C�gC�<�C�TbC�\�C�U�C�(�C��C�WC��>��&>���>���>�~�>�q@>� >�Q�>�[�>��J>���>��>�5@>��>��/>��>��>���>�6">��7>�$�<"9        =��>U!>�c?	�?6�P                                                ��E�� @�^     @�}                             .�s2��b                                                                                            -��\1��                                                                                            *�.N��                                                                                            )4�I-���                                                                    6�t�3�F�6.��3 c�0��,                B��                    A�                @�p�        A�[AצGA��A�	�B�'B��B"�A�a�                                                <*�?ݫ�@��A2߳A�q^A�)�A��B6�<�<�<�<�<�<�<�<�<�<�<�<�C�J$C�V�C�jRC���C��>C��C�uC�pC���C��;C�C�Z�C��TC��~C�)�C�rzC���C��SC�C�0�C�A�C�,�C��C�uC��>���>�I�>��>��>���>���>�f�>�Ǻ>��>�H�>�9@>�j>��y>�s�>�q?>�xX>�;�>�sK>��>�z�                        :H|�>c�                                                �-E�� @�}     @��     5\��                        3X�                                                                    4�Mj                        2<��                                                                    0�                         .��                                                                    0��                        -�h�                                                                    6�L�3��'6+��2���1.�.�8.�P�        B��                    A�                @�p�        A~G�A���A�
�A�&B�}B �B"��A��                                                >_Q?�+�@�N�A2�A�,TA��A��-BL��<�<�<�<�<�<�<�<�<�<�<�<�C��jC��IC���C��<C��2C���C��+C���C���C��?C�TC�/�C�`�C��C��^C�\C�_�C���C�ԯC���C�$�C�*�C��C��C��>���>��T>�ƭ>�ب>��<>�/>�9>�s�>�@s>��!>��Z>�:n>��Y>��/>���>���>�g�>���>��>�t�;�N�                        >��y                                                ��E�� @��     @��     9+8I-�7%�V6��3���    0̤'3p�}                                                                    8%��7~6QR51��3 Ę    0?A2�%�                                                                    4��3ӵX2�_�1��*/��    ,W@&.�\5                                                                    3�[�3��1�C40���.�$�    +��. W                                                                    6}�`3��6(q"2��1��p3J{�3'/        B��                    A�                @�p�        @��A87zA�x�A�V)A�1Bf�B�Q?�n                                                @��A�C$A�A�d�A�ӻA�%-A�vBz��<�<�<�<�<�<�<�<�<�<�<�<�C�wC�dC�͜C��`C��eC���C���C���C���C��/C���C� �C�H�C�w!C���C��aC�&nC�c�C���C�˦C��C�#YC�$C��C��>� >��>��t>��&>��e>��<>�$�>��!>���>��]>���>�r[>�~�>��>��|>��>�=[>���>��i>��<?o�?"�t>�`�>��=�,�    <�Uu?;                                                ��E�� @��     @��     9V4n8�ϙ8O)7o1�6�΂5�2�4�s3�B                                                                    8�I�7ҵc72~�6��5��o4���3���2㱶                                                                    4�i�4/��3���2��R2(!E1!(�0��/=��                                                                    4]�3]��2���2�41T_�0K�/Gi.ow#                                                                    6zg=3��6&?�2�]1�P�6|��6R�        B��                    A�                @�p�    5��r<}�5�S        >��wA�[Aa��                                                    @�W�A�V.A�|�B�BBE�B_:�Bp�B��<�<�<�<�<�<�<�<�<�<�<�<�C��wC�2�C��tC���C�i�C��C�v�C�P�C�=�C�4�C�5FC�?�C�S�C�q#C���C�ʝC�+C�9�C�p�C��VC��C�PC��C�C��>��>λG>�z�>��k>�y>��>���>�M�>�&�>��y>��7>��>�T>��^>��o>�E�>�4�>�L�>�[I>�,�?Lh�?MZk?O�!?S{�?YQ�?;�?!v$?=	8                                                �YE�� @��     @��     9|�8ǭ68+�7���6�7�6�q5�4yV                                                                    8�47�97X�Y6�CM6>75)K44�"36~Q                                                                    5�4R�3�{�3Q|2X�l1��&0�k�/��                                                                    4'�
3��2���2B��1�)0��/�&.��                                                                    6y�N3���6%��2�Ļ1���6�3�6���        B��                    A�                @�p�    4�Ѫ                                                                                @jǇAs%A��Br�B8��Bn&�B�?�B��E<�<�<�<�<�<�<�<�<�<�<�<�C���C�*�C��[C�SkC��YC�t�C��C�tC��C��	C�FpC���C���C�eKC�62C��C� ,C�6�C�[EC���C���C��C��C�ZC��?
�d?��>��W>�b�>�X�>ⶡ>��>�r>���>ø�>��>��y>�w[>���>�*>��Q>��>�$>�)�>��<?Aգ?B:?C��?Fi�?Jr�?Q?[��?o                                                ��F L @��     @�	�    9�Y�8�6�84f7�36�v6q555h4��                                                                    8��^8"�7c�>6��36
_53�P4? 3<'�                                                                    5
7�4\'3���3#1
2f}�1��0�
C/���                                                                    4.�s3�,2���2N"�1���0�E|/��.���                                                                    6|13���6'z:2��1���6���6��p        B��                    A�                @�p�                                                                                        @O4�A^�OA�+|A�4�B.�B`O B�fdB�<�<�<�<�<�<�<�<�<�<�<�<�<�C�8C���C�T�C���C��C�7�C��PC�>CC��^C�`�C���C��vC��C��fC�S7C��C���C��OC���C���C��C��C��C��C��?.i?)�?3j?4�?.�?�~>�ʛ>�U>踮>�6�>٥>�J�>�g.>�ɫ>��
>�wr>���>��>��'>�G�?4VB?4Ǿ?6
�?8u�?<D@?A�4?J6�?V��                                                �!F � @�	�    @�     9wV�8ʋ�83#:7�=�6��a6m5"8<4%$�                                                                    8�6�7��17bGr6�[�6i�5;|84L�3P�                                                                    5"t4U"63�}o3$c]2k��1�!�0��W/���                                                                    4$aj3��X2��2O��1��D0�8/ׄ�.�_y                                                                    6�>G3�l�6*[�2��1�-S5�Zz5�@S        B��                    A�                @�p�    6�#�                                                                                @tM�Ax.�A�`�BƞB9��Bn�]B��9B�B<�<�<�<�<�<�<�<�<�<�<�<�C�C�Z�C�=�C�"�C�(C��C��C�n�C�.�C��C���C�W�C�NC���C�9�C���C�qC�#C��|C���C���C���C�cC��C�?��?��?47?�? ��>��>�	�>�0<>�E>�]>�4�>�m�>�ac>��>� O>�O�>�c�>���>�v�>���?EL7?D��?F��?H>�?KM�?P�O?[N�?jӠ                                                ��FD @�     @�(�    9K#8��8:�7���6�,�6�54$c+                                                                    8�L7�17D66�+$6 ��50��4H��3O��                                                                    4���44��3�T�3�G2V˵1�30�S/�߄                                                                    4 �3d��2�P528L�1��0��/�[�.�]�                                                                    6��3ļ�6-�%2�c�1�M443��4��        B��                    A�                @�p�    7Zw�                                                                                @~l#A�XAغ@B	��BB=)B}&jB�M�Bԁ@<�<�<�<�<�<�<�<�<�<�<�<�C���C��5C�C�#�C�9PC�LfC�ZQC�`nC�[�C�P�C�=BC� C��BC��FC���C�7tC���C���C�OyC�RC���C��C��C��C�>�\>�.z>ۻN>�1>ޔ�>��v>��L>�(�>���>� >ޗ+>܇W>��r>�Sz>��>��C>ǎn>�>�U�>�7K?Jn[?K�w?M��?Q9�?V��?`?o�3?�C                                                ��F� @�(�    @�7�    9 j48��7���7bJ"6��	5�I�5`[4<�                                                                    8J�7�*�7�6��t5���5�4'6Y32g�                                                                    4��4�3�_2�'20[�1sJ�0�: /���                                                                    3�7�37D'2�4�2^41^�0��S/�ݢ.���                                                                    6�)�3�h6S�Z3C"�1Y2�2�UQ2�*�        B��                    A�                @�p�    7��?5<�                                                                            @{GA�G)A�4NB��B?�mBxr(B�2OB��3<�<�<�<�<�<�<�<�<�<�<�<�C�7C�\8C��C���C�	[C�EC���C��AC��OC��C�%�C�;/C�GQC�IKC�=�C�"GC��4C�óC���C�W[C�;C��KC�
�C�vC�>��X>��>�5)>�9c>�c0>���>�x{>���>ǽ >��>�ģ>�#�>��>�)>�B�>�z�>���>�w�>���>���?H{?JG�?Lm.?O �?S��?[��?iE�?{�D                                                �MF8 @�7�    @�G     7�TF7w�;7�6�73�-6g�5��4��3��                                                                    6�w�6��36�`
6c5�xH4���4�H3�                                                                    3�!3y�3P�k2�~1��W1�0^��/m��                                                                    26�*2$��2��D1��A10Go:/��B.�.�                                                                    6�
�3��(68΁3O31#n�,��c,�:        B��                    A�                @�p�    6>ClA�_�A�@;�Ϩ                                                                    ?��A4A��A�\cB4^�Bi��B���B�:,<�<�<�<�<�<�<�<�<�<�<�<�C�sC���C��(C��-C�6�C�|�C�ȸC��C�`JC���C��C�wC�I�C�x�C��&C���C���C���C���C�pC�:�C�bC�	�C�MC�>��>�f{>�u�>�m�>��`>�5;>�5�>���>��>���>���>���>�}>���>�)�>Ď1>��r>�5>�j�>�N�>�!�?x�?<4??ǐ?Dʮ?K�?V�:?f(z                                                ��F� @�G     @�V     76���7�Q7!!�6Fe35C��4(�3���                                                                    6E- 5��6�s�6K�?5z��4w?;3T�r2�o                                                                    2�@�25��3{ 2��?1е�0��}/�7�/@Q�                                                                    1�z1e}	26��1�$1�B0	�.�ڠ.r��                                                                    6�;4�t6w�>3m��0��5                B��                    A�                @�p�        A�1�A�=?(�                                                                    ?��KAi�A�uyA�}DB&��BU�xB��oB�P<�<�<�<�<�<�<�<�<�<�<�<�C��vC��
C��#C��0C��C�PC�OeC��sC���C���C� `C�T�C���C��PC��
C�()C�N'C�eC�kkC�cC�G'C��C�	�C�$C�!>�ô>�)>�IR>��>���>�;�>���>��`>���>�*D>��d>�ҁ>���>��	>�<>���>�B�>���>�r>��L>�a�>�I?+��?.�?1�.?5�?;�?C>�                                                �uF, @�V     @�e�    