CDF      
      time       levdcmp       lndgrid       natpft        levsoi        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Sun Nov 13 13:48:31 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,NPP_NACTIVE,LEAFC_TO_LITTER,QDRAI,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,nbedrock,T_SCALAR,NPP_NNONMYC,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1853.nc -O ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1853.nc
Sun Jan  9 16:23:23 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1853.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1853.nc
created on 12/10/21 16:01:52    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1850-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:23 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1853.nc had following "history" attribute:
created on 12/10/21 16:01:52
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�   7Q��7;$7�/k7�6Er�5B�94 �l3�J�                                                                    6�m�6*�q6�I_6<v5yhK4u�V3K8�2��t                                                                    2���2��R3$��2��1��0���/��/2c_                                                                    2�h1� g2O��1�+�1f70uZ.�b.aU                                                                    6�i�3ѳE66�'3S0���                B��                    A�                @�p�        A�K�A\Dd=�v                                                                    ?��.A<� A�ۇA�Q�B& 2BT��B~�B�4�<�<�<�<�<�<�<�<�<�<�<�<�C��aC��C��)C���C��sC��C�+�C�Y%C��C��sC��EC���C�)�C�Z�C��:C��C��5C�#C�5C�BWC�C/C� �C�#C��C��>���>�D >�U�>��5>�>'>� �>��>�U">���>�v�>��J>��\>�o�>�#x>�:>��>��`>�%G>���>���>��N?�?,UD?.m&?1�?4��?9�w?@t$                                                ��D�� @�     @��     6��6�̦7i׸7
F�6>S�5:*8473��                                                                    6 #Q5�SX6���6.�q5pi�4k'�3@�2��                                                                    2U�2&H�2�k�2��+1ȋ�0�%�/�1n/#r�                                                                    1��1R
�2�q1��0�RZ/��+.�Yc.Nv                                                                     6��x3�$z64S�3�00���                B��                    A�                @�p�        A�rA�W�>�~x                                                                    ?m��@��A��A�6�B"TBJxBr$}B�^O<�<�<�<�<�<�<�<�<�<�<�<�C���C���C��sC��pC���C��C��C�69C�Z�C�{C���C��GC���C��C�K�C�C���C��C�-C��C�.�C�$C�eC��C��>��>�>}>��">� �>�`�>��>���>���>�i@>�t>��m>��4>�*�>���>�O�>�,�>���>���>��f>�*S><eR>���? o�?"�<?%?�?(�d?,�?2�                                                ��D�@ @��     @�     7�'K7e�7���7K(68<�52�e4�23u�                                                                    6���6+h6�8i6,)
5h��4a��36�?2�(�                                                                    3{R2���3"S2��1�!�0�kq/�Q�/oO                                                                    25=u1�9~2L��1�j�0�7�/� �.�g�.#                                                                     6���3�)�61��3��0��S0G�0��        B��                    A�                @�p�        A�wDA�^                                                                        @nA7{DA���Aݯ/B��BFUBk��B���<�<�<�<�<�<�<�<�<�<�<�<�C���C���C��C��<C��-C���C�+C�'9C�F�C�bvC��\C��C���C��C�dC�L�C�}C���C��JC���C�OC�!�C��C��C�>�f$>�n�>�A>�C�>�W>���>��>���>�by>�϶>�i�>�,�>�@>�:�>���>�Z�>��>��|>�G>��:>�
�?�C?F? ?!vs?#}�?&I�?*                                                �aD�  @�     @��     8�:K7���7���7?e6X�5u�R4kA�3���                                                                    7�.�7�6�(y6G�I5�s54�&O3��82��A                                                                    4,H3�53:+j2��k1�f1i�/��$/��                                                                    3YG�2���2k)]1ңy1��0#x(/��.@�c                                                                    6�D�3ɤ�6.tA3 U>19pv2P*.2!�D        B��                    A�                @�p�        A�C@��b;�Q                                                                    @m�,AQ�A���A쵗B%�BS�YB}aWB�|�<�<�<�<�<�<�<�<�<�<�<�<�C��C���C��bC�mqC�X�C�I`C�B�C�E�C�Q�C�`+C�t�C��WC��AC��tC���C�%�C�T�C���C��.C���C��fC�UC�-C��C�>�1�>�|f>��>�-�>��7>�͌>�O�>�g�>��
>���>��L>�9>��}>���>���>�4�>�˗>�`>��8>��?2F?*;?,W�?-��?/��?2�
?7=�?=��                                                ��D�� @��     @��     9A_8�n�8 c�7^��6�*�5�@�4�[�3��                                                                    8tBB7��P7"-6���5�Q5��4��3�.                                                                    4��X4!��3�K�2�-2(�N1b�Z0w�C/u;                                                                    4 ��3L:�2��%2#E1U-�0�=/���.��3                                                                    6�َ3�T6-��3�1n��3�qN3vC�        B��                    A�                @�p�                                                                                        @Z�^Ag[8A�A���B2V,Be��B���B���<�<�<�<�<�<�<�<�<�<�<�<�C��CC�5dC��C���C�_CC��C���C�i�C�"�C��C��1C���C�`OC�FC�9tC�=C�O�C�m�C��YC�� C���C��C�C��C�>��>�J�>�>�D�>ϒa>ʅ^>�P�>�_>��8>�_q>�>�>���>��P>��>�O>�~L>���>�3�>�4�>�=�?:C�?:�t?<�?>l�?B%?G�G?P��?^�h                                                �)D�� @��     @�x     9:k+8���7��L7U�6�"5�<4�y�3��b                                                                    8ky�7�U7�>6�95��n5%�4�3
)�                                                                    4�s�4�3��2�_�2#s@1\tF0o�Y/f�:                                                                    3�&E3CK�2�~2WH1Nv�0�;�/�\8.��?                                                                    6��3�	�63&�2�t�1q�S6��O6Y�|        B��                    A�                @�p�                                                                                        @(�aA@�A�/DA��B"vBM��Bx�iB���<�<�<�<�<�<�<�<�<�<�<�<�C�OC��qC�z#C�4<C��dC��_C�=�C�ڎC��C�=�C��{C���C�l�C�-�C��%C��/C��\C��XC���C��C��qC�C�?C��C�,?8�>���>��L>��L>�Q�>�.9>ޫ>��:>��>�TA>Ȭ�>�;�>�)�>�Y�>��P>�?�>���>���>��>��?b!?�}?��?!� ?&p�?,C�?3-�?:�                                                ��D�� @�x     @��     9H
�8�T�7�c7Y[�6��95�hA4�^�3��                                                                    8|�i7��l7 �V6�GS5�Z�4��4dR2���                                                                    4��^4!�_3�	2�i2" P1UC0a�R/S�D                                                                    4$�3Lx2�N�2��1Lʀ0��/��4.��                                                                    6��3�٬6:A�2��1l-O6�#�6���        B��                    A�                @�p�                                                                                        @$A2w
A�ZRA��B��B<Q�B_��B��.<�<�<�<�<�<�<�<�<�<�<�<�C��,C�.�C�ѧC���C�0�C��cC�m.C��aC���C�L/C���C���C�T�C�C��HC�clC�$aC��)C��fC��"C��C��C��C��C�9?�5?�?	|�?��?~�>�Y�>�=W>��>�k�>߶>��>�?H>��>�s$>�9�>���>��m>�x>�g�>�ί?�?�r?�|?D�?j�?�?5�?7�                                                ��D�` @��     @�l     9 28{�|7�޵79#~6��5�^B4š`3���                                                                    8I ?7�37q�6i��5�U�4��3��l2�e�                                                                    4��X4��3`R�2�_2A�1?�0P6�/G�K                                                                    3�Ҋ3'��2���1�o�12n80qb�/��.|CR                                                                    6��3�<�6BO31o��6�B6Z�d        B��                    A�                @�p�                                                                                        @
�A#MkA��0A�y]B��B3�|BV#�B�5u<�<�<�<�<�<�<�<�<�<�<�<�C�V�C�]C��bC��*C��tC���C�Y_C�rC��C��zC�}[C�=�C��bC��/C�UmC��RC��C�bzC�+�C��C��C��qC��C��C�E?T�?��>���>��>�h#>�@�>��>��>��>�F
>�`>ޓ#>��C>�r>�΢>�,>��>�}*>�.�>�<?}�?��?��?a�?9�?
��?m�?�	                                                �UD�@ @�l     @��     9>j8_�7�n$7(P6uA5��>4���3��	                                                                    8,�7��
6��6T��5���4�z�3��72�i�                                                                    4���3�*3I��2�]�27F1/��0@��/;64                                                                    3�]3p�2~�1�
�1#8Y0]��/s].lz\                                                                    6��-3��P6H��3% �1Vq6
T5�5        B��                    A�                @�p�                                                                                        @L�A!�=A��*AĔUB
�%B1��BSSB�W^<�<�<�<�<�<�<�<�<�<�<�<�C��C��(C��C��cC��1C��C��*C��C��C�uVC�XmC�4*C��C�ԑC��'C�M�C��C��eC�y�C�F�C��C���C�/C��C�Q>�8>���>��>�MC>�p�>�Oj>�>�I>�%>⑐>�~>���>���>�K�>��>�K�>�w�>��'>��>���?��?u�?�_?��?a?7@?6�?6�                                                ¹D�  @��     @�`     8�Wx8-�7�ܪ7�~6ED5�ō4���3���                                                                    8�7[��6���6)5x��4��G3��'2��                                                                    4WY�3�4 3��2��1ϵ+10�[/��                                                                    3��2�j2H�u1�7�1/04��/I�f.HVn                                                                    6���4j6y�*3��o11"X4IB4ٵ        B��                    A�                @�p�        =��                                                                            ?���A��A��A��B��B,K4BL�jBx�p<�<�<�<�<�<�<�<�<�<�<�<�C��C�E�C�r�C���C���C���C�jC�;�C�[�C�q|C���C���C��C��LC�rmC�OQC��C��lC��(C�z7C�9}C�C�RC�bC�\>��|>�R�>��>���>�A�>��,>�|�>�As>�V�>п7>�۟>҇^>ҥ�>�>зA>�fz>�P�>��~>�0J>��">�4�>���>�'X>���>���? n�?�D?�|                                                �D�� @�`     @��     9"8Z<77��<7)r�6vg#5�ޮ4�,�3���                                                                    8#7��06��>6V	�5��t4х3�SO2���                                                                    4�S3���3I�=2���2�@1.Ư0<ǿ/8D�                                                                    3���3@�2~�?1Ᏽ1#�_0\�/nu�.h�x                                                                    6��p3�a6L{^36��0���0�.�0�1�        B��                    A�                @�p�        ?y�7�Nh                                                                        @N+�A[�A�U�A�vB%�PBOb�Bq�	B�ol<�<�<�<�<�<�<�<�<�<�<�<�C�x�C���C�
�C�B�C�z�C���C��,C�2�C�d�C���C��C��C��LC���C��C��C��C���C���C���C�Y0C�7C��C�8C�f>��}>��L>�!�>�&�>�9�>��w>��>���>���>�a>�/�>��>ǭR>��y>ɤ
>ɦM>�Ց>�;�>�	P>6?27�?1ւ?2��?1y�?/�A?-��?+5?+.w                                                ÁD�� @��     @�T     6��?6�ݨ7�W7�26QMO5��[4�л3� �                                                                    5���5��%6+,�67M5�0�4���3���2�                                                                    2(�I2/��2��_2���1܏%1�m080�/<��                                                                    1T߄1^ �1�e�1���1L�0@�9/h��.nH                                                                    6�;=4K26�A�3v\#0��                B��                    A�                @�p�        A#��A��JA$V9                                                                    >��p@�U)A�-�A�m!B$q�BS�B~��B��<�<�<�<�<�<�<�<�<�<�<�<�C�J�C�pQC��0C�ɍC�BC�JFC���C���C�'uC�a?C��C���C�C�?�C�n�C���C���C��C��C��_C�h�C�$�C��C�C�o>�Z>��Z>��>��>�%+>��W>�g�>���>�W�>��a>�ݖ>�-�>�X�>�l�>�CS>>>��D>�>�>ðk>�h3<��=�݆?��?*��?.��?3y?9��?ASb                                                �ED�� @�T     @��     