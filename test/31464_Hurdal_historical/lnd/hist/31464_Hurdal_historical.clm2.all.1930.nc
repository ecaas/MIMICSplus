CDF      
      time       levdcmp       lndgrid       natpft        levsoi        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Sun Nov 13 13:48:39 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,NPP_NACTIVE,LEAFC_TO_LITTER,QDRAI,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,nbedrock,T_SCALAR,NPP_NNONMYC,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1930.nc -O ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1930.nc
Sun Jan  9 16:23:28 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1930.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1930.nc
created on 12/10/21 16:51:14    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:28 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1930.nc had following "history" attribute:
created on 12/10/21 16:51:14
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�   7�>�7>��7�W6z��5�r4»3�۱3���                                                                    6�ȝ6q?+65�F5�8j4ʪ�3���3Gc2��*                                                                    3	��2�4�2��52J1('�0Ld/Y̱/�                                                                    2-ɇ1��\1�|y1%־0Th/��.���.C�$                                                                    6�<G3�^�6@E�3q1�M�                B��                    A�                @�p�    7Zw
A�T@k��                                                                        @%fA�9Aؕ�B	�BB�B|��B��QBֳ�<�<�<�<�<�<�<�<�<�<�<�<�C���C��:C���C���C��VC��C�AC�m�C��|C��,C��C��C�IUC�y�C���C�آC��bC��C�&�C�(C��C��C��=C�#C�
>�7Z>��C>�;>��2>�:�>��>�>�eu>���>��3>��>���>�,@>�ߥ>��Y>�Np>��s>��>��&>���>�
?K}A?Mŷ?Q�?V}�?_�Z?o�1?�                                                  &iF�^ @܄     @܋�    7�n73]�7��7)\6F3'5C�4 �V3�5�                                                                    686b�6�6U�|5z[�4w(�3J�62�C�                                                                    2��_2�3H-2���1Ϲg0�/�N/�q                                                                    1���1��2|�1�:�11�0��.Ԙ�.0O�                                                                    6�d
3׏�6=Л3�1�I�                B��                    A�                @�p�    6��A��AE                                                                        ?���Aq��A�B��B9\�Bo��B��rB�g�<�<�<�<�<�<�<�<�<�<�<�<�C�~SC��OC���C���C��fC��-C�}C�A�C�e�C��hC��&C��JC��C�(C�W�C���C���C���C��C�C�C��XC��gC��C��>�v7>�xm>�f�>���>��|>�x�>�1B>�!�>��>���>��u>���>���>�U5>��$>���>�;e>�l�>��>��>���?Ag�?C�R?F�K?K�|?R��?^5�?sp�                                                &�F� @܋�    @ܒ�    7=�7	U7���7%�(6D�&5@ T4s�3{�                                                                    6D�6-$�6�U6Q��5x�4r��3D\{2�M                                                                    2�� 2��+3C�52��1�h�0�V�/��B/��                                                                    1͏P1��62v��1�� 1],/�R�.���.&��                                                                    6� �3�n�6:�R3	i�1Ɯ/h�j/7�@        B��                    A�                @�p�        A��A��<8e�                                                                    ?��CA:/_A�a�A���B1�BekpB���B��<�<�<�<�<�<�<�<�<�<�<�<�C��C��gC��IC��'C��@C���C�C�$/C�C�C�`�C��C���C���C��vC��C�M�C�| C��~C��C��C���C���C���C�#C��>��>�m?>�%�>�)�>�>�>��>���>���>�:_>���>�e�>�9>�-0>�[q>���>�j0>��>�i�>�t�>� 9>�N ?~�?:Yy?=[a?A�B?G��?P��?^��                                                &�1F�� @ܒ�    @ܚ�    8�=89�7�V�7<��6}v�5�E/4���3���                                                                    8T�7jӺ7	�6n9�5�4�!3�J�2�_                                                                    4r�O3���3d��2ŭ2�1"��0*K�/B��                                                                    3�e?2�)32��1��"1'�0M��/W.u��                                                                    6���3я�67fg3��2�e2X1���        B��                    A�                @�p�    6�Ψ@��@�uw                                                                        @�[�Ap/�A� �B#B6$�BjGB�<�B�e�<�<�<�<�<�<�<�<�<�<�<�<�C�f�C�!�C��C���C���C���C���C�|sC�{�C���C���C��BC���C��tC���C�"7C�N�C�zGC���C��/C��[C��hC���C��C��>�8>���>���>���>��3>��>���>�D�>�/�>�g�>���>���>�0>���>��>��>�xF>���>�!>� (?>Z??�I?A��?C��?F�c?L��?U�m?f��                                                &��F� @ܚ�    @ܢ     9L��8��Z8��7pD�6��5�9�4��4 t�                                                                    8�@�7�Ϛ70\�6���5�S5$�4�;3"B�                                                                    4ևR4,t�3�Y�2�Ք23D�1o*�0��s/��#                                                                    4}�3Y��2��M2�1bq�0�}/�oN.��                                                                    6�N�3�eX68_y3j�2.�3�s�3_�        B��                    A�                @�p�    6w                                                                                @lH�AtoeA��tB�B9rUBo��B�q�B���<�<�<�<�<�<�<�<�<�<�<�<�C��EC�3�C���C���C�duC�LC���C�u�C�1�C��lC���C��C�r&C�VC�E�C�DC�P�C�iC��lC���C�ϔC���C��5C�_C��>��>��>�Tm>Ԭ<>��>�
$>��I>�˳>��<>�R*>�>�>��n>�2>��[>���>�� >��z>���>���>�i�?B�*?C�?D��?Gb�?K�?R9?]�z?r��                                                &��F�N @ܢ     @ܩ�    9V8�ae8�j7t��6�I�5�3v4�^l3���                                                                    8�9�7��J74/C6�lG5��5�043�j                                                                    4�pU41��3���3 "�25��1p00��b/}X�                                                                    4�63`G]2�޽2!�.1eg�0���/��.�2                                                                    6���3��6C�}3QX2�J6�a6��j        B��                    A�                @�p�                                                                                        @DEAU1�A���A��B(�BYbGB��2B�w�<�<�<�<�<�<�<�<�<�<�<�<�C��*C�fC���C�`�C��C��uC�ZNC��C���C�Q&C��C��C�wFC�5�C��JC���C��'C���C��C��-C��IC��!C���C�	C�`?�/?�i>��|>��}>�{�>��>�ރ>٥4>ә'>΢G>ɺ�>��>��g>�ϲ>�;S>�f>���>��->���>�U�?..?."�?/.�?1`�?4�M?:�?AT,?K/�                                                &�]F� @ܩ�    @ܱ@    9Z��8���8˙7x��6�^5��4�~93���                                                                    8�g7��78)�6�9�5݄O5/54�	3�{                                                                    4��44��3�Ҽ3vB27ʯ1p�	0�#�/uk                                                                    4��3dm�2�
02$�F1h(r0�$�/��.��l                                                                    6���3�?�6P�73 ,,1�.M6���6��3        B��                    A�                @�p�                                                                                        @5e AGb�A�&A�"�Bs�BK�Bs�B��<�<�<�<�<�<�<�<�<�<�<�<�C���C��C���C�}:C�6C��C��.C�ZC���C�s�C��C���C�tC��C�C�o�C�)�C��sC�ԀC���C���C��NC��*C�
�C�:?3�?L!?��?�L?��>�y>�!�>�"u>�9�>�x�>܉�>֛�>��n>�K>�nS>�N>��>��!>�:>�5N?$�q?#��?$<G?%T?&�?)ھ?-�;?2z�                                                &��F�� @ܱ@    @ܹ     9V=8�"�8�7US6���5��@5�P3��y                                                                    8�5�7ڲg7;�6�C]5�h�5l4%m3 $                                                                    4�j�45}�3�4'3�2>Wg1|QO0�9�/��6                                                                    4��3e@[2��2)�1png0�[�/�V0.���                                                                    6�Q4 �6]��3Ex�1��L6���6Y��        B��                    A�                @�p�                                                                                        @B�AR�A��DA�ZQB%��BSzRB}*�B��<�<�<�<�<�<�<�<�<�<�<�<�C��LC���C�W�C�0C�5C��7C��C�[�C�9C��9C��"C�]�C�C��yC�e�C��C��6C�eC�)=C� VC��-C��C��C�
UC�?	�?5~?�?\�? �_>�D�>��/>�>��e>�>��!>���>ۛ[>���>���>��C>�S�>���>��>��.?-8?,�;?-g�?.j ?0S?3s/?7��?>;C                                                &�%F� @ܹ     @���    94�`8��D8��7f��6�v5�K4��D3��                                                                    8d�Q7���7&k�6���5�DO5�4�a3x                                                                    4���4"]3�L2��]2.u�1j�0�ܰ/���                                                                    3�3G��2�r2��1\^�0�Of/�	K.�                                                                    6�'U476fֵ3gd1��5�W5�o�        B��                    A�                @�p�                                                                                        @MƬA[}�A�ڵA��.B)��BW��B�>B��J<�<�<�<�<�<�<�<�<�<�<�<�C�WC�s�C���C��$C��C���C���C��C��C���C�i�C�HrC�NC��C���C�_�C��C��C�z�C�BC�HC��C��>C�	�C��>�6�>�R>���>�ߣ>�Q>�-�>�K�>��=>���>�g>�>>�]=>�v^>���>�u�>�x�>�GK>�O�>��8>���?2��?2E?3V�?4%3?5��?8l�?<�D?C�#                                                &��F�B @���    @��@    9)$8���7�{�7bzO6���5��/5H84
�5                                                                    8U��7�JW7 �6�	�5ң�5�a4([93/f�                                                                    4�:4s{3�a2�b�2.ŷ1rz0���/���                                                                    3�ݗ3>"2�z�2�1\��0��/�i^.��@                                                                    6�wc4�J6��+3�ז1�;�2á�2�TW        B��                    A�                @�p�    6��(=:H�                                                                            @w_A{ܝA���B�(B<�Bs}�B�ZWB�D(<�<�<�<�<�<�<�<�<�<�<�<�C�:�C�z�C��9C��sC�C�2$C�Z�C�}�C���C���C��C��=C���C��lC��PC�Z{C�(C���C���C�w�C�0wC��C��7C�	�C��>��;>��.>��h>�Z'>���>̯)>�Cj>ъ�>�s>�p>Ԩ�>��4>�km>�r�>Ѳ�>�!o>��i>�f>�K�>��i?G1o?Gv�?Il�?K�?O�^?Vm?bDd?s��                                                &��F� @��@    @��     8���8-�:7�77@�6��G5���4��3��?                                                                    7��7[��7��6s��5��L4���4D�3�                                                                    4;;3�E�3_��2�5-2�91NR*0yP/�L:                                                                    3l�u2�=k2�(�1�k�192c0�N�/�u�.��.                                                                    6��4��6g�|3n�r1��.���.��        B��                    A�                @�p�    7d @��a@��                                                                        @U��Aq9�A�>�B.BD�
B��nB��YBђ!<�<�<�<�<�<�<�<�<�<�<�<�C���C���C�.�C�r�C���C��C�Y�C��rC��C�2�C�mC��oC���C��nC��C��C� �C���C��_C��+C�Q�C��C��KC�	@C��>�0>��>�G�>��Y>�|�>��>�3p>��>�>���>�<6>�m�>�A>�Ex>ɭ�>�4>�SC>ǜ�>�2�>�t�?'��?<��?O_?Sz?Yg~?c��?rO�?|��                                                &�QF� @��     @�׀    8Xm8Z7�\�7.�6Q�+5]@E4h�T3���                                                                    72��7:�;6�	06[�5�� 4���3��3ԓ                                                                    3�-	3�( 3M��2��1��0��1/��/kQr                                                                    2�+j2���2� �1憾1
�0r/Y.��:                                                                    7f4S3�6�ً3�g�1�pp                B��                    A�                @�p�    6��.AsvH@�)                                                                        @4�Aw��A�PCBO�B;�Bq�2B��&B�d�<�<�<�<�<�<�<�<�<�<�<�<�C���C���C��QC���C��C�O�C���C��C��6C�1�C�jC���C��;C�C�K�C�y�C���C���C��-C��vC�a�C�C���C��C�l>�g>�C>�L�>�@�>�]>>�ۏ>���>��>� >��>��>�C�>���>��J>� d>��X>��>ßf>�R�>��?}.?EGE?G?I��?M�j?T��?`�?vX+                                                &�F�� @�׀    @��@    