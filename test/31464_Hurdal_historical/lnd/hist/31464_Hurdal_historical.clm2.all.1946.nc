CDF      
      time       levdcmp       lndgrid       natpft        levsoi        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Sun Nov 13 13:48:41 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,NPP_NACTIVE,LEAFC_TO_LITTER,QDRAI,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,nbedrock,T_SCALAR,NPP_NNONMYC,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1946.nc -O ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1946.nc
Sun Jan  9 16:23:29 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1946.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1946.nc
created on 12/10/21 16:59:24    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:29 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1946.nc had following "history" attribute:
created on 12/10/21 16:59:24
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�                   /+14X��47@�3Jx                                                                                    .X=�3��3gz~2N�                                                                                    *�+m/�'�/�א.��                                                                                    )�+/.�.d�                                                                    6�U3���6J�W3'�1�                B��                    A�                @�p�        AT��A��A�5�A��A��Ac?ɘ�                                                    >�!@B4$@�. A>w�A��#A�j�B:�IBi8<�<�<�<�<�<�<�<�<�<�<�<�C�^�C�qcC���C�C���C�TFC���C��C�AAC�zC��;C��IC�8�C�}TC��+C�
*C�FJC�tEC���C��C���C�UuC�)�C��C��>��m>�Nw>���>�΁>��~>�x�>�m�>���>��>�
�>�@�>���>�?�>�Q>�/�>�5�>���>���>�M�>���                8��>�>�U�>�                                                (�iG� @�     @��                        3�`�47��3:=                                                                                        3��3hO�2k?�                                                                                        /�Gr/���.��                                                                                        .��.�y�-��,                                                                    6�;3��q6H&93<)1��D                B��                    A�                @�p�        AT��A��A�71A��A��A��?��.                                                    >��m@B3p@�'&A>a�A��mA�l�B:�TBaZU<�<�<�<�<�<�<�<�<�<�<�<�C�'MC�-�C�9�C�I�C�]�C�yC���C���C�
�C�6%C�f�C��6C���C�;C�XC���C��C��C�H�C�e�C�xC�Y�C�,�C� C��>�X�>���>�4S>��>���>�5�>���>���>�V�>���>�E>�х>���>�1>��:>���>��1>�@�>��>��4                    =��>�j�>��[                                                (��G	 @��    @�#`    6��                4�-�45�e3:�                                                                    6{�                3��23eS2k��                                                                    2�_                0�/�.Ĭ�                                                                    1���                //��.�VX-�nv                                                                    6���3�7'6D�|3�1��#1�>1�g�        B��                    A�                @�p�        ATʿA��A�3�A��A��&A��                                                        >�*@B��@�3�A>f�A���A��B9EBa'<�<�<�<�<�<�<�<�<�<�<�<�C�YC�\�C�c�C�m�C�z>C��XC��'C�݊C��C�)XC�PC�z�C���C�ܬC��C�WzC���C��HC��C�0�C�WC�W�C�/{C� �C��>���>��>�4k>���>�D�>�X>��]>�->�(�>��Q>��D>��>�z&>�=>�q/>��>���>��>��>�|�=4G�                >*�,>��m>�bJ                                                (�1G	: @�#`    @�'@    9�t8H�g7w�6�$=5�ep5@4M��3C��                                                                    8E�7}�u6��5��4�<�4r�3��2weI                                                                    4�%13�B3r*2?��1GY�0�1@/���.�6T                                                                    3�Wf3m2$�1r*0{�w/�f�/�.=B                                                                    6�5�3�N�6A�"3l�2ad4-t�4�E        B��                    A�                @�p�        @-��@�D@�h�A�iA>-@=�H                                                        @��An�A�s�A� �A�3}B'��BBmFBd�G<�<�<�<�<�<�<�<�<�<�<�<�C�.[C��SC�`qC�%C���C��GC��C�lC�kC�-C�H�C�i�C���C���C��C�'�C�b�C��*C���C� NC�3pC�O=C�1�C�!C�>�V%>�FM>�[N>��>��H>��>�->��m>� >�k>�~�>�3>�(�>�y�>�:
>�O�>��->��+>���>��6?L2�?)��?�e>͟�>�X>�\$>�D�>���                                                (�G	X @�'@    @�+     9Ei8M�7��7N�6T��5���4�A'3���                                                                    8'O7�ص6��6;V+5�C4��73��2��                                                                    4�5�3�[�35��2��1���1B50&�K/$:m                                                                    3��3��2eU�1�/=1Q�0?C/Ra�.Or9                                                                    6��I3���6?vY2ȼg28�16n�/6G�        B��                    A�                @�p�        :��~                                                                            ?���A�A��A�'�B��B6iBZ��B���<�<�<�<�<�<�<�<�<�<�<�<�C���C�"TC���C�n�C�7C��C�pC�+C���C���C�n�C�J,C�0�C�"�C�#;C�4AC�TC�~3C���C��\C�GC�B�C�2�C�!�C�>���>�O�>֎�>��J>ˑ>�z>�y	>���>��F>���>�Y�>�K|>���>�E>��>�G>��L>�U>�¶>�Wq>��}>�6'?  )?��?�e?mz?�u?;#                                                (��G	w @�+     @�.�    9-t�8��7��78-(6z�E5���4�ɚ3��b                                                                    8[#7���7|�6h��5�N�4� b3׻I2�F                                                                    4��'4�g3k�\2���2�1*�03�/,�
                                                                    3��30��2��1��1&�-0V�e/c8.Z�                                                                    6�Fb3�m6?�2�-[2�n6ɡ>6�\        B��                    A�                @�p�                                                                                        @$k@A5`IA�.A��|B�B6MfBU��B�]�<�<�<�<�<�<�<�<�<�<�<�<�C��'C�OC�SC��OC��/C�4C�ޣC��C�0^C��VC���C�s�C�:�C��C�՞C��&C���C���C���C��C���C�4oC�2mC�".C�? w&>�'>�e>�{T>�H>�>�P>Ѫ�>�p>�e>�dN>���>��>��>�&?>�&�>�#>�%�>�>��\?�v?�?�?� ?J?�?��?*�                                                (�]G	� @�.�    @�2�    9�&8^�i7��?7\�6`�45� �4��3���                                                                    89��7��]6��6F�5��4���3�Qk2�T                                                                    4��K3�=+3D�D2��]1�y1#G�0://:�C                                                                    3�D�3�2xD 1�;h1i�0N?�/k	.l4�                                                                    6�C�3�q6@Y�2���2
G�6��S6��        B��                    A�                @�p�                                                                                        ?�- A�A�r�A�k�B9B)KBN,uB~[`<�<�<�<�<�<�<�<�<�<�<�<�C��3C���C�p_C�&�C�ڱC��*C�$DC���C�\wC�WC���C�tC�'JC��mC��&C�LGC�eC��C��C���C��ZC�'�C�19C�"�C�,?�"?	9�?�?�t>��>���>�>灮>��S>ۊ�>�)>��>���>���>�=|>�(�>�� >��D>��4>��>�gW>�\�>�/]>�׼>�q>�d?�`?�e                                                (��G	� @�2�    @�6�    9 �8zz�7�<�7/��6o��5�.�4��3��                                                                    8J�7�2r7�<6^/�5�Dh4���3��2ŹI                                                                    4�y4Ε3`r�2�'�1�1"I�0+��/$��                                                                    3Ղ�3&~62���1��1:�0L��/X��.P2�                                                                    6��3���6B��2ݫ�2 �s6���6�X        B��                    A�                @�p�                                                                                        @�dA%~�A�d�A��YBXB)PBE�Bk2<�<�<�<�<�<�<�<�<�<�<�<�C�w�C�"�C���C��yC���C�t�C�;�C��C���C���C�L�C�C�ʥC���C�2C��C��pC�[C�.UC��C��C��C�/YC�#C�??	:?�>��7>��g>�_>���>�C>� A>�>���>ߥ�>�7�>֠�>ѷ�>̀�>�d�>�ʸ>�
>�Y>��o?��?��?Zk? �r>���>���>�>��>                                                (�%G	� @�6�    @�:`    9E[w8���8�;7q46���5�dx4�n3�r�                                                                    8yK7̛
71��6�V�5ױ�5?t4*�3 �                                                                    4Ϻ4*y�3�H�2��M23��1pgk0��/�p{                                                                    42$3WV�2�No2 \11c.0�Շ/�y.��                                                                    6�
�3��
6Ec3��2	�64_5�r�        B��                    A�                @�p�                                                                                        @`��Ai8�A�%|A��/B/G�B^B�hB���<�<�<�<�<�<�<�<�<�<�<�<�C�(C���C��LC��C��vC���C�̛C���C��C�y�C�VXC�,�C��2C�ňC��oC�=�C���C���C�y�C�NHC�%BC��C�-sC�#gC�S>>�K�>� >�H>�a�>�a�>���>�
�>���>�ۻ>�V�>�d�>��>�F=>���>�>�>Ȣ�>�k.>��g>�B�?=�?;ώ?=?=:~?=ԝ??A(?A��?G��                                                (�G	� @�:`    @�>     9�8|	�7���7Hz�6��A5�S4���4@                                                                    88-�7�.~7��6}<�5��5��4�;3)W�                                                                    4�x�4��3m��2�
2ƥ1[H�0�J#/�+x                                                                    3��3'��2�82I�1D�0�~�/��.�Q�                                                                    6�R�4(�c6��3��2�:3>��3A�        B��                    A�                @�p�    5�?���                                                                            @gnAs�xAΕ>B�2B9	lBn��B��5B��<�<�<�<�<�<�<�<�<�<�<�<�C�.C�p4C��6C�2C�MaC��$C��C��C�E*C�d�C�|�C���C���C���C�q"C�M5C��C���C���C���C�I�C�&C�,C�#�C�h>�F�>�M$>���>��M>�� >�:>�9�>�"�>��d>���>�~S>�f�>ҝP>��>С">�C>�.�>Ǿ�>�\>�d??.Q?B�?DjO?Gt?K�?Q�?\b�?p2�                                                (��G
 @�>     @�B     8�Y�8<��7͘�75�r6�Z�5��v4�E3�                                                                    7���7n��7��6eS5�r�4�4-3	�                                                                    4Q��3��F3Xj�2�2a�1E�0g�^/}}�                                                                    3��o2�2��1�h11R�0y�/�m�.��                                                                    6��73��x6T��3.�q1���0$�#0        B��                    A�                @�p�    6�(@���@ˊ�9zLq                                                                    @YO	Af\AA�vpB��B=W�Bu~<B�0zB�x�<�<�<�<�<�<�<�<�<�<�<�<�C�)�C�?�C�e�C���C���C���C�5C�y�C��OC���C�'�C�]]C��+C��hC��C���C���C��SC���C���C�h�C�1IC�+�C�#�C�~>�v�>�_�>�4�>�b�>��w>���>�.>���>�}@>��>� >�7�>�!E>ļ�>�̖>���>�b>�@>�<k>�Z?s�?3��?H�??K�?Ph�?Xe?dF@?v�g                                                (�QG
. @�B     @�E�    8}@�8r7�"�72�66�D5��4�3�                                                                    7���7+�6��6a��5�y4�Y�3�!L3��                                                                    4H�3���3KM02��2��1-�\0(nH/tgE                                                                    3(\ 2�U2�f�1��1)��0[QL/T�3.�\,                                                                    6��P44��6�M73�>1��c                B��                    A�                @�p�    7M]A@��A%S_<E��                                                                    @2��AX��A��B	�BB�/B�ϝB��{B�E<�<�<�<�<�<�<�<�<�<�<�<�C��BC���C��C�_C�=C�oyC���C���C�:C�?]C�qC���C��[C�	�C�;�C�hC���C��*C���C��2C�u6C�=_C�,%C�$C��>�s>�uS>�n�>��>��>���>�LX>�H>�8>���>�gA>�H�>�.�>�,�>�+�>��P>��3>��>��>�uB?	�?&�?Lُ?P��?W\�?a��?n��?z%�                                                )G
M @�E�    @�I�    