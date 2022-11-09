CDF      
      time       levdcmp       lndgrid       levsoi        natpft        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Fri Nov  4 12:00:19 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,HR,HR_vr,NPP,NPP_GROWTH,NPP_NACTIVE,NPP_NACTIVE_NO3,NPP_NACTIVE_NH4,NPP_NAM,NPP_NECM,NPP_NFIX,NECM,NAM,NEE,NEP,NACTIVE,LITFALL,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,T_SCALAR,NPP_NNONMYC,LEAFC_TO_LITTER_FUN,GPP,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN,SMIN_NO3_LEACHED 31539_Modum_historical.clm2.all.1872.nc -O 31539_Modum_historical.clm2.all.1872.nc
Sun Jan  9 16:25:45 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1872.nc /nird/home/ecaas/31539_Modum_historical/lnd/hist/31539_Modum_historical.clm2.all.1872.nc
created on 12/13/21 20:32:58       source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31539_Modum_hist_for_decomp    Surface_dataset       !surfdata_31539_Modum_simyr2000.nc      Initial_conditions_dataset        -31539_Modum_Spinup.clm2.r.1201-01-01-00000.nc      #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         9./31539_Modum_hist_for_decomp.clm2.h0.1850-02-01-00000.nc      Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:25:45 2022: Appended file /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1872.nc had following "history" attribute:
created on 12/13/21 20:32:58
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�6G�5��o5e�H5W�m5r�5�4�\�4��N                                                                    5|��4�3�4�'�4�C4Ih�47q4�]3���                                                                    1�14�0���0�S�0��0�#G0�-0&n�                                                                    1	�'0d.�0SN0��/۱�/���/�9�/R;=                                                                    6�&;3�?6(�B5��6X��6P�t66��6)��6&a6G6S[5��                                                6�u3z�7���        1��    7:ڷ:ڶ��v                                    B��@Q�5                @Qе?	�                ?�G�    5�e�*�pA�IQAͼ�A��A���A͎�A�r?>�&�                                                    >��Q@WqA4�A���A�19BF|�B�*�B��R<�<�<�<�<�<�<�<�<�<�<�<�C��C�:$C��rC��C�J[C��HC��
C�2�C�s�C��sC��C�(C�e�C���C��|C�%�C�[�C���C���C���C��C�u�C�\*C�WC�S�>���>���>���>���>��>�c�>��>�cR>��W>��>���>�Y�>��k>�co>�-'>��b>��>��7>�,;>���=~�=��v=��">3��>��?X,?Rd6?`0/                                                ��E�� @�^     @�}                         0�G3��4�7V                                                                                        /���3m�3�E�                                                                                        ,Ab2/�^T0S$                                                                                        +tF$.�w'/7�q                                                                    6��3�Ē6�/4ǐ*                    1��<5'q5�h                                                6�O�3�=7��]        1��    6�0��0��-                                    B��@Q�5                @Qе?	�                ?�G�            A�A�Aժ.A��A�z0B�GB8�A�U�                                                    <-��?�y@��4A<�cA��eA��dB>�B�0`<�<�<�<�<�<�<�<�<�<�<�<�C��'C��#C���C��C�.C�VrC��C���C��NC�1�C�iLC���C���C�&C�nBC���C���C�2�C�^9C�zxC���C�z C�^C�W.C�S�>��#>�f�>�E>�Yu>��>��>�Q
>���>��>�Q�>�*�>�F�>���>�9t>�;�>�V>�:>>��>�9�>��2                    ;P��>�.?%�_                                                �-E�� @�}     @��     7��5��_3=QA            45��4Tv                                                                    6'��4��2o#E            3e�X3�B�                                                                    2�p1��.΀l            /�Oj0E(                                                                    1�?�0G'w.l            .�C//�                                                                    6��3�X7���5O��7)X�6o��4$�            5VV�5��X                                                6��m3u07�o�2$��    1O��0�����6���6���5[��0�*8/N�+(c        *!~0dp�    B��@Q�5                @Qе?	�                ?�G�            A��oA��cA�RA�E4B��B8�
AjXQ                                                    ?�	o@7�rA^�A=0%A���A�"PBl��B�]$<�<�<�<�<�<�<�<�<�<�<�<�C��JC��CC���C��'C��C���C��.C��C��C�C�?*C�o\C�� C�ܱC��C�b�C��DC���C�VC�GC�o�C�w�C�_�C�WaC�T>��>��>�i>��>�,w>�J�>�v�>���>���>��>� Z>�z�>�,�>�=�>��*>���>�ua>�#>�U�>���>I4k=�"�<�            ?',�?%�                                                ��E�� @��     @��     9+85$�7��5�G4�z3���4�N�4ya�                                                                    83r�7d�W61~5q3#�2�
�3�V!3��                                                                    4���3Ŗ2�E1�^&/�D�/06�0K?0�                                                                    3ý2��1���0�v�.�q�.^��/�]�/+Ґ                                                                    6��<3�(�8F��7(s�9�i8��7�J6��5
ȝ4��5�Q5��                                                6��3�7�F<4i`    1�F22+�����7���8�[7�p�3,�1��|.�        -t�P2�N    B��@Q�5                @Qе?	�                ?�G�            @�F@�AA��A�g�B�B$]?	�                                                    A/	�A�m�A�TA�a�A�X�A���By_B��><�<�<�<�<�<�<�<�<�<�<�<�C���C�>�C��C��C��C��.C���C��'C��GC�
�C�.C�U�C���C��4C��.C�)�C�jTC���C��C��C�LXC�o�C�`�C�W�C�T>���>��>��>��`>�u�>�v�>�Z�>�>��q>�W'>�"7>�'|>�c�>���>��>�i�>��>��m>��	>�Ӆ?q??y}�?M��>�"�>�>2�:?3�?9�m                                                ��E�� @��     @��     91��8��H7�V9726�+6 O�5�o_4GK                                                                    8`��7��h7y�6`��5��5"�4�q�3�:�                                                                    4���4}3f�r2�G62r�1���12Kb0>H                                                                    3��R37+�2���1�gl1<��0�Ϋ0a6�//��                                                                    6��3Ց�8�)�7�h�9:{�9"q�8���8��7�\�7�86Ĩ�5��8                                                6�ܵ2ʥ�7���4��Y    1�"2&�𷶟F7���8;8@�6$C�4�|}1b�2        4�i�6�n    B��@Q�5                @Qе?	�                ?�G�    6���*�"�:��            >��@�	w                                                        @���A�;A��?B��BJ��Brv�B��xB�|�<�<�<�<�<�<�<�<�<�<�<�<�C�3RC��.C��C���C��C���C�X�C��C���C�̀C��mC��,C��C�ďC��aC�-C�D6C�}�C��>C���C�)�C�cXC�a0C�W�C�T+>ᒢ>��}>ˁ�>���>���>��>��B>���>�e�>���>�v�>���>��>��[>��>���>��k>�j>�PT>�M:?S?Tnf?V��?YO�?]��?U3�?bc?q{�                                                �YE�� @��     @��     9K��8���8T|7`�6���67�F5��5��                                                                    8��S7�yy7%��6��65�!�5hA05a�4�U                                                                    4�4'�A3�@�2�m�2E:1Ȕ71w�17a                                                                    4={3T �2��G2`!1x�e0�\�0�ki0g�                                                                    6���3�J{8���7�*9V N9<��8�R89J�7���7O7��6��U                                                6���2��57�
4p�-3�1�`�2���n27�n�89��8!�6�`�5E�`2�        5W<G6z�U    B��@Q�5                @Qе?	�                ?�G�    6U��*�@�                                                                                @r�yAy��AҁB�B<D�Bs��B�ُB���<�<�<�<�<�<�<�<�<�<�<�<�C���C�<	C��<C�}�C� �C���C�IC��sC�oC�4C���C�i�C��C�۩C���C��lC�{�C���C��>C���C�jC�TQC�`KC�XC�T>?	�)?X/>��>���>�b >甶>ߌo>�t>Ц�>���>�?�>�#>��#>��]>�j >���>��{>�ć>���>�Ⱦ?E�R?F:�?H&�?J�$?OC�?V�'?cey?|]�                                                ��F L @��     @�	�    9[�}8��8��7q�q6o6D�B55��                                                                    8��7��M72�6��p5���5x��5H�4�|                                                                    4�N450�3�u�3��2T+�1��r1��41DwO                                                                    4e	3dߏ2�:2&�1� �1��0�
]0x*�                                                                    6�?�3��8�H7��l9g�
9LYI8���8Ay7��7S37D�6�"�                                                6��_2�57��4��(.�I1��2A���.J7.�8)�8"��6�r�5<}�22X        54Ȟ6q�    B��@Q�5                @Qе?	�                ?�G�    4b�(�ad                                                                                @\��AiA�A�l�A�M�B3[�BgB�AB�W<�<�<�<�<�<�<�<�<�<�<�<�C�a�C���C���C�,5C��C�h�C��dC�w{C��C���C�;0C��tC�oC�ZC���C�V�C�C���C���C��/C��C�E�C�^TC�XGC�TR?�j?=??M?F�?	;�?��? �>��>�]�>���>ތ>�V�>Е�>�q>��+>��3>��>���>��d>���?;BC?;�?=n�??�6?C��?Ist?R��?a��                                                �!F � @�	�    @�     9K$�8�&�8
zS7p��6��]6M�R6z[5ke                                                                    8�M!7��T7.�[6��5���5�҅5#�4�q                                                                    4ݖ�4,�{3��3^S2W�}1�<K1�@w1 3�                                                                    4�m3Y�02���2%�i1�F�1�X0�la0!�V                                                                    6�|�3�+*8���7�79W.�9Cu8�Я8;Y�7�@i7T`�7	� 6y�#                                                6�-i2�`�7�2�4�.%J�1���2O�~�yro6yv8d}8
R�6B�4���1�X	        4̎�6"UE    B��@Q�5                @Qе?	�                ?�G�    6���*ۯ�                                                                                @}��A��A��"B/�B>�5Bv�9B���Bȶ�<�<�<�<�<�<�<�<�<�<�<�<�C��bC�y�C�a�C�J�C�0QC��C��C���C�m�C�1�C��bC���C�J�C���C���C� �C�C�uJC�>wC��C�zC�<C�[�C�XYC�Tf?�?�B?�C?�!?p�? �s>�T
>���>� >�}�>�@�>刱>߄t>��>�!#>�f�>�kr>���>�Q>�|�?I�C?Iӏ?K�?N%�?R�?Yh�?f6�?r�:                                                ��FD @�     @�(�    9*tn8�q�7�l7Xrq6�[ 6C��5���4�Z�                                                                    8WO�7�Z7�~6��5�X 5v��4�U�3�<�                                                                    4���4�P3��2��2F�51�K�13�0){V                                                                    3���3;��2��2�1{EP1��01�/V                                                                    6�3�8=H�7ǜ�95d�9)V�8���8%^+7���7Dy6���5� f                                                6��]39}7�n�4�M.[�1��24�X5��]���7��7˰�5�Q4t^�0�N�        3��5�ac    B��@Q�5                @Qе?	�                ?�G�    6��+)4�                                                                                @��A��!A���B��BHkB�	�B�=�B�Ph<�<�<�<�<�<�<�<�<�<�<�<�C��C�B�C�^�C�uC��'C���C���C���C��
C��	C���C�g`C�@�C��C��6C���C�1�C��)C���C�j�C�:�C�;QC�X�C�XPC�Tz>ެ	>��m>ᆙ>��Q>�F�>�t~>�D�>慧>��>�*�>㮈>Ꮰ>��D>�J�>��>��3>�{>�s�>�,Y>��g?OOz?P��?SHM?W<@?]ƒ?i?y~�?�                                                  ��F� @�(�    @�7�    9	��8^�;7�``73R6��6#�05�vj4�-                                                                    8-��7��6�l6b��5�45O�4Ⰶ3�8�                                                                    4�p3��u3Z�x2Û�2%�[1��\1C��0&�                                                                    3���3Re2��1��1Qu�0��0wNx/Q��                                                                    6���3�,�7�͈7��9��9
u�8��&8|�7�%7"��6�'�5�l�                                                6���3<�7�K�3�-��o1��71vBF7�!췇!�6�+e7 3��2VBh-�<        0��3��    B��@Q�5                @Qе?	�                ?�G�    6��,7i>R�                                                                            @�+zA��[A�3�Bp�BE��B�=B��B�p�<�<�<�<�<�<�<�<�<�<�<�<�C��C��C�{C�P�C���C��mC���C�#�C�KOC�i�C���C���C��dC��HC��C�p�C�E�C�AC��RC��eC�geC�CFC�V�C�X2C�T�>���>���>��j>��>�Є>��>�rv>˲<>�9>�6i>�ܝ>��>ә�>Ӄ�>ҍ_>М)>��8>�a7>�Ӗ>Õc?NkL?OL�?Q��?U1?Z�B?dǃ?s,?�                                                  �MF8 @�7�    @�G     8qm�7�\�7N$�6��6�x5�'k4�	 4�                                                                    7�{#7�O6�2J5�
�5A�4���4�3��4                                                                    4��3h�2���2F�1���1 ��0x�0#�l                                                                    3&T�2��"2V1z�z0ң�0J�w/�o/N�R                                                                    6���3�t	6�d7#�]8�-8�8|7��}74�6��5�kE5���                                                6�'Y3

e7���/��G+��1A\-K��7�EI��EI���i2�}�.N<�-�W&� �        (�/�.'�%    B��@Q�5                @Qе?	�                ?�G�    7�+��A7=�@tZ�                                                                        @�"A�.�A��B��BO��B��|B�VwB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C���C��|C���C�5'C�u�C��C�
�C�[�C��1C��aC�&sC�dC���C��VC��C��C�pC��uC���C���C���C�PvC�U�C�XC�T�>��&>�_�>�a�>���>��>��E>�ÿ>�=�>�[>��>��z>���>�>� >�dD>�ʧ>��>��>�\�>�&�?6�?S�?X�O?\�?b��?k�:?y;�?�                                                  ��F� @�G     @�V     8`<�7�h�7\��6�]&6�5�r95	B4�+P                                                                    7���7	�U6�hf5�	�5HJ@4�g�4-a3��                                                                    3���3n?2��72V4z1��b17�0�}b0v�                                                                    3}=2�x�21�I�0�~a0A��/��F/Im�                                                                    6�t�3��86�`I7*)�8oe�8��8,��7���7��6��6J�5�J                                                6�_	3��7���        1DB    7�5E��5E� @�                                    B��@Q�5                @Qе?	�                ?�G�    7.z, ��A_u�?f�                                                                        @b�+A��	A�].B�SBUĹB���B�CB�Ph<�<�<�<�<�<�<�<�<�<�<�<�C��UC��<C��C���C�4C�M�C��:C���C��C�$�C�[�C��C��lC�CC�@C�s�C���C���C��zC���C��_C�^C�VC�W�C�T�>���>�WN>�b>�P>�an>�>�i$>�o>�^>>�+n>�2�>�``>���>��>�le>��%>��>�zF>��>�\�?&A?Z��?^(G?cz�?l3�?yQ�?��?�                                                  �uF, @�V     @�e�    