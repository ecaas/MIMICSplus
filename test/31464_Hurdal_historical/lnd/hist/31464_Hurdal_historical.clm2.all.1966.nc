CDF      
      time       levdcmp       lndgrid       natpft        levsoi        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Sun Nov 13 13:48:43 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,NPP_NACTIVE,LEAFC_TO_LITTER,QDRAI,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,nbedrock,T_SCALAR,NPP_NNONMYC,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1966.nc -O ../test/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1966.nc
Sun Jan  9 16:23:31 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1966.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1966.nc
created on 12/10/21 17:09:48    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:31 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1966.nc had following "history" attribute:
created on 12/10/21 17:09:48
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�               4�W�4KV�5-$4�E�3���                                                                                3���3�l�4A|c3�X@2ȓ                                                                                 0
9/�u�0�I�0!R�/(Uq                                                                                /.q�/v/���/K�.T��                                                                    6��J4;6Y {3��2e�o                B��                    A�                @�p�        A;�A�J0A�A�_�B|?�?A                                                        <#�1@ �@�q3AlA�}BO��By\�B�`�<�<�<�<�<�<�<�<�<�<�<�<�C�(C�2�C�CUC�W&C�nC���C��C� �C�5�C�d�C���C��C��C�K?C���C��C��C�C�C�e�C�v�C�x>C�S�C�>�C�7�C�1	>�a2>���>���>��>��!>���>�
�>��>���>���>���>��%>���>�G>� �>�
�>���>���>���>��=            =�py>eN%?/+�?41??:V                                                +��G%� @䬀    @�`            5�p4���424OH4i��3��                                                                            4�%/3ֵ�3`�H3��;3���2�u�                                                                            15I�04�/��g/�ph/���/�q                                                                            0d��/c�q.�i�/
��/Rq.I��                                                                    6���4�L6Va3��2y�                B��                    A�                @�p�        A;�A�J0A�A�#�B�B@W�                                                        <#�@ �~A ;YA|�]A��qBK5Brk�B�>�<�<�<�<�<�<�<�<�<�<�<�<�C�LxC�SdC�]�C�i�C�x�C���C���C��C�RC�<2C�gdC��UC��<C���C�=C�}{C��gC��;C�#�C�D�C�_�C�VKC�@HC�8C�1+>��>�j�>���>�}>�/�>�q*>�c�>��e>��>�ٔ>��>��m>�%�>�>�}/>��>���>��I>���>��        =9�W>��>G!�?)��?-9�?1�                                                +�G%� @�`    @��            6���4�m4���5#��41{33���                                                                            5Ɋ�3��b3���4N��3`/�2�)�                                                                            2)�0@�30&�(0�j,/��/�!                                                                            1U�Q/s��/R�3/��.��.H��                                                                    6�mg4 ��6Sd3G52w�k0�!/�-b        B��                    A�                @�p�        A;�A�@OA��vA�4 A�!? ��                                                        =�u�@!r�A@��A~ҢA��$BJ]?Bq0)B�f4<�<�<�<�<�<�<�<�<�<�<�<�C��cC��RC���C��uC��pC���C��NC��C�
�C�*SC�NC�ueC���C��C�C�BQC�~�C���C���C�YC�AC�R�C�ApC�81C�1M>�
�>�>�'<>�<�>�T�>�8>���>���>�Xo>�� >��~>��n>�.>���>��]>���>�'�>�rc>�k�>��        >7ܶ>�->���?(�?+�Y?0�
                                                +�qG%� @��    @��    7B�|    5�"�4�0�5"^5Fj'4'V�3�"I                                                                    6v�    4դ�4��4An�4z�	3S_�2��                                                                    2�ja    1350a��0�K 0�B"/�QQ/�                                                                    2^"    0b^	/���/� l0�t.���.Fg�                                                                    6���3���6OG73z�2�EN2���2�&        B��                    A�                @�p�        A:QLA��VA�6�A� A_                                                            ?* @$ڑAJ�nA��6B�jBI��Bo�kB�S<�<�<�<�<�<�<�<�<�<�<�<�C��3C��3C��3C��3C��3C���C��C��cC�kC�"#C�A;C�c{C���C���C��GC��C�P�C���C���C��C� �C�I�C�A�C�8]C�1p>�w>�w>�w>�w>�w>��Z>��v>���>�!�>���>�?>��>���>� 2>���>�}>���>���>��M>�YZ=���    >_��>%�?'��?'��?*Ru?.D|                                                +��G%� @��    @什    9-N8K8G7��6���6E�5���4��(3�8z                                                                    8Z�H7�Yp6���6��5x�v4�j�3�G[2���                                                                    4�� 3�GR3ݸ2c�h1��|1 �0,6�/G�h                                                                    3��3�'2I�1���1�0J�/Y��.|M                                                                    6�q�3��6Lr�3�-2��A5���5��        B��                    A�                @�p�    6Wm]@/'W@�0LA�AB�~?�tA                                                            @��AnjA��A��uB.�B^��B���B�S<<�<�<�<�<�<�<�<�<�<�<�<�C�#�C���C�2�C��C��PC��#C�t�C�`9C�X�C�X\C�_�C�n�C���C���C�ˣC���C�.�C�c�C���C���C��C�=mC�A�C�8�C�1�>�P�>�H\>��->�/�>�H6>��R>�!�>��Q>�p>�[�>��H>�s0>��O>�5�>�S�>���>���>��>��e>�D�?J	?�\?"X>�=~?;7=??�?E��?R��                                                +�9G%� @什    @�`    9j<�8�M�8��7v�6�9'5��Y4���3�0�                                                                    8��q7�}m7:��6���5�H15��4�j3b                                                                    4�*�4=��3��K3�25{�1q��0�f/� n                                                                    4��3o�2ſ�2%R$1e>'0���/��-.��                                                                    6��3�h�6K��2���2�v�6���6��#        B��                    A�                @�p�    3�'�                                                                                @P�A`9oA�Q-A���B/�Bb�B�2?B�a�<�<�<�<�<�<�<�<�<�<�<�<�C�f�C��-C�"�C��mC�.�C���C��C���C��C��#C�E�C��nC��&C�r+C�K#C�;ZC�B�C�\�C��*C��
C��,C�/:C�@FC�8�C�1�?�=?��?^D>��>>�	?>���>�u;>�L�>�@�>��>�C>���>��H>��d>�Q]>�i�>��>�=�>�X�>��J?4� ?5�~?7+c?9�=?>C*?D��?Mz�?[!&                                                +��G& @�`    @��     97+`8�97��7K�z6�=&5Ý�4��3���                                                                    8g_57��f76��25�2K4�4��3�                                                                    4��4�I3}`E2��2�.1OS�0c�/a��                                                                    3�%3>T'2��2o�1A :0��/��.���                                                                    6�>Q3�36M�D2�92�a�6��66�<        B��                    A�                @�p�                                                                                        @ةA/��A��Aѽ�B
2B@A�BgCbB� 6<�<�<�<�<�<�<�<�<�<�<�<�C�H�C���C�y�C�<�C��NC���C�S�C��PC���C�8�C��NC���C�%�C��FC�i�C��C��hC��eC���C��yC�ڙC�!C�=�C�8�C�1�?$?�?��?��? 0>� 
>��>�/>�*H>�Ba>�3>���>���>ž�>���>��>�m@>�U�>��>�S�?D�?YF?��?h?�V?u�?!D�?&�6                                                , G&8 @��     @��     9C��8��N8�7_GR6�V.5�^^4�>|3�-                                                                    8w6�7ï�7&X�6��5�zU5;�4��37n                                                                    4�_G4$#3���2�d2%��1]��0p�/nѤ                                                                    4��3OT�2�J�2y1QV90�(�/��.��?                                                                    6�~L3�ו6P�.2���2��j6�2�6h�C        B��                    A�                @�p�                                                                                        @<}MAL�]A���A�	B!^BL|�BrT�B�gG<�<�<�<�<�<�<�<�<�<�<�<�C� �C��rC��\C���C��bC�ZC�(9C��C���C��C�E?C��C��C�rC�%C��C�wSC�5�C�<C���C���C�C�:�C�8�C�1�?7c? >���>��7>��s>�>��,>�Dt>��>�L�>�>ژ�>��7>а,>�'�>Ű�>���>��{>�>��$?)��?(?(�a?(�4?)�h?*��?,�#?0M�                                                , eG&W @��     @���    9,y8��57�'(7Qm�6�	.5�\�4�3�]�                                                                    8Y�P7��7��6�E]5�A}5 pX4��3��                                                                    4�4�3��72���2��1W�0oE�/sG=                                                                    3���3<Zy2���25�1G��0�%/�q.��A                                                                    6��3��?6S�i3v<2���5�a5��        B��                    A�                @�p�                                                                                        @?C~AQ�A���A��RB$QbBQSWBy��B�)�<�<�<�<�<�<�<�<�<�<�<�<�C�k�C�~pC�� C��YC��	C���C��kC��C�jEC�RBC�2�C��C�߳C��C�i�C�"�C��"C���C�[�C�0�C��C�4C�7�C�8�C�2>��t>�l>�^>���>�*B>��>�>�a�>���>�5>���>�'�>��>�j<>�&S>˅�>��>�>��>�~m?+��?+o(?,K�?-4?.x�?0�4?4�U?9��                                                , �G&u @���    @�Π    9�Q8�<�7��h7G��6�f�5�64�O3��                                                                    8DĜ7���7�\6|fI5�f�4�l4Lh3�                                                                    4�94
 �3r��2��<2.1U/g0s�U/��                                                                    3Є3.Q�2�7U2�&1A`�0���/�l.�e�                                                                    74U\6��3���2�p3Q��3.}�        B��                    A�                @�p�    6>.�?V�z                                                                            @e�NAp:HA���B�B4��Bg7B�7hB�k�<�<�<�<�<�<�<�<�<�<�<�<�C���C�.�C�cC���C��:C��C�%C�4C�NnC�`|C�m!C�r�C�p�C�fC�O_C�,0C��VC��6C���C�g�C�38C�TC�5_C�8JC�2>>���>�X�>�&�>�>�J>Ǻ�>�La>̲�>�e�>ϑ�>�d�>��\>С�>��C>�hV>��>�+�>��j>�>��@?>x_?@"�?A��?B��?E!�?IFM?P��?_�                                                ,-G&� @�Π    @�Ҁ    8��28f�7�V�75Z�6��I5�Q74�A3��<                                                                    8��7���7>6e�5���4��R4Q�3;N                                                                    4���3���3X�w2�C2��1K�E0r:�/�
�                                                                    3���3��2��t1��h14<-0��,/���.���                                                                    6��o4j6WzW3Ɏ2pq/?�d/��        B��                    A�                @�p�    7T(^@�l?mP�                                                                        @k�pA��OA�O>B
��BBiYB|��B���B�$	<�<�<�<�<�<�<�<�<�<�<�<�C��C�#C�]kC��C��vC�oC�TC���C��?C��C�FvC�vQC���C��7C��mC���C��%C���C���C���C�S�C�&�C�3�C�7�C�2\>�8>��S>���>�|�>��	>��>��>���>���>��>��*>���>�9S>�I<>ƺu>�J->��>Ś�>÷�>��?3  ?L�}?OB�?R$�?V�?_uE?n^C?�i                                                ,�G&� @�Ҁ    @��@    8v��7�v[7h�{7 %�6U��5��4Q��3�
�                                                                    7��7s,6�v6!��5���4��x3�zd3J(                                                                    4��3cE2���2��_1⟙1f�/�^]/t
                                                                    3%32���2�|1��1!n0C�/q~.�!P                                                                    6�'�4��6z�B3T�2i�;                B��                    A�                @�p�    7��pA@�Avc�A*Fk@��?z�                                                            @A��AI�GA�XB	o�BR�iB���B�sZBԐ�<�<�<�<�<�<�<�<�<�<�<�<�C�+�C�ZC��HC��C���C�4�C�p�C��C��EC�"�C�[C���C�ʀC� �C�4�C�`C�~�C���C���C��NC�c?C�3C�3�C�7�C�2x>��>��q>���>�U*>��9>��!>���>���>��>�Q>�,R>�\�>���>��$>���>�X>�4W>�">� n>�[�>�*>�Z?%4??�?d�,?l��?v��?~��                                                ,$UG&� @��@    @��     