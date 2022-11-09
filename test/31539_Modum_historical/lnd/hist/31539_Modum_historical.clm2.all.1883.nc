CDF      
      time       levdcmp       lndgrid       levsoi        natpft        hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      �Fri Nov  4 12:00:20 2022: ncks -v CWDC_TO_LITR2C_vr,CWDC_TO_LITR3C_vr,CWDN_TO_LITR2N_vr,CWDN_TO_LITR3N_vr,FROOTN_TO_LITTER,FROOTC_TO_LITTER,HR,HR_vr,NPP,NPP_GROWTH,NPP_NACTIVE,NPP_NACTIVE_NO3,NPP_NACTIVE_NH4,NPP_NAM,NPP_NECM,NPP_NFIX,NECM,NAM,NEE,NEP,NACTIVE,LITFALL,LEAFN_TO_LITTER,QDRAI,PCT_NAT_PFT,mcdate,T_SCALAR,NPP_NNONMYC,LEAFC_TO_LITTER_FUN,GPP,TSOI,SOILLIQ,SOILICE,W_SCALAR,NDEP_TO_SMINN,SMIN_NO3_LEACHED 31539_Modum_historical.clm2.all.1883.nc -O 31539_Modum_historical.clm2.all.1883.nc
Sun Jan  9 16:25:45 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1883.nc /nird/home/ecaas/31539_Modum_historical/lnd/hist/31539_Modum_historical.clm2.all.1883.nc
created on 12/13/21 20:45:21       source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31539_Modum_hist_for_decomp    Surface_dataset       !surfdata_31539_Modum_simyr2000.nc      Initial_conditions_dataset        -31539_Modum_Spinup.clm2.r.1201-01-01-00000.nc      #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         9./31539_Modum_hist_for_decomp.clm2.h0.1850-02-01-00000.nc      Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:25:45 2022: Appended file /nird/home/ecaas/all_sites_decomp/31539_Modum_hist_for_decomp/lnd/hist/31539_Modum_hist_for_decomp.clm2.all.1883.nc had following "history" attribute:
created on 12/13/21 20:45:21
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�            2��3�sB3|kJ4�Z�4�
�                                                                                1(�2���2�l/3���3��                                                                                -���/�;/	�@0#N0#�                                                                                ,�A�.D..�/N�/N��                                                                    6�g3���6DĻ5"z.            2�D�4�T�4���5�\ 5�_�                                                6�)�3�7��        1<�)    7	5��	5����                                    B��@Q�5                @Qе?	�                ?�G�            A���A�g�A��WA�Bj�B)Ô?���                                                    <#�
?4�%@��AM��A��eA��CBt[{B���<�<�<�<�<�<�<�<�<�<�<�<�C��C�kC�UC�1�C�J�C�h+C���C���C��C�F�C��C���C��$C�@C��uC���C�gC�D�C�i�C�~&C���C�n*C�_�C�]C�Z�>���>�,v>��H>�ٟ>�7>�g�>��u>��!>��
>�f>�X�>��(>���>��e>���>��e>��/>���>��>�*�            ;|r�=�~f>�p?/'2?4�g                                                SyF<� @ǆ�    @ǖ     5�iu        4Q��3���3�4��4��                                                                    5/        3�t2�a�3�R4�o3�s�                                                                    1�Y�        /�/b>/�{�0~�0!�                                                                    0�.'        /��.1S�.�"�/�~�/C�y                                                                    6���3�[�71X5w��6��^        5.�4~=�5�*6>5���                                                6��3 �7�3�.��    1Eh+��]6te[�te[�6r�1<�W,��\+Bw�(�l�        &'U�,_��    B��@Q�5                @Qе?	�                ?�G�            A�H�A�ĩA�BuA��<B�NB�>-�                                                    >��B?h!�@�9�Ar�A��B��Bq3�B�e6<�<�<�<�<�<�<�<�<�<�<�<�C�u�C�xC�{�C��|C��1C���C���C���C� kC�&�C�RCC��C��{C��C�/�C�s�C��C��IC�'C�L9C�m�C�oC�`�C�]C�[>�>�+�>�X�>���>��Y>�1>���>��5>��>��_>��\>�r�>�&�>�9�>���>��m>�]�>��6>��>�"�=�I�        =�=�s�>��p?+�D?0��                                                S�F=  @ǖ     @Ǥ     8o�S7�Z�6XRJ5?~3�H;4�"4�4�ѱ                                                                    7�j�6�<�5���4(P42˹�3�E�4ͬ3���                                                                    4��3��1�6�0���/0 0 �0���0f�                                                                    3%W2<?1/�/��5.^yM/B��/��/<��                                                                    6��3ׯ�7��[6�+98}% 8+��7@��5��I4�w�5��s6�>5�/                                                6��V37� ~3(    1p�01 ����Έ6�Έ7��6a�1���0k�-Ѥ)        +���1�L#    B��@Q�5                @Qе?	�                ?�G�            A4.@AxQA�!A��B�mA���                                                        @�'TAT]vA|�A��A�%zBB�Bl;�B���<�<�<�<�<�<�<�<�<�<�<�<�C��)C�ʷC��SC���C��aC��3C��GC�ׅC��
C��C�@�C�iC���C�ŔC���C�;C�zC��C��7C��C�M�C�i�C�a�C�]2C�[>���>�O�>�Eq>�c�>�l�>�w>��G>���>��>�>�>�>�'�>�i�>�=>�>�`9>��>�]3>���>�>�?9��>�Y{>��*>vE�=�1�?��?&��?*C�                                                TAF=� @Ǥ     @ǳ�    8㋝8P��7���7I,5�V 5�4�׏4�s[                                                                    8�~7��;6ӝ�62w`4��^4Hp4�
3�[�                                                                    4xw�3�$36�,2�H�1=I�0�J0�0��                                                                    3��Z3�2g*1��a0o�/��E/��Z/<5�                                                                    6�q&3�884~�7��8��8�L�8��-7�26��6:N6B�5�{�                                                6�0k3�J7��4Q��    1��2���D37D@7���7��23�01ɹ�.EG�        -�[m2�V    B��@Q�5                @Qе?	�                ?�G�    6 _+�o�?��=A�]:��>:�:A��.@�a�                                                        AVlA��EB��B.sIB�Bh�B�3eB��;<�<�<�<�<�<�<�<�<�<�<�<�C��pC��C�P�C�yC��~C���C��WC��C���C��C�2xC�T�C�{C���C���C��C�L�C���C��kC���C�,�C�_�C�a�C�]NC�[>�r>�Q:>�dK>�>�5�>�?�>��u>�}�>���>�ۨ>�Z�>��>�>�`�>�&>��>�[>���>��s>��|?v6�?v�?w��?xة?�?F-�?H�F?L��                                                T�F> @ǳ�    @�    97F68�H7�<�7J�6�WG5�x#5F�/4�f`                                                                    8g�7��C7�.6.k5�g5!YZ4{G�3��                                                                    4� $4�3�o2ܚ525��1�~�0�?08*                                                                    3��H3@�2�;d2S�1esP0�40	5H/@F�                                                                    6��3�V�8���7�C9@��9*�8�+=8'�7��[7��6l�55�[�                                                6�%2�ȅ7��4�=�*�]1�21sh�y97y�8"��8 gh5�4��'1$a*        4���5�*`    B��@Q�5                @Qе?	�                ?�G�    7�%,��                                                                                @���A�#�A�˩B��BN��B��B��)B�Ls<�<�<�<�<�<�<�<�<�<�<�<�C�c�C���C���C�C�C���C���C�A�C���C���C�v�C�H�C�%jC��C��C��C��C�=�C�k/C���C�͛C�#C�R�C�`�C�]eC�['>�O�>�1�>��>�I�>��$>�*>��>���>��>�ف>�D6>�F�>���>�S�>���>��p>��n>�{>�ܝ>��?S�(?UfM?Xr�?]'�?d�1?qn?{��?�                                                  U	F>� @�    @��     9G��8�7=7���7Yn96���64�u5�>O5���                                                                    8|*�7�h7 ��6�S5�u�5d�5	փ4���                                                                    4���4%�3��2�nj2A/$1�"�1nV=1$�                                                                    4	�3P�R2���2��1t�0��0��O0P,�                                                                    6��)3�M8��I7�?9R�99�v8�"8/#\7�x�7H,�6��>6�                                                6���2���7��!4wq-�\1��2.��m��7m�)8/,j84�6���5D��2a�        5L�U6z7�    B��@Q�5                @Qе?	�                ?�G�    5�T�)��R                                                                                @ffAp��A�o�B�EB7܄Bm��B�P�B���<�<�<�<�<�<�<�<�<�<�<�<�C���C�L C���C���C�M�C��TC��FC��C��zC�\C��C��ZC�_�C�GC��}C��7C���C���C���C��C��oC�C�C�^�C�]nC�[2?
5_?�>��]>�J�>��>�u�>��?>���>�)�>�[!>ɢ7>�?u>�j�>��>�W>���>�[�>�Gn>�D�>��??�/?@�?BVp?E?I�n?P�n?[�L?oS9                                                UmF? @��     @��     9G�I8�-�8 ��7]�6�L 68w5ݱi5��F                                                                    8|j7�9�7"��6��*5�y5h��5B4Ѣ�                                                                    4��?4&*�3��K2�h�2D�C1��1r�159L                                                                    4	��3Q�72���2x1x�G0��g0��T0d�                                                                    6���3��8�H7�ʺ9S9;r�8�8�8-��7���7Cߍ6�� 6���                                                6���2��7�tR4��;.8�b1���25V�GA�7GCL8%608I6���5F�Z2A?�        5B�6~�    B��@Q�5                @Qе?	�                ?�G�    2��&P��                                                                                @^@'Ai��A��/A�w]B3H�Bf�?B���B�б<�<�<�<�<�<�<�<�<�<�<�<�C�C��iC�t<C�6�C��>C���C�U�C��C��C�Y�C�
cC��hC�fFC�MC���C�j�C�)C��#C��yC���C���C�6�C�\)C�]`C�[=?��?��?T^?��>���>���>�Bh>�:>�@>��<>��>�k>��C>�j�>��>>�w>��>�R�>��<>��?;��?<-5?=��??�+?Cc�?IL?Q�'?`�                                                U�F?� @��     @���    9=�38�07�D�7[)86�[�6:�5��5�+�                                                                    8o��7���7��6�j�5���5k�5c�4�7�                                                                    4�'J4 ��3�G2�Q�2D��1�1�1w�1&*P                                                                    4�r3K�2���2&91x�e1 U/0��H0Q�                                                                    6��_3���8|,7ޏ�9IU�95�8�
u8(�07�B�7?�i6��6�ϧ                                                6���2��,7���4��v.M7B1��N2B�׶��q6���8Xf8�U6J��5:T1�~j        4ԃ�6)7�    B��@Q�5                @Qе?	�                ?�G�    5�;6*\n�                                                                                @k*KAr_?A�}�B��B6�#Bj�:B�H�B���<�<�<�<�<�<�<�<�<�<�<�<�C�(�C�1C��*C�յC���C���C�j�C�4�C��GC�ǛC��|C�G5C���C���C�[C�CC��[C�n�C�;�C�C��C�/�C�X�C�]7C�[H?O�?y>��4>�m�>�%�>�Q�>��>�'>��>��>��>�C)>�A�>��>�.J>ɏ�>�v
>�>U>�)�>�Q�?A�1?A|x?C%F?D�O?G�S?M2_?V@�?eq�                                                V5F@  @���    @�      9.�8���7�1u7Q�R6�6568o�5�j5R(                                                                    8\��7��a7�6�j�5�Q�5h��5@C4Eu@                                                                    4���4V�3���2��92?Ud1�j21|�0��>                                                                    3�ȏ3=��2��s2��1q�<0�k10��4/ם�                                                                    6�ߙ3؄�8@�y7�[w99�J9*��8��^8g�7�Z`79yi6�o�6$!l                                                6�z3��7���4�X�.2��1���27�Y6"/-�"+�7�7�O5�ɷ4t�'0���        3�1D5��    B��@Q�5                @Qе?	�                ?�G�    6� �+�c                                                                                @v�A|hGA�o�BIwB>,�Bv��B���B�.�<�<�<�<�<�<�<�<�<�<�<�<�C�رC��JC��<C���C��)C�ǡC���C���C��C���C�c8C�?C�VC��C��C�[{C��C��C���C�^�C�2%C�11C�U�C�\�C�[Q>�>�7�>� >�o >��>�
>���>�->��>�a�>�AV>ޮ�>۰�>��>��v>�2N>�j>��=>�>�?k?F�2?G��?I��?L�]?Q�s?Y�6?g��?~�                                                V�F@x @�      @�     9��8h@f7���71`�6���6
n�5�,�5�;                                                                    84�7��i7 5)6`L5��5.��4��q4CeX                                                                    4���3���3]��2��t2"H:1�-@1,��0���                                                                    3ě�3 +.2���1��A1L��0���0Z1�/�]l                                                                    6��3�I�7��7�ȟ93�9^�8���8��7�|;7	�6�� 6!z�                                                6�[-3>'�7���3r-#-4��1��*1��7��h���M4�6��c34�j1�[�-�DI        0>��3T�    B��@Q�5                @Qе?	�                ?�G�    7�!,#rJ?V��6��S                                                                        @��QA�A��@B�BRCB���B�m}B�e)<�<�<�<�<�<�<�<�<�<�<�<�C�r�C���C��C��C�Q�C��BC�ÚC���C�(C�K�C�ibC�C��~C���C�{�C�])C�0�C��C���C���C�Z�C�9�C�S�C�\�C�[X>��>��=>��>���>��>�E�>��e>�La>�>�\2>�FW>ѩg>�a5>�Y�>�Y.>�O%>�i�>��>ř>S?T�?Y?\d�?_�U?e.�?m]�?v�w?u�                                                V�F@� @�     @��    8�M48C��7��l7]�6i	#5��5��4�9�                                                                    8#I7w<v6���6FǏ5�.14��4B"f3д�                                                                    4p��3ո�3J�2���1�z�1E$�0��04k�                                                                    3��?3�M2�41��1 �N0y�/�;/c�y                                                                    6���3ܣR7	$�7�Y8��8��8�'�7��,7R06�`6"��5��                                                6���3	�+7���22�;+=B1��0�L7�+2��+/���[5n�*0ԉ/���'��        +�:�0��B    B��@Q�5                @Qе?	�                ?�G�    7`&�+N
g@��?�9                                                                        @�t�A���B�*B.9�Bo��B���B�-%B�Ph<�<�<�<�<�<�<�<�<�<�<�<�C�HUC���C�ŗC� �C�<4C�{DC��[C���C�+�C�W C��C���C��}C���C���C���C��BC���C��oC��gC�{C�F�C�R�C�\9C�[^>��>�A=>�rt>��>��8>�:�>���>�S�>�O�>�ސ>�Ru>ÓC>ŋ1>�AR>Ȑ�>�7�>�
�>�h>�D?>�$
?6Л?k��?t�(?w��?{N�?~)?�?�                                                  WaFAl @��    @�-�    7�6g��7��6�ީ6E5��[5	>4���                                                                    6G�B5�D65a�6 5x�4İs4-�3͹�                                                                    2���1��2���2^��1�#1*D0���01�l                                                                    1�	c1�O1�"1���1�G0V� /�~/`��                                                                    6�'�3�
!6l�6�Wo72�7�7��7��g73�6��6�?5��                                                6�ؙ39\7�i5        1Gڗ    7�?D��?D�R�                                    B��@Q�5                @Qе?	�                ?�G�    6YM*�}�A��A�4�A�?��                                                                @�@A�A�A�B��B9��Bp�kB��uBĞ<�<�<�<�<�<�<�<�<�<�<�<�C�8*C�U�C���C���C��qC�3�C�t�C���C���C�4�C�q>C���C��C�"DC�X�C���C��1C��yC��QC���C���C�TwC�R�C�[�C�[a>�@->��'>�;>��8>�r>�n�>��>�r�>��>��>�k�>��M>�B>���>��>��>Ô>>ą�>Ă�>ø�>�ZK>��?��?B��?K�S?S��?`_�?t~:                                                z%FA� @�-�    @�=     