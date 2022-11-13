CDF      
      time       bnds      lndgrid       levsoi        levdcmp       cft       glc_nec    
   ltype      	   natpft        levlak     
   nvegwcs       string_length         levgrnd       hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      Sun Jan  9 16:23:29 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1939.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1939.nc
created on 12/10/21 16:55:49    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:29 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1939.nc had following "history" attribute:
created on 12/10/21 16:55:49
     NCO       4.6.9        5   time                standard_name         time   	long_name         time   bounds        time_bounds    units         days since 1850-01-01 00:00:00     calendar      noleap     axis      T              L@   	time_bnds                                LD   levsoi                 	long_name         Dcoordinate soil levels (equivalent to top nlevsoi levels of levgrnd)   units         m      axis      Y         P     J�   levdcmp                	long_name         2coordinate levels for soil decomposition variables     units         m      axis      Y         d     J�   levlak        	         	long_name         coordinate lake levels     units         m      axis      Y         (     K`   mcdate                  	long_name         current date (YYYYMMDD)            LT   mcsec                   	long_name         current seconds of current date    units         s              LX   mdcur                   	long_name         current day (from base day)            L\   mscur                   	long_name         current seconds of current day             L`   nstep                   	long_name         	time step              Ld   lon                	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��           K�   lat                	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��           K�   area               	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��           K�   landfrac               	long_name         land fraction      
_FillValue        {@��   missing_value         {@��           K�   landmask               	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����           K�   pftmask                	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����           K�   nbedrock               	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����           K�   ACTUAL_IMMOB                   	long_name         actual N immobilization    units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Lh   AGNPP                      	long_name         aboveground NPP    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Ll   ALT                    	long_name         current active layer thickness     units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Lp   ALTMAX                     	long_name         %maximum annual active layer thickness      units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Lt   AR                     	long_name         !autotrophic respiration (MR + GR)      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Lx   ATM_TOPO                   	long_name         atmospheric surface height     units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L|   BAF_CROP                   	long_name         fractional area burned for crop    units         s-1    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   	BAF_PEATF                      	long_name         "fractional area burned in peatland     units         s-1    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   BCDEP                      	long_name         -total BC deposition (dry+wet) from atmosphere      units         kg/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   BGNPP                      	long_name         belowground NPP    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   BTRAN2                     	long_name         root zone soil wetness factor      units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            L�   BTRANMN                    	long_name         *daily minimum of transpiration beta factor     units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            L�   CH4PROD                    	long_name          Gridcell total production of CH4   units         gC/m2/s    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   CH4_SURF_AERE_SAT                      	long_name         :aerenchyma surface CH4 flux for inundated area; (+ to atm)     units         mol/m2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   CH4_SURF_AERE_UNSAT                    	long_name         >aerenchyma surface CH4 flux for non-inundated area; (+ to atm)     units         mol/m2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   CH4_SURF_DIFF_SAT                      	long_name         @diffusive surface CH4 flux for inundated / lake area; (+ to atm)   units         mol/m2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   CH4_SURF_DIFF_UNSAT                    	long_name         =diffusive surface CH4 flux for non-inundated area; (+ to atm)      units         mol/m2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   CH4_SURF_EBUL_SAT                      	long_name         Aebullition surface CH4 flux for inundated / lake area; (+ to atm)      units         mol/m2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   CH4_SURF_EBUL_UNSAT                    	long_name         >ebullition surface CH4 flux for non-inundated area; (+ to atm)     units         mol/m2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   COL_FIRE_CLOSS                     	long_name         Ttotal column-level fire C loss for non-peat fires outside land-type converted region   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   COL_FIRE_NLOSS                     	long_name         total column-level fire N loss     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            L�   CONC_O2_SAT                       	long_name         /O2 soil Concentration for inundated / lake area    units         mol/m3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg_plus_lake         P     L�   CONC_O2_UNSAT                         	long_name         ,O2 soil Concentration for non-inundated area   units         mol/m3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       P     M   COST_NACTIVE                   	long_name         Cost of active uptake      units         gN/gC      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            M\   	COST_NFIX                      	long_name         Cost of fixation   units         gN/gC      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            M`   COST_NRETRANS                      	long_name         Cost of retranslocation    units         gN/gC      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Md   CPOOL                      	long_name         temporary photosynthate C pool     units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Mh   
CROPPROD1C                     	long_name         #1-yr crop product (grain+biofuel) C    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Ml   CROPPROD1C_LOSS                    	long_name          loss from 1-yr crop product pool   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Mp   
CROPPROD1N                     	long_name         #1-yr crop product (grain+biofuel) N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Mt   CROPPROD1N_LOSS                    	long_name          loss from 1-yr crop product pool   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Mx   CWDC                   	long_name         CWD C      units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            M|   	CWDC_LOSS                      	long_name         coarse woody debris C loss     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            M�   CWDC_vr                       	long_name         CWD C (vertically resolved)    units         gC/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     M�   CWDN                   	long_name         CWD N      units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            M�   CWDN_vr                       	long_name         CWD N (vertically resolved)    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     M�   
DEADCROOTC                     	long_name         dead coarse root C     units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N<   
DEADCROOTN                     	long_name         dead coarse root N     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N@   	DEADSTEMC                      	long_name         dead stem C    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ND   	DEADSTEMN                      	long_name         dead stem N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            NH   DENIT                      	long_name         total rate of denitrification      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            NL   DISPVEGC                   	long_name         1displayed veg carbon, excluding storage and cpool      units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            NP   DISPVEGN                   	long_name         displayed vegetation nitrogen      units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            NT   DSL                    	long_name         dry surface layer thickness    units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            NX   DSTDEP                     	long_name         /total dust deposition (dry+wet) from atmosphere    units         kg/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N\   DSTFLXT                    	long_name         total surface dust emission    units         kg/m2/s    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N`   DWT_CONV_CFLUX                     	long_name         Xconversion C flux (immediate loss to atm) (0 at all times except first timestep of year)   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Nd   DWT_CONV_CFLUX_DRIBBLED                    	long_name         Gconversion C flux (immediate loss to atm), dribbled throughout the year    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Nh   DWT_CONV_NFLUX                     	long_name         Xconversion N flux (immediate loss to atm) (0 at all times except first timestep of year)   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Nl   DWT_CROPPROD1C_GAIN                    	long_name         <landcover change-driven addition to 1-year crop product pool   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Np   DWT_CROPPROD1N_GAIN                    	long_name         <landcover change-driven addition to 1-year crop product pool   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Nt   DWT_SEEDN_TO_DEADSTEM                      	long_name         #seed source to patch-level deadstem    units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Nx   DWT_SEEDN_TO_LEAF                      	long_name         seed source to patch-level leaf    units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N|   DWT_SLASH_CFLUX                    	long_name         Wslash C flux (to litter diagnostic only) (0 at all times except first timestep of year)    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   DWT_WOODPRODC_GAIN                     	long_name         6landcover change-driven addition to wood product pools     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   DWT_WOODPRODN_GAIN                     	long_name         6landcover change-driven addition to wood product pools     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   	EFLXBUILD                      	long_name         Cbuilding heat flux from change in interior building air temperature    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   EFLX_DYNBAL                    	long_name         0dynamic land cover change conversion energy flux   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   EFLX_GRND_LAKE                     	long_name         Bnet heat flux into lake/snow surface, excluding light transmission     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   EFLX_LH_TOT                    	long_name         !total latent heat flux [+ to atm]      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   EFLX_LH_TOT_R                      	long_name         Rural total evaporation    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   ELAI                   	long_name         !exposed one-sided leaf area index      units         m^2/m^2    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   ER                     	long_name         8total ecosystem respiration, autotrophic + heterotrophic   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   ERRH2O                     	long_name         total water conservation error     units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   	ERRH2OSNO                      	long_name         &imbalance in snow depth (liquid water)     units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   ERRSEB                     	long_name         !surface energy conservation error      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   ERRSOI                     	long_name         #soil/lake energy conservation error    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   ERRSOL                     	long_name         "solar radiation conservation error     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   ESAI                   	long_name         !exposed one-sided stem area index      units         m^2/m^2    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FAREA_BURNED                   	long_name         timestep fractional area burned    units         s-1    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FCEV                   	long_name         canopy evaporation     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FCH4                   	long_name         2Gridcell surface CH4 flux to atmosphere (+ to atm)     units         kgC/m2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   	FCH4TOCO2                      	long_name          Gridcell oxidation of CH4 to CO2   units         gC/m2/s    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   
FCH4_DFSAT                     	long_name         SCH4 additional flux due to changing fsat, natural vegetated and crop landunits only    units         kgC/m2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FCOV                   	long_name         fractional impermeable area    units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            N�   FCTR                   	long_name         canopy transpiration   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FFIX_TO_SMINN                      	long_name         )free living  N fixation to soil mineral N      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FGEV                   	long_name         ground evaporation     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FGR                    	long_name         Oheat flux into soil/snow including snow melt and lake / snow light transmission    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FGR12                      	long_name         %heat flux between soil layers 1 and 2      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FH2OSFC                    	long_name         +fraction of ground covered by surface water    units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   
FINUNDATED                     	long_name         .fractional inundated area of vegetated columns     units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            N�   FIRA                   	long_name         !net infrared (longwave) radiation      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FIRA_R                     	long_name         'Rural net infrared (longwave) radiation    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FIRE                   	long_name         %emitted infrared (longwave) radiation      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            N�   FIRE_R                     	long_name         +Rural emitted infrared (longwave) radiation    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O    FLDS                   	long_name         Iatmospheric longwave radiation (downscaled to columns in glacier regions)      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O   FPI                    	long_name         $fraction of potential immobilization   units         
proportion     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O   FPSN                   	long_name         photosynthesis     units         umol m-2 s-1   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O   FREE_RETRANSN_TO_NPOOL                     	long_name         deployment of retranslocated N     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O   FROOTC                     	long_name         fine root C    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O   FROOTC_ALLOC                   	long_name         fine root C allocation     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O   FROOTC_LOSS                    	long_name         fine root C loss   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O   FROOTN                     	long_name         fine root N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O    FSA                    	long_name         absorbed solar radiation   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O$   FSAT                   	long_name         +fractional area with water table at surface    units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            O(   FSDS                   	long_name         $atmospheric incident solar radiation   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O,   FSDSND                     	long_name         #direct nir incident solar radiation    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O0   FSDSNDLN                   	long_name         1direct nir incident solar radiation at local noon      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O4   FSDSNI                     	long_name         $diffuse nir incident solar radiation   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O8   FSDSVD                     	long_name         #direct vis incident solar radiation    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O<   FSDSVDLN                   	long_name         1direct vis incident solar radiation at local noon      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O@   FSDSVI                     	long_name         $diffuse vis incident solar radiation   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            OD   FSDSVILN                   	long_name         2diffuse vis incident solar radiation at local noon     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            OH   FSH                    	long_name         Ssensible heat not including correction for land use change and rain/snow conversion    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            OL   FSH_G                      	long_name         sensible heat from ground      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            OP   FSH_PRECIP_CONVERSION                      	long_name         ;Sensible heat flux from conversion of rain/snow atm forcing    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            OT   FSH_R                      	long_name         Rural sensible heat    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            OX   FSH_RUNOFF_ICE_TO_LIQ                      	long_name         Dsensible heat flux generated from conversion of ice runoff to liquid   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O\   FSH_TO_COUPLER                     	long_name         �sensible heat sent to coupler (includes corrections for land use change, rain/snow conversion and conversion of ice runoff to liquid)      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O`   FSH_V                      	long_name         sensible heat from veg     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Od   FSM                    	long_name         snow melt heat flux    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Oh   FSNO                   	long_name         "fraction of ground covered by snow     units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Ol   FSNO_EFF                   	long_name         ,effective fraction of ground covered by snow   units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Op   FSR                    	long_name         reflected solar radiation      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Ot   FSRND                      	long_name         $direct nir reflected solar radiation   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Ox   FSRNDLN                    	long_name         2direct nir reflected solar radiation at local noon     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O|   FSRNI                      	long_name         %diffuse nir reflected solar radiation      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   FSRVD                      	long_name         $direct vis reflected solar radiation   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   FSRVDLN                    	long_name         2direct vis reflected solar radiation at local noon     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   FSRVI                      	long_name         %diffuse vis reflected solar radiation      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   FUELC                      	long_name         	fuel load      units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   F_DENIT                    	long_name         denitrification flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   F_N2O_DENIT                    	long_name         denitrification N2O flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   	F_N2O_NIT                      	long_name         nitrification N2O flux     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   F_NIT                      	long_name         nitrification flux     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   GPP                    	long_name         gross primary production   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   GR                     	long_name         total growth respiration   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   
GROSS_NMIN                     	long_name         gross rate of N mineralization     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   GSSHA                      	long_name          shaded leaf stomatal conductance   units         umol H20/m2/s      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   GSSHALN                    	long_name         .shaded leaf stomatal conductance at local noon     units         umol H20/m2/s      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   GSSUN                      	long_name          sunlit leaf stomatal conductance   units         umol H20/m2/s      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   GSSUNLN                    	long_name         .sunlit leaf stomatal conductance at local noon     units         umol H20/m2/s      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   H2OCAN                     	long_name         intercepted water      units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   H2OSFC                     	long_name         surface water depth    units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   H2OSNO                     	long_name         snow depth (liquid water)      units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   
H2OSNO_TOP                     	long_name         mass of snow in top snow layer     units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            O�   H2OSOI                        	long_name         Avolumetric soil water (natural vegetated and crop landunits only)      units         mm3/mm3    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       P     O�   HEAT_CONTENT1                      	long_name         #initial gridcell total heat content    units         J/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P    HEAT_FROM_AC                   	long_name         Lsensible heat flux put into canyon due to heat removed from air conditioning   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P$   HIA                    	long_name         2 m NWS Heat Index     units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P(   HIA_R                      	long_name         Rural 2 m NWS Heat Index   units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P,   HIA_U                      	long_name         Urban 2 m NWS Heat Index   units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P0   HR                     	long_name         total heterotrophic respiration    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P4   HR_vr                         	long_name         3total vertically resolved heterotrophic respiration    units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     P8   HTOP                   	long_name         
canopy top     units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   HUMIDEX                    	long_name         2 m Humidex    units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   	HUMIDEX_R                      	long_name         Rural 2 m Humidex      units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   	HUMIDEX_U                      	long_name         Urban 2 m Humidex      units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   ICE_CONTENT1                   	long_name         "initial gridcell total ice content     units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   JMX25T                     	long_name         canopy profile of jmax     units         	umol/m2/s      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   Jmx25Z                     	long_name         Bmaximum rate of electron transport at 25 Celcius for canopy layers     units         umol electrons/m2/s    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LAISHA                     	long_name          shaded projected leaf area index   units         m^2/m^2    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LAISUN                     	long_name          sunlit projected leaf area index   units         m^2/m^2    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LAKEICEFRAC_SURF                   	long_name         $surface lake layer ice mass fraction   units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LAKEICETHICK                   	long_name         @thickness of lake ice (including physical expansion on freezing)   units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LAND_USE_FLUX                      	long_name         |total C emitted from land cover conversion (smoothed over the year) and wood and grain product pools (NOTE: not a net value)   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LEAFC                      	long_name         leaf C     units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LEAFCN                     	long_name         "Leaf CN ratio used for flexible CN     units         gC/gN      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LEAFC_ALLOC                    	long_name         leaf C allocation      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LEAFC_CHANGE                   	long_name         C change in leaf   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   
LEAFC_LOSS                     	long_name         leaf C loss    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LEAFC_TO_LITTER_FUN                    	long_name         leaf C litterfall used by FUN      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LEAFN                      	long_name         leaf N     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LEAFN_TO_LITTER                    	long_name         leaf N litterfall      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LEAF_MR                    	long_name         leaf maintenance respiration   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LFC2                   	long_name         3conversion area fraction of BET and BDT that burned    units         per sec    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LIQCAN                     	long_name         intercepted liquid water   units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LIQUID_CONTENT1                    	long_name         "initial gridcell total liq content     units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LITFALL                    	long_name         "litterfall (leaves and fine roots)     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   LITR1C                     	long_name         LITR1 C    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            P�   	LITR1C_vr                         	long_name         LITR1 C (vertically resolved)      units         gC/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     P�   LITR1N                     	long_name         LITR1 N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Q@   	LITR1N_vr                         	long_name         LITR1 N (vertically resolved)      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     QD   LITR2C                     	long_name         LITR2 C    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Q�   	LITR2C_vr                         	long_name         LITR2 C (vertically resolved)      units         gC/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     Q�   LITR2N                     	long_name         LITR2 N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Q�   	LITR2N_vr                         	long_name         LITR2 N (vertically resolved)      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     R    LITR3C                     	long_name         LITR3 C    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Rd   	LITR3C_vr                         	long_name         LITR3 C (vertically resolved)      units         gC/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     Rh   LITR3N                     	long_name         LITR3 N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            R�   	LITR3N_vr                         	long_name         LITR3 N (vertically resolved)      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     R�   
LITTERC_HR                     	long_name         "litter C heterotrophic respiration     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S    LITTERC_LOSS                   	long_name         litter C loss      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S$   
LIVECROOTC                     	long_name         live coarse root C     units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S(   
LIVECROOTN                     	long_name         live coarse root N     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S,   	LIVESTEMC                      	long_name         live stem C    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S0   	LIVESTEMN                      	long_name         live stem N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S4   LNC                    	long_name         leaf N concentration   units         gN leaf/m^2    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S8   MEG_acetaldehyde                   	long_name         
MEGAN flux     units         	kg/m2/sec      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S<   MEG_acetic_acid                    	long_name         
MEGAN flux     units         	kg/m2/sec      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S@   MEG_acetone                    	long_name         
MEGAN flux     units         	kg/m2/sec      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            SD   MEG_carene_3                   	long_name         
MEGAN flux     units         	kg/m2/sec      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            SH   MEG_ethanol                    	long_name         
MEGAN flux     units         	kg/m2/sec      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            SL   MEG_formaldehyde                   	long_name         
MEGAN flux     units         	kg/m2/sec      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            SP   MEG_isoprene                   	long_name         
MEGAN flux     units         	kg/m2/sec      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ST   MEG_methanol                   	long_name         
MEGAN flux     units         	kg/m2/sec      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            SX   MEG_pinene_a                   	long_name         
MEGAN flux     units         	kg/m2/sec      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S\   MEG_thujene_a                      	long_name         
MEGAN flux     units         	kg/m2/sec      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S`   MR                     	long_name         maintenance respiration    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Sd   NACTIVE                    	long_name         Mycorrhizal N uptake flux      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Sh   NACTIVE_NH4                    	long_name         Mycorrhizal N uptake flux      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Sl   NACTIVE_NO3                    	long_name         Mycorrhizal N uptake flux      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Sp   NAM                    	long_name         AM-associated N uptake flux    units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            St   NAM_NH4                    	long_name         AM-associated N uptake flux    units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Sx   NAM_NO3                    	long_name         AM-associated N uptake flux    units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S|   NBP                    	long_name         �net biome production, includes fire, landuse, harvest and hrv_xsmrpool flux (latter smoothed over the year), positive for sink (same as net carbon exchange between land and atmosphere)   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NDEPLOY                    	long_name         total N deployed in new growth     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NDEP_TO_SMINN                      	long_name         *atmospheric N deposition to soil mineral N     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NECM                   	long_name         ECM-associated N uptake flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NECM_NH4                   	long_name         ECM-associated N uptake flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NECM_NO3                   	long_name         ECM-associated N uptake flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NEE                    	long_name         �net ecosystem exchange of carbon, includes fire and hrv_xsmrpool (latter smoothed over the year), excludes landuse and harvest flux, positive for source   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NEM                    	long_name         TGridcell net adjustment to net carbon exchange passed to atm. for methane production   units         gC/m2/s    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NEP                    	long_name         Unet ecosystem production, excludes fire, landuse, and harvest flux, positive for sink      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NET_NMIN                   	long_name         net rate of N mineralization   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NFIRE                      	long_name         fire counts valid only in Reg.C    units         counts/km2/sec     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NFIX                   	long_name         Symbiotic BNF uptake flux      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NNONMYC                    	long_name         Non-mycorrhizal N uptake flux      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NNONMYC_NH4                    	long_name         Non-mycorrhizal N uptake flux      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NNONMYC_NO3                    	long_name         Non-mycorrhizal N uptake flux      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPASSIVE                   	long_name         Passive N uptake flux      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPOOL                      	long_name         temporary plant N pool     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP                    	long_name         net primary production     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   
NPP_GROWTH                     	long_name         Total C used for growth in FUN     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NACTIVE                    	long_name         Mycorrhizal N uptake used C    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NACTIVE_NH4                    	long_name         Mycorrhizal N uptake use C     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NACTIVE_NO3                    	long_name         Mycorrhizal N uptake used C    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NAM                    	long_name         AM-associated N uptake used C      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NAM_NH4                    	long_name         AM-associated N uptake use C   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NAM_NO3                    	long_name         AM-associated N uptake use C   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NECM                   	long_name         ECM-associated N uptake used C     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NECM_NH4                   	long_name         ECM-associated N uptake use C      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NECM_NO3                   	long_name         ECM-associated N uptake used C     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NFIX                   	long_name         Symbiotic BNF uptake used C    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NNONMYC                    	long_name         Non-mycorrhizal N uptake used C    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NNONMYC_NH4                    	long_name         Non-mycorrhizal N uptake use C     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NNONMYC_NO3                    	long_name         Non-mycorrhizal N uptake use C     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            S�   NPP_NRETRANS                   	long_name         Retranslocated N uptake flux   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T    NPP_NUPTAKE                    	long_name         Total C used by N uptake in FUN    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T   NRETRANS                   	long_name         Retranslocated N uptake flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T   NRETRANS_REG                   	long_name         Retranslocated N uptake flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T   NRETRANS_SEASON                    	long_name         Retranslocated N uptake flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T   NRETRANS_STRESS                    	long_name         Retranslocated N uptake flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T   NUPTAKE                    	long_name         Total N uptake of FUN      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T   NUPTAKE_NPP_FRACTION                   	long_name         frac of NPP used in N uptake   units         -      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T   OCDEP                      	long_name         -total OC deposition (dry+wet) from atmosphere      units         kg/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T    O_SCALAR                      	long_name         8fraction by which decomposition is reduced due to anoxia   units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     T$   PARVEGLN                   	long_name         (absorbed par by vegetation at local noon   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Tt   PBOT                   	long_name         Jatmospheric pressure at surface (downscaled to columns in glacier regions)     units         Pa     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Tx   PCH4                   	long_name         #atmospheric partial pressure of CH4    units         Pa     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T|   PCO2                   	long_name         #atmospheric partial pressure of CO2    units         Pa     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T�   PCT_CFT                       	long_name         #% of each crop on the crop landunit    units         %      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            T�   PCT_GLC_MEC                       	long_name         5% of each GLC elevation class on the glacier landunit      units         %      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       (     T�   PCT_LANDUNIT                      	long_name         % of each landunit on grid cell    units         %      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       $     T�   PCT_NAT_PFT                       	long_name         =% of each PFT on the natural vegetation (i.e., soil) landunit      units         %      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       <     T�   PFT_FIRE_CLOSS                     	long_name         Stotal patch-level fire C loss for non-peat fires outside land-type converted region    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U   PFT_FIRE_NLOSS                     	long_name         total patch-level fire N loss      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U   PLANT_NDEMAND                      	long_name         &N flux required to support initial GPP     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U   POTENTIAL_IMMOB                    	long_name         potential N immobilization     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U    POT_F_DENIT                    	long_name         potential denitrification flux     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U$   	POT_F_NIT                      	long_name         potential nitrification flux   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U(   PSNSHA                     	long_name         shaded leaf photosynthesis     units         umolCO2/m^2/s      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U,   PSNSHADE_TO_CPOOL                      	long_name         C fixation from shaded canopy      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U0   PSNSUN                     	long_name         sunlit leaf photosynthesis     units         umolCO2/m^2/s      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U4   PSNSUN_TO_CPOOL                    	long_name         C fixation from sunlit canopy      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U8   Q2M                    	long_name         2m specific humidity   units         kg/kg      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U<   QBOT                   	long_name         Hatmospheric specific humidity (downscaled to columns in glacier regions)   units         kg/kg      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U@   QDRAI                      	long_name         sub-surface drainage   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            UD   QDRAI_PERCH                    	long_name         perched wt drainage    units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            UH   QDRAI_XS                   	long_name         saturation excess drainage     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            UL   QFLOOD                     	long_name         runoff from river flooding     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            UP   QFLX_EVAP_TOT                      	long_name         -qflx_evap_soi + qflx_evap_can + qflx_tran_veg      units         
kg m-2 s-1     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            UT   QFLX_ICE_DYNBAL                    	long_name         4ice dynamic land cover change conversion runoff flux   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            UX   QFLX_LIQDEW_TO_TOP_LAYER                   	long_name         >rate of liquid water deposited on top soil or snow layer (dew)     units         mm H2O/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U\   QFLX_LIQEVAP_FROM_TOP_LAYER                    	long_name         ;rate of liquid water evaporated from top soil or snow layer    units         mm H2O/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U`   QFLX_LIQ_DYNBAL                    	long_name         4liq dynamic land cover change conversion runoff flux   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Ud   QFLX_SNOW_DRAIN                    	long_name         drainage from snow pack    units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Uh   QFLX_SNOW_DRAIN_ICE                    	long_name         1drainage from snow pack melt (ice landunits only)      units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         ice            Ul   QFLX_SOLIDDEW_TO_TOP_LAYER                     	long_name         ?rate of solid water deposited on top soil or snow layer (frost)    units         mm H2O/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Up   QFLX_SOLIDEVAP_FROM_TOP_LAYER                      	long_name         zrate of ice evaporated from top soil or snow layer (sublimation) (also includes bare ice sublimation from glacier columns)     units         mm H2O/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Ut   QH2OSFC                    	long_name         surface water runoff   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Ux   QHR                    	long_name         hydraulic redistribution   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            U|   QICE                   	long_name         ice growth/melt    units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         ice            U�   QICE_FRZ                   	long_name         
ice growth     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         ice            U�   	QICE_MELT                      	long_name         ice melt   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         ice            U�   QINFL                      	long_name         infiltration   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QINTR                      	long_name         interception   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QIRRIG_FROM_GW_CONFINED                    	long_name         3water added through confined groundwater irrigation    units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QIRRIG_FROM_GW_UNCONFINED                      	long_name         5water added through unconfined groundwater irrigation      units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QIRRIG_FROM_SURFACE                    	long_name         ,water added through surface water irrigation   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QOVER                      	long_name         'total surface runoff (includes QH2OSFC)    units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QRGWL                      	long_name         isurface runoff at glaciers (liquid only), wetlands, lakes; also includes melted ice runoff from QSNWCPICE      units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QRUNOFF                    	long_name         @total liquid runoff not including correction for land use change   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QRUNOFF_ICE                    	long_name         Btotal liquid runoff not incl corret for LULCC (ice landunits only)     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         ice            U�   QRUNOFF_ICE_TO_COUPLER                     	long_name         Ktotal ice runoff sent to coupler (includes corrections for land use change)    units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QRUNOFF_TO_COUPLER                     	long_name         Ntotal liquid runoff sent to coupler (includes corrections for land use change)     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   	QSNOCPLIQ                      	long_name         Rexcess liquid h2o due to snow capping not including correction for land use change     units         mm H2O/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QSNOEVAP                   	long_name         Nevaporation from snow (only when snl<0, otherwise it is equal to qflx_ev_soil)     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QSNOFRZ                    	long_name         $column-integrated snow freezing rate   units         kg/m2/s    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QSNOFRZ_ICE                    	long_name         9column-integrated snow freezing rate (ice landunits only)      units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         ice            U�   QSNOMELT                   	long_name         snow melt rate     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QSNOMELT_ICE                   	long_name         snow melt (ice landunits only)     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         ice            U�   
QSNOUNLOAD                     	long_name         canopy snow unloading      units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QSNO_TEMPUNLOAD                    	long_name         canopy snow temp unloading     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QSNO_WINDUNLOAD                    	long_name         canopy snow wind unloading     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   	QSNWCPICE                      	long_name         Qexcess solid h2o due to snow capping not including correction for land use change      units         mm H2O/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QSOIL                      	long_name         HGround evaporation (soil/snow evaporation + soil/snow sublimation - dew)   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   	QSOIL_ICE                      	long_name         'Ground evaporation (ice landunits only)    units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         ice            U�   QVEGE                      	long_name         canopy evaporation     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   QVEGT                      	long_name         canopy transpiration   units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   RAIN                   	long_name         Eatmospheric rain, after rain/snow repartitioning based on temperature      units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   RAIN_FROM_ATM                      	long_name         >atmospheric rain received from atmosphere (pre-repartitioning)     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   RETRANSN                   	long_name         plant pool of retranslocated N     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   RETRANSN_TO_NPOOL                      	long_name         deployment of retranslocated N     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            U�   RH2M                   	long_name         2m relative humidity   units         %      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V    RR                     	long_name         /root respiration (fine root MR + total root GR)    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V   RSSHA                      	long_name         shaded leaf stomatal resistance    units         s/m    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            V   RSSUN                      	long_name         sunlit leaf stomatal resistance    units         s/m    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            V   SABG                   	long_name         solar rad absorbed by ground   units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V   SABG_PEN                   	long_name         2Rural solar rad penetrating top soil or snow layer     units         watt/m^2   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V   SABV                   	long_name         solar rad absorbed by veg      units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V   SEEDC                      	long_name         /pool for seeding new PFTs via dynamic landcover    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V   SEEDN                      	long_name         /pool for seeding new PFTs via dynamic landcover    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V    SLASH_HARVESTC                     	long_name          slash harvest carbon (to litter)   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V$   SMINN                      	long_name         soil mineral N     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V(   SMINN_TO_NPOOL                     	long_name         #deployment of soil mineral N uptake    units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V,   SMINN_TO_PLANT                     	long_name         plant uptake of soil mineral N     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V0   SMINN_TO_PLANT_FUN                     	long_name         Total soil N uptake of FUN     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V4   SMINN_vr                      	long_name         soil mineral N     units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     V8   SMIN_NH4                   	long_name         soil mineral NH4   units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V�   SMIN_NH4_vr                       	long_name         soil mineral NH4 (vert. res.)      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     V�   SMIN_NO3                   	long_name         soil mineral NO3   units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V�   SMIN_NO3_LEACHED                   	long_name         soil NO3 pool loss to leaching     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V�   SMIN_NO3_RUNOFF                    	long_name         soil NO3 pool loss to runoff   units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            V�   SMIN_NO3_vr                       	long_name         soil mineral NO3 (vert. res.)      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     V�   SMP                       	long_name         Asoil matric potential (natural vegetated and crop landunits only)      units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       d     W8   SNOBCMCL                   	long_name         mass of BC in snow column      units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOBCMSL                   	long_name         mass of BC in top snow layer   units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOCAN                     	long_name         intercepted snow   units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   	SNODSTMCL                      	long_name         mass of dust in snow column    units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   	SNODSTMSL                      	long_name         mass of dust in top snow layer     units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOFSRND                   	long_name         .direct nir reflected solar radiation from snow     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOFSRNI                   	long_name         /diffuse nir reflected solar radiation from snow    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOFSRVD                   	long_name         .direct vis reflected solar radiation from snow     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOFSRVI                   	long_name         /diffuse vis reflected solar radiation from snow    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   	SNOINTABS                      	long_name         8Fraction of incoming solar absorbed by lower snow layers   units         -      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOOCMCL                   	long_name         mass of OC in snow column      units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOOCMSL                   	long_name         mass of OC in top snow layer   units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   	SNOTXMASS                      	long_name         ksnow temperature times layer mass, layer sum; to get mass-weighted temperature, divide by (SNOWICE+SNOWLIQ)    units         K kg/m2    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOW                   	long_name         Eatmospheric snow, after rain/snow repartitioning based on temperature      units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOWDP                     	long_name         gridcell mean snow height      units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOWICE                    	long_name         snow ice   units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOWLIQ                    	long_name         snow liquid water      units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   
SNOW_DEPTH                     	long_name          snow height of snow covered area   units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOW_FROM_ATM                      	long_name         >atmospheric snow received from atmosphere (pre-repartitioning)     units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOW_PERSISTENCE                   	long_name         BLength of time of continuous snow cover (nat. veg. landunits only)     units         seconds    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         natveg             W�   
SNOW_SINKS                     	long_name         snow sinks (liquid water)      units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SNOW_SOURCES                   	long_name         snow sources (liquid water)    units         mm/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   SOIL1C                     	long_name         SOIL1 C    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            W�   	SOIL1C_vr                         	long_name         SOIL1 C (vertically resolved)      units         gC/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     W�   SOIL1N                     	long_name         SOIL1 N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            XH   	SOIL1N_vr                         	long_name         SOIL1 N (vertically resolved)      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     XL   SOIL2C                     	long_name         SOIL2 C    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            X�   	SOIL2C_vr                         	long_name         SOIL2 C (vertically resolved)      units         gC/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     X�   SOIL2N                     	long_name         SOIL2 N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Y   	SOIL2N_vr                         	long_name         SOIL2 N (vertically resolved)      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     Y   SOIL3C                     	long_name         SOIL3 C    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Yl   	SOIL3C_vr                         	long_name         SOIL3 C (vertically resolved)      units         gC/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     Yp   SOIL3N                     	long_name         SOIL3 N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Y�   	SOIL3N_vr                         	long_name         SOIL3 N (vertically resolved)      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     Y�   SOILC_CHANGE                   	long_name         C change in soil   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Z(   SOILC_HR                   	long_name          soil C heterotrophic respiration   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            Z,   SOILC_vr                      	long_name         SOIL C (vertically resolved)   units         gC/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     Z0   SOILICE                       	long_name         4soil ice (natural vegetated and crop landunits only)   units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       P     Z�   SOILLIQ                       	long_name         =soil liquid water (natural vegetated and crop landunits only)      units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       P     Z�   SOILN_vr                      	long_name         SOIL N (vertically resolved)   units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     [    	SOILRESIS                      	long_name         soil resistance to evaporation     units         s/m    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   SOILWATER_10CM                     	long_name         @soil liquid water + ice in top 10cm of soil (veg landunits only)   units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            [�   	SOMC_FIRE                      	long_name         C loss due to peat burning     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   SOM_C_LEACHED                      	long_name         .total flux of C from SOM pools due to leaching     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   SR                     	long_name         'total soil respiration (HR + root resp)    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   STORVEGC                   	long_name         )stored vegetation carbon, excluding cpool      units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   STORVEGN                   	long_name         stored vegetation nitrogen     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   SUPPLEMENT_TO_SMINN                    	long_name         supplemental N supply      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   SWBGT                      	long_name         !2 m Simplified Wetbulb Globe Temp      units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   SWBGT_R                    	long_name         'Rural 2 m Simplified Wetbulb Globe Temp    units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   SWBGT_U                    	long_name         'Urban 2 m Simplified Wetbulb Globe Temp    units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   TAUX                   	long_name         zonal surface stress   units         kg/m/s^2   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   TAUY                   	long_name         meridional surface stress      units         kg/m/s^2   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   TBOT                   	long_name         Fatmospheric air temperature (downscaled to columns in glacier regions)     units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   TBUILD                     	long_name         'internal urban building air temperature    units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   TG                     	long_name         ground temperature     units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   TH2OSFC                    	long_name         surface water temperature      units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   THBOT                      	long_name         Patmospheric air potential temperature (downscaled to columns in glacier regions)   units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   TKE1                   	long_name         (top lake level eddy thermal conductivity   units         W/(mK)     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   TLAI                   	long_name         total projected leaf area index    units         m^2/m^2    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   TLAKE             	            	long_name         lake temperature   units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       (     [�   TOTCOLC                    	long_name         >total column carbon, incl veg and cpool but excl product pools     units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            [�   	TOTCOLCH4                      	long_name         \total belowground CH4 (0 for non-lake special landunits in the absence of dynamic landunits)   units         gC/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \    TOTCOLN                    	long_name         -total column-level N, excluding product pools      units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \   
TOTECOSYSC                     	long_name         Atotal ecosystem carbon, incl veg but excl cpool and product pools      units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \   
TOTECOSYSN                     	long_name         *total ecosystem N, excluding product pools     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \   TOTLITC                    	long_name         total litter carbon    units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \   
TOTLITC_1m                     	long_name         $total litter carbon to 1 meter depth   units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \   TOTLITN                    	long_name         total litter N     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \   
TOTLITN_1m                     	long_name         total litter N to 1 meter      units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \   TOTPFTC                    	long_name         )total patch-level carbon, including cpool      units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \    TOTPFTN                    	long_name         total patch-level nitrogen     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \$   
TOTSOILICE                     	long_name         /vertically summed soil cie (veg landunits only)    units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            \(   
TOTSOILLIQ                     	long_name         8vertically summed soil liquid water (veg landunits only)   units         kg/m2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            \,   TOTSOMC                    	long_name          total soil organic matter carbon   units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \0   
TOTSOMC_1m                     	long_name         1total soil organic matter carbon to 1 meter depth      units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \4   TOTSOMN                    	long_name         total soil organic matter N    units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \8   
TOTSOMN_1m                     	long_name         &total soil organic matter N to 1 meter     units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \<   TOTVEGC                    	long_name         (total vegetation carbon, excluding cpool   units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \@   TOTVEGN                    	long_name         total vegetation nitrogen      units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \D   TOT_WOODPRODC                      	long_name         total wood product C   units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \H   TOT_WOODPRODC_LOSS                     	long_name         "total loss from wood product pools     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \L   TOT_WOODPRODN                      	long_name         total wood product N   units         gN/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \P   TOT_WOODPRODN_LOSS                     	long_name         "total loss from wood product pools     units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \T   TPU25T                     	long_name         canopy profile of tpu      units         	umol/m2/s      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \X   TREFMNAV                   	long_name         (daily minimum of average 2-m temperature   units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \\   TREFMXAV                   	long_name         (daily maximum of average 2-m temperature   units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \`   TSA                    	long_name         2m air temperature     units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \d   TSAI                   	long_name         total projected stem area index    units         m^2/m^2    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \h   TSKIN                      	long_name         skin temperature   units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \l   TSL                    	long_name         Rtemperature of near-surface soil layer (natural vegetated and crop landunits only)     units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            \p   TSOI                      	long_name         <soil temperature (natural vegetated and crop landunits only)   units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg       d     \t   	TSOI_10CM                      	long_name         $soil temperature in top 10cm of soil   units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            \�   TSOI_ICE                      	long_name         %soil temperature (ice landunits only)      units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         ice       d     \�   TV                     	long_name         vegetation temperature     units         K      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]@   TWS                    	long_name         total water storage    units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]D   T_SCALAR                      	long_name         'temperature inhibition of decomposition    units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     ]H   U10                    	long_name         	10-m wind      units         m/s    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   U10_DUST                   	long_name         10-m wind for dust model   units         m/s    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   URBAN_AC                   	long_name         urban air conditioning flux    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   
URBAN_HEAT                     	long_name         urban heating flux     units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   VCMX25T                    	long_name         canopy profile of vcmax25      units         	umol/m2/s      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   VEGWP             
            	long_name         Fvegetation water matric potential for sun/sha canopy,xyl,root segments     units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   VEGWPLN           
            	long_name         Kvegetation water matric potential for sun/sha canopy,xyl,root at local noon    units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   VEGWPPD           
            	long_name         Epredawn vegetation water matric potential for sun/sha canopy,xyl,root      units         mm     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   VOLR                   	long_name         !river channel total water storage      units         m3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   VOLRMCH                    	long_name         (river channel main channel water storage   units         m3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   VPD_CAN                    	long_name         canopy vapor pressure deficit      units         kPa    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   Vcmx25Z                    	long_name         1canopy profile of vcmax25 predicted by LUNA model      units         	umol/m2/s      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   	WASTEHEAT                      	long_name         Csensible heat flux from heating/cooling sources of urban waste heat    units         W/m^2      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   WBT                    	long_name         2 m Stull Wet Bulb     units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   WBT_R                      	long_name         Rural 2 m Stull Wet Bulb   units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   WBT_U                      	long_name         Urban 2 m Stull Wet Bulb   units         C      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   WIND                   	long_name         #atmospheric wind velocity magnitude    units         m/s    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ]�   WOODC                      	long_name         wood C     units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ^    WOODC_ALLOC                    	long_name         wood C eallocation     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ^   
WOODC_LOSS                     	long_name         wood C loss    units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ^   WOOD_HARVESTC                      	long_name         &wood harvest carbon (to product pools)     units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ^   WOOD_HARVESTN                      	long_name         !wood harvest N (to product pools)      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ^   WTGQ                   	long_name         surface tracer conductance     units         m/s    
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ^   W_SCALAR                      	long_name         .Moisture (dryness) inhibition of decomposition     units         unitless   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     ^   XSMRPOOL                   	long_name         temporary photosynthate C pool     units         gC/m^2     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ^h   XSMRPOOL_RECOVER                   	long_name         0C flux assigned to recovery of negative xsmrpool   units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ^l   ZBOT                   	long_name         atmospheric reference height   units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ^p   ZWT                    	long_name         =water table depth (natural vegetated and crop landunits only)      units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            ^t   ZWT_CH4_UNSAT                      	long_name         Fdepth of water table for methane production used in non-inundated area     units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ^x   	ZWT_PERCH                      	long_name         Eperched water table depth (natural vegetated and crop landunits only)      units         m      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         veg            ^|   mcdate_2                	long_name         current date (YYYYMMDD)            ^�   mcsec_2                 	long_name         current seconds of current date    units         s              ^�   mdcur_2                 	long_name         current day (from base day)            ^�   mscur_2                 	long_name         current seconds of current day             ^�   nstep_2                 	long_name         	time step              ^�   lon_2                  	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��           K�   lat_2                  	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��           K�   area_2                 	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��           K�   
landfrac_2                 	long_name         land fraction      
_FillValue        {@��   missing_value         {@��           K�   
landmask_2                 	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����           K�   	pftmask_2                  	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����           K�   
nbedrock_2                 	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����           K�   CWDC_TO_LITR2C_vr                         	long_name         .decomp. of coarse woody debris C to litter 2 C     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     ^�   CWDC_TO_LITR3C_vr                         	long_name         .decomp. of coarse woody debris C to litter 3 C     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     ^�   	CWDC_vr_2                         	long_name         CWD C (vertically resolved)    units         gC/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       P     _\   CWDN_TO_LITR2N_vr                         	long_name         .decomp. of coarse woody debris N to litter 2 N     units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     _�   CWDN_TO_LITR3N_vr                         	long_name         .decomp. of coarse woody debris N to litter 3 N     units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     `   GROSS_NMIN_vr                         	long_name         gross rate of N mineralization     units         gN/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     `t   LITR1C_TO_SOIL1C_vr                       	long_name         !decomp. of litter 1 C to soil 1 C      units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     `�   LITR1N_TO_SOIL1N_vr                       	long_name         !decomp. of litter 1 N to soil 1 N      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     a<   LITR1_HR_vr                       	long_name         Het. Resp. from litter 1   units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     a�   LITR2C_TO_SOIL1C_vr                       	long_name         !decomp. of litter 2 C to soil 1 C      units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     b   LITR2N_TO_SOIL1N_vr                       	long_name         !decomp. of litter 2 N to soil 1 N      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     bh   LITR2_HR_vr                       	long_name         Het. Resp. from litter 2   units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     b�   LITR3C_TO_SOIL2C_vr                       	long_name         !decomp. of litter 3 C to soil 2 C      units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     c0   LITR3N_TO_SOIL2N_vr                       	long_name         !decomp. of litter 3 N to soil 2 N      units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     c�   LITR3_HR_vr                       	long_name         Het. Resp. from litter 3   units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     c�   NET_NMIN_vr                       	long_name         net rate of N mineralization   units         gN/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     d\   SMINN_TO_PLANT_vr                         	long_name         plant uptake of soil mineral N     units         gN/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     d�   SMINN_TO_SOIL1N_L1_vr                         	long_name         +mineral N flux for decomp. of LITR1to SOIL1    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     e$   SMINN_TO_SOIL1N_L2_vr                         	long_name         +mineral N flux for decomp. of LITR2to SOIL1    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     e�   SMINN_TO_SOIL1N_S2_vr                         	long_name         +mineral N flux for decomp. of SOIL2to SOIL1    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     e�   SMINN_TO_SOIL1N_S3_vr                         	long_name         +mineral N flux for decomp. of SOIL3to SOIL1    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     fP   SMINN_TO_SOIL2N_L3_vr                         	long_name         +mineral N flux for decomp. of LITR3to SOIL2    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     f�   SMINN_TO_SOIL2N_S1_vr                         	long_name         +mineral N flux for decomp. of SOIL1to SOIL2    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     g   SMINN_TO_SOIL3N_S1_vr                         	long_name         +mineral N flux for decomp. of SOIL1to SOIL3    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     g|   SMINN_TO_SOIL3N_S2_vr                         	long_name         +mineral N flux for decomp. of SOIL2to SOIL3    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     g�   SMIN_NO3_LEACHED_vr                       	long_name         soil NO3 pool loss to leaching     units         gN/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     hD   SMIN_NO3_RUNOFF_vr                        	long_name         soil NO3 pool loss to runoff   units         gN/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     h�   SOIL1C_TO_SOIL2C_vr                       	long_name         decomp. of soil 1 C to soil 2 C    units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     i   SOIL1C_TO_SOIL3C_vr                       	long_name         decomp. of soil 1 C to soil 3 C    units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     ip   SOIL1N_TO_SOIL2N_vr                       	long_name         decomp. of soil 1 N to soil 2 N    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     i�   SOIL1N_TO_SOIL3N_vr                       	long_name         decomp. of soil 1 N to soil 3 N    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     j8   SOIL1_HR_S2_vr                        	long_name         Het. Resp. from soil 1     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     j�   SOIL1_HR_S3_vr                        	long_name         Het. Resp. from soil 1     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     k    SOIL2C_TO_SOIL1C_vr                       	long_name         decomp. of soil 2 C to soil 1 C    units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     kd   SOIL2C_TO_SOIL3C_vr                       	long_name         decomp. of soil 2 C to soil 3 C    units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     k�   SOIL2N_TO_SOIL1N_vr                       	long_name         decomp. of soil 2 N to soil 1 N    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     l,   SOIL2N_TO_SOIL3N_vr                       	long_name         decomp. of soil 2 N to soil 3 N    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     l�   SOIL2_HR_S1_vr                        	long_name         Het. Resp. from soil 2     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     l�   SOIL2_HR_S3_vr                        	long_name         Het. Resp. from soil 2     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     mX   SOIL3C_TO_SOIL1C_vr                       	long_name         decomp. of soil 3 C to soil 1 C    units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     m�   SOIL3N_TO_SOIL1N_vr                       	long_name         decomp. of soil 3 N to soil 1 N    units         gN/m^3     
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     n    SOIL3_HR_vr                       	long_name         Het. Resp. from soil 3     units         gC/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     n�   SUPPLEMENT_TO_SMINN_vr                        	long_name         supplemental N supply      units         gN/m^3/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown       d     n�   mcdate_3                	long_name         current date (YYYYMMDD)            oL   mcsec_3                 	long_name         current seconds of current date    units         s              oP   mdcur_3                 	long_name         current day (from base day)            oT   mscur_3                 	long_name         current seconds of current day             oX   nstep_3                 	long_name         	time step              o\   lon_3                  	long_name         coordinate longitude   units         degrees_east   
_FillValue        {@��   missing_value         {@��           K�   lat_3                  	long_name         coordinate latitude    units         degrees_north      
_FillValue        {@��   missing_value         {@��           K�   area_3                 	long_name         grid cell areas    units         km^2   
_FillValue        {@��   missing_value         {@��           K�   
landfrac_3                 	long_name         land fraction      
_FillValue        {@��   missing_value         {@��           K�   
landmask_3                 	long_name         &land/ocean mask (0.=ocean and 1.=land)     
_FillValue        ����   missing_value         ����           K�   	pftmask_3                  	long_name         (pft real/fake mask (0.=fake and 1.=real)   
_FillValue        ����   missing_value         ����           K�   
nbedrock_3                 	long_name         !index of shallowest bedrock layer      
_FillValue        ����   missing_value         ����           K�   LEAFC_TO_LITTER                    	long_name         leaf C litterfall      units         gC/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            o`   NFIX_TO_SMINN                      	long_name         1symbiotic/asymbiotic N fixation to soil mineral N      units         gN/m^2/s   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            od   OFFSET_COUNTER                     	long_name         offset days counter    units         days   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            oh   
OFFSET_FDD                     	long_name         #offset freezing degree days counter    units         C degree-days      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ol   OFFSET_FLAG                    	long_name         offset flag    units         none   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            op   ONSET_COUNTER                      	long_name         onset days counter     units         days   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ot   	ONSET_FDD                      	long_name         "onset freezing degree days counter     units         C degree-days      
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            ox   
ONSET_FLAG                     	long_name         
onset flag     units         none   
_FillValue        {@��   missing_value         {@��   cell_methods      
time: mean     landunit_mask         unknown            o|   FROOTC_TO_LITTER                   	long_name         fine root C litterfall     units         gC/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            o�   FROOTN_TO_LITTER                   	long_name         fine root N litterfall     units         gN/m^2/s   cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown            o�   
FROOT_PROF                        	long_name         1profile for litter C and N inputs from fine roots      units         1/m    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       d     o�   	LEAF_PROF                         	long_name         -profile for litter C and N inputs from leaves      units         1/m    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       d     o�   	NDEP_PROF                         	long_name         %profile for atmospheric N  deposition      units         1/m    cell_methods      
time: mean     
_FillValue        {@��   missing_value         {@��   landunit_mask         unknown       d     pP   date_written                             p�   levgrnd                	long_name         coordinate ground levels   units         m         d     K�   time_bounds                    	long_name         history time interval endpoints            p�   time_written                             p�<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�<#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�f=L��?��@ff@�33A��AI��A���A���B	L�B3�A0I{Bq��?
?�           A0I{Bq��?
?�           A0I{Bq��?
?�           <#�
=#�
=�Q�>#�
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�fF� @߹@    @��     '��           ��4�6�b:B'�fB'�f7@��D"g�    R+�K6���?�  ?��4ě�1�A2��
/Ƽ�.0'@/j�m0�$��� y�h8��`1�b�6b�<��:=ũ=I��80�-2?i�:(�<(Pj<��w=H#�=��=��_=�Ն=ص|=��=�BT>��>W�8��1ح%5�;>$�8��n1���0� 8l҇@``9@n\?@d4�@d
@c�K@c��@cQ�@b�@b�H@b!@a��@a�{@��{@��{@��                    E�}5�G�	�F��FR�NE�RE�D1�eCC�'BA��                                                A^.CL9B)RA�v�A?!r@�(m?�SQ>�U=�R�                                                                    EEbT@ʆ�F#kIA���1�ZNFa�.B4��<�D/,�k�                                                    {@�ο�%)��%)@o��7y 1%��i�l�fg/=�C$�ӈ>�&R���/��\3�r    >�y�=���1�1[>q]`�_4?� �l<#�
    Aӿ.Aӿ.C�_�C�_�C�#�?�  >�u3��C�e�6D�c6�V�A6��AF�>�y�AT<a@G]�A���@a�@w�A7R@� eA�W�=�޽�	�    �=��    �=���<S�<'0�?}p�?}p�@_�^?<�:@}��?6�?5��@x�7?�T�F�H�1�ZN/e,�|"2*z�6P�E    4a�dF9<F�F��F��@)�@��.C��@rr�?�I�?O��?+	??�>��>���>��K>�u8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�8�    ���@���@{@��6bU�        5S�7.�7_�B6��5x��4ְL                                                A�l�����{@��CV��A �B7gV@l��=5�{@��{@��    C���BP*�6$B    6Tm7� @�C�3��6�^�    >��C�j7��GBr�Dl��Co��B|�0A�-AL{�@�5#@�i>�                                                 ?,�sAr�@��?�P-?0>`>��>�#=M�J<?��                                                                    DA��F�C�E��D�,�C�j�C=;xB��YA��@�#&                                                @���B�ѹB:�A$p�@���@�o?���>���=�E|                                                                    C�$�E�z�EwYC�\yC<!�B��LA��Aq@W�                                                @+l�BWfpA��@�� @Ք?��i>�m�>/�=4��                                                                    5��6?��A$%\>Z
�B�3?5��?�A7)��(J�+�3i+2��)��(�	)�Hy+�۟+��(��7@��                        �D�`3��1���            7D�`��c�D�`3�P�                        >lG��h�                                                                2�Dz            2�Dz{@��,�=�>L��>L��>L��>��:?~s�?7�b?��?5��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A���G���>*�MA��j                                                B�                                          B��                    A�                @�p�    $1� ~Y    4�2Q�2*z�=bQF6=!�=��t4�+�;%�w;(`                �-
        ,��1        {@��4!4?f�6�+��`�{@��{@��{@��178�             6�6    6�6{@��    6�6    3'C3a�{@��3V�{@��8!cL7B�7�o�    3�/{@�ε3_�3 �3"��3"��<�lB    B�9�6��Fu�ND�E�>���=�|zA�            ?�=    :YJ�    >���>72�@�C�@�N�@@�]?���>�'> D
����
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/?�!>���>72�@��@�M}@@�H?���>�'> C�����
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/;�5    (�H4��3���=�q�9:6�ـ2�p1z�5�+�                                                �!���C����������ef@�A����û̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6NS4	��@�7u�5`f�@4�@�@ '�@�Y�=�*C7��<56�)G6w8Ay�?ArC��9e?� 8Ay�J#�04=}�8T�!CF)D��FDƀ�D8W�C�7�C��BL�A��@��                                                A�)CY�FCF��B�W�B&7�A���@��@�?�                                                                    E��iG� GK�F�miF9-E�(bDɋ�C���C N�                                                C�<EOz.EA�E~D��iC�1C�BB-�A:�                                                                    E�F�I�F��OF{-	F]�F.�E�6�E���EG�Z                                                D)��D�S�D�[D��dD��
D~n
D0��C�{SC�>                                                                        6 �fGX�oGM�iGcWF͸>Fz��F��E��EP�                                                A��A��A��@��f                                                                >�V�@d~A!jB �B6��Bk��B�S�B���<�<�<�<�<�<�<�<�<�<�<�<�E�`�E�jEc�ER�D���DW �C�$'C�`�                                                {@��{@��{@��{@��{@��A�/�B�e#�M�f�6�?sC��@�S    @��@��{@�ξ�ܾ��C��-{@��C��]C��bC��-{@��@o��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G)>�@bD�B$Gp�D�A�D�@�D�@�A��A��Fb$B=Y�B���C�*0FK��FK��D��D��Fc��B=T�                ?��C��;C��C��F? �C�z�C�|�C�|�C���C��yC���C���C���C�+�C�l�C���C���C��C�[�C���C��^C��C�Q�C��C��JC���C��C��~C�[�C�2C��C�?C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�w�DrZ>�c�>���>��L>��X>���>�q>�>�X�>���>���>�Ǒ>�7#>��U>�\�>�>�~>�CV>�%�>�h>��d@�?偛        @m��ƾ�ƾ��ƾ��K���?���9��5��]��ƾ4ƾ4ƾ4�F}        =NNA�6�    ��p{��p{{@��@e5�FUu�6��7[�        0��        =.G5?A�i?HS�?N��?XsJ?i�                                                ��V�5UܖA�  ?k�    <#�
'��           ��        4�v6M:�6h�5kP�4A��3�d�                                                                            3�R5��l5���4���3t��2�!                                                                    G�	�F��FR�NE�RE�D1�eCC�'BA��                                                        0	�e1�l�1��0��'/�@�/��                                                                            /.N�1�1�L0�/ ^�.>�                                                                            3ID-5,\	5a��4{p:3u��2�7                                                                            3;�5)g�5O8N4�N�3��2��                                                                            1�3�j3)��2q�&1��-0�                                                                            3e*5O�5}D�4���3�Y�2�К                                                                            4w��65��6a�q5��+4p	3��                                                                            1��3p�3��2�Y�1퇜172V                                                                            4w��65��6a�q5��+4p	3��                                                                            3��5�1�5�%t4�w�3� �3o�                                                                            1�2�"3)b?2a�k1h��0���                                                                            3��5�1�5�%t4�w�3� �3o�                                                                            2U m4��4��3�D2�r.2~��                                                                    9V��8�I�;� �;9�H:�X�:��9L٧8�I�                                                                            1qY�3Pݣ3t��2���1�!0��{                                                                            2��j4���4���3ɉ�2��"2G�                                                                            �J�d�=iK�mX��r�#�U���~�                                                                            �[s��mD�Y�����ӱU��O�                                                                            1ӕ�3��,3�?�2Ȝ�1�U1=;                                                                            �T���u�oش/���,��9g                                                                            ���q�q}���~��������,nU                                                                            ��?���$q�èt���̰�V��2`                                                                                                                                                                        *�K)�(m                                                                                                    4h�?6;@G6g��5�L4j�a3϶�                                                                            0�A�2��2�t�2W1��0kId                                                                            3c�~59�5o�4��%3��j2���                                                                            /�v�1���2 �1��0�/���                                                                            4_3�68��6w��5�/4�$3�̘                                                                            0᫞2�'�3U 2�21�$0�X�                                                                            4$��6��6@��5D��4-�>3���                                                                            2<'�4/��4\�3`��2Fal1���                                                                            3	�4�h�5�4�3D�2b8�                                                                            1c3�31��25��1 N�0�E                                                                            4I8�6;��6kf�5p�!4T(03��                                                                            2e��4V�4���3�z\2rw1ÅD                                                                            11�'3��40��3�>p3,�b3(�                                                                            0��2d��3��2�m2��2�p                                                                            1YH�3��A4W��3���3SRx3MkB                                                                                                                                                                        '��           ��6G=`                    >�    6�r�3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @߹@    @��     16:55:49        F�@ @��     @��     '�]            � 4o��6��B'�fB'�f7R��D"g�     	\+��i6W�'?�  ?��5'��1f�w2��`/Э�.�p=        &/� ""\n8ˬx2�jl8�q>B�=�_b=I�j87��2G(�:1V�<(�<��=H,H=� �=��T=���=ضm=���=�C>�>W�8�ө2-��6�.�<ʒ8p�D11W�0�ӛ8e��@\R1@n"&@d3@d�@c�?@c��@cQ�@b�@b�H@b!@a�@a�	{@��{@��{@��                    E��6
!�G���F�?�FS�E�Y�E�wD1�CD"WBB�                                                A�uCTB�1�AބyA?i�@�$�?���>�֮=˶�                                                                    EED�@�j*F#SFA��2��FaGB3��<�D/-�e�                                                    {@��@��@��@ll�7���&nT�%2$z�o�q/6�%��>쏨 	\>�`/�@|4�C>    >�y�?@K	1�[?t�转�:�3�<#�
    B*��B*��C��>C��>C�2??�  ?8
�3qC�6�6�r:A4l�A�.�>�y�B@�A6��BO�y@��-@�|�BlA4��B2c8�u�W�]�_!6F��u�W    �u�W�>N�@tM�?}p�?}p�A��@ոA-F?˗z@�SA @7^F�WU2��0aS-n�82�A.7�W    4���F�q�FӢ�G9qGC�?u�-@���CY�@�k�?�;	?O��?&r�>��8>��q>�C_>��>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��/    ��i��i{@��6��6%̊5�lC7/�8�;7^��6}D=5v!�4��                                                A�c������{@��C��Ay�>BGތ@e��=ڱH{@��{@��    C��BP't5�"�    6Q�27�@�.�3�@6���    >]�C��G7��IB%�Dx��C��B�8+BBAO��@�\s@�?J�                                                ?=\�A~_�@�#`?�[?8�>��s>��=e_�<R�E                                                                    DC�{F��JEΦjD��5C�^hC=��B��A�>�@���                                                @�.VB�r�B?�dA/��@�y@ݰ?��>��#=���                                                                    C��E�j�E�ZC�E�C?�B�t�A�D�A@�(                                                @/�BY3lA�0�@�|�@��?��?�>6`�=:E                                                                    6-d6���A��>S��B�g?0��?��*�U(�d�,�+�n�*�U(�ۻ*�Rw,G�+�^")��7R��                        �"�J3q1�             7"�J��\m�"�J4	�                        >�/�����                                                                2�e            2�e{@��,�bJ>L��>L��>L��?�  ?F?6�?	�?4�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Bl��G�;�>+U�A��9                                                B�                                          B��                    A�                @�p�    %�-�!ȕ�    4o��3uh2�A.>�66�i�>���5�K�;P}�;S*                5_�g    4�$�2^��        {@��4�56l6����)��{@��{@��{@��1�6~7�             6�Е    6�Е{@��    6�Е    4�Z�7-H{@��7?�{@��7���7��7^[�    4��M{@��4N^�4�>55<�    B���6��FFO�yD�?�� >Ď`A�c,            ?�V�    :�$�    ?+��>�k�@�ܑ@�l@|��?�Y�>���>D�Ë���
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/?�J?+�/>�k`@�O:@��@|��?�Y�>���>D�K����
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/;��    *��R7A�68�'=�U�9��h6���2��1��%5�=                                                �!���C���ĉrE�w�i�S9�%���m9̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6�E 4�l�?O��80�6���@�"Q@�GF@�V�A'u�> �86�6
O�Gfչ7��?@�CT�B@�(�?B�V7��J���587�RCFC�DٜD�7uD8��C���C�BME�A�@�?u                                                A�C�CY�CF7uB���B&��A��@�E�@??u                                                                    E���G��G.�F�}F9-]E��DɆQC���C Nd                                                C�
EOK�EA��E�wD���C���C�;B-MA:��                                                                    E��F�DEF���F{(vF]�F.�E�7yE���EG�R                                                D)��D�LD�|D��D��D~nD0�C�|�C�>�                                                                        6V�GXr�GM��Gk�Fͺ�Fz��F��E�	JEP�                                                A�E%A���A��}                                                                    ?�z@�3A���B �CB4�ZBh�7B�B�o<�<�<�<�<�<�<�<�<�<�<�<�E�E�E�Q Ec�/EUD���DW!!C�&�C�al                                                {@��{@��{@��{@��{@��A�-9B�A�%�-��f�\7'O'B߼p?��    @��@��{@�ξ��k���kC�|v{@��C�6`C��"C�|v{@��@ll�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G �>���D�A�Ge#D�AwD�	�D�	�A�A�Fa0B;��Bv�C�%oFK�jFK�jD���D���FczB;��                ?�Y\C�_�C�8�C�{�?   C�^CC��VC��VC���C��3C��C��C��0C�#C�E�C�rNC���C�ɸC��*C�/ C�g�C��NC��RC��C�L�C�n�C���C��C�a�C�4�C� 
C�sC�� {@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�aD)�>�S�>�c�>�w>��G>�#>��o>���>�W5>���>��\>�9�>���>��>��D>�X�>���>�L�>�('>�=>�ft@; �@'J        @����2��������r0������ͫ�ܳ�Ī4�ƾ��ƾ��ƾ���T:n        =�=�AΏ,    �1o>�1o>{@��@�:�FUQS6�[�7B�        0�	0<���<�y�?�r?A�-?Ex�?KlL?T�	?c��                                                ��$�6�ƕA�  ?k�    <#�
'�]            � 5
ȵ4�s:6U��7)6f+65fb`4<13�>�                                                                    4/N^3�N/5��)6G�i5�^�4���3m�72�\                                                                    G���F�?�FS�E�Y�E�wD1�CD"WBB�                                                0��0�1�Wd2��1�B0�͸/ŀ�/�                                                                    /�]//�i1�u1љ1�0��.�z<.5�3                                                                    3#K�3+vr5#,~6��5_�l4v@�3n�q2��                                                                    4�3��58��6��5QM4�'3�؇2إ�                                                                    2�1_�#3K�3�ח3*��2{,�1�s�0��                                                                    4��+3��F5a�6")5z�4�03�]�3eX                                                                    5��5.�6OS7�6`��5���4n�L3�4�                                                                    2�EF2��3n�J4<63���2縧1�14                                                                    5��5.�6OS7�6`��5���4n�L3�4�                                                                    5��4c�5�z6^{�5��u4�Vx3˼�3��                                                                    2
�1���2�3�.�3)EC2b�1k�0��l                                                                    5��4c�5�z6^{�5��u4�Vx3˼�3��                                                                    �:-x�@f)4ME5E0
4��63�ր2��2n2�                                                                    9�9E3�;�+�;a7D;��::��9t��8ߏ�                                                                    2Í1�3n`	4#��3w5�2���1�= 0�a,                                                                    4":�3���4�5
5k4�4��3�S�2��20�                                                                    ��={��5�$AZ����kM*�m'�O[[���                                                                    �����M�1M��=9�W�ر��@�N�'�E��                                                                    3!O�2��p3��94i�63�p�2Ʒ�1��1_�                                                                    � B³�v�������4r�,~9�'D����b                                                                    ���ܯ�e��kP�:k˱����⯽x��$ߕ                                                                    �IWA�W���g��pQ����ޱÀA��׋                                                                                                                                                                        ,��,(S                                                                                            4P8D4Yf{6<��7��6eĴ5{{�4d=�3Ɯ&                                                                    0҅�0���2��U3�#F2�~�2HJ1D�0`�_                                                                    3K�!3T�u58�j6�J5m��4���3zC�2��*                                                                    /�;�/�S^1���2�.1��1�0�%/v��                                                                    4G��4P�n65�7{�6u��5�3)4�$�3��.                                                                    0��n0���2�h3�3H|2��17S0�1�                                                                    3�1�3�o�6J6�V6>��5@r@4(C�3��-                                                                    1�)1���4T�5-�4Z7�3[�I2@M�1��|                                                                    2�(Q2���4�j�5�)�5K�4�3��2X�                                                                    0�	�0�0�2�0�3�x�30V�21�j1ep0v�B                                                                    3�<�4ҋ6"�7�6i_k5k6N4M�?3�oS                                                                    2
k�2^V4:.5%84�[3�h-2k	l1��_                                                                    0��0�Z*3�"4Y�4/'�3�,
3'�3 �H                                                                    /h��/� u1膮30%3��2��2o�2��                                                                    0��0ǧ3/�4�#�4V�3���3L��3D>�                                                                                                                                                                        '�]            � 6Dߛ                            6��<3�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��     @��     16:55:52        F�~ @��     @���    '��      ?     ��4�u;6u�B'�fB'�f7M*�D"g�    !�;�+�Ɖ6#j?�  ?��5�־2S�2�1+1��1ʪ0�'�1w(M/#�'9�d�8�Z�> �>F9=�&�=W��:4�[2��:8;<(��<���=H4u=�#�=��U=��9=طr=��=�C�>�n>X?9��*8 �8=x�b=�^9V.z2� 0��&8`ҁ@XoT@m�j@d1@d�@c�5@c��@cQ�@b�@b�J@b!@a�@a�A�=�;�?�{@��                    E���6lRG�n�F�ZWFT��E�V4E��D2�CD��BBw7                                                Ax�C�TB��A�7MA?g�@�%�?��>�\=�#c                                                                    EE�@�D�F#2A�~e3%�F`�B2l�<�D/.�ݧ                                                    {@��@���@���@h$=7�:�$�=��"�h�Z @�.o�&2�>��!�;�?���/�25'%�    >�y�@m%1�b�@,� @�T����<$��    Bu�Bu�C�IC�ICv�?N�W?Ɔ�3�tC�$5�a�6��A1ڧB���>�y�B��mBGBC �AIįA�/2B�HiA�A�Bt��@�,��.K"Y��@�,    @�,A<-/A6G�?}m�?}m�A��~@ؚ�A�ӕ@"��@��fA��|@�UQF�Xx3%�0�-�>�3%f�7�Ew3� �5��FǠ�G �uGO��GpT	>���A-FC��@���?���?Mh]>���>ƀg>�,U>��1>�'	>�^X8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�[8    ������{@��7J 7��i7���8QW�8r�7VB�6{X'5t*f4��                                                A�U���4"��4"{@��C���A�c�Bo��@ZD�>^\{@��{@��    C��$BO�X5��    6N�F7.�@���3��6�    =6wGC��37���B/�D��C�<B�a�B ��AR�@���@��?&&                                                ?H}A��Y@�J�?�2?0�/>�Ԡ>�=yZ�<b�                                                                    DD��F�w�E��D��YC�)C>�?B���A���@��                                                @�L�B��BD�A6C�@���@�?�!D>���=��'                                                                    C�8E�32EU�C�i�C?HSB�9CA�B�A"��@#%C                                                @2�lBZ�%A� �@�.>@��?�*�?�>=\x=?Im                                                                    6���7��A��>L�%B Q,?*�n?�B#*��)��,,�+�6.*��)T�e+��,E�,��)3��7K<�2Q�s0Z��#Ժ�            �/�w3M]�1�a�0Z��0Z��#Ժ�6/�w����/�w4y2    %�<�2D�2D�&&+�    ><Y6��5��x/�z�.�7p(K�                        )�*/�k/�c()}��    /��2�l
            3�b9Ԓ}-��>L��>L��?:[d?�  ?s�?<%�?�9?::)?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B�Y%G�G�>,D�A���                                                B�                                          B��                    A�                @�p�    '���#��|2h��5>*3A�13%h�>��@7jϤ?(�6�v�;G�_;H��                6H��    4��12�8    6�x�{@��2��5�`73֎�a�8{@��{@��{@��2��77	T�            7L��    7L��{@��    7L��    5}��7���{@��81k{@��7_z6d�6�Lq    5�o{@��5O5��5{WU5{WU<��    B��6��>F3�rDǔ@z}�?���B���            @��2Q�s:�h�2Q�s=��=�P3@�@���@��?�t9>�kF>e�ً���
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/@~b=�9=�Nk@�2@��@���?�t8>�kF>e�����
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/9b�    *��5�+�6c�@:���:[�7�k^3 V2�}6</�                                                �!���OM�Ġi�ĭ(�ĝN�Ĉ�.�`���'��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7�t4�ډ>��9 �K7�A�@���A��A��i>.\�8P��6|=G���7�?N�SC�:�A�?Q�7�J�c6�y�7:b�CF�CD�G\D�ZbD9�pC�b�C@BM��A��@���                                                A��CCZG\CFZbB��pB'b�A�@@���@�?��                                                                    E��2Gg�G�F��qF9 7E�_DɂEC��C O>                                                C���EO"{EA�vE��D��C��sC�JB-(A:��                                                                    E��F�?;F��#F{%�F]#F.��E�8VE��EG�l                                                D)�'D�D�D�D��D��_D~n9D0�'C�~�C�?�                                                                    /��6�٦GXY�GM��Gr�FͶ�Fz�hF�nE�DEP
3                                                A��OA���?�9�                                                                    ?�9$A�A���A� �B.��B`j=B�B��r<�<�<�<�<�<�<�<�<�<�<�<�E�4�E�=�Ec�<ER�D��jDW"oC�)�C�bv                                                {@��{@��{@��{@��{@��A�,�B�>�'\v��e�)7]�?B��?�o    @�2 @�2 {@�ξbۊ�bۊC���{@��C�o�C��C���{@��@h$^{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G �5>���D�@:GSMD�@D�I�D�I�A��A��F`n�B9�7B��C���FK��FK��D��VD��VFb_LB9�
                @�C�kjC���C��_?   C��tC��3C��3C��3C��&C��yC��hC��C�C�4�C�YTC�{sC���C��C��bC�(�C�_�C��C��mC��C�36C�R�C�mBC�a�C�7C� �C��C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��5D ��>�w>�w>�=>�-�>�f�>�غ>���>�u�>�W>>��>��>�K�>��P>�]j>�m&>���>���>�S>��>��+@)�h?��m        AC]���A6����������o� �� ��� ~v�2�R������������ĉ�         =�F&A�k�    �&2��&2�{@��@z��FU�6^��7!Y        0��>�]�>�W_?6Ù?9`�?<�?B1�?I��?U�                                                ��07B�A�  ?k�    <#�
'��      ?     ��7TE
6�ʼ7{<?7޾6\�5b��47�(3�U                                                                    6��5�>6���6?��5��z4�!~3h
�2��                                                                    G�n�F�ZWFT��E�V4E��D2�CD��BBw7                                                2��#2CY3�Q2�lt1��0��O/�˝/
�                                                                    2�#1v[c2&�1�`�1w�0:�.��./!Y                                                                    5yr�5o�J6@�6 �5V�(4rz*3iM�2��b                                                                    5�o�5���6gtk5�lA5KFY4�^�3�YG2�Se                                                                    3kA�3="B4$�3��B3%��2��1���0�                                                                    5��u5�ț6�q�6^�5xrl4�:�3��W3z                                                                    7��6�7w{}7��6X��5�4n�3�A�                                                                    3�e3�4�͕44��3�%f2���1��11�                                                                    7��6�7w{}7��6X��5�4n�3�A�                                                                    6G�65��6���6U�U5�"�4��@3��F3�                                                                    3M5�3R4�4	�3�T}3#��2c;�1m��0�y@                                                                    6G�65��6���6U�U5�"�4��@3��F3�                                                                    �� ��L�+5%,�5?��4� �3��2���2`�R                                                                    7�xd8  �;�UX;��;0+�:^S9���9�1                                                                    4�;3̈�4�i�4�3p��2���1� 0��                                                                    5oP5T��5�H5b&�4�I_3��|2��N2 �!                                                                    �:�6�h�@����aȴ�h��IꞱ�B�                                                                    �KԱ�S�OҼ���
�O0�����I�>�                                                                    4n,d4S�4�;04`��3���2Řc1��217�                                                                    �D��;jQ�	���eٵ��*��#�7��$                                                                    ��=���|Ⲋ�f�3�N��B!��F���R���                                                                    ��Y����鳞��f�䲺 [��ٰ�sx�v�                                                                                                                                                                        ,L�,�\�                                                                                            6�'16� ;7^o�7`�6\��5w�4_>�3�:o                                                                    3 �3��3��@3��o2���2T�0��&0X��                                                                    5�� 5��6Y�@6
a�5d`4��3t�P2Ѯ                                                                    2�`2��2�P 2��	1��1�V0
��/m��                                                                    6���6��*7U^�7	b�6lf5�8�4�)�3�!�                                                                    3[3�F3׻3�ˤ2�4�2��1�$0�5                                                                    6�-6c�7C�6�;�677�5=�4#٪3�l�                                                                    4,��4)��52��5��4Qdf3X�2;A�1��p                                                                    4��<4�Ң5���5��M54��3g|2O�(                                                                    3h�3	
�4Pp3�ڡ3)4�2.��1Q�0m4�                                                                    68�65]�7>��7
�E6_��5g=4HB�3���                                                                    4R�+4OF�5ZFv5��4�}3��2d��1�c                                                                    2۔?2�Qz4(��4QTG4(!?3�3#2�34�                                                                    1�p
1��3G
3)'�3��2�9�2�1�8�                                                                    3/�3� 4N�4��4M}�3��3Gw3<y.                                                                                                                                                                        '��      ?     ��6A�w%�<�                        6�l�3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��     @���    16:55:54        F�� @���    @��@    '�%      ]     �p5&�X6���B'�fB'�f7��[D"g�    &rJ+��&6:�?�  ?~'�5k�)3�
2���1�ï2�Em1��<    ,�7�(�?=B�=H&>1��> ǣ=�n�=�9��*2��:=�<(��<��=H<�=�&�=��l=�ٳ=ظ�=��=�D�>��>X�?�7�?]�D?V�>�)>t�>;[>D>u�@Wmq@j�@`�N@ae>@a��@b@�@b��@b��@bpX@b!�@aӨ@a�HAƖ�< ^?{@��                    E��6�k�G��@F�H�FT��E�Z<E��D2JgCE$ BB�b                                                A��C^/B�EA�A:A?m�@�0?�#�>��o=̝                                                                    ED��@�$�F#�A�d3��&F`z�B1p�@d�.�-R                                                    {@��A�j�A�j�@cm�8 ��^{�&�ާ*�_.���&)X>�&rJ@�1-/�G�5'��    >�t�A��)1���>8
�BZg?�<=M6�    B�.hB�.hC��C��C��{?O�@���3�2C�E�5���6�e�A/�7C��>�t�C![B|�jCCl9A�V�BE��C.A�"nB�2�@�A��M����@�A    @�AA�<�BV�?I6?I%�A��A%��A��=@VN�@ū\A��@��NF�bx3��&1I�*.��3㵂8ZV6�/5U�XG
�GAA�G�r�G̈�>]�jA��JC��@�b�?Q�?��>�@>��>�,>>��>�5m>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M̏�    @T]�@T]�{@��7f�!8�M�8���8���8��7V��6��o5�'�4��                                                A�F�?��d?��d{@��Cv�B'1MB��u@L�>�|n{@��{@��    C�)�BO�5���    6K+�7Z�@���3.B7i�    =�$�C�~g7��?B4��D�^C���B�g�A��kAV�H@��@)�S?��                                                ?MN�A���@���?ǃ�?+��>�	>&G�=�W�<o�e                                                                    DE�;F��E�ZD�d�C��C?��B���A��@�b�                                                @�7B��aBHu�A5V@��@nY?�#T>��7=�-�                                                                    C��E���E �.C�)C?�[B�;;A�7UA&�[@&��                                                @4��B[�iA���@��@RJ?���?�>D4�=D��                                                                    6�N�7a�A�>G	A��[?%�?���+Z��*df,�e�,�r+Z��*/0�-ܫ,�5,��s)��7q7|4�2�2��,�            7��14��2	�[2�g2��,𮷴�1��o77��44=    (ﳀ4$�4~.W�     =5��8�b7��{22�X1�-}�                        +��82L2�.���    22��2��            4��9�ؓ-n�>˩�>��?��?�0?v��?D�???F.I?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C?ȩG�)�>+FA��                                                B�                                          B��                    A�                @�p�    ,J��(2��4�?5�3ոF3��?Y��8��?��d7���;w�L;u&�6[            7(    5��J4�*Q    8�W {@��1X 5M�|8���Wz�{@��{@��{@��7��7rWG            8��    8��{@��    8��    ��Б7CF�{@��8��Z{@��6�I>6%�)5�~�    2���{@��6(�f6��7G:�7G:�<�|�    B��a6��F^eC���A0�K@'��B��d            @Ñ4�2�:�U4�2�<�T;ׇ�@���@��}@���?�Up?�>������
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/@L�<��;ׂ1@�b#@�
@�M?��6?&�>{�\����
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/;�T/ss)B�4��A5,�_< �c;��<<e�<D��;��;�d                                                ������qĤ�[Ğ��Ĕ�'ą��b��//̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7��5��;=��9�Z�8�>�A���A�\B>A���>R8X�7�G]L6�7�>�?TC~�@��I>�ħ6�7�K��8��7k��CGP�D�rD�d.D;uC��1CB#BN��A��+@�?�                                                A�P�CZrCFd.B�uB'�1A�B#@���@�+??�                                                                    E��GPG��F���F9�E� �D��C�OC Q`                                                C�ޡEO 
EAr�E�\D���C���C�}B-QA:��                                                                    E��F�;9F���F{&kF]�F.��E�9NE��8EG��                                                D)�D�>�D� �D���D�ǥD~nrD0��C߀QC�@�                                                                    22��6�sGXA4GM�Gv|FͲ�Fz�jF�}E��EP�                                                A+�@��                                                                        @�h�AioA��A��B1YzBc�B���B�C�<�<�<�<�<�<�<�<�<�<�<�<�E�#ME�,�Ec��EP�D���DW$uC�-�C�c�                                                {@��{@��{@��{@��{@��C[BQ�q,��f7�v4B�4?�H    A��A��{@�ξA���A��C�V {@��C��
C���C�V {@��@cm�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G ��>�dnD�=NGOD�=$D��+D��+A��A��F`�9B98XA�_kC��FK{0FK{0D���D���Fb�B93*                @J�?C�M�C�d�C�;�?   C�#YC��;C��;C���C��C���C��C� <C�)C�1"C�N=C�j6C��5C���C��C��AC�0�C�e�C��C��zC��WC�%�C�NC�[�C�9C�!!C��C��7{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�HjD ��>��>���>���>�ü>�&>��@>��X>�I�>�ş>�6'>���>���>�Ƴ>�9>��9>��h>���>���>���>��2@!��?�d        A���7ޠ�7�76��\]1Ǵ�TǴu�Ǵ
���z������������đ	�        >��B0�    ?��9?��9{@��@q9�FT�f6~�7��        :&�(?YVH?;p�?<�_?>�?@��?D��?K�?U��                                                ��b37�(A�  ?]s-?A�)?|4'�%      ]     �p8xa07��57�k7D�6\�5t��4b�%3��                                                                    7��6���6��6G��5�Y~4��3�X2�t�                                                                    G��@F�H�FT��E�Z<E��D2JgCE$ BB�b                                                4`*3C"^3A��2�$�1皯1 ��/�4�/ͧ                                                                    3$�^2v|A2t�1���1F�0"`�/rM.>}	                                                                    6�{�6o�6���6�d5V�4�?c3�i�2�-y                                                                    6�Zz6�Ƣ6�]�5�\N5N�4���3�?
3 �e                                                                    4[24'ѭ4Z�43���3)$E2���1�΀0�Q�                                                                    6�R#6�H6��J6*05|�4��_4m�3�                                                                    7�~�7�To7���7� 6Y\�5�t4�=O3���                                                                    4�34�8.4���4<�3���3 v;2��1FV�                                                                    7�~�7�To7���7� 6Y\�5�t4�=O3���                                                                    73m�7DA7�6_`�5�<4�V4@3%V                                                                    4:/b41ZJ4HI�3��3%(C2z�P1��M0���                                                                    73m�7DA7�6_`�5�<4�V4@3%V                                                                    �H�����5�U5I�4���3��Q2랱2p�;                                                                    /�d1���;�R�;��;3v]:n*9���9)�                                                                    5
��4��l4�V4��3t��2��
1���1��                                                                    6W8Y61�i6�>5k�4���3���2��2�"                                                                    �X�ӵ6.���XA��g�aY6�{YI�xʐ���]                                                                    ���:V������g\�N�5����w�8�N2�                                                                    5Vp�51#�5��4ja�3�r2�ؾ1��1�n                                                                    �d���:څ�J����Kj�״8N�J￲���                                                                    ��L��������;���-����'�����,��                                                                    ���^��.߳�	j�p{����r��3l�����v                                                                    (/�(�̊/L�C/@�/�/� �/x��/(B                                                                    *���+~��                                                                                            7��z7���7�o{7�B6\��5�Z�4�t�3�[                                                                    4;�+4X�4&@�3�D�2��2�
1��0l�                                                                    6���6��[6��6��5dB�4��3���2�u�                                                                    37��38�3"��2�3�1��1Ժ0+��/�d�                                                                    7� 67�}37���7�)6k�[5���4�,�3���                                                                    44{44z�3�"e2�w2 j1;0�Ǹ                                                                    7/�7�$7eeD6칣66�`5K��4I�S3�[1                                                                    5IJ5(��5�K5E�4P� 3i�2f��1�C�                                                                    6"�5��z69^�5�K5��4$ѡ3#$N2a8�                                                                    4"p�4��4S�3ڞ�3(�"2<]K1:r�0���                                                                    7V��74�e7�/�7�U6_�Y5yI�4v��3�S                                                                    5u�!5N�s5�6y5%T�4n3�sZ2� K1¨                                                                    3���3�<4w�N4Z!t4'��3���3I�3'S
                                                                    2���2���3H�30Dh3��2�y�2"{�26-                                                                    4o44�N�4�MU4M)&3�!�3u�G3L��                                                                                                                                                                        '�%      ]     �p6>�V(ﳀ                        6���3��v@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @���    @��@    16:55:57        F�� @��@    @��     '��      |     �@5���7��B'�fB'�f7ț5D"g�    (��+y.�6ʀ?�  ?s�B1�-3�\�)��/��=��J        .*�o>�5U>Fθ=�e�=��)<ڱm8��=1$��2X�:>��<(��<� =HB=�)4=���=���=عm=��E=�E&>�>X�@֮@���@u�?���?���?~�?Q��?�P@e�^@e�v@[�H@\d�@\�m@]��@^TJ@_�@_��@`@`Vz@`}�A@� =@�{@��                    E��!7cN�G���F��FTY$E�#E��D21�CEBB�|                                                A��C6�B�,�A��A?4�@�#?�_>��F=̠�                                                                    EDӑ@�EF"�A�[�4H��F`=�B10A8f�/8��                                                    {@��Buo�Buo�@`��8T��b�    ����/JQ4&3��>���(��@�Et��M�2��    >���B60K1�ɾA7Z�A"'tA!u�        B�3�B�3�C��C��C�R�>�g@�D�3��C�W6���6�߀A.�uCM��>���Cb�B�ECxd%A��HB���CSr�B�zB��NBbuv�al�/�Bbuv    BbuvBc�;���6彮76|JA���AA��A���@j�R@0w�@㈴?�b�F��74H��1�|l/0�4v$8��)6��35�j�GMDG=�oG�@�G��=?I    6�1�    >6O�>���>؛�>�3>�n�>��o>��q>��8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�M�    A��A��{@��7�N�9��t9d6�8���84�:7���6�V6T53��                                                A�> @��@��{@��<:
�BkA�B���@B��>�h�{@��{@��    C���BNt68�    6I$#7�@���3
Ǣ7C��    =<��C��7�7VB+PD��C�BB�\cA�VzAL	@ǧ|@��?�                                                ?A�cA���@�X�?���?�>�<{>��=n��<]��                                                                    DD��F��jE�oD��C��C>��B���A�y�@�);                                                @��B�|�BFq�A/ܑ@�j@�?��@>�A�=�C�                                                                    C�"E�~WEs�C�.C=�mB�I�A��A%R@%Z�                                                @2�BZ�A�� @��9@��?�S4?N�>A�=B��                                                                    7e0�7�?�A�	>D�IA�Qq?#��?�J�,J^;*�׌-$x,��,J^;+!�.#�4-�,�-LND*L;7�k�5��2�&0�=�            8
�K5%�2,�A2��2�&0�=��
�K1�_8
�i4�:�    ,��K5��5rX2�q�    =:(8{0�8O�3��O2!�0��                        .���3S�532�T2K�    3��72�Je            5 �`:�a�-*R?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Cq|
G�<=>-�A�@�                                                B�                                          B��                    A�                @�p�    .�!,*��5,��6J�4H��4vP�?�;�8Q�H@[��8�
;� �;�ٙ6'�h            7���    /�B6���    2��{@��    1[��    �9��{@��{@��{@��52`C7	#            6��    6��{@��    6��    6��7    {@��2��{@��2��Y2���1�r&    6��7{@��5鐦7��7�7�<�+    B�\�7n�E�WC䢶A���    C16G            ?�y�5�:'zt5�<�Ц<��Y@�b�@^į@)?�2�>�<U>��5����
�<���1��1����⊤�ˊ��<�p6�S�O�=0�*�]��,?�eB<�˰<���@��S@]cH@	v-?h;�>��>*7����
�<���1��1����⊤�ˊ��<�p6�S�O�=0�*�]��,=�E0��(�A6��6��w<ێ�<���>9�>x�L>��=ʺ�                                                ė�rĔ?tč˒Ă��fAk�=���_uêAb̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         :)��        {@��{@��{@��{@��                1l�4$�        :���1l�            CG��Dښ�DƴD<C�CaBO�sA�?@�u                                                A���CZ��CF�B�B(A�a@��s@??u                                                                    E���G[rG?F�u�F8��E��	D�o^C���C N                                                C���EO�EA}sE�D���Cյ�C�B-A:�                                                                    E�{F�=fF��5F{)F]-F.��E�8�E���EG�c                                                D)��D�B	D�6D��zD��!D~nAD0�C��C�@H                                                                    3��77]l�GXN�GM��GuEFͨ�Fz��F�`E�EP�                                                </Yf                                                                            @c$�An��A��B@B77�Bm�B���B�Ri<�<�<�<�<�<�<�<�<�<�<�<�E�-�E�5'Ec�<EI�D��DW#rC�/C�c�                                                {@��{@��{@��{@��{@��D\CB�.Gꫧf?�8C*�@.��    AB�AB�{@�ξH�HC��{@��C�k]C�k�C��{@��@`��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G8?<�3�D�0�GwwD�0�D��0D��0A>�A>�Fa��B;�9</YfC���FKz�FKz�D��8D��8Fb��B;�                @���C�A�C��HC�<�?   C�sC�k]C�k]C���C���C�?�C��C��(C�AFC���C��'C���C�`aC�E�C�7C�4�C�@8C�Z(C�LC���C���C��(C�/yC�Q�C�:9C�!�C�C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C�Y>��>�P7>���>�&�>ȿ�>�9>��h>�i>�ip>��>��i>�<>�6�>�
�>��x>�B>�/1>���>�3�>���@X�?���        A�b�ǚ�Iǚ\�ǚP�����L�Z���(�w����5��5F��5ĔG        ?�qB;�)    @�Gu@�Gu{@��@]- FT��7
<�7�        ;���?>D�??L?@��?C�?H�??Oڇ?[
?nI�                                                �|��7���A�  ?\�@�9%@��-'��      |     �@9@,�8���7�>7U�^6�m�5�J.4�3�6�                                                                    8r�7�ɥ7 ~�6�$�5�G_5�4c83�                                                                    G���F��FTY$E�#E��D21�CEBB�|                                                4��&4d#3�co2�02��1W|$0o�~/u��                                                                    3���3IV2�}�2߲1H��0��/���.�J�                                                                    7a��7C�n6�265Z5��24��D4��39�                                                                    7��7K�56�#�6��5�̟5��4H��3A�                                                                    5"�5 ��4!�3�z�3\��2�k�2/L�1(��                                                                    7���7x�A6�+�6:�p5��55 l4uB�3k�k                                                                    8�Ww8�dD7��I7>
�6��B5�!�5��4�                                                                    5���5��O5?�4{tg3�33T2��1��D                                                                    8�Ww8�dD7��I7>
�6��B5�!�5��4�                                                                    8��7���7>y|6��~5�LY5B�t4��r3�m�                                                                    5k�5!K4���3���3`�2��*2��1$�                                                                    8��7���7>y|6��~5�LY5B�t4��r3�m�                                                                    �|S���C5�x^5���4���43H3z�V2Ȋs                                                                    1w6%1��6;8�:�N�:���:'j�9��S9n|                                                                    5�5�SL4�C4<�o3�La3�2b	1Yq�                                                                    7#Df70�6R�Q5��q4��43��3_92f"                                                                    �(����̵EX������ҹ;�������                                                                    ��w@��f��Ҍ?��沈���DH#�����*                                                                    6"�?6�}5R�4�m(3�%�33�2_��1g�                                                                    �1a�����|s���P�K�紛^:��y��7{                                                                    ��W����	�
)�~�B�ل,�(de�iᅯ��6                                                                    �����u��� ���������-��N���f��                                                                    (l.�(ʤG/��/t�0�9�0���0��/�>�                                                                    *�v�+[�                                                                                            8��+8x{=7�	7Ex6���5�95�u4*��                                                                    5��4�:�4f��3�^�3"E2u�71���0�_�                                                                    7�7sk6�k_6D�5��4���4v+3;/�                                                                    4�k3�13a��2��x2'��1��T0��!/��                                                                    8�8nZ�7�ŵ7B�_6�7�5��}5(�4K��                                                                    5�4���4]13�z!3-R2��1�\�0�X                                                                    8Z�7��7�?�7 $V6{�5��;4�l	3�2�                                                                    6ՙ6
;�5��Z57�4�~3�l�2�{x2��                                                                    6�_06�{N6��6hZ5J�$4�-�3�a�2��                                                                    4�ڀ4�hY4�%�4��3g�*2��W1��O0��h                                                                    8&��8�@7�j�7C��6�u5��4��}4
��                                                                    6>v�6(�5�5_��4�a,3��Z3�2��                                                                    4�`�4�'"4�ڹ4���4f��4G!3ʣ43�i*                                                                    3�N33�m33�l3n��3:]�3 ��2��[2\v?                                                                    4�vF4�*4��q4�q�4��Z4B�E3��3��l                                                                                                                                                                        '��      |     �@6<��,��K            F��    ;��6���3ಲ@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��@    @��     16:55:59        F�4 @��     @��    '��      �     ��5��;7��B'�fB'�f8��D"g�    )Y +F,7R
?�  ?n��    4�)A\a/��}�b�        /���+�7r>�1>)��=�(�=F��<9��7��^0]چ2fRs:4��<(f�<�
>=HA8=�(�=���=���=ع�=��=�E�>�J>Y@��@���@#^}?�R�?�c�?���?���?��@h��@h��@^��@^�`@^�g@^uL@^i�@^c@^b(@^f�@^m�@^su?X=���?�o�                    E�K�7u�G�LF�E�FS¼E��_EK�D1�!CD��BB�(                                                A<|C�-B���A�P�A>�<@�� ?��}>�s�=�k�                                                                    ED��@�P�F#IA���2�AGF`��B2��A@R�.ЖH                                                    {@��B�5oB�5o@f��8�zy#*�    ���/3(�&�.�?   )Y AF�	��Fb2��    >�y�Bf�l1���Ap��A�gA�X        B��FB��FC�1C�1C�O>��Y@���3��C���7�g6��_A1��CS�l>�y�Ch��B�|
CrܕA��B�UCM�_B%B�B�i6B/o��S+��rB/o�    B/o�B7� ;��    6��A��AIG�A�'�@�6@)��@��p?��F���2�AG0K8a-o��2��x8ı6�95�̋F�kfG-�6Gݭ�G��=��y            >f�>�4w>Ť�>��>�W�>��n>��>�b�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�Y�    Ag�Ag�{@��7�t=9��Z9u{�8گ8@Z�7���6�p6Y�5/&                                                A�V3A_�A_�{@��9i�aB�/�B�|9@Cy}?��{@��{@��    C��BJ	�6��3    6JM�7 �@�B�2���7_[-    =�d�C���7�	B:=Dy�C���BZ��A�v6A:G@�2@ �>�g�                                                ?/��A}�@��?��?
٢>�
�>A=JXN<Gt�                                                                    DB�F���E�u�D�f�C�*$C<6B��DA��@��                                                @��@B���BAEA(�@�i�@b?�,2>�]�=�{�                                                                    C��E���EcgC�`�C:B��uA�K8A �@"B�                                                @.!�BW��A��g@��@��?��
?�:>:�==��                                                                    7q�97�$A��>J�SA���?)%\?��C,���+d�-��-�p,���+�.�4.f�-�;�*���7�\"5�2��z1
�-���-���(A;7�٦5w�2��2�2\2��1
_��٦2��7��24�Q>    3,�4�	�4�'3?�    <���8r�88Y(M6_f�5�o3�                        5�E6:3j6.�~44��3/�7�k2ǻ�            5�==��a-/��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Cr]�G���>+�4A��                                                B�                                          B��                    A�                @�p�    /9Um+*��5S��6X#62���2�Z?Ɗ>8\��@��98,�L;ɍ�;�                8�n    3J��6�ip    2L{@��            ���F{@��{@��{@��7�.�8K�e            7P�    7P�{@��    7P�    6��
    {@��2L{@��2EBR2.��0�6    6��
{@��6���7�M�8T��8T��<���1�+�BqN�7�DE�W�C�g�A�c�    C:��            =�k�5�`8l�|5�`<�[{<�5�>��>d1C>��=��=���=n-͋D�,{��k��ޭ���o?���Y��,{�p���T�=�a�+Jq�C�=�Q<�U<�(>���>a�s>#=��=,Y�<�\�D�,{��k��ޭ���o?���Y��,{�p���T�=�a�+Jq�C�<�j�    *��h6L6>6�LL;y�\;�:��=I��=&�=�                                                �~y����M�������Ľ��Ġv�ă�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         9i�a        {@��{@��{@��{@��                1��                1��            CG��Dک?D� AD;�~C���C)BPJ�A���@�]�                                                A���CZ�?CG AB��~B'��A�)@�J�@��?]�                                                                    E���Gz�G*F�c�F8�TE���D�S�C���C G�                                                C��REO>>EA��EwD�p�CՋlCk~B-A:��                                                                    E��F�CF���F{-$F]F.�_E�7�E���EG�G                                                D)�eD�JKD�
D��xD���D~m�D0�C�~^C�?z                                                                    7�k7n�AGXq�GM��GmF͗(Fzs�F��E�JEP
o                                                                                                                                @D�nAT��A�J�A��B("TBW܉B� B���<�<�<�<�<�<�<�<�<�<�<�<�E�F�E�K_Ec�E<SD��QDWC�,�C�b�                                                {@��{@��{@��{@��{@��D��A��.�=v�f�=8�/Cy$�@a�    A|�A|�{@�ξ�`���`�C�n�{@��C��rC��zC�n�{@��@f��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��<2v�D�*�G�^D�*�D���D���AϷAϷFc�=B@��    C��iFKz�FKz�D��:D��:FdvNB@��                @���C��vC��;C���?   C��VC��rC��rC�"C���C�s�C�$jC���C�f�C��[C��JC�P-C��C���C�rC�2�C��pC���C���C��;C��SC���C��C�E�C�:[C�"WC�WC��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C���?	%e?#�>�A�>���>︵>���>�>�E>ӿ,>ΐ6>�~>ľ�>��>��>�mF>�(V>�l>�B>���>�y�@T��@��        Bܠǵ��Ǵ�jǴx��\�/��.F��-��ƎG��!b��!���!a�s0        ?KY�B3��    A?A?{@��@��FT��7�(�79        ;�Y?-"4?-�?.$y?0I?3��?855?>c�?F�e                                                �'�7;K�A�  ?k�A	��A z�'��      �     ��9O�O8��c8	�7fg}6�Hc5ه�4�3��                                                                    8�OZ7��7-"�6���5��5	c!4f�3��                                                                    G�LF�E�FS¼E��_EK�D1�!CD��BB�(                                                4�F&4+��3���2��2*\1dW0xIE/r��                                                                    4	�v3Y2��"2�1W0�0�6�/���.�>~                                                                    7u�7T�6�"c6C1�5�94��Y4yN37�'                                                                    7�g7F�A6���6�5�#@4�40�3+]4                                                                    5'_4�k�4b�3ޛ�3X �2���2�1�o                                                                    7���7sl6Ă�62�5��k5Au4W$�3Qq�                                                                    8�,Y8�i8�X7I�6�I�5��5н4�:                                                                    5�4�5��@5�g4�fE3�`3Z�2��D1�-�                                                                    8�,Y8�i8�X7I�6�I�5��5н4�:                                                                    8jj8q�7H6�Ww5��@5J@�4��n3�O_                                                                    5�%57�4��3�^�3k�2�_C2�1`@                                                                    8jj8q�7H6�Ww5��@5J@�4��n3�O_                                                                    �S�+���5�]�5��L5
��4G� 3�e�2�>�                                                                    4�nT5t��91[�9 Lf8��8�-�8A08��                                                                    5�{B5�BK4�C�44�`3�F$3�h2F'}1@��                                                                    71��7s�6_�5��4�͵4<p3b��2`E{                                                                    �6���!\���ŵT���
d�ߠW�毲
B=                                                                    ��Q��1ٲ�Ȁ��`����&�Ph��q���o                                                                    60��6��5]��4�y
3�vP3;�!2ca1`�v                                                                    �@�a�%Ѹ���I��&�ZG���5��8"���;                                                                    ���ҳ��&��I��6��8l�3k�s�ȯ�@!                                                                    ������?�-�~��w<��f�8Y�V+d�c��                                                                                                                                                                        ,;�,�y-                                                                                            8�d8���7�� 7T�t6���5�ܔ5�4(��                                                                    5�5�4yO13��3-�f2���1�S�0�`%                                                                    7�4>7�ڌ6��6SE
5�FN4��"4!/39@5                                                                    4��4O�3t:�2؃l23��1��0�`~/��t                                                                    8�v8�t7술7Q��6���6\�5/,L4I��                                                                    5�(5��4o&?3��Y39�h2�uQ1�m0�N�                                                                    83f8�*7�  7,�`6���5�w�4�ҷ3�c                                                                    6)_O6��5Òm5Ef$4��r3�d#2���2 8�                                                                    6�26ӟ�6�H�6�65Y�4���3�\�2�R�                                                                    5��4���4�	�4��3x��2���1³0�9�                                                                    85"`8 
O7�'C7Su6���5��5 �4	                                                                     6O�66�65�M5qC�4�4�3�zd3=�2��                                                                    4�v4�gf4��)4�[4w�O4)@3��3��]                                                                    3�83��)3�]�3�ż3H&:3�A2��(2Y�                                                                    5��4�)4��4�Ĭ4�\�4N�24 a�3���                                                                                                                                                                        '��      �     ��6=�m3,�            H!H�    =�aD6�A�3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��     @��    16:56:02        F�r @��    @��@    '�Q      �     �5�\�7��gB'�fB'�f8)N�D"g�    (�у+;�7LUG?�  ?u��2$y�4@��)�I/�.��O�        /�3*�R%>���>� =�wp<��s9��{2��|+�72`>::%��<'s�<���=H9�=�&#=��5=���=ع%=��l=�E�>�_>Y.@� �@qX�@ծ?��U?�C�?m�s?G<�?)D2@h�o@h�1@^�S@^��@^�_@^��@^�L@^��@^�@^�_@^��@^��>X�=���?j:$                    E��<7��G��1F��FR��E��E�*D1AMCDFBBH&                                                A��C�BB�%�A�g�A>�@�S?�)>��w=���                                                                    EE8�@ʘ�F#KKA��F10�Fa%�B4��A5�.Ķ                                                    {@��B�D�B�D�@p��8�B2��?    �t[�/.X��S�?   (�уA������2��X    >�|B@-�2� A/�d@��!@�
�        BO��BO��C�יC�יC���>�Op@�'h3C�X�7�6�s�A5��C>8�>�|CQ�B���CXm�A�dMB��HC5��BK$B�ȉB[
@��T    B[
@    B[
@Bal�            A��A74iA�L�@sp@�|@��?�)�F���10�.���+�;1C o8��6�l�6=�G~�GEP�H�G�/>$�            >@��>�N0>�N�>��>�${>��>��n>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M㤞    A}:A}:{@��85�9��J9��:9M�8e��7���7�G6=]�5`6�                                                A�nJA��1A��1{@�Π; :B�}(B��@MN(?f {@��{@��    C�7BI��6�ƹ    6M��7�@� �2�Q�7w�F    >$�C�~7��pB�eDkECm�pB2��A��A%�3@��O?׏�>��o                                                ?5�AmG�@��?_=g>�Lv>r,�=�E=)��</j�                                                                    D?�*F���Eʺ�D��_C���C8�B���A� d@�V�                                                @��B�i�B9�fAJ@�.9@��?~:1>�; =�*                                                                    C���E� �E_NC��C4ѼB���A�L�A�
@�*                                                @(?_BTC�A��O@�$�@{�?�ѳ>���>1�+=6��                                                                    7�H&8�0A�>SE�By?0�?��P,�a+m�-�
-�,�a+�i.�̍.��-�*� 7��4�2��,�K�/m:j/m,�)Y��7m5
��1��Z2��2�s,����m1��7m�4�z�    3���4��4��.֧�    :�##8Jx�8C�66��X5O�F2lG�                        5D��6�3L6�հ3�f�3t�7@��2�7D            5�>-?-C!?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Cbq�G�(>+8�A��                                                B�                                          B��                    A�                @�p�    .�;8*���5Q��6}s?2���1X��?�X@8O�@C�8#�I<nf<6n�            7���    /�G�6�m�        {@��            ��[1{@��{@��{@��8( 8�m�            7���    7��{@��    7��    6�k     {@��    {@��                6�k {@��6޽=7�%�8��8��<���1�H�B�2�7EE�& D+�QA�Lt    C-            <Jl4���6�2 4���=�<���=f�-=��<��(;�/Y;f>_:⌯�5ËO�������ފʡ���O��O�p�m�TG��=�q�+t��j�<I~�=<�m=e֛=��<��;�ȴ;d�y:��ߋ5ËO�������ފʡ���O��O�p�m�TG��=�q�+t��j�8m)u+z�*�Դ6pgt7�9���9#.�8��83R57�s[7�s�                                                Ą�]ă�1�{���lo�T���3�G��CèB�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �; :        {@��{@��{@��{@��                                                    CG/�D�v�D�߰D:��C��CrBOo}A��@� &                                                A�/�CZv�CF߰B���B&�A�r@�o}@�? &                                                                    E��EG�~G@�F�IF8�)E��D�-kC���C <�                                                C���EOyqEA�WEc�D�M�C�QCO�B,��A:��                                                                    E�F�J{F��9F{2�F][F.�E�5�E��2EG�O                                                D)�jD�UD�<D���D�ɶD~l�D0� C�{�C�>                                                                    7@��7�#/GX�sGM�6G[�F�{�Fz\lF��E�]EP�                                                                                                                                @p��Av��A�i�B��B9�3Bol�B�J�B�<<�<�<�<�<�<�<�<�<�<�<�<�E�e�E�d�Ec�!E'_D��nDW�C�'C�`�                                                {@��{@��{@��{@��{@��Do8B/!.V��fZr86��C�B@r �    A���A���{@�ξ��T���TC�%{@��C�EC�EC�%{@��@p��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�/;�iD�1G��D�0�D���D���A)�A)�FeF�BC��    C�\jFKv�FKv�D�ұD�ұFe� BC��                @�p�C�s�C��WC�G�?   C�tOC�EC�EC��)C��C�\FC��C��;C�o�C��C���C�m$C��C���C�v�C�#�C��YC���C�G�C��C�]C�,C��C�9�C�9�C�"�C��C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�w`C��!?�8?
��?#?B�?Z>�<>�U0>��>��>�4>�A�>֋�>�Z>˘s>�J�>���>��>�[&>�	?>��Z@?�Q@�I        B�Ǚ3Aǘ�
ǘ,v���a��U��1�uP�c����w�������tĈ�        ?�B#?C    AD�.AD�.{@��@�FUAG7���78S        ;�ۯ?D3?D+�?E��?H+�?K��?R>b?\�5?m�D                                                ��Q^7�A�  ?V'�@��A z�'�Q      �     �9t��8�֋8%��7���6���6�5�*4]�                                                                    8�z	7�_�7Q:f6��%5��t5)�4;��3=�i                                                                    G��1F��FR��E��E�*D1AMCDFBBH&                                                5 es4M��3���3��2O��1�w�0���/��                                                                    4"/M3��q2۪A29S�1��0�n�/� v.�iZ                                                                    7��	7~ߋ7 �6l��5��/5p�4?�f3oI&                                                                    7�g47T�6��
6v�5���5 � 49�N3Diz                                                                    59��5�r4]��3��`3ie2��n2"`1+�r                                                                    7�(�7��6�O�6<��5���5)4c
3p#                                                                    8�c8ե(8��7pZ6�3�6�59E%4?�e                                                                    5��5��5&�M4��&4��3��2�KJ1�^�                                                                    8�c8ե(8��7pZ6�3�6�59E%4?�e                                                                    8.�8 ��7i�6�l�6�5q��4���3��                                                                    5/��55�t4�:�4��3���2��24�1=��                                                                    8.�8 ��7i�6�l�6�5q��4���3��                                                                    �B�2��6�P5�5-��4��3���3��                                                                    5"B 5���7���7��/7U86���5��s5z��                                                                    5��5�QO4ж34@��3�~�3HI2Q@�1](�                                                                    7R��7=,%6�ڏ5�8�5�84c�3�2S2�0�                                                                    �X+��B$j���뵁r������	��#��4@|                                                                    ��ǳ0L�
6���/���X�����#�؍&                                                                    6Q0b6<�5��[4���4�.3aź2���1�4�                                                                    �c�Y�G.&���"������YH�!�                                                                    �����b�93��H��۱\���P���w�                                                                    ��3���
մR�X��l2�'/��c_ٱ�����{                                                                    )�M)�d�,4ј+ٰ2+`z�*�}X*_)�t�                                                                    ,h��,���                                                                                            8��t8���8��7��
6�ɱ6�58��4[��                                                                    5:�r5#v�4�004�Z3SSo2�oG1�3�0��l                                                                    7���7�a�7��6��5̪"5[�4J��3q�                                                                    46�!4 !�3� �3
�2Z�`1�a�0�c�0��                                                                    8�6"8�A8}�7}�U6ӊ�6"�$5\T4��                                                                    53+�5�L4�U4G3bQ2�Ty1��Q1�                                                                    8/j�8��7ψ
7R�6���5��Y5Ň4E                                                                    6Hy�6405�-�5p�4��3���3�-2'*k                                                                    7�=6���6���6)Ĝ5��:4�ݲ3֔�2�e\                                                                    5" E5~64���4BE3���2δ91�<1Y                                                                    8Ve�8@�7��F7�c6�%h6Ǩ5"F�42�@                                                                    6u�6\�6�L5��b4��4Q�39up2LPJ                                                                    4��G4� 94�PX4�^4�lv4P�44^T3���                                                                    3��3èB3�C�3��{3t��3(��2���2� :                                                                    5�V5�?5	�4�%s4��4 #4!Ƞ3�ƿ                                                                                                                                                                        '�Q      �     �6@��3���            F�/�    =-��6���3濤@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��    @��@    16:56:05        F�� @��@    @��     '�      �     ��5���7���B'�fB'�f8�mD"g�    (���+7��7X,�?�  ?yvr2T;�4P��)��*/��n��        .�Uh*��V>��@>	j�=�Խ<��k9G	�2Vb�+y��2S�E:��<&au<�Ǜ=H.�=�!�=��@=���=ط�=���=�ER>�^>Y?@ͮ�@p��@��?�m9?���?\�
?-�>�)�@gdP@gۧ@]��@]И@]�Q@^p@^-�@^T�@^zf@^�u@^��@^é>���=��\?�~�                    E��7��G��	F�jFQ��E�"
E>
D0k�CC"�BAu�                                                A�UC�2B~�;A�9�A=�@���?�4>�Շ=��                                                                    EE��@���F#�uA��1�Q6Fa�B6k�A55v-l�n                                                    {@��Bo��Bo��@wK�8��^#Ϫ�    ����/(�Ϧ���?   (���AE��ϔ�2�    >���B�Z2oA��@�kV@�JU        B���B���C�<�C�<�C��>�J@�=?3"��C颭7�6�HA9��C��>���C*�B��CD�_A��6BQxC �IB�B���A�-��S<�    A�-�    A�-�Aؕ?            A��A�A͋�@Q��?۵L@�&U?�2�F���1�Q6//>�,[��1���8�,N6�$6�KG	
�G%�H
)�G�q/=���            >@5>�:`>�>�>Շ�>��>�M�>�w�>՗r8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M뗱    A�E+A�E+{@��8�P9��\9���9y�8l�37�<7��6Fcw5rl                                                A���A���A���{@�Π!�QBD�B���@Z�i>�{@��{@��    C���BJNs6�k[    6R�O7f@�ϕ39�7��    =���C�+�7�O�A��D[�	CI�+B��A��XAԛ@~�j?��]>�k                                                ?
BJA\j?@`�v?7�>>ͷ�>TWz=��<=�0<��                                                                    D<�AF�O|EľXD���C��C2�^B��BA�mr@�C                                                @��UB�dB1EGA n@�̾@R@?o�>�NC=��                                                                    C��]E��EZC�1iC.��B��;A�>kAU|@O�                                                @!uWBPdA�P�@�Yx@�?"�>�?>%�{=+b                                                                    7�d8�A"�>[��B�?70?��,�o+~ �-���-*[,�o+�k.�ޱ.i�-�q@*ʱ�7���4�D2���+�)/�
�/�(��6|5bb1�22���2���+�7��|1�ܔ6}�4�Vw    2��I4��4�|.MD�    :���8�08"��6�4�&C1���                        4�\�6��6'�2Ǚ�2���6��82���            5]Z=��-�N?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CQ�cG��[>,��A���                                                B�                                          B��                    A�                @�p�    .�܏*��75�q6��'2�1��`?�D�818�@M6�7�?�<��<�&6�KJ            7�<�    0��46oJ        {@��            �(5�{@��{@��{@��7*b�7Ť            6��m    7C�[{@��    7C�[    6o8G    {@��    {@��                6o8G{@��6��&7��7��7��<�6I10/�B��?7 V0F��C��A?��    C��            <�pj4�I7%��4�I=~�<�X=�L�=ws�<�N<p��;��;u碋 ��<}�������.�ʆ���R��<}�p�U�T+i�=�-�+^$�U�<�×=z�<�H�=Ĉ�=v�<�7<o�@;릥;s �� ��<}�������.�ʆ���R��<}�p�U�T+i�=�-�+^$�U�8�Ң+���*�?6tV�6���9��9k�S8�8�]�8<�81��                                                Ā�)�|'��o��]�l�D�� U����P�pU̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �!�Q        {@��{@��{@��{@��                                                    CE�-D��DƸD8�C�DC(BM$�A�M�@���                                                A��-CZ�CF�B��B%DA�(@�$�@M�?��                                                                    E��G�!Gj�F�#lF8k�E�sQD��|C�}C )                                                C��hEO�0EBhEHOD��C��C(ZB,��A:jG                                                                    E�vF�SfF�� F{:F]�F.�E�2~E���EG��                                                D)�mD�b	D��D���D���D~j�D0��C�x7C�;�                                                                    6��87�F�GX΃GN�G@}F�U�Fz;vF�TE�LEP9                                                                                                                                @pBcAw��A�J�Bt�B;�Br��B�]�B�=�<�<�<�<�<�<�<�<�<�<�<�<�E��JE��Ecc�E	�D���DWHC�C�\�                                                {@��{@��{@��{@��{@��DݟB��.J���f:+8;��C�ql@�u�    A���A���{@�ξ0�[�0�[C��Z{@��C��C��C��Z{@��@wK�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�;��1D�6QG��D�6(D��D��A�%A�%FfbMBG�    C��FKn�FKn�D��ID��IFf�BGa                @���C��gC��dC���?   C��oC��C��C�[^C��C��C��PC�f�C�OC��C�|�C�4�C��C��+C�?�C���C���C�(�C���C��JC�YC�8aC�"�C�3C�8C�#bC��C�Q�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C�)?U�?�T?�?
�?�?��?�>�b�>�J�>�Y>���>�ʦ>��e>��>��>��>��>��>���>��8@��?�MH        A��Ǆ�6ǄNhǃ�Ū����S����:�F����^���_��^�Ă��        ?�EB
v    A[q*A[q*{@��@Xd�FU�7��+7ui        ;Z��?DT�?E�?F��?Iς?Nb�?U�U?b}?x�E                                                ����6���A�  ?F9@��eA z�'�      �     ��9z�;8ʱ�8,4X7�C�6Ћ	6a�5��4%�                                                                    8�r
8 �7Y�}6��+6� 53�z4JV3Qx�                                                                    G��	F�jFQ��E�"
E>
D0k�CC"�BAu�                                                5�L4T� 3��[3��2Z��1�x0��/�g                                                                    4&[�3�g�2�a�2A�1�GC0�͑/��.��?                                                                    7�F7��U7w�6v�e5ʥ�5+�4M�3��                                                                     7��7DV76�:�6��5�NQ4�ؗ4&o�38r�                                                                    53�4��4?aN3��3YK
2�:2"f1 �}                                                                    7�Zh7o��6�+&6-	�5�&�5�\4Kl<3ao�                                                                    8�C8�8��7u�46��16[W5@@4K�                                                                    5���5�35"MF4�k�4m3���2�t�1�s�                                                                    8�C8�8��7u�46��16[W5@@4K�                                                                    83�8(�7jy6�}>6��5wx4�=�3��                                                                    53BN5<Y�4��Y4#�3�Fx2�3�26�1EhC                                                                    83�8(�7jy6�}>6��5wx4�=�3��                                                                    �N~B��%�6-��5�ep5<@�4�/�3λ�3k,                                                                    4��5"�8I�8�7��V7��6�8�6	�N                                                                    5��>5��N4��416B3��3��2;�C1P-                                                                    7Z�H7H�6��P5�"<5�4k�3�t2���                                                                    �_^[�J[Q��5���U�֒��En�0�<�G�                                                                    �"d�Y�x���~$�ű�{a�0f���s                                                                    6Xɿ6F��5�4�Y4��3iڕ2��"1��                                                                    �jw.�O2����ݶ(�%��9�ֵ���5�0�9                                                                    ����}��>�����̲��h�0�����D8                                                                    ��"����_�\+6��,�0�K�r�̱����u�                                                                    )�gb*!�,���,f5�+�6�+t�Q+0*��q                                                                    +��g,&��                                                                                            8�Q8�/A8$l7�R6�s&6�5E��4q>�                                                                    5@k�5*`4���4	Y�3]�22���1��m1�m                                                                    7�p[7���7�6�5֩�5$E4X�3�B�                                                                    4<�)4&�`3��K3^�2e[1�
20�~�0щ                                                                    8���8�T�8��7��6��6,�5ḱ4��                                                                    58�n5#_4���4c03m�2�q}2��1# �                                                                    85B8$57غ�7\�76��5��5z4!�                                                                    6O&�6;�k5���5|2�4���4�S3#�%29�                                                                    7x�7�c6�"_62RF5��4�$3��3�)                                                                    5'eJ5�(4�'H4K˾3�͹2ܻ�2��1�x                                                                    8]��8H��8q�7��6��C6g5/\:4E��                                                                    6}/�6e^�6]�5��4�7(4&�,3Hig2b"#                                                                    5�/4�B4�xC4��4�l�4_ _4%Q3                                                                    3Ԉa3˼+3�xt3���3���34M�2�X�2�?`                                                                    5 �V55I~4�nk4�4�Z�4.��3��B                                                                                                                                                                        '�      �     ��6Enx2��I                        6��3뫅@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��@    @��     16:56:08        F�� @��     @���    '�      �     � 5�+F7���B'�fB'�f7��D"g�    )�w+j�7TZ?�  ?}l�    4[Y(���/����i�        /(ښ+"��>�f'>5�=�8�=X|�<U�7��R0���2A��:
�><%k<��!=H"r=��=���=�դ=ض{=���=�D�>�X>YQ@׃.@��@��?�|w?���?�6?u��?M�o@f�@f��@\��@\��@\�@\ӣ@\��@] �@]K;@]s�@]��@]��?��P=�O{@��                    E�8�7q��G��AF�x�FP��E�ZzE ��D/��CBH�B@�i                                                AC(�B}�)A�92A<P�@� ?�m>��}=�\�                                                                    EE�@�:`F#��A�M�2�'QFbk0B8<�A;�.��                                                    {@��B�PB�P@}��8]6��O�    ����/;o�&��>���)�w@1�P��<1��    >�y�A�c2�@�4��3��2,�        B���B���C�9C�9C�P�>�+r@�G�3'g�C��7rP6��yA=�B�a�>�y�B��BBF C'W�AsB�yCj�A֦WB�Z@юi��    @юi    @юi@���:��6��)6���AMXA@�\�A�5�@+�?�,@��}?�&�F��|2�'Q0EJR-r�2ų>8I�6?�c5�{F̑�GE�G�эG� =FG7    8&�    >+	�>��>�|>�Le>���>�W�>��:>�S�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M簥    A5*_A5*_{@��7�Q29�'>9sRG8Ү58B�`7���6��u6)�5M%�                                                A���A/~IA/~I{@��;�n�B	-�Bm%@j>�>�}{@��{@��    C�N*BJ�\6���    6W��78@���3�7Uij    =D�\C�B�7���A�K�DQ|C2�UB]A���A�@o0?�%�>�O�                                                ? �gAQ@Fy�?,Nz>Ŭ�>L`,=�S\<���<
�                                                                    D:bF�PE�o�D���C�$C/^AB���A���@���                                                @��B��B*0�A��@�9?�Ԕ?dj>���=�ſ                                                                    C�FE�@�E�)C�^�C*Y�B�;�A��$A�$@z                                                @S�BM`�A�zy@���?�!�?uYM>ݪB>��=#Cu                                                                    7gfE7�F�A)5H>d�@B[?>��?�M,>7*��-H
v,��e,>7+�.hZ-��-[,�*h7�3�4���2Q�-��'/l��/l�*��Y��$4�d 2�12N��2NA-�.�6�$1�ܶ��4�?!    0nq 4�n4�9/�ڿ    :��(7��>7��4�!3kk_0�7j                        2%TE4�(~4�l�1��{    4�I�2��d            4��<S,��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C98JG��>,�]A��                                                B�                                          B��                    A�                @�p�    .�g�*��24��j6G��3��2Ƿ
?^�_8��@��7��O;�q7;ƍ�                7k��    0���6H}8    1���{@��    1,y�    �bC,{@��{@��{@�ζք6v1            4��p    4��p{@��    4��p    6H��    {@��1���{@��1���0�ݺ1-F�    6H��{@��5��7&�6��6��<��f    B�g^7	w�F$��C�Ю@���    B�%            =�xb4�w8>Z14�w<͉�<�J2>�G�>�� >	Z�=�U�=2m<�:�� !� ��U��ޙȊ�]B��㱊� �p�@�T ,�=r��+;8�6=��<̈́"<�>`>�ie>�?>� =�h�=�)<�>V� !� ��U��ޙȊ�]B��㱊� �p�@�T ,�=r��+;8�6;�J    '�06-r6��:ޙ:���:qh9���;7�C;��1                                                ĳ�İ3ĩ�Ğ��ďM�vi�D�j��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         9�mZ        {@��{@��{@��{@��                .6R�5�!�        ;�M�.6R�            CDjDٺjDƈ{D5IFC�)PC��BJT�A�X�@��[                                                A�jCY�jCF�{B�IFB#)PA���@�T�@X�?�[                                                                    E��	G�vG��F��F8.�E�@XD��CC�`C "                                                C��EO��EBDE,sD��CԺ�C 1B,��A:J`                                                                    E��F�Z�F���F{@VF]?F.�E�/�E��OEG��                                                D)�?D�l�D�(`D��mD���D~i�D0��C�t�C�9?                                                                    4�I�7s< GX��GN(�G#�F�0@Fz�Fz7E���EO��                                                ;�{~                                                                            @UP�Ad'�A���A��~B0�kBc�B���B��w<�<�<�<�<�<�<�<�<�<�<�<�E���E��hEc3�E�\D�u�DV�C��C�Xe                                                {@��{@��{@��{@��{@��D�B��.���fM�8|C��3@���    Aa.zAa.z{@�ξ:�Y�:�YC�4t{@��C�JC�*C�4t{@��@}��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��Gد;���D�7�G��D�7�D�srD�srAs�As�Ff�`BIy�;�{~C�&�FKdxFKdxD��?D��?Fg��BItf                @/F�C�%=C�srC� ?   C�
TC�JC�JC�;�C�UC�g�C�v�C��"C��C�r�C�Y�C�:�C��C��*C���C�V�C�C���C�K�C��}C��tC�zC�E�C�3�C�6�C�#�C�C�>�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�	mC�@�>�{I>�ֆ>��m>���>��g>�$1>��~>���>��>�4*>���>��>�/>�Y�>�l�>�4u>�*>Ⱦ
>�F>��@?ٸ�        A�1��N���Ng��M�XŃ����[��&���x��Ec�Äh�Ä��Ähĳ)�        >�ޒA���    A&A&{@��@`��FU��7��37�t        ;[x�?7�-?8��?:
 ?<��?@W�?F�?N�|?[��                                                �Zv6�SA�  ?k�A	��@���'�      �     � 9@DU8��&8�r7r�{6��.5�!�5
��4pl                                                                    8r�7��70��6�`N5�;5��4/)�33�R                                                                    G��AF�x�FP��E�ZzE ��D/��CBH�B@�i                                                4���4(�3��_2���29a91���0���/���                                                                    3� �3U,Z2���2!
51j)�0��/���.��                                                                    7e�.7SDJ6��6L2i5��!5.41pL3a(U                                                                    7m��7��6Ttt5�T�5[��4�I4�3��                                                                    5'c4�a�4/�3�gg31�h2�T�1�D1R�                                                                    7�[�7/��6��c6
P85�>F4�04"jq35��                                                                    8��8��W7���7I�6��-5�C5!�}4*��                                                                    5�Z5��4���4rp-3��3U�2��1���                                                                    8��8��W7���7I�6��-5�C5!�}4*��                                                                    8G�8$E78�-6���6 SN5O#/4���3�                                                                    5	f5�J4l�3�U3b,L2��2�+1".                                                                    8G�8$E78�-6���6 SN5O#/4���3�                                                                    �l�\�� �6��5��5 ��4w��3�2�3��                                                                    2�[j3��	9W�9J�8��18i�7�J�7R�                                                                    5��G5F��4�\{4�3��2�=/2#n1'��                                                                    7)��7 �6Um5���5A4G�3v�2�:                                                                    �,9��!lF���"�bnݴ�-��Y��j�+�                                                                    ���B���U����Ѭ粨��m�Y��D�ν�                                                                    6(=K6X5RV4�Ѐ4�3EYj2t��1��\                                                                    �4Lq�$�C����%�j���U���O���                                                                    ��Jݳ��F�b���t]��$�Fuٰ�ѻ����                                                                    ���ڴ��3�A���U�.C�Q�^�}������                                                                                                                                                                        )t�*n�                                                                                            8�Yh8��Q7�H�7\��6�b=6{P5)�E4Mx�                                                                    5��5be4y�3�f3:Q�2���1�&�0��                                                                    7�^7�,�6�C�6[V�5�r_54: �3aL                                                                    4��4�03s��2��.2@�A1��}0ұM/�4$                                                                    8�b�8�r[7�?*7Y�?6���6��5J^�4u�                                                                    5�5��4n��3�)�3GF�2�-1�;�1
ԯ                                                                    8��8��7��77��6�Պ5�e�4��4}                                                                    6�6�'5��i5Q�4��3���3�X2jP                                                                    6���6ӳ�6��6{5n�4��3ɿ�2�o�                                                                    5�4��4�6�4)�63��%2���1�1 �                                                                    8*К8 �7�"7`��6���5�Ca5��4*|]                                                                    6C7�66��5�ܹ5�Tj4��4&�3.^k2B�E                                                                    4ʝ�4���4�_74�%}4�b�4@�Q3�J3��#                                                                    3���3�b�3���3�}�3\k�3�?2�r/2���                                                                    4��R4���4��'4���4��P4k��4W�3�+                                                                                                                                                                        '�      �     � 6J20nq                         6���3�ً@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��     @���    16:56:10        G  @���    @��    '�}      �     �5v�o7��\B'�fB'�f7�t�D"g�    (?�+�d�78r�?�  ?��    3�=q'��0�е��x�        .y*o�v>l+Z>���>0{�=�σ=�;�<x04��+2*�:Ԇ<%E<�=H]=�=��=�� =صl=��o=�D�>�r>Y�@א(@��=@#�?�t?���?�#�?��w?~�d@j��@j��@`�@_̸@_�o@_0�@^�M@^��@^;�@]�j@]��@]��@]�b<��_{@��                    E��m7%�G���F� �FPj�E��E z�D/S�CA��B@l                                                A�dC�@B}�A��OA<�@��Z?��>˄�=��                                                                    EFW@ˆF$.A��O3Jd�Fb�5B9U�A5'N.	5�                                                    {@��@���@���@�+8��V$�!��=/8@M&9K�?��(?舿C�㬡 1��!    >�y�@�52nv@-G�¤�R2        B��B��C�"�C�"�Ct��>�b?�,3n.PC�@|6���6�%A@��B_�>�y�B��A�y1B���A&CyAz�1B�JA�KUBb�4���f@?Gd"~{����f    ���f���S?W۔<�`�<��AG�@l�Aoq�?ؾ�?�le@�Lh?��fF��
3Jd�0�E". 1�3�]�7��v5_�]5�n�F���G *�Gj?9G}��>+��    =���:��~>GaS>�!x>��J>��v>��>�o�>�(>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�͞    ?��?��{@��7��39/*�9&�8�:�8S7nX�6���5�)5�                                                A�� ���G���G{@��?B�A���Bl��@u�*>&�S{@��{@��    C�|�BK��6�}�    6���6���@���3�Q>7�9    =��OC��M7���A�#DQ�C5V�B%&�A���A'��@�r+?��.>�q�                                                ?��AS|�@Lz�?Q��>��>u�t=��M=7< ��                                                                    D9WFF��9E�=�D��Cދ�C1�B��A���@��                                                @���B��5B(O�A�@�w@O�?f��>��=��(                                                                    C��E���EHOCͤyC+�B��A�boA�j@KT                                                @\�BMA>A��i@�e|@]�?{(�>�Bd>��=%>�                                                                    7!�M7�_A.�>kϰB��?D��?�@*��)n��,qau+ѻq*��)�},7?J,��U,DhK)y��7vN3���1|Lp-���.a��.`|�*���t�4��1�{h1zn�1x�~-�,g7t�1��!�t[4���    ,��3�:3��.0!�    < ��6�6��^3��1��S/��D                        /v!39O2�#o0���    3#b3<            4��;��,���?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B�fG���>,��A��                                                B�                                          B��                    A�                @�p�    .'* T�3�~�6BZ3L��3��&>�_�7�B?J�6�و;lJ�;m")                6<ʗ    45`U�    5�ד{@��/�wx3�ց    �)��{@��{@��{@��6���7�            6��    6��{@��    6��    5V�    {@��6)�&{@��6 �k5�t95WU;    5vp�{@�δ���6��7%/7%/<���    B��o6�țF@�:C��@H��?�BRz�            >���3��9$�3��<�=�<���?�\|?d)m>�[>i��>0=�2�� � ��Us�ޙ���]0��㢊� �p�+�T �=ry�+;(�5�>��H<�:�<�z
?��?b��>���>Mj�=��+=��� � ��Us�ޙ���]0��㢊� �p�+�T �=ry�+;(�5�<��R    (�4>5��16N�<'�m;���;v�\<��= �j<��                                                đZ�ľ�Mĸ��İ)>Ĥ$�ē2O�z���D."̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ,sp,sp=�4W.���.���@���@���@�Y�@�1�=�3-̱�-̱�>럺5ƈ�;��:�)V9}YE=�s�5ƈ�@hB    3r�fCC��Dَ�D�['D3��C�7C��BHu\A��@�]&                                                AÄ�CY��CF['B���B"7A���@�u\@�?]&                                                                    E���G<G��F���F8
�E��Dș@C�7�C �                                                C�EP�EBZ�E�D��UCԋ.C��B,��A:2�                                                                    E�F�^F���F{CZF]F.�.E�.\E��EG�%                                                D)��D�q�D�-	D���D��_D~h�D0��C�r�C�8                                                                    3#b7))GY�GN8�G�F�FzuFr�E��{EO��                                                ?��6��y                                                                        @N?�A_��A��A�[B-J�B]�aB���B�w�<�<�<�<�<�<�<�<�<�<�<�<�E���E���EcE��D�e�DV��C��C�U�                                                {@��{@��{@��{@��{@��D"!�B�E-�&ʧf"67��ZC�M�@��    @�l3@�l3{@�ν�����C�k�{@��C�y�C�wC�k�{@��@�#{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�;�D�6�G��D�6xD�%�D�%�A�A�Ff�\BIY�?��C��FK\FK\D��(D��(Fgy�BIT�                ?�
aC�C�C�(uC�6�?��C��C�y�C�y�C��@C�?C���C���C�6C��TC���C�/AC�e C���C���C��]C��4C�ʼC���C�vpC�7?C��>C���C�p�C�<2C�5�C�$!C�WC��w{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C�9]>��>��>���>�]@>��#>� �>�>�>̘>��>� �>Ձ$>��E>�y�>֨Z>�j>���>�Ԝ>ȗ�>�̈@��?��        A�q���p��m���<%����HH�G���Gz�Łď�Ħw�Ħ{�Ħw���&        >;A�_�    ���q���q{@��@QrAFVIM7{�Q7��        ;ay?3x�?5�?77?8�)?;D�??M�?EVA?N�s                                                �+F�6���A�  ?k�A	��@�1�'�}      �     �8���8[�e7ï�7,�63�5�M�4�� 3�/&                                                                    8 CC7��u6�.�6Z&�5�.!4���46�3	+@                                                                    G���F� �FPj�E��E z�D/S�CA��B@l                                                4�9_3��3MzH2�W�2�C1=L�0ZC/de                                                                    3�Hx3��2��1��1)<0o�/���.�	�                                                                    7��7��6���6�5v!�4���4*�3*�]                                                                    7z�6ɜ�60`5��5=Xb4�C3�ߓ3�                                                                    4�Yj4|�m3���3���32���1�0�S�                                                                    7A�66�ji6W��5��H5gl4�T4A�3!A�                                                                    8z��8w�7�F�7�6l�5�Sw4�,4                                                                    5G�%5\z4��{40�3��3862g�71��                                                                    8z��8w�7�F�7�6l�5�Sw4�,4                                                                    7��=7��7��6^�o5���5� 4LL�3]�}                                                                    4�u�4��4&��3�y�3'd2��1��0�c>                                                                    7��=7��7��6^�o5���5� 4LL�3]�}                                                                    �MF=�˔�5�.t5t��4�%�40OU3��2�c"                                                                    1	�B1�#:Fm�:D�9�n9��8���8�                                                                    5dH^5
|F4dP�3�4�3`��2�q�2 �16U                                                                    6�6�\c65t:[4���4�j398�2FGL                                                                    ���r�������!1<�����QM���N�23                                                                    ��ײ�����@���bG�s\-�/,�梫���                                                                    5�W75�8�5�i4p�3���3�27ρ1EIE                                                                    ��0'�ᡡ�R���Ģ�(>����"���(���*                                                                    �p�t�d b��ṲI�W����ධQ#K��!�                                                                    �;͢�6Ce��uH���>��CN�j��>I�XM�                                                                                                                                                                        *g��+�                                                                                            8AV�87%d7��<7A_6z��5��4��4��                                                                    4�y�4�+�4,�3�"33�2S8�1���0�3�                                                                    7=f 73j6�l�6#l5���4̊54
�3*�:                                                                    3�~O3�e�3)F�2��'2
��1]��0�t�/�4]                                                                    89u�8/��7��:7x6��5�1�5Fq49��                                                                    4���4��4%�W3��3<�2h&C1�9#0�4�                                                                    7��k7�jW7w�s7͹6U�@5� �4�P�3���                                                                    5�GV5��5�q�5}f4tk�3���2�)1�W�                                                                    6�c�6��m6H5�fW5,ҡ4u��3�]�2�y                                                                    4���4���4d�*3�3E�2�`�1��a0ģ�                                                                    7��g7�I7�D�7��6��}5��=4���4�                                                                    6_5���5��56��4�^E3�R�3��2��                                                                    4���4��4�-4rqU4E{4&�3�'�3��                                                                    3X�3^Te3X�3C��3��2��2�<�2O
                                                                    4��p4�#4���4�(�4q]�4-�O3侶3���                                                                                                                                                                        '�}      �     �6��[,��F�P�    =%�T            6�ƴ4#[�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @���    @��    16:56:13        G 3 @��    @�`    '��      �3     	�5}��7U�,B'�fB'�f7{��D"g�    '6�+��W7?�  ?��3��~3q��/
ǂ1�ڞ0(��1-�    -k>�)`��>w>P`�>1��=�/=��`<��%5) �2(�/:Fj<%*�<�?=H�=�)=��<=��{=ص =��:=�D�>�>Y�@�m�@V��?���?�1?�N�?b�#?<�I?�+@l0�@lM�@a�x@a��@a�\@aK�@a@`��@`de@`4@_�6@_��@��2<_��{@��                    E��7(�/G�� F��gFPBE���E eDD/-	CA��B@L�                                                A�bC��B|��Aڰ�A;�@��?���>�M�=���                                                                    EF@�@˹F$)SA��33��/Fc(,B9�N@��$-�׿                                                    {@��@I͸@I͸@��8�%:��"�� ��a�/,mJ�ist?T�'6�?��L-2�3�+-    >���?(�2�>�M7������=���    B-v*B-v*C��UC��UC�?�>�u/?$�\3-�uC�/6��.6�o>AA�xA�ş>���A��p@���A��:@���@M��A�k�@��oA�A���x�?�{"�����x�    ��x��զ�@Q�[>n��>n�@A.�?lbg@��?@�?��@>��?P�aF�i�3��/1z��.�4�3�ٖ7��3��5��ZF�t\F��GB�GaK�?H9?��@�,�?�K�>���>�R�>�G�>֎Y>�/�>�jP>���>��8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��$    ?.��?.��{@��7�U9. 9-8���8h@7t�<6�He6�o5(w�                                                A��(���׿���{@��AlA(��BA��@{��=�0{@��{@��    C�UBM�H6��t    6`�B6�s@�53 <�7�    >�YC���7�s�A��5DT��C;�_B;�A�ULA>�@���?�S�>�                                                ?��AW�7@Uƅ?p[E?
1>�*�=�P="b,<27�                                                                    D9,YF��E��gD��&C�+�C3�B�XhA��@�R�                                                @���B���B'��A��@�;�@��?lD�>��=���                                                                    C���E�͕EcQC�n�C.��B��A�ѱA�u@�x                                                @m�BM��A�3�@�U	@�?��">��> [�=)5                                                                    7'��7���A1�>o5[B�?GV�?���*_Q)� ,R��+�=�*_Q)3d+*֠x,4-�,��)QzT7{�1��/<�
,�Y            ��.H3>�X1�/>Ű/<�
,�Y7�.H����.A4��    '��81� �1�x�.DV    <�]��O�4�V�07�.���-MR                        *��R/� �/�1�.'v�    0Jm2�`            2��(:��d,���?z��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B��G�R>*��A��                                                B�                                          B��                    A�                @�p�    -a)�g1�!�6	Q+3�Hg3��>�\6�/>��5�Y�;V;�7�6�#_            5��    5B��4�9    7�W{@��1kƻ4R�5���|#�{@��{@��{@��7�X�8G�            7U��    7�+�{@��    7�+�    ����5�/�{@��7$��{@��7H��6�66��    4H�{@��5Vy74��b8q8q<���    B���6�7�FX��C��*?�A�>��Aqc            ?,��1��9�$�1��<�22<�\@c��@�<?��a>�	]>�=�=��� � ��Ux�ޙ���]4��㦊� �p�0�T �=r}�+;,�5�?Mx<�0:<��@b�@[??�f>�W3>�4=l��� � ��Ux�ޙ���]4��㦊� �p�0�T �=r}�+;,�5�=[o0�o*)ѧ�5|j�6/��<��y;��];u��=�ȩ=�=}Z;                                                ��/S�\��g�x�Z���F��)Z�� �×�z̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 3�tL3:\�>�t5���5=i�@:~(@4��@ƙ@�
�=�8�54�4�|#Dr�o7d�G<�wT@U�>a��=�&7d�GG�x 6��E7<�eCC-�D�r�D�=XD2�zC��,Ct�BG�A�d�@�3                                                A�-�CYr�CF=XB��zB!�,A�t�@��@d�?3                                                                    E�y�G�G�VF��zF7�yE�D�{�C��B���                                                C�EPqEBg�EpD���C�h�C�KB,njA: �                                                                    E�xF�_�F�� F{E3F]�F.�E�-�E��QEG�=                                                D)��D�tD�0 D���D���D~h�D0�jC�q�C�7r                                                                    0Jm7,��GYGNA�G�F�{Fy��Fn*E��2EO��                                                @��?!o%                                                                        @y�Az��A�3B�B;��Bq��B�߂BľV<�<�<�<�<�<�<�<�<�<�<�<�E���E���EcWE��D�\rDV޵C��C�T{                                                {@��{@��{@��{@��{@��C�9xB)�,Ð§fm;7��CCj�b@X�~    @�]�@�]�{@�ξ�v2��v2C��{@��C�	�C��C��{@��@�.�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��;�`�D�1oG��D�1FD�<BD�<BA9DA9DFe�=BGg4@�@�C�FKVeFKVeD���D���Ff�}BGb                ?LY�C���C�/.C��?l!C��,C�QC�QC�<�C�q1C���C��$C��C�O�C��C�ҽC�	�C�AmC�x C���C���C�!C��C�$�C��C���C���C��aC�JC�5�C�$pC��C�I{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�ƘC��]>�!�>�(4>��A>�=S>��>�`>���>�s>��u>�0e>���>��)>��>���>�po>�!N>ˢw>���>���>�{�@[�l@)F�        @����0s��)n��#�ć6���q5��RF��6�Į��=���=���=��|�+        =��}A��    =nz�=nz�{@��@��BFVyY7A�7�        ;lOa?=D?F�??H�??J�>?NdH?T��?_�q?s                                                �_@06ivA�  ?Fa@��4@A�'��      �3     	�8�v8c��7��57/�K6�(�5�X�4�|�3�                                                                    89�7���6���6^"y5�&14�T�43�3�A                                                                    G�� F��gFPBE���E eDD/-	CA��B@L�                                                4�^�3���3SZ2��f2�t1?t�0_�/p�q                                                                    3�43�2�F1�C1+P&0q�/��Y.���                                                                    7��7��6�s�6@�5y�4�Ԝ4�?34�                                                                    7�U6�Q�6M��5�T5Z�j4���4�23�L                                                                    4��X4���4�3�132�%2�}1���1��                                                                    7C/L7�S6{�^6	K5���4�S�4��3=%                                                                    8xrd8�(�7�(7�	6r�h5�	E4���4F                                                                    5F��5d��4���49�3���3%-F2q�'1�                                                                    8xrd8�(�7�(7�	6r�h5�	E4���4F                                                                    7�k!7�i7lD6f}l5�'e5��4S��3m׃                                                                    4��e4ӛ�4/e�3�\�3/�2���1�Z1=�                                                                    7�k!7�i7lD6f}l5�'e5��4S��3m׃                                                                    �AVﳽol5�15q�V4��}4,�3�.�2�I                                                                    4(�.��;@�:�3�:MC9��]9�68�7�                                                                    5el~5��4�{~4
��3�L�2րw2�y1.��                                                                    6ߠE6�6ԉ5z�04�R|4��3?D2S��                                                                    ����k����u�$'C��cس�~��<��
l�                                                                    ������j���@��6`�v~ɲ1P{��'ᱦ��                                                                    5�i�5�U5��4w�3O3Hd2=�R1R��                                                                    ������ǶWw���ʵ*5����E��:����                                                                    �oB�l�����L􂱵��3�U6��c                                                                    �:�-�=��(y��R!�����6E�B�*�d8�                                                                    (^sc)8$}/�ma..8��0�J�1B;0���                                                                    +Y�,��                                                                                            8@�8=��7���7��6}�i5ĩ�5 k�4#Ͻ                                                                    4�5�4��M40�3�ݞ3|�2U&$1�w�0��n                                                                    7<+�7:M6�T�6�5�7`4�g�4�33�,                                                                    3�@�3��3-9�2���22�1_�?0��D/�v                                                                     88A�86#�7��M7�16��5�%�54Z4Cl�                                                                    4�K�4�'y4)�>3���3�2jD[1���0�]�                                                                    7��L7�+97~G�74�6X|F5��l4��13��                                                                    5�35�ï5�M�5<e4wiu3��2��2 _�                                                                    6���6�p�6Mz�5�Hz5.�4x��3��2��L                                                                    4���4��J4j�93�	�3G��2��1�	�0�y{                                                                    7ॖ7�)7�d�7"��6�K�5��K4�L�4	J0                                                                    6 ^�6�5���5:�4�293�ʟ3�2�                                                                    4�*�4��e4�и4w�4H\4��3��v3�6�                                                                    3W7�3f�3^�X3G�q3!��2�2��32Z��                                                                    4��A4�n_4�p�4���4ty�4/�i3�8W3�B�                                                                                                                                                                        '��      �3     	�6R��'��8                        6��3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��    @�`    16:56:16        G R @�`    @�
@    (�      �R     `5:��7]�B'�fB'�f76L�D"g�    $�p�+�Cg6�>�?�  ?��5/ 2�~91��r2'L�/Q��        +6'&:ya�=���>E7{>ڳ=�J�=	۟6�x2-�:�1<%k�<���=H#0=��=��%=���=صI=��~=�E<>��>Z?ɣ�@�u?��?SI�?*P�?�g?�>�D@i�v@k$�@`�@a *@a	�@a�@a O@`�M@`��@`��@`_�@`C`{@��{@��{@��                    E��e6��HG���F�~FP�E��XE O,D/�CAcVB@�                                                A��C�|B|��AڊA;�m@��?��^>�=��                                                                    EFM�@��VF$4rA��j3�ګFc+�B9��@n0,�v                                                    {@����a���a�@|m�7݄t%$҇�β���/Z��%R�>�ܚ$�p���"�.��5��    >���<�DV2�=��������T9A��    BX��BX��C��.C��.C`�%?#�)>I��3WE&C��T6���6���AAf�@�I>���A��@�
Af@7��?���A��@gT�A�Z����O>�I��B����O    ���O��v>��?|/A?|/A@:�?c�@Y<�?��>�-@@e�?p9:F�1�3�ګ1�4 .�H�3��6��    5wtgFHg�FR�~F��(F�Ax?���:���Ba��@z�g?JbM>��N>�P�>���>ȶf>��>���>��^8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�u�    ��L���L�{@��7�^8�[P8��8��8	s7i��6���5���5#��                                                A��o�:S��:S�{@��B��@���B5z�@z��<�8�{@��{@��    C���BM��6B_�    6��)7\�@��3X�96��w    >���C��Y7�m�B�_DYx�CB��BIz�A�2`AM�@���?�L>�ȓ                                                ?NEA]pD@_A<?�Ui?�>��0=��I=+��<:�'                                                                    D9
F��E�D�bHC�cC5��B��~A��z@��5                                                @���B�z9B'q�Ab�@��f@
sA?p�)>�֯=���                                                                    C��>EϦ:E�QC���C0�B��6A�A��@�                                                @}�BN6�A�@���@	7�?��>�}�>"N=+U2                                                                    6���7w�A0�>nk�B�3?F��?�w)�C(v�+��+��)�C(L��)=�+�+o`x(��76L�                        �ɭ 3WE&1�($            7ɭ �5�ӷɬ�4r�                        >L����                                                                3 -'            3 -'{@��,�̯?"��?'c�?�  ?�  ?*?{;6?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A��G��>+ԶA��                                                B�                                          B��                    A�                @�p�    *��N&��    5��3��3��=)-�6��=��,4.��;Y;�;6΢]            �+��    1k�0�         {@��4M�c4B4�    �ğ{@��{@��{@��1�W7X|�            2z:�    6�{{@��    6�{    �h�n66%�{@��5�{�{@��7l��6�\y7�F    2�v,{@�ζ-�a1�L5� 5� <��    B�&�6��F~
�D�h�>��=�"@�nO            ?��]    :e�    <	3<&R\@�[�@_J�?Ď?!W>�*�>*m�� � ��Ux�ޙ���]4��㥊� �p�0�T �=r}�+;,�5�?r[<n<&LU@��S@^�?�-�?�@>L��=�h�� � ��Ux�ޙ���]4��㥊� �p�0�T �=r}�+;,�5�=��0�� $�4�~m5��|<�B�;��`;@��=���>)��=�rK                                                ��b��.��n��Z��@���~��)�^�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6��4�@�?P�7�?66
?�p?�o"?�P8@Q��=�{H7E�@6#7�FmD7p��>f�TBa�=LF�>h�]7p��I�&`4@y�7���CCD�_�D�#�D2c9C��WC��BG��A�!@��Y                                                A�CY_�CF#�B�c9B!�WA���@ǡ�@!?�Y                                                                    E�n�GqG�F��2F7�lE���D�]�C��B��K                                                C�EP�EBpE $D���C�F�C�~B,U�A:�                                                                    E��F�`VF�҉F{G	F]F.��E�,�E��yEG�	                                                D)��D�t�D�2<D��LD��&D~h�D0��C�p�C�6�                                                                        7u�GY�GNGqG�hF� Fy�|Fi�E���EO��                                                AI+o?�/�                                                                        @hcAw��A�[�B��B<+ Bs�^B�`UB�J�<�<�<�<�<�<�<�<�<�<�<�<�E��E���Eb��E��D�T8DV�HC���C�R�                                                {@��{@��{@��{@��{@��C�$BE��*w��eԻ7��IC:6^@8}�    ?@�?@�{@�ξ*�+�*�+C���{@��C��gC�|)C���{@��@}�Y{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G@�=ƽPD�(HG��D�(D�TzD�TzAiAiFd�BE1�Af�aC�0)FKPNFKPND��D��Ff�BE,�                >�6C�{nC�C�C�� ?<�C�]BC���C���C��C���C���C�%iC�Y�C���C���C��C�B�C�{'C���C���C�'fC�_�C��C��C��%C���C���C���C�X@C�7#C�$�C��C��G{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�N	C�M�>���>�r>�;�>�d�>��	>�a�>�c8>���>��n>��>��>�B>���>��6>�U�>�X>ġu>���>�,�>�l�@�L?�Z�        @2��ƿ�9ƿ��ƿ���m���������������q :ƿ�]ƿ�]ƿ�]�mlH        =L��A��g    ������{@��@l~�FV��7G7�        8�yW>��?D��?G�?J>?O�?V�{?d�?y��                                                ��5	��A�  ?F�?��p<#�
(�      �R     `8nb�8z7���7%�6t)�5��4���3�i�                                                                    7��272�"6�jC6Q�5�5E4װ�4 º3��                                                                    G���F�~FP�E��XE O,D/�CAcVB@�                                                3�WF3���3F6o2�<2 2�13L0V�/h��                                                                    3,2��h2z_�1�1!�P0b{/�1.��                                                                    6��s6�ea6�U6
�5k�4��4W]3.H                                                                    6�6�/|6PQK5�C5^m�4��\43�]                                                                    4/41J4��3�dZ36��2��1�L`1%Z                                                                    6��:6��6~�@6��5���4��T4 %3>V                                                                    7�!A8 f7�Z�7�6h~N5���4�4��                                                                    4�Pe5�*4�J�434�3���3��2j��1�I7                                                                    7�!A8 f7�Z�7�6h~N5���4�4��                                                                    70�17o	76�y6\e�5���5�4L��3h�                                                                    4.�4�UI4'É3��23*��2��u1��1��                                                                    70�17o	76�y6\e�5���5�4L��3h�                                                                    �2 5���5���5_x�4˿�4�3zj�2��L                                                                    �u�1�� ;K��:��!:_[u9��9T�8�ȼ                                                                    4ި#4���4�vP4w3��2�!�2h1.�a                                                                    6W756�&�6�5n�4�G47/38Qb2M��                                                                    �V@V���ӵ�6��㷴|,���&���`�ݢ                                                                    �6��Z����k����i$%�&1ı����3�                                                                    5U�5���5!Q4k�3��}3�27�1L�q                                                                    �_�:��t��I�㵼Ű�!E�|����[2��h                                                                    ��5�����9��Au����|�LL��n�                                                                    ������X���{�^B���$�UU�:��\��                                                                    (c��)6{�/�/c\`.��x0�V�1&��0�E�                                                                    &z��''%                                                                                            7��d7�"7��*7�6p�5�~4��4S                                                                    47��4n��4%�S3���3 ;m2G��1�bs0�W�                                                                    6��6�R�6���6�5xaM4�>�4�z3-�n                                                                    33߸3i�3"d�2��[2��1Qqz0��~/ĥ�                                                                    7�4�7��7�F�7܈6�\�5�]�5�	4<�]                                                                    40!�4e>4�3��3	& 2[T91�H�0���                                                                    7-�7g�47n�6�`�6L��5��4�3�A�                                                                    5F�_5�a�5���5�4i�H3���2�\�1�KX                                                                    6~6;56A(5�"@5%[�4h͒3���2���                                                                    4 �$4U�r4\��3�'%3<�82��1��I0ȤR                                                                    7T~�7��K7��7��6z�5��4��h4ļ                                                                    5r�5���5��u5/��4���3�5T3 ��2�D                                                                    3��V41�4���4i7�4=074��3��3���                                                                    2ˎ�3}�3Q]�3<u`3�%2���2�ě2Si�                                                                    4��4Y�4�Uo4���4g:�4$�3��3��                                                                                                                                                                        (�      �R     `6y�    F�%     =�            6�l�4��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�`    @�
@    16:56:18        