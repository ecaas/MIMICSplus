CDF      
      time       bnds      lndgrid       levsoi        levdcmp       cft       glc_nec    
   ltype      	   natpft        levlak     
   nvegwcs       string_length         levgrnd       hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      Sun Jan  9 16:23:29 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1938.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1938.nc
created on 12/10/21 16:55:19    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:29 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1938.nc had following "history" attribute:
created on 12/10/21 16:55:19
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�fF�. @�^     @�e�    '��      }�     �P4�2/6�FeB'�fB'�f7<��D"g�     z�+�'�6�'?�  ?��5)-�1��a20�0���0|`        &���"�ݎ9d?�3?S�8n>�,=�T=G�A8;g�2Dq_:,�k<+�n<��d=H �=� 	=��,=�ף=ض=���=�?�>�>U?���>���8G�=�5�=¤�>oZ�?"�^?kr$@Uq@WG�@NӨ@P!"@Ql�@R�?@T2�@U�@V�i@X�@X��@Ys�{@��{@��{@��                    E���60��G�A�F���FT��E�EG�D2��CD�gBBWV                                                A{	C�B�S�A�PAA �@��6?��_>Α�=���                                                                    ED�u@ɪ�F"�8A���3��nF`kRB2@2�2,�Y�                                                    {@������@gH77�o��L,ަ ��q��/IZ�����>�@ z��0$/��4Џ�    >�y�=��1�v�H� ��p��Kg<�u    B��B��C�:ZC�:ZC~��?C��>�{�3HWC�ɖ6yI6���A1XA"�>�y�A[u�@M��A�#@h��@��A@��@�"A�*���`�/)О�2��`    ��`���=N��?}��?}��@b��?B �@��B?:̍?9U�@�?��F��v3��n1�i�-��l3M��6])�    4�&OFc��F��F�G �:?�G?��C&�@��T??J�?I�4?)�$>�n>�P>���>��i>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���Mƙ�    ��ۏ��ۏ{@��6�37�μ7tg=7+O�7�u7G�@6���5��	4���                                                A�-�����{@��C^�A1�B:�l@dy�=3�{@��{@��    C��BO��6&0�    6N�7*�@��I3,6� �    >��ZC�|�7��B'ѰD�8C��fB�T�A��AWf2@�iS?�Φ?�w                                                ?Ak�A�G[@���?É_?�I>�4�>	��=9�Z<Oe�                                                                    DI	0F�8�E���D�w�C��IC@i�B���A�)�@��^                                                @���B�!BI��A2d&@���@��?|��>���=�V�                                                                    C��Eٰ�E"��C�GC?�B�f�A��A�8@�|                                                @5�Ba�ZA��@�I�@ǧ?��?>�X3>)��=4�$                                                                    6M�Z6ʳ�A/�>Q�B&�?.�?��1)���(I+��V+/߿)���(�^)��}+�}+���(��O7<��                        �s��3HW1��             7s�洁ː�s��3�Ё                        >�÷{`                                                                2�KG            2�KG{@��,�G?�  >�"�>L��?�  ?}��?aR�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A�ӼG�ˈ>*�A�M�                                                B�                                          B��                    A�                @�p�    &Q��"81    4���3��y3M�X=|��6HH>�4���;)��;,�k                ��l�    0�\�.���        {@��4`�H4';3    ��q�{@��{@��{@��2�{�7�            .�"c    .�"c{@��    .�"c    �i��6	+{@��4"D�{@��7�d7 7�&d    ���-{@�ε�{�3�24='g4='g<�\N    B��c6�y�FqĽDC-D>��B>��Ak�            ?��)    :2]    >�Χ=m�@��@��@Fy?�Ζ?��>l���5֋���p0��K^�ŊV���u����j鲊N�r�8�q�'&H�|�?{I�>�ƴ=lel@��@���?㔤?�L>DnX=��$�5֋���p0��K^�ŊV���u����j鲊N�r�8�q�'&H�|�>x#    $���8~oh9$�>1(�:��>��7?��>�\�>(Bi                                                ����a5e��[����@��#���.Ĵ.ė�8̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6��4�`�?�%�80� 5���@u(@�@�@�=�=��i7��+5��OG/�]8uL?%C%�=ߙv?l8uLJ�И4%�y8!8�CH��D�dQD��7D<�rC��<C	�mBN&�A�@�H                                                Aȅ�CZdQCF�7B��rB*�<A��m@�&�@?H                                                                    E�pG/�G�&F�հF:*�E��D�wVC��nC �5                                                C�w�EN��EAjfE��D�d�C�qPC?�B-�9A;�                                                                    E��F�3�F��,F{F]
-F.�ME�:�E��MEG�                                                D)�D�46D��oD���D���D~ppD0�C߁�C�B�                                                                        6jͿGX�GM�G��F�F[F{#rF��E�EPt                                                AOwA�X`A�`�                                                                    ?M�@�6�A���A��B%�BS�B|ʴB��[<�<�<�<�<�<�<�<�<�<�<�<�E��E�+�Ec�&E�%D�;�DWQ?C�DYC�i$                                                {@��{@��{@��{@��{@��C0�dB��&P��f��73VOCQ@Vp    @�@�{@�ξW�E�W�EC�Ԏ{@��C�˝C��bC�Ԏ{@��@huC{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G[>-DD�YKG{D�Y!D���D���AW�AW�F`�OB;�B[�_C���FK�4FK�4D���D���FbsXB;
�                ? 15C�|�C�p�C���?C"C��NC�~|C�~|C��C��xC��	C�ьC��RC�/C�h^C���C��C�%C�8�C�rCC���C��C�0C�f�C��C���C���C��C�S�C�'*C��C�C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��oD ��>�x�>���>�a�>��f>���>���>�/)>�>��>�g�>�5j>�?�>�u->���>���>�t�>��C>�6^>Ü�>���@(��@Q�        @o������5��jĪ���Ɩ*�Ǝ��ƊkĲ8���]���]���]�ĩ�@        =M��A��    ��b��b{@��@��~FT�6��7�Y        3��=��>c�?e?-N�?/��?2�?7�@?=�                                                �τ.5|�lA�  ?k�=�<#�
'��      }�     �P7��n6�(�6O��7��6N}�5|��4���3��j                                                                    6�c5���5�@64u15�ja4���3�o2���                                                                    G�A�F���FT��E�EG�D2��CD�gBBWV                                                3"N�2)�1�ǝ2��(1ؤ%1��0)�[/,��                                                                    2Mw1<j�1	��1�W[1�v0'��/Vr�.Z -                                                                    5�$�58�a5 }G5�#5Hy{4���3�s�3�                                                                    6I5Xx5(h5���5?̹4��m3�J�2�q                                                                    3��3
l�2�T�3�3�3�2v+E1�oG0ު�                                                                    6,��5�X5M�u5���5jk�4���4f�3m�                                                                    78�k6�~6J�L6�8�6I�f5�ym4��3�                                                                    4�3�$�3jM,4&^
3��E2�,2?�,1Q�:                                                                    78�k6�~6J�L6�8�6I�f5�ym4��3�                                                                    6�u�62S5�3 6G c5���4�84'�33s�                                                                    3��^3(��2�\�3�F3��2iޜ1���0ͮ#                                                                    6�u�62S5�3 6G c5���4�84'�33s�                                                                    ���>��4�f5>j4��H3��3C#K2�	                                                                    9��^8��;+ˆ;W�:�N:"f�9��9��                                                                    4I�3���3W%�3���3`�2�71�&v1��                                                                    5�-�5+\�4�o�5S�K4��3��O3�y2��                                                                    ��fn�T��dg�4 �Q�ʳ�k��M\��г                                                                    �E4۰�ı+�زpX��?F���5��좱l��                                                                    4�C�4*`3�4R j3���2�12�1.$                                                                    ��Iݵ_״�龵����
oH�= ����@��Z                                                                    �߂�����ht��+���������$.��E~�                                                                    ��<]�g]2��eZ�XP�|��p4���#�                                                                                                                                                                        !��O'2�}                                                                                            6蝱6j`�6:�7�N6NE?5��g4��3��m                                                                    3k082��{2��o3��K2�c�2m	1`0���                                                                    5��L5e�_56Ѥ6��5Uq04���3�܏3n�                                                                    2feg1�$^1��22��y1�=1�0u�/���                                                                    6�"�6`��63+7Ʒ6\�"5��?4��4��                                                                    3a��2�PB2���3��2��2$;1���0��b                                                                    6[�k5㿾6W�6��6*^25Ss4��3���                                                                    4{#�4$l4�4�[4B��3q>:2�[�1�t6                                                                    51��4�
"4�	�5��5	��4*�o3i��2��                                                                    3J�2�T�2���3Ħ�3V�2B�1��L0�>�                                                                    6�JP6.6�7 �6P:=5��4��83�1�                                                                    4�y�444�\5��4m�k3�m$2��1�8�                                                                    3 82�L3Q4C	4753ß�3���3@P                                                                    2P�1��t1�(s3��2�څ2�U2h,2g�                                                                    3C�2�@y3*F�4n`i4=�A3�Z3�yb3k�                                                                                                                                                                        '��      }�     �P6ADP        >�xL    G��>��=$�6���3�p@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�^     @�e�    16:55:19        F�f @�e�    @�l�    '�M      }�     ��4��6�_B'�fB'�f7EƯD"g�    !�F+���6bl�?u�?ux�5[~t2q�E2�^2?t1.&O1�(�1FŨ'ŝ�#�%�:��7�ڳ=i��>,[=�Hq=UqE8W�2M�:5G�<+��<��=H(�=�"�=��=���=ض�=��=�@$>��>U\>�q8�3[<��<�ux7�Gy6k�=�/{?H��@R	@W=@N�(@P �@Ql�@R�@T2�@U��@V�?@X�@X�@YsY{@��{@��{@��                    E�_6�`�G��tF���FUq�E��EV�D3 CD��BB��                                                A�C�4B��A��RAA0�@��?��^>Χ�=��                                                                    ED��@ɏF"��A���3U^�F`-XB1��?��-�&P                                                    {@������@d�7��l�%���kg�#��/H���N�x>��!�F�id/��'4�z    >�y�?5�i1�����<�>٩��D7K;��@    B~�BB~�BC�k�C�k�Cm??Pj?B1�3K+�C�� 6 �Z6��gA0�B�C>�y�B5"�Ad��B{��A�[A"��B5��AG��B<�r��Zd��y�xI���Zd    ��Zd���@���?~�[?~�[A""�@Ae�AK�;?�X{@$��A;`~@7��F�]3U^�1O -���3}�7�_    5��F޵�F�G?�GN�u?	�@��C�9�@�6d?l0?O��>�v)>�cb>� �>��W>��>��8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���MŸ�    �>k��>k�{@��7l)8�P27��7�#7�`�7KI6c '5�s4�s                                                A�#���l���l�{@��C���Az8�BMDN@]9=��`{@��{@��    C�H9BO��5� s    6}'�7��@�]3W�!6�J    >b�
C�q�7��B2wD�:�C��.B��hA�v�A_�@�Y�?��_?-]                                                ?M�A�J�@��?�R�?$'�>�҈>�\=@��<V�                                                                    DJ|%F�ڸE��D��C�CB+B�T�A���@��>                                                @�SoB��iBNkA;��@�L!@6�?�}x>�t=���                                                                    C�uE�g�E%V�C�ޔC@7B��A��{AD@�%                                                @8xzBbFA�_�@���@~�?�]'>�t�>,n�=88�                                                                    6��7�A�>LbA��?*P7?��*<�(���,
c+h�*<�(�9+&�+�&�+�e)	�07EƯ                        �:Vx3K+�1���            7:Vx���n�:Vx4�B                        >� ��3�>                                                                3 ��            3 ��{@��,�6�?NN>L��>�.S?�  ?|{�?7:E?^$z?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B���G���>,��A�8                                                 B�                                          B��                    A�                @�p�    '���#]k�    5@�|3_�3��>*!M6��2>�AV5Ǻ2;Fm;Km                �ݝ3    5�s1��        {@��4�^4��I5\���ڭ{@��{@��{@��23D7;<�            5]�L    5]�L{@��    5]�L    ��w7M�{@��7���{@��75��6�Q�6�͟    ���{@�εôF4�2m6O��6O��<���    B�}z6���F@[�C�=n?�i}?I�B&�            ?��w    :Lq]    ;ӱ�:��@���@�=l@%�>?m��?+�>����5֋���p0��K^�ŊV���u����j鲊N�r�8�q�'&H�|�?��;ӯ�:ɾ�@���@�9�@%�3?LQ�>U3=���5֋���p0��K^�ŊV���u����j鲊N�r�8�q�'&H�|�>'�    %j �4g=j3D;Z9��9�
�6+8�>p3>쎎>T8/                                                �=��C��ă����R��χ�Ľt�Ħ�1ċ��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7�C5wS>���8���7»A]�@�d�Ag:A4�N>�%8<�6\��G�=75p�?S��C��?Q��?T�#75p�J�4���7��	CH�D��D��FD=:C���C	�BNWGA���@���                                                AȹCZ�CF�FB�:B*��A��@�WG@��?��                                                                    E�<G �G�GF���F:"IE��D�swC�ߢC ��                                                C�iEN�EAF9E��D�^�C�d+C<�B-��A;�                                                                    E��F�13F�� F{TF]
�F.�E�;�E��eEG�y                                                D)D�0JD��]D��&D���D~p�D0��C߂C�CA                                                                        6�ĨGX�GMk�G�0F�BF{�F��E�TEP�                                                AgU�A��@��[                                                                    @+�D@��A��NA��:B&��BUj-B��B���<�<�<�<�<�<�<�<�<�<�<�<�E� YE��Ed
9E��D�9DWQ�C�C.C�is                                                {@��{@��{@��{@��{@��BA�vB��'%���f��7RH{B�0k?�ޜ    @���@���{@�ξu��u�C�/�{@��C���C���C�/�{@��@h �{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G �f>�SlD�T�GhD�T�D�/D�/A�WA�WF_�B9E�B&3C���FK�gFK�gD���D���FaɹB9@�                ?���C���C�t�C�[?pRC��+C���C���C��1C���C��C��~C��C��C�HeC�sbC��4C��*C��C�$�C�Z9C��9C�ѵC�C�<�C�bC�xlC��RC�Z?C�*<C��C�:C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��D�o>�f�>�j�>��T>��>�t]>�)�>�%>�x�>���>��\>���>�x�>�"�>�W>�n>��>�B>�4&>�t�>��1@8bc@,�        @��v��:���2���,oĩ<W�����������ґ�L�p�L�p�L�pĘݱ        =v�bA�7%    �[2��[2�{@��@�؃FT�g6��7��        7N�V?�r>�f�?/C�?0U?2b?5�"?:q?Ao                                                ��M36�"&A�  ?k�;���<#�
'�M      }�     ��8E�6���7�t7y	6P��5P@�4���3���                                                                    7y��5�*6:�67�O5��4��M3�&�2�wq                                                                    G��tF���FUq�E��EV�D3 CD��BB��                                                3�X2,�H2�B�2��y1�#�0�{ 0ǋ/-�`                                                                    3�?1ZP�1��1��~1
g�0	��/;�.[�D                                                                    6hr5T��5�~&5��P5J�54^%Y3���3 �                                                                    6�gn5hs!6��5�6G5I�p4��r3�n3t�                                                                    4+3ͤ3ʖ�3�`�3&ޠ2^`�1���0���                                                                    6�a�5�x6+�>6�d5vq44�q63�#�3 ��                                                                    7�[6�K�7��7�}6M�
5f:%4�x3�Ǚ                                                                    4��3��H40/;4+)3�F
2��2)ܼ1V��                                                                    7�[6�K�7��7�}6M�
5f:%4�x3�Ǚ                                                                    7�6�6g�6K��5�x�4�Y�4 �36�                                                                    4�L3'U�3���3��(3!�2F�)1�d|0��                                                                    7�6�6g�6K��5�x�4�Y�4 �36�                                                                    �	���8h4�5?0
4�LE3�I�3(�j2��?                                                                    N�/�K;1T�;N�:�E:L�9� 	9��                                                                    4�A�3�O43�t4�3lf@2��1݈�1��                                                                    6+�Z5'�[5}#�5X<�4���3�@d3 �2                                                                     �,�e�!�g��[�����T7��U�ʲ�*ȱ��                                                                    ���Q��.��+�t���Ad<��9������n��                                                                    5*��4'@l4}s�4V��3��W2��!2�1 �a                                                                    �6�õ&{����u��➵�����}ϔ���                                                                     ���s��Rв%�>�.���s�����T�G                                                                    ��lĲ�[b�:�o�\*S��񛱰N�� ���$�                                                                                                                                                                        '!�'��                                                                                            7�n�6�"u7�u7�k6P�45b�`4�*t3��                                                                    4*3�3�r�3�Y-2��1��Z1D'0�і                                                                    6�h�5�a�6D�6^5W�4n-Q3���3x�                                                                    3U2�R2���2�[�1���1F0WV/��                                                                    7�bi6���6�L7R6_?�5ykB4Εt4
�                                                                    4�3�3��D3�^�2��2)�1j�0��                                                                    736E�6��
6ع�6,5�5-��4}s�3��^                                                                    5 :]4�4��4���4D��3FX�2��X1���                                                                    5��4�(�5�nj5�!�5(�4>�3L�12���                                                                    4z'2�w�3���3�&z3
02 G�1j\0�.�                                                                    7+Z�6 q�6���7qv6Rz�5T4��3�o                                                                    5CՎ47]�5 I}5]4p�53rl�2��1�o                                                                    3�8_2�4�3��U4F��4�3���3{�g3A��                                                                    2��1�f2��53 p 2���2���2KD�2f�                                                                    3��;2�$3�4r��4?�n3ę�3���3l�;                                                                                                                                                                        '�M      }�     ��6pr�    F��>?^��= �F��o    =
�6�\A4^_@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�e�    @�l�    16:55:21        F�� @�l�    @�t�    '��      }�     �`5�t6_�B'�fB'�f7q�D"g�    #�Ϯ+�$6)y�?�  ?��5��02�;s2��2i�Y27e�2!�|    )���%��c:DƊ;�'�>2�>a*=�]�=/o�8^Te2T��:;��<+��<��=H1!=�%�=��=��C=ظ=��q=�@�>�S>U�>�u�>�
?>]w�=���=�v�=�BW=�x�?)�d@M�@T��@L�}@NsF@P"k@Q�@S�
@U6�@V�d@W�@X�v@YT�A��#;ց�{@��                    E�_�6�bTG�-�F��QFU��E�+�Ef�D3Q�CEf�BB�i                                                AL C}�B�`�A�2:AA?�@��?�!�>�3=�W�                                                                    EDj�@�n�F"�8A�ƍ3f@F_�B0�>��n.�H�                                                    {@��Az�ZAz�Z@_�d7��%X`��'�></(�0&YZ>��=#�Ϯ@l��/�`�5?R    >�y�A�1�E@N�kA�>w�9Rk=��    B|�B|�C�WcC�WcC���?V�@��3��C�h=5�(�6��A,�,B���>�y�B���B?%C7OAM��A��%B���A�I�B������4��u �^����4    ���4���qB�I?d?dAmo2@���A� �@"�@JQ�AI��@/MF�{�3f@1@C�.kU3��Y7�	5�>5=DGc�G)�G�{�G���=�n�A��Cn
A��?pŝ?��>�K�>�%>���>�� >���>��^8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�Z}    @Y�@Y�{@��7L]O8���8P@�8��37��7Jb86q��5���4���                                                A��h���h��{@��C�8A��B�]d@R#O>\�{@��{@��    C�?7BS�5�ŭ    6M417d9@���3-i7�    =u4�C��7��TB<�D�SVC�"B�ǦA���Ao�S@� �@??��                                                ?Z�A�'{@�%~?߿H?/!q>�h>*OO=b %<f$�                                                                    DK��F��E���D�JUC���CD�PB�>�A�g�@�ap                                                @�)B���BS�A?N@�+B@�P?�Ƶ>�w=���                                                                    C�ρEن�E(VC�Z CB��B��BA��JAa@#~]                                                @<+�BdS?Aō�@���@^p?�X?��>5��=>_r                                                                    6ɪ�7F��A�>F�PA�?%w�?�x�*��)�2,y;a+�l*��)��,`H,��6,V6)�! 7g��3�+�1s0�.w_R            5�F%3��1׊{1w1s0�.w_R��F%��d�5�F&4-.6    'R3v�J3sp0r�u    =��67e�6��1MՄ0�P.Zo�                        *�� 1$��1�	/V.z    1M�t2�p�            3�Ǻ: -m>�Y>���?�:?{?uZ?<SB?	�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C��G���>*��A�x                                                B�                                          B��                    A�                @�p�    )f�r%GD�3���5��#3z�3���>��7��i?��6��};f��;c�g                6��q    5��3��|    8�1{@��1{�6+t�8�R����{@��{@��{@��7S��6�J            8��Q    8��Q{@��    8��Q    5F3�7iD	{@��8���{@��5���5��5H    5�.P{@��5�R�6j�N6j�6j�<�
7    B���6�=kF��C��@���?�/B�s            ?�Y3�+�:ps3�+�;�P�;���@���@�B3@Q
�?f?*K�>�л�5֋���p0��K^�ŊV���u����j鲊N�r�8�q�'&H�|�?�fC;�P;���@���@�,�@P��?~U~>�&�=�̋5֋���p0��K^�ŊV���u����j鲊N�r�8�q�'&H�|�>ȩ    (�c�3�u~4���;[;*O:��(;1��>�pK>�                                                �����d���֎��ѯY��j�ļ�ĩbď�j̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7&�N6"�<�P�9?��8�:gA���@�U@A��A�ٰ=�tI8dj37\��G��5���>��[CW`@��O>���5���K0�8�Ӄ7;CIAD��zD��D>�C�:C
;�BOSA�&&@�                                                 A�ACZ�zCF�B��B+:A�;�@�S@&&?                                                                     E���G�G��F��7F:.E��D�r�C��C ��                                                C�V�EN�lEA%�E��D�X�C�X-C<9B-�1A;/                                                                    E�fF�-�F��GF{ F]
�F.��E�=E��TEG��                                                D)~�D�+D���D���D��D~qD0�C߃bC�C�                                                                    1M�t6��GW��GMR�G�F�>�F{�F��E�zEPf                                                ANmA#�{                                                                        @���AZ&sA��-A�c.B(y�BVb�B��B�qX<�<�<�<�<�<�<�<�<�<�<�<�E��"E�\Ed�E��D�7�DWS�C�D�C�i�                                                {@��{@��{@��{@��{@��B6��Bd�m)�fW�7��B��?��H    @쐷@쐷{@�ξ�ܾ��C�о{@��C�L�C�V�C�о{@��@_��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G �`>��7D�RXGS{D�R/D�z�D�z�A!.�A!.�F_h�B7\�A�"�C�u:FK��FK��D���D���Fa%%B7W�                ?�?C��C��#C���?�C�DC��3C��3C���C��lC��JC�دC��ZC�/C�<<C�_�C��CC��<C��vC��dC�$.C�X�C���C�ȚC��IC�'�C�H_C�d�C�Y�C�-'C�EC�`C��d{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C���>�w>���>���>���>�ۮ>�G�>��>���>���>�[>�=K>�Qn>���>��>�>�-�>�_v>�b>��7>��@d�@+h�        A2���
-��
	`�	�����oe~�n�9�nu�Ŧ���e���{���dİjh        >`��B	"/    >V.>V.{@��@�U5FTQ�6gbg7��        8p�?W�a?1��?2T!?3�?4X�?6�+?:�?@�-                                                ��697���A�  ?k�?�(<#�
'��      }�     �`8FA�7�:7�h`7J'6M��5Y`4HI3�                                                                    7zm�6�V�6؃�6;P15��f4�3|�2�o�                                                                    G�-�F��QFU��E�+�Ef�D3Q�CEf�BB�i                                                3�
�3��33ܥ2��)1���0��/��Z/-!�                                                                    3e	25�2c1�1ċ�1H�0թ/��.Z�S                                                                    6h860��6��I5��E5G�4g�3~��3�                                                                    6��U6F 76�5�b�5TL�4�-�3�Gi3�                                                                    44I�3��54c:�3���30T�2��1���0�-1                                                                    6�W�6r C6��6-�5��'4���3��3*��                                                                    7�V�7�O�7���7P�6Mu&5v�4~�3ث                                                                    4���4��/4��42�H3���2ސ�1���1\�[                                                                    7�V�7�O�7���7P�6Mu&5v�4~�3ث                                                                    7�G6�#�7�6R�{5��)4͖U3��53:�                                                                    4��4�4C�m3�G�3��2YC�1x�0تF                                                                    7�G6�#�7�6R�{5��)4͖U3��53:�                                                                    �"#��H�5c��5? �4��I3���2඲2�o�                                                                    .��1[�;VO-;t� :���:�9���98�                                                                    4�S�4�Z�4�j�4 3xD�2��E1���10�                                                                    6+�C6}�6��5]�)4��3�b~2��;2!�                                                                    �,vz����
ȵ��P�F�^�/�[�����4                                                                    ��]�βL��$�yV��>U���=ޱZ)�m                                                                    5+5�53�5/�4\��3��>2�gq1�S1"P�                                                                    �6�ٶ
�>����@ߵ
1�"i�2�=��@r                                                                    ���#���������1�����z����ʬr�F��                                                                    ��+��]ѳ���`S責#ұ�nİ�4"�#s                                                                                                                                                                        *X �+5�                                                                                            7�Q�7`%17��7	��6M�v5l�84t%3�"�                                                                    4��3⟽4�@3�^2� 2 b�1
F�0�p�                                                                    6�L 6[��6��p6��5U.4x��3���3�                                                                    3�u2��3u2�+1㦨1��0��/��'                                                                    7�Fr7W�7��?7��6\9�5�14��\4�o                                                                    4�93�cv4C�3��2�L�2�1$�70���                                                                    7�6٦f7T��6�ҙ6)sC54�42_Q3��                                                                    5�,4��+5s�4�^e4A�M3N\�2Kڦ1�~�                                                                    5�.D5���6+��5�qJ5��4�3#�2���                                                                    4?3�4Dh�3��03}�2&��1$�0�G                                                                    7+�7"7���7�z6O5\�C4Z�3�=                                                                    5C|5p5���59�4l�A3|8M2y'�1�E�                                                                    3��[3���4d��4JT�4s�3�\,31�3@b�                                                                    2�ԧ2���39�3#�2���2�=�2�2v�                                                                    3�Np3� �4��4wKb4<�Q3̍3X`3k#�                                                                                                                                                                        '��      }�     �`6@�Y'RBiϰ    :�w�            6��3�k)@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�l�    @�t�    16:55:23        F�� @�t�    @�|     '�      }�     � 5��W6��QB'�fB'�f7��D"g�    (C+�p�6�Z�?�  ?z�50�М3���)h̦1��Y��?0A�    .K��*>��>o�>K�L>�z=��=>�;J�l3�:�2M�V:@h�<,<�,�=H8v=�(�=��=��x=ظ�=��*=�A�>�>U�@��@�-@n�?�U7?���?}�?O��?��@VQ�@VU�@M	�@M9Q@M�@@N-@N�a@O��@P�+@Q�P@R�7@S)AgK�<�Yf{@��                    E�`�7=�LG�7F���FU��E�eE^D3Z�CE|,BB�K                                                AM�C��B�o�A��AA);@��?�,M>�,�=�b�                                                                    ED?i@�T�F"t�A��!4M{F_ccB/n�A-HI.�.�                                                    {@��BV�BV�@]J%8+�����������)@��k�_�j?_~(C@l(�����2	�    >�k�A�"1���@�@��@��D        B��PB��PC�vRC�vRCu�}>�
:@��3>C�k6>~�6��A+��C&�>�k�C9ˌB��Cd��A��YBmC>DB~�B���A�E��#6��A�E    A�EA�`'@J�=��-=�sbA�{~A*�sA�+s@DO"@M�A?!@�F�k4M{1��/M4|�x8r��6N�	5���G�rGE�G���G���> ��    ?���?��^>`�m>�J�>�H>��X>ŀ>�^�>�xC>�t18��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��=    @�4�@�4�{@��7��9S%�9>�8�4�8�7�	�6ʶ�6��5%�                                                A�?�8?�8{@��@*�2B'��B�O?@F�M>��{@��{@��    C�AABNt%5��    6Er7��@��36�7�    =9C��7�~�B5�<D���C��DB���A�GqAd�@Ж�@	E�?�                                                ?Q׼A���@�\?�ы?%a>�m>!x=Y�D<W�                                                                    DKAF�a�Eݡ�D��KC���CD�kB�1�A���@���                                                @���B��KBS!�A:��@��@��?��>�'=��<                                                                    C�T�E�:E'�{C���CA�MB���A��A\�@"�                                                @;�Bc�6A�^@�Β@	5?���?��>6�+==@�                                                                    7Cx�7��A�q>A�MA�?!�r?��+��r*6��,��),.�g+��r*s�%-[<-��,���)Ϛ	7{�U4�:�2���0}Q0n�0^�(.t��7��34â�2)�2�i�2��Y0nT���31��7��M4p��    +,T�4��4���2|c    =@x�8&�/7갹2�"�1��]0�=k                        -��2��2�%�1�v    2�B�2��;            4�ň:���->y?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C`��G�V>,Q~A��U                                                B�                                          B��                    A�                @�p�    .et)�J�4��W6,�A4MOu4|�d?z�8�@ �7�׽;Xr
;Q\X6l�            7��    3���6>Q&    72��{@��    4��    �WM0{@��{@��{@��6�h>7n�            6lpF    6�<�{@��    6�<�    6JG5I5{@��7m�{@��6h�T5�B5��f    6K��{@��5�h7:C�6��u6��u<���    B�O�6�d2F uKC��A���?p�C��            ?�Jo4�:�:LϠ4�:�<��W<�j�@�+@�[�@+d7?��P?+�>����5֋���p0��K_�ŊV���u����j鲊N�r�8�q�'&I�|�?�}�<��<�d�@��@�-z@[�?U�>{C=����5֋���p0��K_�ŊV���u����j鲊N�r�8�q�'&I�|�>Vf1-v(�t�5�v�6<H0=,D=/�>P��>6��>��w>�W�                                                �^��ć7Ă��q��V{�2,��rÚj�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 5��5��=�p�8e1�8?��A���A�fA�B�A�L�=Ej6��6��dC�G�6��;���?��=`<ݦ^6��J ��7�:6���CI�4D��D��D?��C�&gC
�@BP$UA��L@���                                                Aɟ4CZ�CG�B���B+&gA��@@�$U@�L?��                                                                    E��yG.G�$F��eF:rE��MD�gCC��BC ��                                                C�O;EN�BEA&E�JD�N$C�D�C3�B-�A;Z                                                                    E��F�.0F���F{�F]aF.��E�==E��tEG��                                                D)"D�+�D��TD��!D��uD~q.D0�C߃�C�Cj                                                                    2�B�7;��GW�GMT�G��F�8F{�F��E��EP8                                                ?6�b=#�v                                                                        @gZzAsj�A�Y�B��B9(Bo6�B�=�B���<�<�<�<�<�<�<�<�<�<�<�<�E���E�TEd�E��D�4?DWS�C�EwC�i�                                                {@��{@��{@��{@��{@��D�BB-�Ё�f��7�B�H�?��    A	A	{@�ξ���C���{@��C��=C��#C���{@��@]J�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G �|=�XD�G�GH�D�G�D��D��A�TA�TF_��B7TZ?A�C��fFK�!FK�!D���D���Fa�B7O.                @KYC�i�C�]C��?�C���C��VC��VC�v�C�PC�0C��C��8C��C���C�ŐC��JC���C�ܗC��IC��C�6�C�eC��C�� C���C�pC�EC�S�C�/�C��C��C�r {@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��4C�l�>���>�3�>�S>�K�>���>��7>�t�>�v\>�S>��>�y>�?Z>�g�>��>�$r>��>��o>�d�>��>�G�@G�@��        A�4��`Ҿ�`@��_�$Ŕ����z��緕����4����pm��p���pmċ�        >�p�B6�t    ?���?���{@��@���FT 6�"r7q�        ;��P??ͣ?BX#?D#�?F�?K2?R�?]P�?r
�                                                ¦�7���A�  ?Vi�@�ȴ@���'�      }�     � 9�8|��7�@�79�o6��5���4ΰ�3��                                                                    8C܉7��H7��6j{5��O4��4��3��                                                                    G�7F���FU��E�eE^D3Z�CE|,BB�K                                                4���4�s3b�b2���2<�1?�0X�]/cu�                                                                    3͈[3'��2�Q�1�"11$�0rh�/��.���                                                                    75�J7"��6���6��5�.�4��4��3+Ն                                                                    7j�b7.��6�dM6̏5�BU4� 4#հ3,�!                                                                    5�
4�%4w�3���3YwO2��2�n1.d                                                                    7��7Ua�6�z�6%�5�m�51>4H>3R��                                                                    8��^8���7ל�7&G�6�MS5��5<4��                                                                    5ya�5�O�4�K�4]��3��N3;\32���1��                                                                    8��^8���7ל�7&G�6�MS5��5<4��                                                                    7�N�7ͅ�7'��6�.Z5��5-	�4`�3s6�                                                                    4���4�Q�4q�3ԩ^3LPW2��
2 I�1tX                                                                    7�N�7ͅ�7'��6�.Z5��5-	�4`�3s6�                                                                    �H�j��}�5���5s<4�7�4�V3k��2�X                                                                    0�q1+��;9[7;6Wd:��:��9�R=9Hޜ                                                                    5�]5m�4�B�4&853��3
r827*�1@��                                                                    7ME6�/|67�e5��4��4 A�3E�2SZ�                                                                    ��X��2���l��*Zմ���ox���&��                                                                    �Īk��{ɲ�$ಜ#�wu�-���Pܱ���                                                                    6Í5�d57��4���3��
3 >�2E��1TO`                                                                    � ����6�r����#h�4(n�������!��	W                                                                    ���ϳ��2��/%�^���}A�.`�R;���?                                                                    �^~��K�o�^���nڲ߸߲�v�;��V8p                                                                    (s�~(�
7/��s/�!|0�[Z0��y1���16D�                                                                    *�`�+��                                                                                            8h'\8Nǁ7��&7,�6�8/5Ȭ�4�?�4��                                                                    4긔4��4Ge3��]3h2Y1�n�0��@                                                                    7cl`7J��6��^6+PB5���4ҝ�4
ן3,�a                                                                    3��33��V3B�C2���2dv1dEO0�E�/���                                                                    8^�d8FZB7�ї7*~6��x5܎u5v4<                                                                    4�'�4ȋ�4>�!3�M)3`�2o�1��0�U                                                                    7��7ȘR7�=,7
=6\9�5�h4�0�3�ߎ                                                                    5�Nr5�@]5�j{5��4{��3���2Ҁ�1��~                                                                    6���6��6X�5�j�51��4u�$3��2�g                                                                    4�D�4�@�4w��3�U/3Kb?2�w�1�W0¾�                                                                    8�+7�+�7��7(�^6��65��*4��4 ��                                                                    6�65��%5Ak4��3�t�3 ��2F�                                                                    4���4���4��<4}he4H�.4�3�֙3|�[                                                                    3��3yЂ3i��3L�#3"DY2���2���2L4�                                                                    4��4��4��4��>4um�4,Wr3�x3�n7                                                                                                                                                                        '�      }�     � 69*�+,T�            G�ȑ    =�VC6���3�-+@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�t�    @�|     16:55:26        F� @�|     @߃�    '�y      ~     ��5�%f7J�B'�fB'�f7�<zD"g�    )�6+~�7�?�  ?s��    3��\)*�/������        /Nv+CD;>�h�>M�=�T�=�4<�,18Ϙ�1R�_2\�':>��<+ϋ<�4d=H<�=�*=���=��g=ع�=���=�B$>��>VB@ܡ@�]�@)�<?���?�o�?�2�?��?�m�@]�@]�8@S�@SQ/@R��@R�@R�@Q�l@Qbg@Q!�@P��@P�L@Y�f=&C\?���                    E� �7E�YG���F�6�FU-lE��`E7D3&�CEA�BB�                                                ANC6�B�/4A߳�A@�g@���?���>���=�N�                                                                    ED2@�W�F"jA��{3�l�F_IYB/�XA@�v/*��                                                    {@��B~nkB~nk@]|�8J�$���    �h��/J�&&=u?   )�6A 8���@�2.b4    >�y�B0o1�AA7Ġ@��,@�0�        B��B��C��C��C��>��N@�v�3�&Cפl6��_6�/�A+�CG�+>�y�C[��B��Cp��A���B��iCM��B��B��;Bb�?�V�    Bb�    Bb�B]~;�3��/5���A��DA@��A��@ib|@!1h@�O�?Ӧ�F���3�l�1-�F.C��3�0�8�c�6��G5�-�G��GW�_G���G�k=t�>    5j��    >V�>��X>å�>���>��>���>��P>��v8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�!    A�EA�E{@��7��9`�9F.�8�ڷ8  �7�ܧ6��6>J5��                                                A��@�`�@�`�{@��=�)Bxd�B��}@>��>�#�{@��{@��    C���BK��6v�    6D7�7b�@��!3s,79�    =s]�C��e7�y�B)��D�M1C�?B}�Aٺ�ARJ�@�#?���>���                                                ?C�A�z@���?�|�?�U>�z>`t=C�<H��                                                                    DI�F�!rE�1(D��/C�)�CBϊB�W�A�)�@���                                                @���B��BOP�A4��@��@@*?�!_>�_R=��                                                                    C��3E�iE%�C��C?c B��A��DA�e@ ֎                                                @8uBa��A�H@��A@�?�?0?�c>2��=:��                                                                    7Hu7ś�A�>A��A鈝?!I�?���,��*���-n,�B�,��*�K-�Y^-�sr-�*#W]7�X5A�2��1��0���0qE�/_�8
�5�K22 e2�q�2��l1���
�2.b48
�4� �    1���5۽4�~\3��v    <*��8l��8G�5j�4 �2�                        3���5CM5!I4�)2��6	�2��            5l�=
-).?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CpmG��3>+��A���                                                B�                                          B��                    A�                @�p�    /�*�n�50�S63�3���3���?��_8G�D@[(�8(�;�;�;�                7�X�        6���    1��{@��    3��    �3a{@��{@��{@��6ͬ�7�Y�            6���    6���{@��    6���    6�    {@��1��{@��1��1��0��    6�{@��6W�7��7�b7�b<�n(0y�Bp6���E��C�<�A��v    C/�            >�`5`�9|� 5`�<���<�R�?�I�?��M?3��>�H�>���>W+�/\��q��f��BN�ł����q�j��N�Њ8溊'P�vc>�F�<���<�J�?�q�?�q�?�">��=�VQ=A.��/\��q��f��BN�ł����q�j��N�Њ8溊'P�vc>�l    )��|6}6��#<k�@<8_U>��>8ϊ>^��>&�E                                                �������������O��՘9ļj�ğ�Ā̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         9�Cd        {@��{@��{@��{@��                    2ULl        ;�s�                CI��D�SD�f�D@,�C��=C
��BP͵A�8@��5                                                A���C[SCGf�B�,�B*�=A���@�͵@8?�5                                                                    E���G"iG��F��F9�1E��D�R(C��~C ��                                                C�N�EN��EA>oE�/D�=�C�(C$zB-��A:��                                                                    E��F�1�F���F{cF]F.��E�<�E���EG�                                                D)�D�1D���D��D���D~p�D0�7C߂�C�B�                                                                    6	�7BGX8GMi�G�F�,zF{^F�XE��EPT                                                =�mp                                                                            @>gAQ��A�k>A��B&�BV�SB��~B�/<�<�<�<�<�<�<�<�<�<�<�<�E�E�bEd�E�MD�-nDWP�C�C=C�iB                                                {@��{@��{@��{@��{@��D%K8A��P.����f��8ƜC2��@26�    A1h�A1h�{@�ξ�����C�o�{@��C��tC��_C�o�{@��@]|�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G@W<��D�>?Gw<D�>D�pLD�pLA�nA�nFa8�B:��=�mpC��FK�OFK�OD���D���FbUB:��                @��C��C��C��Q?   C��C��tC��tC���C�H�C��C��C��/C�v�C�<�C�"C��VC��@C��C���C��DC�~�C���C��$C��FC�ؾC��NC�'�C�I�C�0�C�OC��C�}{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C���>܈L>�J>Ο�>�*�>��u>�`�>���>�P�>���>�U�>�W�>���>�\m>�{�>�+H>���>��i>�P�>�N>�R�@I��@�A        A�L�ǖ��ǖ &Ǖ����������"k��JƂ2��ǈ�ǈ8�ǈ��N        ?��BH�I    @�R�@�R�{@��@�`PFTN77�7gv        ;��{?*��?+��?,в?.�|?2!{?6�5?=r�?F!S                                                �[��7f��A�  ?k�A	��@��'�y      ~     ��9$��8��W7�q7>�6��)5��;4�	�3�H�                                                                    8O��7�<87��6p�5�E�4��e3��42��b                                                                    G���F�6�FU-lE��`E7D3&�CEA�BB�                                                4���4
3j^2�ou2EE1?��0R�p/QI                                                                    3�6�3.q2��I1��13��0r5j/�4}.��                                                                    7Au�7*o6�S�6!tU5�F�4�d�4 �G3)	                                                                    7s�7*�6��6 e�5v�4���4��3��                                                                    5�4�se4X�m3ƕ3Ks�2��S1�@1%�                                                                    7��7Pd�6��&6�`5�a�5&v4/�h34qJ                                                                    8�^]8���7��67(�'6�,�5̤G4� "4	                                                                    5���5��4�EW4_43��#37��2y�?1���                                                                    8�^]8���7��67(�'6�,�5̤G4� "4	                                                                    7��7׭�7)��6��H5قu5*��4V��3]�                                                                    4�G4���4ql�3ձ�3L5v2�V1��Y1 :�                                                                    7��7׭�7)��6��H5قu5*��4V��3]�                                                                    ��&3��.�5�5<5~��4���4#Z�3mo�2�!�                                                                    3�o�4/��:<Wz:��9�t89�I�9D�Q8��                                                                    5�S�5gȃ4���4�Y3�Z�2��_2 �^1%�                                                                    7�b6���6;$�5��@4ס�4�3=��2@��                                                                    ��S� ����غ�.�ʹ�����it������                                                                    ��?��ƫܲ�̲��{R��-���f���1}                                                                    6�5�Ϡ5:��4���3�E�3y2=��1Ag�                                                                    �<w��{S��ZQ�7 ���ų�%��Կ�                                                                    ���p��������c��és����M1�p�k                                                                    �l�D�T���֢��樲�ò~��6$��D��                                                                                                                                                                        +,�H+��                                                                                            8w$F8X�7˾�70��6�n�5�[�4�.%4&�                                                                    4���4�g+4M��3�x3�%2Z<�1���0�k(                                                                    7r7S��6Ǘ�6/ms5�,�4�Ur4��3(                                                                    3�Ȁ3���3I�
2��2֙1e�0��/�H�                                                                    8m�8O67�q37.*6��.5�O5p�4-)�                                                                    4�4р�4E�}3�|�3�2o�E1�=0�&                                                                    7�"�7�M7��G7�6_�\5��4�L3��
                                                                    68�5�3�5�qv5!�}4|�3��P2��*1ݘ                                                                    6�d�6�!�6`�5���54��4u�@3���2���                                                                    4�Nl4�KD4��4Ѱ3Nt�2�s 1���0��                                                                    8x�7��7�r;7-!�6��5��24�$3��~                                                                    6"�O6-�5��5E�'4�!�3�m�2�rP2k$                                                                    4���4�7�4�Tq4��4K�F4\3�	�3heM                                                                    3�6_3�F�3qWW3Q�J3$�P2��U2�ޕ2;�e                                                                    4ψ�4�
�4���4�4yCV4,c83ٚB3�                                                                                                                                                                        '�y      ~     ��67�C1���            G���    =+�-6�03��'@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�|     @߃�    16:55:28        F�Z @߃�    @ߋ@    '��      ~-     �p5�Z�7���B'�fB'�f81�D"g�    )z�+D&7*�5?�  ?snj    3�7�)3v;/�q��
��        /��
+�,>���>0��=Ȣ�=R��<Q+75H/�K�2c� :6�<+<�<�-�=H;�=�)�=���=��l=ع�=��$=�B>�>Vv@�A�@���@+~A?�R?ʹ?�
�?��%?���@b��@bf�@X@�@W��@Wkf@V��@Vf�@U�P@UC�@T�Z@TZ�@T=�%�=}��?���                    E�ɑ7]��G�z�F���FT��E��>E�D2�!CE7BB��                                                A�C��B��aA�;A@��@���?��X>ί�=�D�                                                                    ED[j@ɏ�F"��A��0��F_�B1(AC��.��!                                                    {@��B�e�B�e�@cg98�k��
D    �x8Q/G<%%��?   )z�AU	^���2�I    >�y�BJ�1� �A;��@�0x@�p1        Bz:�Bz:�C��LC��LC�y�>�<T@�l�3K�Cڮ�6��6�SEA.�|CKT>�y�C_@�B�HDCi��A��B��CD��B M�B��Be�)=;�    Be�)    Be�)Bebb;T�w    5���A�
�AA�A��@|�'@Vs@�L'?���F��0��.���+�%y1
w)8�96�}5ʶ*G^�G[s�G��G�H�=Ȩ�            >�:>�^>�G>���>��7>�
N>��D>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�d�    AG��AG��{@��7�)z9~k�9\g@8�+�80W�7�wy6���6
~�5y�                                                A�AB��AB��{@��94?)B��oB���@?/?��{@��{@��    C���BH�6�    6E��7�$@�N�2��7ZP.    =�N�C�7�TbB�9D~�C�yGBZ&�A��AC=�@��?�?">�ZG                                                ?4&�A�<�@�ы?�C?
��>���>��=4��<@��                                                                    DH�F��E�5�D���C�-�C@��B�Q�A��@�X(                                                @��B�EHBJ�UA-��@���@�2?�4�>���=�+�                                                                    C�M�EׯE#�C�dC<w&B���A��MA_�@ya                                                @4%�B^�!A��@�д@�}?�%�>��\>/63=8��                                                                    7[0�7إ4A��>HG{A���?&�?��+,.f+�{-M�m,�@�,.f+M+O.S�@-�?�-��*}�7��4��2B�.-mя.��.���)���7�~�4�92P-2AԚ2A�|-ly���~�2�I7�84�nU    3��4��4���/��    :te}8XTY8Bߞ6�ذ5k@�36                        5z��6�NV6��47g3��L7��W2��            4��>g��-,?�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Cns'G���>+sA��v                                                B�                                          B��                    A�                @�p�    /T�7+A��5^�(6GxF1��1#??گ�8Y|3@��80�?;�#;ú�                7���    1��6���    2&�A{@��            ��Xt{@��{@��{@��70^o7�LW            6�0    6�0{@��    6�0    6�oP    {@��2&�A{@��2&�A2*0��    6�oP{@��6���7��07ۣf7ۣf<��2/{B�g�7�JE��D$1A�WE    C5�k            <	,�4ēE6e�4ēE<ε�<�V�=4v<�2�<5`E;��;\2:���"��ފ�Ê����H����׊��ފj�q�N���8�׊&I�<8�<ΰ<�H�= <�n�<3H�;��;EF:�[�"��ފ�Ê����H����׊��ފj�q�N���8�׊&I�8�1�    *��6;�U6�2�9�9�9a߹9��8ʧ�9n�8�%                                                �YR��{�;���H�=l��`���ĤA�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         94?)        {@��{@��{@��{@��                                                    CI�RD�(�D�q�D?��C�m�C
�QBP�cA���@���                                                AɪRC[(�CGq�B���B*m�A��Q@��c@��?��                                                                    E���G<�G��F��F9��E�ԆD�:{C�C �M                                                C�NEN�EA\�E��D�*C��CBB-��A:��                                                                    E��F�6MF��rF{�F]�F.�hE�;�E��EG��                                                D)��D�7�D��_D���D�×D~p�D0�C߁�C�B�                                                                    7��W7Y"\GX1.GM�1G�F�F{�F��E��EPs                                                                                                                                @7'�AJ��A�b�A��B"�$BP��Bz��B��}<�<�<�<�<�<�<�<�<�<�<�<�E�4E�(aEd]E��D�$&DWK�C�?�C�h�                                                {@��{@��{@��{@��{@��D#�A�!G/Sާf��8��CmY�@K�e    Ak�Ak�{@�ξ�L���L�C�y�{@��C��[C��`C�y�{@��@cg9{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��<C$}D�@�G�wD�@uD��.D��.A��A��FbűB=߫    C��CFK��FK��D��tD��tFcYoB=�                @�5�C�L+C��NC��_?   C���C��[C��[C�g~C�|C��9C��C�K�C���C��C�Y�C��C��>C��DC�n C�;eC��C���C��@C�ױC��C��NC�cC�=�C�1dC��C��C�_v{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C��? ı>���>�*�>�	>�!^>߸�>�Z>��>�/@>�1�>�E�>Òy>�:>�(p>��>���>�t�>�?>���>�ږ@j��@+F        B��ǣ��Ǣ�\Ǣr�����"��"�!t�Ɗ ���[���\#��[���         ?4�B7�    A��A��{@��@���FT@r7i<�7��        ;��`?&d�?&-?'	�?(�s?+�?0�?5��?<(�                                                ���7(`9A�  ?k�A	��A z�'��      ~-     �p9:�48��&7���7S!#6���5�M^4ٰ�3�5�                                                                    8k��7�n7@26�XL5�7x4�j4	}>3�                                                                    G�z�F���FT��E��>E�D2�!CE7BB��                                                4���4�^3���2ݍa2.^1R)B0db�/^�                                                                    3�O�3C�i2��2�1F�i0���/�>p.���                                                                    7[��7?E�6©V63B�5�3�4ׂ�4sd3(k�                                                                    7���7/jR6��.6A5|�94ڧ�4��3Z�                                                                    5��4�+o4MK=3�T\3P�2���1�{c1GO                                                                    7�v�7Ve�6��6!ul5���5�U4/�{38�'                                                                    8�Q�8���7�C�79��6�̓5ݒ�5֧4mo                                                                    5�u5���5�{4r3�<3EY�2��1��                                                                    8�Q�8���7�C�79��6�̓5ݒ�5֧4mo                                                                    8��7�}780�6���5�h58<�4e"�3i��                                                                    5	*�5
��4��03�y3]F�2�td2	$1/�                                                                    8��7�}780�6���5�h58<�4e"�3i��                                                                    �H�23�s�5�l�5�St4�wC48�;3��2��                                                                    5A��5�RT7��7;.M6�i6+M5�J5,a�                                                                    5���5o>�4��4#
U3��w2�cL2!/1)n7                                                                    7�7$�6LUN5�T�4�o�4,<x3KT�2K��                                                                    �#`|����B%t��Q+�͘���w��t                                                                    ������g���ߖ��"}��	�>����Ǳ���                                                                    6�46��5K��4���3�
3+�2K@�1Lw;                                                                    �-k��Ƕ�_����J\洗����C�⊌                                                                    ���S��Li���|�~��6ϱ$F��^Pn��Nx                                                                    ���D�o64��j����� ��)|�Eh8�Q                                                                                                                                                                        +�G:,<��                                                                                            8�w�8r��7�_�7D66��,5���5��4��                                                                    5U4���4e�3���3!�2o��1��D0��                                                                    7���7m��6޽�6B�y5�4��+4Ҕ3)yy                                                                    4 �3�Y3a4z2�u2&��1{^B0�O�/���                                                                    8��<8ij7��7A9�6�>�5�ߕ5�348ce                                                                    5;�4�4\�N3�x3,Ha2���1��0��U                                                                    8��7�y[7�]7��6w3�5��'4�R�3�{�                                                                    6�6�}5�3�54&4�BH3���2�?1���                                                                    6�C�6�H56z��5��<5G4��A3�D2���                                                                    4���4�w4�1�4~�3dL2��1�v0���                                                                    8"	�8�7��7@�6��5��L4�d3�^                                                                    69/�6$uD5ؔ�5\4��X3�
�3��25�                                                                    4��4�I�4��4��L4a�74�U3��3w��                                                                    3��03�~�3���3i��36X�2�-�2��2H�                                                                    4���4ݒ�4�/�4���4��04=2�3��3�C�                                                                                                                                                                        '��      ~-     �p69T3��            G�t    =A76��3���@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @߃�    @ߋ@    16:55:31        F�� @ߋ@    @ߓ     '�A      ~L     �@5��,7��9B'�fB'�f8+xD"g�    )��+<�7@k�?�  ?uo�    4��)ys/���6�        /�5"+�f>�Wy>!=�%b=*��;�9�6���/[2`�):))<*]�<�=H5�=�'H=��3=��m=عR=���=�B�>�2>V�@�uK@��@4�c@n#?�|`?�O?���?�"@f@e�@[�@[+�@Z��@Z:"@Y�s@Y�@X��@X�@W��@W�j=���=�,{?z1:                    E�j�7g��G��F��NFT;QE�)�E�vD2��CD��BB��                                                AN�CGB�}\A޹$A@6@�J�?�m>�j�=�9t                                                                    ED�D@��F"�xA��0��,F`qB2E3AG��-���                                                    {@��B��_B��_@h-�8��A�l�    ����/Vh�&�?   )��AV�ˬ��2g�    >�y�B5�z2�A&!DAMiA��        BW��BW��CƧYCƧYC���>��@�]�3��C�~�7˂6�ߖA13�C0~I>�y�CB4B��mCP��A�7BtyeC.��B�B�zXB38��c    B38    B38B9�Q            A��WA'ήA�d�@h�&?� �@���?�MF��0��,.��V+�Hn19��8�	m6���5ғ�G��G<$�H��G���>
�h            >v�>���>��b>�7�>��i>���>�Ŝ>��'8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�dB    A���A���{@��7��9���9g�8�ܦ83��7��6՗�6�5*                                                 A�/�A���A���{@�ΠG�xB�8�B���@E�q?
�q{@��{@��    C�N�BIn6��    6I"�7�?@��2օp7w�    >
�hC��7��tBL
Dt_C���BC�\A���A9�t@��x?��a>��                                                ?'��Ay��@��v?u�?��>�j@=�V<=+B�<<��                                                                    DFTF���E�5D��kC�ZmC>��B�q�A�m�@�6�                                                @��B�nBE��A'6�@�� @0?�zk>���=��l                                                                    C���E���E z�C�"C9��B��A�GGA"�@Y�                                                @093B[��A���@�+|@�5?�.�>�s�>+�n=7�                                                                    7bk�7���A�N>M�A��y?+��?���,�v+h	-��1-�4,�v+���.�d�.��-���*��7��W4�\u2GW�-���.�A�.��v*H�7�'�4�j�1�;�2F>2E�	-�i���'�2g�7�(�4���    3��c4�74���02�A    :�߈83b8/�6�{5H��30�6                        5V#A6�Z6~v4_��3�<7U�2�.�            4� ;>F��-��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CXKVG��S>+��A�5                                                B�                                          B��                    A�                @�p�    /��_+�Q5BW6P%�1� 1N�%?�+�8?T�@{(�8��<	h�<��                7�S�    1p��6�k�        {@��            ���{@��{@��{@��7���8��            7��    7��{@��    7��    6�M�    {@��    {@��                6�M�{@��6�*7�u�8v�8v�<��x1�_lB��c7�ME�w�D%ԵA��    C܊            <&j�4�=N6�mo4�=N<���<¨�=>f�<냲<b�f;���;?W�:����ЋM��P'�����A����Ŋ�M�kÜ�O�k�9���'�V�
4<"�3<��<p=<;�<�\1<^X;���;0�e:�����ЋM��P'�����A����Ŋ�M�kÜ�O�k�9���'�V�
49]��    *m�Z6U�6��":
��9��19���9���9g�y9d�                                                �E?��H���FP3�C��<��2{X�%zN��"̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �G�x        {@��{@��{@��{@��                                                    CILD� �D�&�D>�6C���C
 YBP>yA��Y@��o                                                A�LC[ �CG&�B��6B)��A� Y@�>y@�Y?�o                                                                    E��|GY0G�F��PF9��E���D� �CC �G                                                C�LWEO]EA|�E��D�UC��C �B-��A:�
                                                                    E��F�;QF���F{ iF]�F.�E�:�E��LEG�                                                D)�pD�?D��SD��5D��9D~p.D0��C߀oC�B3                                                                    7U�7a�pGXN�GM��G��F�pFz�F��E�EPs                                                                                                                                @(�A<m$A�gLA�ŷB��BBh�Bf��B�v�<�<�<�<�<�<�<�<�<�<�<�<�E�.sE�7LEc�rE��D�HDWE�C�;�C�g�                                                {@��{@��{@��{@��{@��D$�A�~/04V�fz�8y�C���@]��    A�A�{@�ξ�}��}C�JV{@��C�e-C�e-C�JV{@��@h-�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��Gɤ;Į�D�H-G�D�HD�?�D�?�A�A�Fd!�B@%�    C��=FK�pFK�pD��D��Fd�)B@ �                @�vYC�^�C�CeC�zJ?   C���C�e-C�e-C��C�t�C�)C��CC��C�%�C��#C�g�C� tC��1C�� C�GAC� C���C�zC�FQC�!�C��C��C� C�3�C�0�C�]C�C��
{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��OC��?�??	�=?h{?*�>��G>���>�T>�_>ᩨ>ܚ�>׈�>җ�>���>�P>�� >��>��:>���>�j�>��@2T�@ ֵ        B��ǙMǘ�1ǘI8�	7���|����ƃ��������:����G�        ?*YB">    AKɁAKɁ{@��@��FT�7�Z�7�e        ;���?<�?p�?޴?�l?�?{i?!$�?$�                                                ���6���A�  ?k�A	��A z�'�A      ~L     �@9Fٽ8�h�8�n7X̤6�i95��4�k_3�+r                                                                    8{-�7�
�7#��6��	5�A�4�i�4(�3 QU                                                                    G��F��NFT;QE�)�E�vD2��CD��BB��                                                4Я�4"
E3��2わ2�}1Pֿ0^��/U1�                                                                    4�E3L��2��"2��1H�
0���/��x.��D                                                                    7k�7H{�6�(67��5��=4��4�3!%0                                                                    7��7+�o6���6 ��5s�S4̉j4�3e                                                                    5c`4�24<�-3� 3F�2���1�N�0�߈                                                                    7�s7R6�6��t6Ar5���4��I4#*3.	�                                                                    8��8���7��=7<��6��*5٤�5�4�                                                                    5���5���5�l4r~43��z3@;�2~��1�׍                                                                    8��8���7��=7<��6��*5٤�5�4�                                                                    8�7��7<r6��[5�)54��4\��3^                                                                    5~5��4�&�3繼3[V�2�]V1�T11 _G                                                                    8�7��7<r6��[5�)54��4\��3^                                                                    �>yR2Ⱦ�5�x5�5B�4:��3�d�2��a                                                                    5�?5�{�7���7pk�6��46Z̕5Ч�5L��                                                                    5��R5k��4�%�4 C�3�L�2�9�2ץ1�f                                                                    7+A7m�6R��5�S�4�I�4)�3D�>2BFu                                                                    �.�X�V`��E��G���ѳ̂���䤱�^                                                                    ��cH��p:���#��Ip���9�=ث��l��t                                                                    6* ]6��5Q�4���3�J|3(�-2D41B��                                                                    �8�&��L��/7���[�L?���rU��4��ح                                                                    �������"��:��p2��:��#�X�+�up=                                                                    ��$�{)�$Ea���²���(�ѱ@Ϫ�H֐                                                                                                                                                                        ,��,�Կ                                                                                            8�:8~�7뱞7H�R6�*�5�Wa5s�4Ԇ                                                                    5�)5 �G4nL�3Ͷ�3"�"2m�e1���0�t"                                                                    7�?7yTg6��+6GK�5�u+4�54
83"0                                                                    4��3�?3iq�2�>I2(<:1y�60�/���                                                                    8��E8t$�7��7Eܢ6���5��5��40[�                                                                    5��4���4d��3���3-�R2��1�I30���                                                                    8�7�<r7��L7!�~6z!:5��v4��l3ų5                                                                    6"(�6G5��59"�4��j3���2��1��a                                                                    6�P�6��i6���6�5J 4��3�__2���                                                                    5	�4�S�4�VX4��3g 2�CJ1�Hl0���                                                                    8-kq8�7ţy7E�)6��j5��s4���3�%                                                                    6F1�6,�?5��e5bG
4���3��_3��2
�                                                                    4�m�4�=�4�K�4���4dy�4?3���3m	?                                                                    3���3���3��{3p`P38��2���2�hh2?�Y                                                                    4�M|4�J4�n4���4���4<JM3�Q3���                                                                                                                                                                        '�A      ~L     �@6<��3��c                        6�g]3�$j@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ߋ@    @ߓ     16:55:33        F�� @ߓ     @ߚ�    '��      ~k     �5��7�,�B'�fB'�f8Y�D"g�    *!�+< 7J$�?�  ?v�}    3���)C�/��l�%�N        0Q�,H1�>�>?>1~�=�v=k�!<pw�7�C�05��2VX�:}a<)Og<���=H,=�#O=��=�٧=ظ==��]=�BM>�5>V�@�@�XN@K�X@+�?�)?݊�?�(`?�j�@iD�@iv@^�v@^A�@]ԙ@]_0@\��@\[?@[�=@[nG@[�@Z��=��=��J?���                    E��7>b8G��OF�7�FSĘE�ҏE�HD2[CDzBB��                                                A�CҒB�%�A�>�A?ܒ@�[?�%�>�%	=�0                                                                    ED��@��\F"�MA�?�0�Q#F`�0B3�QAQ��-j��                                                    {@��B|�;B|�;@m��8o4q%Ƃ4    �p��/c�z��u?   *!�@��֬��z2
&4    >�y�B/��2�:AE.�@�Ŭ@��        B�ݬB�ݬC��C��C�қ>�{@��3��C�H7
�D6�ܨA4�Cr�>�y�C0�B��mCJ� A�0�B\��C(0�B�FB�ϟA�'��Iv�    A�'�    A�'�A�Vb            A���A ��A�@�@R��?�@�ތ?��F��0�Q#.j�+�18�(8��l6h�n5�b�F��^G_�G�#�G�|'=��            =�_>~}�>���>�b�>�i�>���>�	�>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��I    A��`A��`{@��7���9W}k9>�f8��8��7}m�6���5��%5O�                                                A�DA�y,A�y,{@�Ο�x|BE��B���@Q�f>�Ĕ{@��{@��    C��BJf�6��    6Nt�6ЁY@���2���7o�4    =��C )7��B��Dk�Cn�LB8�dA�/�A6l�@��?�.T>��                                                ?��Ao7�@��&?c`>�V_>���=�=Y=%�_<;��                                                                    DD��F�[TE��,D�T	C��C<�TB���A��6@�t�                                                @�t�B��KB@�EA!��@�y@��?|t>���=��b                                                                    C�&bE�8E�C��C7��B� gA�;�AIV@�T                                                @,��BYR�A�6�@��~@��?��L>��@>)=�=6�                                                                    7;.?7�.�A��>T �B��?0��?�HZ,���+��(-�K�-1�R,���+��5.��.&z�-�"�*��7�Y4�2,�.W~�.��.�5�*��7.W4�<�1���2+g�2*� .U�ڷ.W2
&471�4��    3R<�4���4�k0�H�    :�7�ͽ8_-6f75137��                        5o�6@"i61��4h��3[��7wZ2ңn            4���>1�o-��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CU@tG�4�>,3�A���                                                B�                                          B��                    A�                @�p�    0
�, ��5��6'�X1D�1G{?��8$��@@��7�^;���;�V�                7�}    0�GA6�f�        {@��            �S��{@��{@��{@��7
�)7���            6���    6���{@��    6���    6�W�    {@��    {@��                6�W�{@��69�w7�x�7�M�7�M�<�x�1�
B�?7(�F�C��lAa�    Ca4            <1k.4�b�6�U�4�b�<�I�<��`=J<�l�<~��;�};jK�:�g���
����\��Z���#׊��芊���pa��S��=<Ɋ+
��	�<*�<�E<��=D��<���<s�;⾄;U�U:�����
����\��Z���#׊��芊���pa��S��=<Ɋ+
��	�9�h4    )���6J@6���:�F$:z`�:1�9��9���9L��                                                �
�<������?���6ŴΡŜ�Ņ�%�l֏̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         ��x|        {@��{@��{@��{@��                                                    CHg(D��D�D=&�C���C	��BOh�A�5@�S�                                                A�g(CZ�CGB�&�B(��A���@�h�@5?S�                                                                    E��nGs�G(F��dF9�\E��{D��C�xhC {�                                                C�I[EO4DEA��E�2D��Cֽ�C��B-n�A:��                                                                    E�F�@F���F{#�F]�F.��E�9�E���EG��                                                D)�CD�E�D�D���D���D~o�D0�/C�IC�A�                                                                    7wZ7<;�GXjmGM�ZG�HF���Fz�[F�E�EPh                                                                                                                                @ �sA�A���A�{yB#.B3�!BW�kB��<�<�<�<�<�<�<�<�<�<�<�<�E�BE�G�Ec�JE�aD��DW>:C�7;C�g                                                {@��{@��{@��{@��{@��D,A�!�/��6�fm�8d�C��@jL$    A��oA��o{@�ξA4��A4�C���{@��C�Q�C�Q�C���{@��@m��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��;=�8D�N�G�D�N�D�ԒD�ԒA�wA�wFd�BB8@    C~��FK�bFK�bD��D��Fep,BB3                @��C��mC�-MC��?   C��cC�Q�C�Q�C�ʰC��XC�e�C�6%C��{C��(C�r.C�+�C��,C���C�_nC�C���C�f�C�`C�«C���C�S�C�6_C�"5C�-�C�/�C��C�@C�ό{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C~��?}z?	}�?��?�L?��? Xp>�w�>�u/>��>�dD>�T�>��>ۍ>��y>��s>�W�>�j�>�u�>��i>��T@�?�Zn        A�
sǝ%ǜ��ǜ#O�<�L�������Ɲ�����J��
��        ?P�jB3-    A>�3A>�3{@��@c�FT��7�ʀ7�s        ;l">��? �<?QM?;�?#�?
��?%�?��                                                ��:�6�޹A�  ?k�A	��A z�'��      ~k     �9g�8};�7�5D78�,6��s5��}4��3�X�                                                                    8F�47���7�6i��5���4�%3�u�2�}�                                                                    G��OF�7�FSĘE�ҏE�HD2[CDzBB��                                                4�2�4�3a��2���2c`1=h�0P�X/K�                                                                    3Ы�3'��2��*1�10�0o@�/��.���                                                                    7:�
7$�|6�Z�6f5�a;4�ټ3�oF3y                                                                    7S��7�%6U�K5��+5Rm4�4�3��3��                                                                    4�^�4��C4y�3�J3)�2�1b1�!�0�a\                                                                    7�Pa7$�6���6��5�c�4�@{4,�3&0%                                                                    8�g�8��m7�7�6�Z5á4�3�w�                                                                    5zB5�M\4ܛ4K�3�4M3+��2k8�1x��                                                                    8�g�8��m7�7�6�Z5á4�3�w�                                                                    7�)�7� )7K�6yB�5�I�5!��4L�3Sol                                                                    4�o�4�i�4PT:3��3=�2�1�1�Z�0�B�                                                                    7�)�7� )7K�6yB�5�I�5!��4L�3Sol                                                                    �K�f��A�5�5�d�4��s4*�y3v
2�X                                                                    4�!56l�7�x�7���7�f6�265�C                                                                    5�59�4��S4f�3zS�2�8!2	��1A                                                                    7��6�K�6-��5�E�4��4��36ې298�                                                                    �
�ڵ��D��˯�*j���U���� ��e���                                                                    �������޲��޲��f�we��,n������w                                                                    6�e5��5,Y4��.3�9w3�26Y�19iL                                                                    ���� �&�o�ֵ֤��2���ೲ�n���                                                                    ��Mֳ�C��3�[�a���s��԰Jqïj|p                                                                    �e�N�������|(�� ����4�t�@�                                                                                                                                                                        +���,�k                                                                                            8n�8Q)*7��$7*�_6�-�5�q�4��i4;�                                                                    4��4�yK4D�3��3K�2W1��0��3                                                                    7iE�7L�6��6)X<5��Y4�F*4�33�z                                                                    3��=3�*.3@�D2��;2> 1a��0�r�/�k                                                                    8dk�8H�7���7( 6�p�5��5w14(}1                                                                    4��_4��4<��3�LY30�2lcC1��s0���                                                                    7�t�7�f(7��7
J76[�5���4�	3�                                                                    6 �5�t�5��J5�4{�3�-~2�x1�#�                                                                    6�/�6�\�6Yo�5��51|�4s{�3���2��Q                                                                    4�6�4���4x�3�m_3J׾2�"-1�kR0��\                                                                    8	�[7��N7�o�7)`6�9�5�"m4ٙ;3�&                                                                    6v6v5��w5A*�4�f[3�p}2��2�                                                                    4��4�n44��4~$4H��4�13��p3bʓ                                                                    3�s�3|�#3jv3MU�3":Q2�#72��27D                                                                    4�U&4�1\4�O�4�H�4u^�4+;3�P�3��=                                                                                                                                                                        '��      ~k     �6A��3R<�                        6�0?3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ߓ     @ߚ�    16:55:36        F� @ߚ�    @ߢ@    '�	      ~�     ��5��r7���B'�fB'�f7���D"g�    )�A+a�B7@��?�  ?}��    3�+](���/��аw        /Ȕ�+��>>��R>8`=ى0=n�.<�E7�"0u��2E��:`�<(m�<��=H!�=��=��m=�ׁ=ض�=��=�A�>�->V�@�?�@���@3%j@J�?ߨ�?ǔ?��?�sV@l��@l��@b8�@a��@acg@`�@`e�@_��@_P�@^՗@^u*@^=|?	Q�=��F?ۯU                    E��&7Q&+G�=RF�8FSe�E���Ek�D2,�CDQjBB�                                                A�C{�B�NA��&A?��@��g?��N>���=�B�                                                                    EE+@�7�F#)�A�p�2
�Fa2�B4�LAE�o.2Ɗ                                                    {@��B
��B
��@s�8P����ц    ���9/,d����?   )�A@�*���1�O�    >�y�A���2
�@�O=�|�=���        Bv$�Bv$�C���C���C��1>۪�@k|3H*C��m7M�6�,�A6�[B̩)>�y�B�VB*�KC�}Ab�/B҆B銚Aç�B��3@½�����    @½�    @½�@�<q `��    �i9A7�m@ؕA�{@��?��@�" ?��~F�a2
�/�B�,�|�2AD�89\�6%�5�F�ɮG?�G���G�r�=�z�            >-\>���>�N�>���>��R>��%>���>�1�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�Y�    A={@A={@{@��7΀�9j	�9R�8�Қ8)y|7��6͝26�5�                                                A�W�A=D�A=D�{@�Ο�i�B�TBhw�@`H�>�̚{@��{@��    C��@BK�6�g�    6T4�6�@�@�O�3�o7Y��    =�z�C��7�v�BXDf��Cb�bB5�WA��lA9)�@���?ֵf>��                                                ?��Ah�@|��?^��>��|>�(�=�v='�<A�                                                                    DCZ�F�JE�W�D��`C���C<D�B�H!A��_@��                                                @��B��zB=)�A��@��w@ܟ?z��>�Q=�}                                                                    C��iEՓ�E�8C��CC6f�B��SA�SA��@#�                                                @*Q�BWe�A��@��@/Q?�L>�@t>(~�=6��                                                                    7M��7ˀ�A"��>Y�WB�/?5��?�%�,�V*��T-:2,��a,�V*��p-�wy-��-C	S*T�L7�S4�ϟ2(��.	��.���.�R�*�]���4���2
Z�2'��2'!T.�86�1�Oٶ���4��    2��4���4��0\�    :�B�7��7���5�?4/��2HXV                        3���5m{75]�%3~A2E�L5�hB2�I�            4��=}�8,�z�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C$]�G�Rh>,NIA���                                                B�                                          B��                    A�                @�p�    /��+w�w4�D66u�2ID�2D��?Rox7��N?�[7tp�;ԑ�;��                7hK�    1���67�5    0��{@��            �x�0{@��{@��{@��7�S�8*-8            735o    735o{@��    735o    67�G    {@��0��{@��0��
�.k�    67�G{@��6&*x7�8.ު8.ު<�'�0�|OB��7 �F(s[C�I@�8h    B���            =,�4��g7��4��g<�W�<���>I��>�J=�1�=o<��/;��i��B�
�Z��׽��(����U���ӊ��Z�p+ڊS�Y�=.�*� ��='r�<�S0<���>G�>@�=���=��<p��;�k��B�
�Z��׽��(����U���ӊ��Z�p+ڊS�Y�=.�*� ��:���    *M�R6�P6��; h�:���:�|�:�7�;ڻ:���                                                �DuN�D�5�By:�A���A�'�Bb�D��D7�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         ��i�        {@��{@��{@��{@��                                                    CG�:Dڤ�D��D;�C���C	VBN��A��@�1�                                                A��:CZ��CG�B��B'��A�V@Ω�@�?1�                                                                    E��G�G-*F���F9�E��)D���C�a]C w                                                C�E"EOR�EA�E��D��9C֝SC�uB-^A:۾                                                                    E�(F�C�F���F{&�F]iF.�E�9E���EG�\                                                D)�D�K|D�	�D���D��dD~o�D0�C�~�C�A�                                                                    5�hB7OUeGX�GM��G�%F��Fz�+F��E��EP�                                                                                                                                @0x�AB��A���A�B}�BB�nBdlB��<�<�<�<�<�<�<�<�<�<�<�<�E�Q}E�U�Ec��Ew�D��DW7�C�3qC�f�                                                {@��{@��{@��{@��{@��D&g�A�	/%vާf��8�oC�]�@j)�    Ai��Ai��{@�ξO	վO	�C�ts{@��C�JmC�JmC�ts{@��@s�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�:��D�Q�GoD�Q�D��dD��dA�A�FeF}BCv    C�óFK�-FK�-D�� D�� Fe�oBCp�                @1)C��=C�J�C�j\?   C�Z�C�JmC�JmC�N.C�L�C�H�C�AAC�4eC� �C��C��C��2C���C�wsC�E�C�yC��-C�y�C�+VC��C��1C�qvC�ABC�.�C�.�C�7C�rC�L�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�Z�C��>�r1>��>��>�>��>��">�L�>�{>�A>��>��)>�f>�+�>�>�Q�>�0O>�M>�YZ>�cH>�d�@ �?�K        A����H�S�He��G�zŴ�o������B�ɤG�-�k��>l��>���>l�I��        >� FA�#s    A�mA�m{@��@rFU�7���7"         ;j��?!m?�l? 8K?�??�?s#?�}?��                                                ��6���A�  ?k�A	��A z�'�	      ~�     ��9*�8���7�"�7M��6�~*5���4�vV3�ҩ                                                                    8W�=7�W�7K�6���5���4�4e�2�%                                                                    G�=RF�8FSe�E���Ek�D2,�CDQjBB�                                                4�9�4��3}
2׾�2��1Jv�0X�n/P�                                                                     3�c�39V2��!2B�1@ܼ0��/�ѿ.���                                                                    7J�g76x6��$6-�Y5��,4���4�@3k�                                                                    7a�7��6lV\5�/5i��4�̽4 Ѯ3�                                                                    4��34�:4!:�3�3<I2���1���0��S                                                                    7��f7/��6�m�6�r5��4��4r3/�v                                                                    8��8��N7��70�K6��]5�e&4�=3���                                                                    5�F�5�T�4��4_��3��36V]2s	.1��                                                                    8��8��N7��70�K6��]5�e&4�=3���                                                                    7��C7蘤7*�i6��5�R45,a�4S��3X��                                                                    4��54)4dL�3Շ[3N��2��]1��)0��                                                                    7��C7蘤7*�i6��5�R45,a�4S��3X��                                                                    �a�R���.5�i�5�U�4��W46-�3}�N2�͞                                                                    3��4)��8�9J8��68f-7��{7�k6��>                                                                    5�i5F�;4���4�3��>2�٢2�a1"�R                                                                    7�:7	?�6@�V5��t4�x�4"ώ3=>�2=�x                                                                    �!��
������=����R۳Ǝ�㤵��P�                                                                    �ۥ����ȭ��^����8���c/����                                                                    6�s6\'5?K�4�@�3�^3!��2<��1>�                                                                    �^A�T.�������C]഑7����Ӆ�                                                                    ��!l����=��sڞ�н2�D�Q^&�o��                                                                    �y+ܴdn)�[���]���#���;���Du                                                                                                                                                                        ,8,�I2                                                                                            8�\$8g�7�&f7=S6���5ӌ�4�5w4P?                                                                    5�[4��4[�$3�ƽ3�@2eH�1��0�xY                                                                    7}r�7bY�6Թ�6;�L5��V4��4
E�3=                                                                    4  3��^3W�2�d2 �j1p�J0���/�>Y                                                                    8x-)8]�P7�L�7:aE6���5�i5p�4,)�                                                                    4��4��4R��3��3&S�2| 1�i0�Y                                                                    7�G-7��p7��%7�6p��5� 4��3�cG                                                                    6(�6 }�5��5/�4��03�$�2��1��                                                                    6�4X6���6s�k5�5B��4�3�3�F<2�E�                                                                    4��4ϩ�4�j�4&
3^O�2��e1��i0��                                                                    8�b8	jD7��-7<�6��5��:4��3�\�                                                                    6+N'6�5���5V��4��3�B32�                                                                    4�<�4�ܹ4���4�~�4\(;4�K3��53g�6                                                                    3��3���3��;3d��31�2��2�sf2;m�                                                                    4��,4�F�4��4��34��]46�x3��$3��                                                                                                                                                                        '�	      ~�     ��6G�2��                        6�N�3�{@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ߚ�    @ߢ@    16:55:39        F�P @ߢ@    @ߪ     '�m      ~�     ��5�.<7}��B'�fB'�f7���D"g�    (���+��N7)I?�  ?�82Z;�3�`c)��n/�Vt���^        .��3*�Y�>�l�>N_=�a�=���<�4?7��'0+��28J0:�<'�]<�Ň=H�=��=���=�ՠ=ص�=���=�A�>�+>V�@��3@v�9@
�?��!?�נ?y�^?S$?�L@n/�@n�r@d0-@c��@c�q@cV�@b��@b��@b�@a��@ad@a2�@u��=d�{@��                    E�u]7U�{G�ݛF�FR�=E�(�E%�D1αCC�CBBe�                                                AQxC�B�A�V�A?.`@��j?���>͓�=��                                                                    EEL�@�kF#W�A��3��Fa�@B5i&A3F�. VN                                                    {@��Af��Af��@u�86�$���    ���/=���A?�Z(���@��|����2�l`    >�@�@���2�p@hJ��z���        BU�DBU�DC��_C��_C�/7>��3?�&y3e�CC�V�6�6�t�A8T�B9>�@�BQ��A�~�B��ZAL^A>�aBXiAdaBK��{QJ?�� `B��{QJ    �{QJ����>3�;7�;;f�@ĝ/@9�A@�m?��3?E�8@Z��?^�F��3��1�@.#��3�%�7��a5g�J5Ƚ=F���G3G�Y�G�6i>�u!    ;��*    >F�H>�-l>�ѯ>ի8>��>���>��>��8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�&    @��@��{@��7�D9exr9U2[8���84Ś7�^�6��6 �5K��                                                A�h
@�m"@�m"{@��=BFgA���BX*@jܤ>#�{@��{@��    C�4�BL�e6�b    6���6�~=@��I3���7,��    >~��C�p�7�U�B	E�Dcg�CZC?B89�A�y�A=��@��?�5�>�'�                                                ?��Af�1@u�[?g��?��>���=�o�=(��<A%�                                                                    DB&�F���ẺhD��C�J�C<�B��UA�S@��2                                                @���B�<B:�A��@��X@�u?y�>�8N=�m�                                                                    C��NE�;�Ee�C�A�C5�B�Y�A�,wA@	�                                                @(�`BV�SA���@��@#?�*�>�T�>'��=6�n                                                                    7TL7��A'>^~�B?�?9h�?�J+"�C)ó�,��[,�P+"�C*w�,q��,�;,���)��'7�m)3�>W1���-;�.��W.��*�y����4��1��1��"1��)-�7���1�:s����4�<    -h��3��83�n�/\��    ;��76�>7�3"�o1�֧/�1	                        /���3�83�0�c�    3#�O3m�            4
>7;���,���?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B�SEG�l�>*�eA�k�                                                B�                                          B��                    A�                @�p�    .j?�*Yǯ3�66��3�Ǝ3��">�p7��B?u�e6�X~;��y;���6�I�            6�W�    3V1#5��U    4��{@��1"D1�qi    �?̮{@��{@��{@��8�8g�            7x3�    7��Q{@��    7��Q    5t�X    {@��5�9{@��5n�4��a4��    5~${@��6=d�6Ư8m�8m�<�X    B���6�$�F@�VC�Ȩ@Cz�=�+�B,�k            >�S3�E�9+{S3�E�<��<���?�H?�M>��u>Z"�=��=F�d�� �
�@���9��2 ���銤�ъ��@�p6�S�U�=7�*�c��1>���<��e<��a?��-?�>�\�>X��=���=ҋ� �
�@���9��2 ���銤�ъ��@�p6�S�U�=7�*�c��1;۱�/*IW�5�_6p�3;�t;�i:���:�gV<z��<B�I                                                �z�Ā@��u*�k�3�\�@�C��� AG����̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         <��        @�<@��@6�X@%d�<k��            4_�-:=�        <��
4_�-            CG2�D�qD��D:i,C�)�C��BM׷A� @�t                                                A�2�CZqCF�B�i,B')�A���@�׷@ ?t                                                                    E��7G�KGB�F���F9^�E�n�D���C�@C m�                                                C�=gEOt�EA�;E�5D��pC�rCC�~B-E�A:��                                                                    E�
F�H$F���F{*�F]VF.�BE�7�E���EG�                                                D)��D�Q�D��D���D��D~n�D0�C�|�C�@�                                                                    3#�O7X)<GX��GMאGx#F��Fz�F�TE��EP                                                 <�                                                                            @vfAz8�A�ĔB�B:xBo��B��GB��<�<�<�<�<�<�<�<�<�<�<�<�E�bXE�fEc��Eg�D���DW/MC�-�C�eJ                                                {@��{@��{@��{@��{@��D��BQ�."�f08T6C��h@\sN    A&q�A&q�{@�ξ�r��rC�k&{@��C��`C��_C�k&{@��@u�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��Gŷ:X��D�P�G�qD�P�D���D���A��A��Fe�BC5�<�C�D�FK�lFK�lD���D���Fe�kBC0[                ?�L�C���C��1C�e�?�(C�GdC��`C��`C��mC�:C�0�C�]�C��jC��AC��^C���C��C�C���C���C���C��(C��C�W>C�7C���C���C�gJC�6iC�.)C��C��C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�C�C�}�>Á>���>�%>���>��D>ҷ�>�n�>׽q>�*�>��>�j�>�=�>�w�>��>յn>Ҝ">��>���>��>ú@6�@-�        A�\��ww��:������f��6���6�G�62�R`����������ĉu        >�A�\�    @��@��{@��@���FU_�7d�7M        ;�3�?F�?F>?H�?I��?L�R?R?[�9?m�B                                                �6�6�=A�  ?BV%@��@�!'�m      ~�     ��9'M�8�"7��r7Z��6�|�5�E�5 �4��                                                                    8ST�7�;�7Ι6�N�5ʸ�5�4!�p3+��                                                                    G�ݛF�FR�=E�(�E%�D1αCC�CBBe�                                                4��J4�3�2�Ձ2(m1jL�0�iR/��"                                                                    3��3;X2���2(�1T��0���/��h.��                                                                    7G�78|[6�	w68��5�o4�&i4#I�3WI�                                                                    7Zn.7`?6x�h6�{5���4��e4!f�3E1                                                                    4�&g4��4-��3�^3V�2� X2+;1+<r                                                                    7�|97.�6���6#��5��5�4ED�3q                                                                     8��{8�.7�F77<V16��*5���5��4-��                                                                    5�nf5��4�%�4n��3�\3R��2�c1�n                                                                    8��{8�.7�F77<V16��*5���5��4-��                                                                    7�#B7룇7/A�6�Ӑ5���5G,f4��b3�\�                                                                    4��5�&4i��3��!3d�2�5�2@�1+�k                                                                    7�#B7룇7/A�6�Ӑ5���5G,f4��b3�\�                                                                    �F�ճ���5���5��F5t4R13�DG2�dC                                                                    15y�1Ό�:[o1:�f9��`8�6�8oJl7៶                                                                    5��5C� 4���4&��3���3w:26��1_%�                                                                    7>�7<-6FA{5�=4��4<7[3j�$2��C                                                                    ����}(��#��JY2�����N�s�"��                                                                    ���ǲ�䮲�q���'M���U�ƱÉ�                                                                    6�6
8�5Dlu4���3�w�3;*2i��1���                                                                    �c�/��ߐ��ՙ�W1����%������                                                                    ������>�eಁ�
���u�5��������                                                                    �t㖴g�.�7?��Π��=��i5���\�                                                                    )*�)�f�.��U.X�.-��,- ��/�s�/��~                                                                    +�Ѝ,{e5                                                                                            8}�8j�7��(7H�6�Rc5�Wo53�4EVP                                                                    5 X�4�k4bL�3��-3+K�2�iz1��0߈I                                                                    7x�K7eL�6�C�6Gy�5��\5 94+E�3X`�                                                                    3�vc3���3]�X2�m`21@G1���0�6/��                                                                    8s��8`��7ִ7F
A6�xT6Fg5:Xe4kkv                                                                    4�;f4�u4Y�3���374�2���1��1U�                                                                    7�7��7��b7$3f6��5��e4吴4CV                                                                    6	�6J5�'5;�t4��h3�c�3.2(c                                                                    6��O6�?k6}LB6��5V��4���3���2��)                                                                    4ޠZ4ґ�4���4��3u�2�o�1� 0�K�                                                                    8P�8V\7��x7H��6�zV5�5|5J54!��                                                                    6(\w6> 5���5e\4��b4g�3 T�28��                                                                    4��4�1-4�%�4�4sXd4-��3�j�3��y                                                                    3�{�3���3��93t"�3D��3H�2��*2�8�                                                                    4��	4��4μ:4���4��4T-�4�s3��>                                                                                                                                                                        '�m      ~�     ��6���-h��F�P�    =%�T            6�{{4s�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ߢ@    @ߪ     16:55:41        F�� @ߪ     @߱�    '��      ~�     � 5��7D�[B'�fB'�f7��D"g�    '��U+�X7�7?�  ?�j3H�3��~,�0ɏꯪ��.���    -��)��1>��>fP�>��=�ӥ=3̉8���1pns2/4�:�H<'�t<���=H=�\=��=�Ԑ=ش�=��g=�A�>�<>V�@�>�@m]@ ��?��?��?>S�?%�>@@Y@i�@n�N@d5�@d>@c�s@c�#@cQ�@b�@b�T@b!@a��@a��A��<���{@��                    E�1"7A�G���F�dFR{;E��eE�SD1i[CCfBA��                                                AlCԯB~�'A���A>��@�_K?�.}>�
\=ˊ�                                                                    EEl(@ʌ#F#r�A��4� Fa��B5��A(�-�#�                                                    {@��@`1�@`1�@tɱ8"�y���!!Gr��}?/3դ%	G�?��'��U@����]�3X:`    >��?WF2�>|�P�/.���l'        B%�1B%�1C��C��C�&j>�!"?4�3"��C�-�6���6���A8ةA�.�>��A��S@�sA��@�w2@fmYA�a�@�Y�A�E��>�T|#�/���>    ��>��t@(��<�L�<�#@<�+?��8@�� ?=�? ��@9�O?2�dF���4� 1�kP.���4�a7fM3�XH5�`F��?GH%Gh	 G�݇?,�    >E,Q>#��>k�>�R8>��>�\�>�;�>γY>�[u>�m68��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    @�k@�k{@��7�P%9M�y9@�/8��8$�x7�9�6ޯX6(�5?^�                                                A�q#?4�;?4�;{@��?G��A*R�BA��@p�E=��B{@��{@��    C��BO�6r;�    6[�k6��@�#�3#ö7!�`    >��PC��R7���Bq�Da��CU�QB<"JA�1�A@@�0�?Ўs>�                                                ?�Af�W@r�?pT,?	�T>��=�2=#r<5*                                                                    DA-<F��XE�@�D� �C�"�C;��B�"�A���@���                                                @��oB��VB7��A@@�\R@8�?x3H>��=�;                                                                    C�#VE��dE�Cڬ�C5�NB�H'A���AvJ@�4                                                @'�PBVD�A�x"@�06@��?�k">��>%9�=3��                                                                    7@T�7�FPA(�3>`CjB��?:��?�w[*��t).Xd,s��+���*��t)hu�+�,S�_,2��)u�&7���1���/�j*7"_+�Hp+�4�&���%�3=�1�b/��/���*6�E7�%�1����%�4�{    (F��1��1��<,�Ֆ    <Ӊ��V�5خ0@��/�6,                        +��0i�0_
-B��    0@�z2��3            2�V�:�+m,�ߏ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B+^G�\4>*��A�`                                                B�                                          B��                    A�                @�p�    -��?)p�1�D6#��4�b4�>��6��J>��U5�5D;���;�>i7v�&            5�3�    4�|4���    7�{@��1pe�3d'9    ���D{@��{@��{@��8�88f�            7e�P    7�7�{@��    7�7�    2�2�    {@��7�{@��6��Y6~܃6>�/    3�!�{@��5���4�R�8' �8' �<��    B�,�6��_FXg�C|
?��V>RA}�p            ?<�1���9��K1���<�|�<��4@��B@&�3?�sD>�ψ>Vr�=������
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/?9b�<�y�<��R@�`�@&��?�+>�q>5��=��p����
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/<UW�0��*J5�RR6<Fp<H�K;���;l�:<�(=rF<�JL                                                �)f��?���4��%̮�5��ᲜÖh���<̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 1W)x1W)x>sj3aT�3aT�@%@QQ?��Z@!ǐ= �s2�x2�xB.yE6���;j�>�1<���=O�6���E1� 5���6��CF��D�EQD��yD9B�C��HCP�BM,�A���@��2                                                AƵ�CZEQCF�yB�B�B&�HA�P�@�,�@��?�2                                                                    E��OG��GS�F�t�F9=�E�Q�Dɩ�C��C `A                                                C�3\EO�EA�gE�hD���C�H\C�,B-)iA:��                                                                    E�zF�KjF��F{-�F]/F.��E�6RE��dEG�                                                 D)��D�VkD�|D���D�ƮD~n^D0�C�z�C�?F                                                                    0@�z7DK�GX��GM��Gj�F��Fz�F��E�	EP
�                                                >�׊;@�\                                                                        @���A���A�.LB:BE+B�0B�WsB�b<�<�<�<�<�<�<�<�<�<�<�<�E�n�E�sVEc��EZ.D��DW&�C�'C�c                                                {@��{@��{@��{@��{@��D�	BW�-%c��e�87븽CV�W@<�<    A �AA �A{@�ξ�s��sC��t{@��C��PC��GC��t{@��@t��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�o:PLD�J+G��D�JD�# D�# A߻A߻Fd1�BAb;>�X�C��!FK�<FK�<D��D��FedBA]                ?NeUC�0�C�T�C���?�TC���C���C���C� 0C�1�C�^�C���C��zC��^C�:^C�o�C���C��:C���C�XC�4&C�EiC�G�C�9TC��C��$C���C��kC�B�C�.~C��C��C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��GC�UE>��>��8>��>�}�>��>�"�>�~�>�>�X>�9�>��>ȪA>��>̩�>��
>��:>��,>��>ȅ�>źs@J�@8�        @�|���/��% ����f|��������\�Ħ�
ƾ|1ƾ|9ƾ|1�M�0        =|-Aȗ+    ?ᔪ?ᔪ{@��@�QFU�p719C7e�        ;���?MmM?O\�?Q��?T��?Z0�?c� ?r\�?�                                                  �l�k6���A�  ?E]?��@���'��      ~�     � 9�8~��7�g7G�6���5�#�4�Xx4 U                                                                    8=:87��G7��6{vF5�Z�54�3"�                                                                    G���F�dFR{;E��eE�SD1i[CCfBA��                                                4�=#4��3l�K2��n2Շ1[zF0���/���                                                                    3ƞ3(�2���2��1C�t0��/���.�"�                                                                    72l�7&�/6��/6'�?5�p�4���4U$3K�&                                                                    7C8I6�.�6f/t5��5v�:4Կ�4aS3/�                                                                    4�p74�,�4#[�3��3IO2��2��1�                                                                    7n�!7��6��86��5���5�46�e3U�                                                                    8�@�8�A�7�&w7+n�6�=�5�f5D�4"h�                                                                    5m�"5�>�4�R�4Z0T3�W3D}�2�T31���                                                                    8�@�8�A�7�&w7+n�6�=�5�f5D�4"h�                                                                    7؉�7�<�7{6��<5�l59��4x�v3��o                                                                    4��74��(4R�U3�I�3S<2�\�2R�1�'                                                                    7؉�7�<�7{6��<5�l59��4x�v3��o                                                                    �@����E5���5���4�bL4E�)3�2�\�                                                                    .Th�.�@q;D�:��0:.A9�J�8��8\�                                                                    5��<5.�p4��y4�3�S�2���2(��1F�                                                                    7c�6�3m62<&5�(�4��4/��3__2sm#                                                                    �N��?������8+M���<��ݏ��a�F                                                                    ������E�����ᲊ��H�v�"B��%1                                                                    6>�5�6�50l�4�bN3�t3.ȯ2^i�1spb                                                                    �-t�Z��xu��
�E�X���%�������                                                                    ���J��˯�����k;z��K��)���xG���V                                                                    �[�q�Q�H�����)���_�1�2�_�k�~ZF                                                                    )fe�)�V{/��./lPO.�g-�?�0�qp0Y�5                                                                    +��,"N�                                                                                            8c��8S�>7�O�76e�6�U#5��5fJ4:��                                                                    4�4��S4K�j3���3j�2w�z1�;:0�W�                                                                    7^��7ONa6�5j65 5�t4��y4#и3L�i                                                                    3�d�3љA3Gc�2���2"��1�
�0���/缀                                                                    8ZJ�8J��7�R73�U6���5�K
52;&4^�.                                                                    4ܴ�4�=/4C=�3�@�3(]2�-�1���0�!                                                                    7�ZK7�P�7��v7r�6t�
5�+e4�II3�`�                                                                    5�B�5��15�5*�n4��3�1�2��x2�                                                                    6�Ԛ6��V6e)�5�5E��4��3�U2�Sr                                                                    4��g4���4��P4
�3a�(2���1�pa0�:�                                                                    87.7�*7�N76��6���5��5��4N                                                                    6~65�5P��4��3��e3��2.�                                                                    4�v4�u�4�@4���4`"�4"��3�P�3�=�                                                                    3}�83�H3wp>3^\n35�3��2�9�2r�@                                                                    4��W4��-4� 4�)4���4G0
4�3���                                                                                                                                                                        '��      ~�     � 6N�V(F��                        6��j3�Ԑ@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ߪ     @߱�    16:55:44        F�� @߱�    @߹@    'ޕ      ~�     ��5@��7��B'�fB'�f7U/MD"g�    &Lt+��6�A?�  ?��4�:3-�D2��00/�\/�b6        ,���(uΊ>�H=�@|=���=�s�=wH�<V@n5H�20�/:�<'��<�Ŗ=H�=��=��=��o=ش�=��^=�A�>�a>W,@N��?�??U4�>�K=���<*� <�kV=6g>@d�@n��@d5p@d!@c�[@c�@cQ�@b�@b�K@b!
@a��@a��{@��{@��{@��                    E�6��HG�p�F�-�FR+�E��;E�*D1I!CC1zBA��                                                AڈC��B~2�Aܚ�A>�R@�.W?��>�Օ=�7*                                                                    EEq�@ʔpF#w�A��l3c��Fa��B5JI@��U,�/                                                    {@��?��\?��\@r�(7�^�����"���&#��/c>~#�!�>��`&Lt��@#/[�#4Eb    >�C=���2[�@L_���g���b";�G    B�vB�vC�9�C�9�C��?-MF>���3E�C�x�6��6�f�A8=�@�SO>�C@��?��A"��@�?w�@���@2B�A��r�[7dA�"�2�[7d    �[7d����>�'�>�'(>�<�?��%>���@
l#>��>zE�?��?��F�D�3c��0�@,.{�3Z�}6Z    5v��F\�F���G
keG&�X?Ƌ)>�9[A3D?��?#A? �e?�,?�>�P'>�0�>�->⋒8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�UG    �>E�>E{@��7��	8�8�]8l�87��*7t��6�Uz5�45-�z                                                A�r-��R��R{@��B5i.@���B5/�@q�<���{@��{@��    C�KBO�\63�z    6p�)7F�@� �3;��6��g    >�,�C��7���BVrDa.CRv�B?�A�,xAB�Z@�q�?��0>⺌                                                ?�hAfV�@pi�?v�?�!>��E> R=,++<2].                                                                    D@W�F���E�O�D�u�C��C;�~B���A� q@�S�                                                @���B�J�B5��AW�@��G@��?{H>��=�N+                                                                    C�x�EԋE?CُC5��B�`A�m%A��@��                                                @')�BU�A���@���@�J?��>���>&b2=1�                                                                    7 77}YA'��>^�B��?9}�?��*~�(���,�:+��*~�(�-)�$+��N+�Ȓ).�7U/M                        ��03E�1ǵ�            7�0�����.4X`�                        =�V��I                                                                2�m�            2�m�{@��,пf?E��?D9h?1�T?9�=?`��?��?CK?z��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A���G��>,��A��                                                B�                                          B��                    A�                @�p�    ,.	�(,�    5���3r@73['(=s؊6G�k>W:4�u;R�;S��7�            5b�    3�L�3�A    5qȻ{@��2�;�4�3    ��3�{@��{@��{@��7�ΰ8@h�            7.̒    8��{@��    8��    4���0���{@��5���{@��7�J+5���7s��    5�<D{@�ε<
.2ޅ�8��8��<��{    B6�)Fv�;D<�>�>9�q@�7�            ?�^�    :.>/    =DHF<�/@�iF@�E@L�?@�M>�b=��+����
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/?��~=DG8<�*@���@�;�@J�?@��>�V�=�k]����
�>���5��1����劤�Ί��>�p6
�S�R�=3�*�`��/;^D|/N^)�#�5��5�'�<�`�:��8���6�{�:��<nk                                                �DO��p�Ëê����eîV��~xZ�#�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 4�@A3��Q?� 5�!4�}A?�1$?�E>?t�@:n=�F6u�4�Q�E;�h7���=t�-A3G6\[�=ʆB7���H�"�49��7�b(CF@�D��D��=D8>�C�'�C�BL�zA�@�\B                                                A�@�CZ�CF�=B�>�B&'�A��@̡z@?\B                                                                    E��>G�'Ga�F�cF9�E�8�DɒpC���C TB                                                C�)EO�9EB?Ev�D���C�#�C�B-A:�                                                                    E�F�M�F���F{0�F]�F.�E�5�E���EG�i                                                D)�lD�Y�D�eD���D��<D~m�D0�TC�y�C�>                                                                        7��GX��GM�G^�Fͱ�Fz�F��E��EP\                                                A
"oA2�FAn\@b�p                                                                @7#�AN�A��KB6BY�$B�gB�s*B�b�<�<�<�<�<�<�<�<�<�<�<�<�E�w9E�}cEc�EM�D��DW �C�"yC�a                                                {@��{@��{@��{@��{@��C��KBg�)+۠��fg�7��OC)�@ -�    @�5�@�5�{@�ξz�W�z�WC�!{@��C���C���C�!{@��@r��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G.O=ߌ�D�B�G��D�B�D���D���ASASFc�B?ROB�[C�p
FK�-FK�-D��~D��~Fdf�B?M"                >��C��qC��C�?��C�
@C��C��C��KC�7�C�nJC���C��C�!�C�c�C��C�ϳC�ZC�3HC�`�C��&C��~C��5C��C��RC�ݔC���C���C�P�C�/�C�5C�
C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���Cۚ�>�g�>��>���>���>�xt>���>�!�>��>��:>�ք>��z>��f>�gW>�
�>�x�>�ck>ǊZ>�Ɯ>�U>Ū�@6b@��        @3��ƽ	�ƽ	ƽ{�����I���D6��@��D�ƼpEƼpGƼpE���        =M!�A���    �T��T�{@��@���FU��7<Q7iA        ;	�C>�]�?	�?�!?Y��?h�w?pK?v��?~�                                                ¨��5���A�  >�g=�g@b<'ޕ      ~�     ��8��p8#�+7���7�*6I5|�F4��3�Җ                                                                    7�7N��6��C6=s(5�;�4���3���3�l                                                                    G�p�F�-�FR+�E��;E�*D1I!CC1zBA��                                                48|�3��Z3 (�2�l2��1�@0P�/uk�                                                                    3i	I2�/2JNZ1��f1)90'�//�.� �                                                                    6щ6ֿO6j�5��5x*G4�=�3��39fv                                                                    6�?q6��6 �5��_5X�B4�3�3���3'�                                                                    4�z4L�,3ۡ-3�.A31�G2e��1���1y                                                                    7
ߩ6�o�6<:5�5�]4�?k3ʷ`3>�!                                                                    8+��87G�7���7O6vwd5��44�x�4�                                                                    5b?5&64���4$2�3���2��2Zl1��                                                                    8+��87G�7���7O6vwd5��44�x�4�                                                                    7~;k7�D�6ӉY6I"�5�|	4�5�4�Z3z�`                                                                    4Ug4��&4	�3���37h2i�L1��11�                                                                    7~;k7�D�6ӉY6I"�5�|	4�5�4�Z3z�`                                                                    ��t�n�}5�v�5M��4��3��3#��2���                                                                    7���7/�;pGX;u$:��9��9 \�8�i                                                                    5"��4ܻ�4F1O3�}3�<2��O1���1/�?                                                                    6�_�6��5��A5X��4�
�3�w2���2\��                                                                    ��������l���
♴������ֲ� M���                                                                    �cqw�{q������ٲoO���]��jL����                                                                    5���5�*M4�+�4V0�3� +2�1�1\d�                                                                    ������ƺ�'u����*�v�=H4�lbL���C                                                                    �&gt�)�������1	����ǰ�&(�ᕯ��                                                                    �)����� ?�d�J��<Ա�=V��\Ӱh�                                                                    )y�)�y�0�.�9$,��*�=.�]L05�                                                                    +B�+� a                                                                                            8��8/P7���7	FE6~�R5���4�F�4)��                                                                    4�/4���4	%�3��3s2�16�80�=�                                                                    7�7h�6��V6K5��4���3��}3:�                                                                    3�Q�3��3Z�2��2��1�0HO�/�ʰ                                                                    8 &A8��7��7O�6�@�5���4�fJ4Jw)                                                                    4��
4�`4�3���3�z2$]1Y�y0�Wg                                                                    7~J&7��$7?�>6�g]6S��5S�`4n�_3�c�                                                                    5�N�5�ۗ5[`G5 �Z4r$n3r% 2�T62�                                                                    6M|�6V��6I5�$�5+644+6�3@��2���                                                                    4jׁ4um<41F
3�*D3C��2C�1\Tb0��#                                                                    7�f47�gH7j�K7	�+6�z�5�z�4���4�A                                                                    5��`5���5�+5l�4���3��+2���2��                                                                    48��4L
�4OW�4O��4B1�3�5m3m��3�i                                                                    3$�3$��3'��3'�W3��2�\N2@2]��                                                                    4a�54ybX4}kp4}��4mYp3�i3�73��d                                                                                                                                                                        'ޕ      ~�     ��6c4\    F�%     =�    >	�i    6�l4	}�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @߱�    @߹@    16:55:47        