CDF      
      time       bnds      lndgrid       levsoi        levdcmp       cft       glc_nec    
   ltype      	   natpft        levlak     
   nvegwcs       string_length         levgrnd       hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      Sun Jan  9 16:23:28 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1927.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1927.nc
created on 12/10/21 16:49:42    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:28 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1927.nc had following "history" attribute:
created on 12/10/21 16:49:42
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�fF�� @�r@    @�z     &
9      m�     ��3=N?6��B'�fB'�f7ItD"g�    
�P{+� �6��?�  ?��4V�0�G 1��/�_).��m        ���&�8�ӹ2�+<��/���6�g=��=35I�:)��<#"�<�-=G�}=��A=���=��=ؕ�=�_=�>�>/8��2Q�+<^�+!��+S.���6#=�v@Wn�@`T3@V�'@WT�@W�P@XO�@X��@Y�@Y/h@YGz@YR�@YWu{@��{@��{@��                    E��}4�5�G��G F^t�E�)�E&D8�&CIc�BH                                                A�EC P�B�3_A��AJY@���?���>��=�z�                                                                    EA�f@�BRF O*A��y1��{F]1�B/�<�D/-$                                                    {@�ξ��3���3@f�q7]���oq��o|�7=�/`������>��s
�P{�L
�.��U2f�    >�y�=ﷷ1�D|>�)^�I�� ��<`$M    A��A��C�P%C�P%C�7?�  >�Q�3��C�/6:^�6�T�A05AL�>�y�A<��@&�qAc�*@RŔ?�'BA �@��A������������9�H����    �������@6G%?|p%?|p%@Ej�?@Oo�?&�K?�>@O��?�^F�nF1��{/>�,ݑ�24P%6w��    3��F|�F�ZGzHGq�@ӧ@&��B�p@�e?�=�?O��?+	?>���>�&�>�Ri>�Fe8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�?�    �i���i��{@��5�;�            4h�6#u=6���5z�4���                                                A����o��o{@��Cp
"AHB5�g@c��=3��{@��{@��    C���BR��5�~j    6M�J7��V@�
V3�6��7    >˚C��u7�w�BH��D�4eC��B���B-�A�T@�מ@"��?�C                                                ?r*�A��@�V:?�?x�=>��>*��=�u<e                                                                    Dl�	F�SvF)�D�k�C�%�CM_�B�:�A��@��                                                @�Z�C̲Bu��A+�*@���@!�N?�؁>���=��C                                                                    C��E��kEṔC�ѧCF`ZB�RCA��iA#��@�                                                @H��Bt��A��@��>@��?���?w�>:�U=4Cy                                                                    5�X5��A�d>U�B�?1��?�)�	w(�l+��+Q�)�	w(�:�)�դ+�u+��I(�_�7It                        � �3��1�7�            7 ����{� �2���                        >��)��(                                                                2���            2���{@��,�B�>L��>L��>L��>L��>L��?K�?X�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A�"�G�#�>*\�A�-                                                B�                                          B��                    A�                @�p�    �����    3=N?2�24P%=�w�6_�>��4�;�;A@J;B՛                �C�j    4���1:N�    6�S6{@��3h�C4��n5̍Ӵ�/�{@��{@��{@��3F�8M^�            6Uz�    6Uz�{@��    6Uz�    2qo17
�{@��71{@��8(>�7m��7ٕq    3���{@�δ�)3Ic7/�7/�<��o    B�6��-Fo@yD+��>��>��A��            ?��~    :E�    >�_>Mh�@���@�Ow@�?� �>�l�>K���"� �����������C����􊀾��^�V�D.��/O��tV��?�`�>�_>Mh�@�\�@�V�@�?� �>�l�=���"� �����������C����􊀾��^�V�D.��/O��tV��=\    $��w/z^c1��/<�?�=|zi8�0�\�4柀=�<�                                                �!���C�������)��$�RĳN	Ĝ%�Ā`�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6��46�?���8%m>5��?�Ӿ@9�?��J@y%�=�7��S5yu�G J�8U8>�B�v�@&��>ڮL8U8J��(6��?8�@�CB�MD�E&D�i�D6��C�C
1}BO_�A��`@��                                                A«MCQE&C?i�B���B'A�1}@�_�@�`?�                                                                    E��=G.�G�aF���F>K�E�I�D�V�C�C7*                                                C�n�EJr2E>p�E�WD�e�C�C�B-��A;�%                                                                    E�ՀF��gF�$BFz~�F\�*F.�}E�<E���EG�                                                D)UtD�W�D�4�D�-�D���D~Z�D0��Cߓ�C�iw                                                                        54��GT��GJ��G�}F�#�F|�F�:E�#yEPQ                                                A�msA��-A��A�mLA�|�                                                            ?��@e��A�<A_��B"�BW0�B�(bB���<�<�<�<�<�<�<�<�<�<�<�<�E�Y9E�@�Ee�&E!D�Y�DW��C�McC���                                                {@��{@��{@��{@��{@��A�1?B�D�B�t�fn6�C �/@C^f    @k��@k��{@�ξ�Qu��QuC��T{@��C��C��`C��T{@��@f�({@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�Z>���D���Ga�D���D��+D��+A-��A-��F^:B;�B�w�C�	FKӋFKӋD� D� F_�,B;��                ?bKC�0�C��C���?_wC���C�nVC�nVC�s`C�{:C��	C���C��UC��C�"�C�W}C��C���C��.C�++C�g�C��!C��dC�<C�G�C�c3C�nSC�iC�8'C�
C��C��C�u�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��"C���>��l>���>�Q�>��>�t$>�D�>���>���>�?�>��>�r>�f�>���>��s>���>��>�Uo>��>���>�3�@>4@�        @j�ƿ�3ƿ�}ƿ�uĀ�&�ŏ��Ň��ł}Ğr�ƾJ]ƾJ_ƾJ\�r�        =Mg7A��    ��n���n�{@��@���FQy>6���7��        2�LG            <ɫR?*��?7͕?<�?D?�                                                ��y5��3A�  ?k�    <#�
&
9      m�     ��            3��(5'7!5m��4?� 3�@(                                                                                2�A�4S8*4�u3rAk2�Q
                                                                    G��G F^t�E�)�E&D8�&CIc�BH                                                            /��0�%�0��%/��A/3N�                                                                                .6��/�=W0 �.��.b~2                                                                                2]�4�74v�-3n�^3M�                                                                                2�f�4Sd4��3��q3��                                                                                0cD�20�b2��1�
�0���                                                                                2��4��4��83��3,�                                                                                3s2+5%��5�d�4u 3֛�                                                                                0���2��)2��{1��i1Spi                                                                                3s2+5%��5�d�4u 3֛�                                                                                2���4�[�4�;�3�J�36�O                                                                                0��1�%�2b�1o$�0�>�                                                                                2���4�[�4�;�3�J�36�O                                                                                1�:p3dy$3���2�6�2��p                                                                    9`i�8�};\�];1׮:���:�9M9R8甧                                                                                0�*�2u�g2���1Ɗ:1�S                                                                                1�{3�2�3̗�2�*�2!��                                                                                �v�K�%��nw��NNͱ��4                                                                                ���3�/h���`�M��o��                                                                                0�p�2�$�2��^1��D1!�<                                                                                ��վ4�,�0�'l�����                                                                                ���װd_���K�����D��                                                                                ��lu��J�ĕ��2�%�`                                                                                                                                                                        $ƦG'5�                                                                                                        3n�v5=�5{��4dts3�I                                                                                /�{?1�#�2�B1c�0�?�                                                                                2l��4$��4�6y3z�3�                                                                                .��0�G1K�0�/�3�                                                                                3k&R5*PO5�t4�Ev4cB                                                                                /��o1���241\G0�'�                                                                                3H>�5��5A��4'i�3��                                                                                1d�3�3]'�2?Tx1�H                                                                                2!�u3�{�4_C3H�2���                                                                                08�<1�h�22�1�-0���                                                                                3t�}5#�55l�H4L��3�7�                                                                                1�ڑ3;�3�&r2i�1�d                                                                                0�3M2���3���3&i�3B�D                                                                                /���1��2���2y�29                                                                                0��3��3�.�3Kd�3m��                                                                                                                                                                        &
9      m�     ��6@�                    :Jg�    6��63�Z@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�r@    @�z     16:49:42        F� @�z     @ہ     &
�      n     ��3��6�*B'�fB'�f7.�uD"g�    
��+���6L��?�  ?�g4u�W1u�,2�zj/���.��    .���֣�ϡ�8ۓ�2�0��6���=�E~=�|=��5H�:1�9<#gb<�h=G��=��=���=�­=ؗ*=�Y=��>�>/j8�e�2l, ?�2X2>9T�1޼�5=~x=��@S�@`.<@V�@WT�@W�)@XO�@X��@Y{@Y/U@YGn@YR�@YWm{@��{@��{@��                    E�k�5/��G���G ��F_�9E��E�CD9�CI�qBH1�                                                AgC!�B��
A�)<AKT@�;�?��>�]v=ўY                                                                    EA��@�&�F 6IA���1�28F\��B.��<�D/-$��                                                    {@���-Vl�-Vl@ci7RJ���5�����z/JHp%h�>�њ
���O��/oA3S=9    >�y�>���1�Bv>e>���ֿ�Ǔ<#�
    BE�BE�C�YlC�YlCo9�?�  >��n3�kC�Z)6&6�szA./�A��">�y�B+AQ̨Bn�qAaaA��B*� AB;'B>;
�.�o��:�    �.�o    �.�o��<��?}p�?}p�A.��@5ڐAG��?�L@<�"AS�@]�`F��1�28/���,��)27�6�1S    4^�F�CZF�?$F�p�Gm?���@� Cl�,@�ow?�;�?O��?+	?>�B�>���>���>�Fl8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M���    ������{@��6            5F��7��6�5�V�5��                                                A�	��j_�j_{@��C��Av50BE��@\��=��{@��{@��    C�'�BRs�5�ň    6K(e7��3@��3+�6�E�    >v��C�F7���BZZ�D�C�z[B�U�BQ�fA� S@��S@,P�?W�                                                ?���A�Ÿ@�R�@�}?�T�? ��>-�=���<a��                                                                    Dn�F��F��D���D,�CSC1B��Aë�@��Z                                                @ݝ�C��Bz�WA7��@�0�@*>�?�m�>��[=�j                                                                    C�/NE��ES��C�b@COzvB�2�B �A'�N@ #�                                                @N�BvK�A��/@��@,m?�!�?K>AP�=6lj                                                                    5���6 �A��>O'B Mc?,�!?���)��8(Oav+�WF+'�)��8(�@�*���+���+���(Ƃ�7.�u                        ��d�3�k1��            6�d��@���d�3H1                        >�4趔W�                                                                2��            2��{@��,�:�>L��>L��>L��>L��?+F�?If�?-�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Byo�G���>,ʞA��                                                B�                                          B��                    A�                @�p�    �{	x1    3��2�27�=�&�6���>c'�5y��;ў;$�_                ���e        /if    2�ʕ{@��4B�V4[�Z6Ca�|e0{@��{@��{@��2?��8 ^            6��    6��{@��    6��    2��6���{@��3��`{@��8��6˖Q7�}�    3�q{@�ε�44��/�K�/�K�<��B    B±�6��FR�D2�?�y?>�^�A�Y�            ?��(    :UM�    ?=�o>� @��@�?@*()?�ʫ>�ϭ>n���"� �����������C����􊀾��^�V�D.��/O��tV��?�5?=�o>� @ą�@��@*'�?�ʫ>�Ϧ=�,T�"� �����������C����􊀾��^�V�D.��/O��tV��=MK    "���,�P�/��<�?�>%	N7�h�2���4f=~>�V                                                �!���C�������)�āy�ď[U�r�Z�<��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7  n4���?��c8��26FųAi@��AQ�A3�G>T83H�6$GzO�8$�?E�pCjk"@	��?G�8$�JҢ�4kEV832�CB�TD��D�!}D7�C�M�C
�/BP��A�q�@��t                                                A��TCQ�C?!}B��B'M�A��/@Ѐ�@q�?�t                                                                    E��GGҷF���F>n&E�X�D�R�CC2�                                                C�jHEJH�E>I�E�D�~�C�#�C�pB-�A;�                                                                    E��:F��WF��Fzy�F\�KF.�tE�<�E���EG�N                                                D)T�D�PD�.D�*D��ND~\LD0�xCߕjC�i�                                                                        5��4GTe�GJۘG�KF�5qF|�JF�SE�%�EPQo                                                A�N�A��=A�TXA㜃?�=�                                                            ?�@f�+AdAn�jB-�B^nDB�*B�"<�<�<�<�<�<�<�<�<�<�<�<�E�@�E�)rEe�YE#3D�b�DW��C�Q/C��	                                                {@��{@��{@��{@��{@��A�.�B�F�.��f{6���C
 �@6��    @
�V@
�V{@�ξ$�T�$�TC��{@��C���C���C��{@��@ci�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�>�yD���G]�D��mD��1D��1A2��A2��F]+KB:)�B��iC� FK��FK��D�D�F_�B:$h                ?�-C��KC�u�C��?�JC��yC���C���C���C���C���C��mC��FC��2C�C�3�C�YZC��C��iC��~C��C�PKC��C�ǪC��cC�$6C�?bC�R<C�;�C��C�C��C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���D��>���>���>�,>�MP>��>�i�>�:�>�]J>�lZ>�W�>���>���>�q�>�Q>���>��>�Q�>�T�>��>�]O@��?�5�        @�QU�±�¬��©�Ė����x���e���Y����ƿ"�ƿ"�ƿ"�ć�g        =b�A�]    ��?���?�{@��@YGsFQS�6���7�T        0�H�            =���?<1�?@/?F6?O��                                                ��I6#3�A�  ?k�    <#�
&
�      n     ��            4d1V6�%5u�4D�'3��}                                                                                3�)5C��4�&$3x�h2�K�                                                                    G���G ��F_�9E��E�CD9�CI�qBH1�                                                            /�[1�B�1 �-/��/<j.                                                                                /��0���0"~�/�.m�c                                                                                36$�5J14�H3uS�3	V�                                                                                3� �5Y��4�j{3��x3j                                                                                1b34�=2���1���1 }Q                                                                                3��-5��4�IA3�9�33��                                                                                4O��6� 5�=!4�%3��c                                                                                1�\�3~u�2��2 �C1`��                                                                                4O��6� 5�=!4�%3��c                                                                                3�r]5��-4�{43ۂW3AWs                                                                                1��2���2o�1}�0�?                                                                                3�r]5��-4�{43ۂW3AWs                                                                                2TR:4G6�3��2��,2�1�                                                                    9��f9[�;`-�;4
X:���:A#9p��9��                                                                                1��:3~F�2�8W1��@1%��                                                                                2�+�4{��3�5G2��(2*�                                                                                 �K7�u3�vT��S0���y                                                                                ���Y��d��7{�Q�S�{��                                                                                1�ي3z��2���1��L1+�                                                                                ��d�Ɨg�3~�,q���l�                                                                                �}��T/:����UůO��                                                                                ���z�{[��O���-ڌ                                                                                                                                                                        "�4%�'                                                                                                        4D�36��5��"4kN3�I�                                                                                0�wJ2�_2�1E0��y                                                                                3C.q5�4�Sx3�:3	8                                                                                /��1��c1�u0!/�o                                                                                4AƯ6=�5���4�[r4KW                                                                                0ƕ�2�f2��1�0��                                                                                4$�d5�ny5G�4+_�3�#!                                                                                2<v)4c�3drY2C�1Õ�                                                                                3As4���4!�'3
{�2�J�                                                                                1J�2�28�v1Db0�                                                                                4I��65J5tO�4Qt�3�*�                                                                                2fW�4,�03��62oa1��                                                                                1�v3܋�3���3*V�3LOz                                                                                0nRK2�82�V�2	��2%p                                                                                1�;4�=3�\�3P0�3y�y                                                                                                                                                                        &
�      n     ��6>��                            6���3�9@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�z     @ہ     16:49:44        F�F @ہ     @ۈ�    &      n#     ��3��r6jR�B'�fB'�f7T��D"g�    ث+�\�6�L?�  ?��4��0˰X2� 6/��@.�3-        1��,x�8��C2	��0�M�7aU�=�U=�;�=.�5ە0:8G<#�;<��=G�=��=��5=��R=ؘs=�i=��>��>/�8�$U2�F,�O�3V9�:C2K�1�ot9��H@Oۆ@_��@V�{@WT�@W�@XOT@X�c@Y_@Y/D@YGc@YR�@YWgA�ax;řK{@��                    E��5|�mG�jG9�F`��E���E�LD97HCJJjBHe�                                                A`C!�JB�fuA�B�AL<�@�C�?��2>��'=�֐                                                                    EA��@��F �A���2�d�F\�lB-�q<�D/-��                                                    {@��@'�3@'�3@_A�7��%�S*%4Ě�-�,//�8%��>�Mث?���/c"�3��X    >�y�@+c21��ٿ���A&�ɽ�6</��    B2�B2�C�>#C�>#C���?�  ?¸�3IrC��5�H6�$A+�KB�M>�y�B�V�A�K�B�PPA@�A���B���A��Bcכ@�^��h*�L*�@�^    @�^A]�A��J?}A-?}A-A��/@�فA��@�@�?kA�#=@�	NF���2�d�0�@�,�2=>u7�F�3��48�<F�3G�GqưG�K?�V�Af��C���@��/?�4�?O��?+�?�*>���>�\�>��O>�6�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�p    �T��T�{@��6D��4���    5�6�s�7^��6���5{�h4ϋD                                                A��7�$L��$L�{@��C�,�A�X�Bh�@@Quq>\�W{@��{@��    C��LBR2�5��    6H7�P @��A3
*7sZ    >�L�C��"7�a&Bi>D���C�'B�aMBp��A��1@�z@3�?�n                                                ?�"�A���@��@"@V?�kf>���>.V?=�K0<cLY                                                                    Dp��F�p�FpD�-�Db{CS��B�i�A�n�@�Ɍ                                                @�=	CCB�AC��@�\�@*�.?���>��H=���                                                                    C� �E��EV��Dj�CW��B��B�<A+m@!�z                                                @R�HBw�A��?@���@7*6?���?[�>G=�=9P{                                                                    5���63	A�>H[A��?&�5?�-*��t)'��,9v+�L!*��t)_��+�,J��,��)@�.7S�11��/��<+}�            62�34U1��/�{/��<+}��2���r62�3�\    %|%�1�wi1��-�ō    >���6�X�5)A
/߮�.��Q,	�                        )#K�/���/�R�-�H    /߳�2�A�            2յ�:/-2 >L��>L��>L��>L��?��?G�N?\T?O��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B��G�a�>+w�A�6�                                                B�                                          B��                    A�                @�p�    ���պ2�&3��r2�K�2=>u>���7f[z?*z+6�c�;n�l;r�_                5�c�    5�M�2�    8F�{@��2��5UO8eĵ gn{@��{@��{@��4
�G8Da�            8F�U    8F�U{@��    8F�U    �)�z7�aD{@��8ONH{@��7��n7SU�7n{Z    �"}e{@��5�5��w7ɗ[7ɗ[<�Y2    B�]q6�MhF0D�DL@O!s?���Bw>�            ?�K1��:i>�1��?���?��@��@��-@S�?Ǣ�>��M>~{ɋ"� �����������C����􊀾��^�V�D.��/O��tV��?��1?���?��@� C@�9[@S��?Ǣ�>��M=���"� �����������C����􊀾��^�V�D.��/O��tV��=A�E    *�F�6"�-_Of=ح>>7��z3Jd�2f7�>�}                                                �!���C������Q^ľ,ī#zēbJ�n��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7.�K4u�u?�Ç8�]�6e1�A�5W@�^�A�1A��>C�8~q�5�e<G���8
��?t��C�zA��U?w9�8
��K5�8H$`8r��CC�D��D�ٚD7pgC�ÃC�CBQ��A�a�@��A                                                A��CP�C>ٚB�pgB'ÃA��C@Ѩ�@a�?�A                                                                    E���G
��G�,F��4F>��E�Q�D�PUCC/�                                                C�aGEJ{E>#)E��D��tC��C��B-�*A;�                                                                    E���F��IF�Fzt�F\��F.߶E�=�E���EG��                                                D)S�D�I$D�'aD�&GD���D~\�D0�CߗC�j                                                                    /߳�5�zHGTE!GJ��G�F�F�F|��F��E�(HEPR4                                                A��A��wA�ޏA�V�                                                                ?i�@n��AL6A�*�B)@�BX�	B�MFB�s <�<�<�<�<�<�<�<�<�<�<�<�E�(_E�9Ee�oE0&D�gDW�HC�UYC���                                                {@��{@��{@��{@��{@��A�,�B�`�����fB/6�hB��@-�h    @©)@©){@�ξ�����C�z�{@��C��C���C�z�{@��@_B�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�V>���D���GW�D��bD��D��A7K�A7K�F\�2B8�gB���C��FK�VFK�VD��D��F^nB8|G                @ ��C�A3C�v�C�x�?tC�j�C���C���C���C���C��3C��AC��cC��ZC�/C�'C�FC�h�C���C���C��C��C�QmC���C���C��<C�eC�5C�9iC��C�WC��C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��D3�>�n�>�p�>�t>�w>�o�>���>�Ik>�>���>�Z�>�#�>��>�D�>���>��.>���>��>�п>��3>��4@�^?έw        A@�������۬b�ە����y����_��7��v6'ƿkrƿksƿkrČ�y        =���A�}�    �KVB�KVB{@��@S�FQ"6TO�7i~        0x+�:��    <c?Sg?5��?9��??T?G
�                                                ��s�7E!�A�  ?k�    <#�
&      n#     ��3�@�    4�X5��6`o	5l'�4;¥3��7                                                                    2��    3EǗ4���5��d4�&s3m+�2��m                                                                    G�jG9�F`��E���E�LD97HCJJjBHe�                                                /
�t    /��1F�1�0�W#/Ġ�/w                                                                    ./	�    .�6�0z!�1z�07.�_5.9!                                                                    1�J�    2Ⴎ4��p5R��4vЦ3j��2��i                                                                    3ZZ    3Iqg55��I4�x�3��K2�z�                                                                    0�1�    1��2�`�3q[P2��01�B�0��&                                                                    3<�R    3v5E5$bv5�'�4�>3��3�_                                                                    4>�`    4 �25���6de5��4x��3��\                                                                    1�4    1:t�2��3�82�N�1���11X                                                                    4>�`    4 �25���6de5��4x��3��\                                                                    3���    3y5�5�@#4��3շ3��                                                                    0���    0��2o��32��2i�h1xc�0��"                                                                    3���    3y5�5�@#4��3շ3��                                                                    ��xW    1/�3��m4��k3�?32��x2dR�                                                                    :�x9���;bp�;.�N:�:c�9���9��                                                                    1X��    1}�3#M�3���2�	1��1�                                                                    2�mY    2���4�4�W3�U�2�dr2��                                                                    �W��    ��!���z�\غ�l���I���R                                                                    ���    ��O��޸�D���>��G�p�C�I                                                                    1�(�    1���3�3��2�)E1���1�8                                                                    �^��    ���2�L���dش-T�%[U���                                                                    ��2(    �!R�Ѻ��X���۫��N��"\                                                                    ��Ă    �>�_�
?���d���k���̰��                                                                                                                                                                        -���$~	                                                                                            2��    4��5��6X��5|��4a�;3�B0                                                                    /6��    0��2&�J2�w�2�@0���0]-�                                                                    1�    2���4�uj5`,4��%3wg�2�k                                                                    .3i    /�Gu1%wM1�1��0�/r�                                                                    2�l�    3�y`5�K�6g�J5��4���3��                                                                    //W�    0{8�2$FQ2��U2��1tK0��>                                                                    2.�p    3��=5��636*5?�m4#3�3��                                                                    0G�    1�Yk3���4L�03[_X2:�&1�ߵ                                                                    1g�    2��c4[��5�D4�3�Q2V�g                                                                    /!�_    0�6(2{\r3%�r21EJ1�]0us�                                                                    2U�l    3�<5�T�6[	P5j�M4Gw�3�k�                                                                    0tn3    2��3��4zS�3��2c��1���                                                                    /�Y    0�@�2�;L4��3�I�3"3�3��                                                                    -���    /�451č3*2�t�2c2 6|                                                                    /�4    0�N�3�=4CE�3�v�3F?3A�                                                                                                                                                                        &      n#     ��6;�i%|%�                        6���3�ݍ@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ہ     @ۈ�    16:49:46        F܂ @ۈ�    @ې@    &e      nA     �04A��6b��B'�fB'�f7r;D"g�    ���+�?6��?�  ?1�4�NX1��12΃�08S2/�4v0���1+%��!�*8�%2��7<�g=n.�=�s;=�B8<�*�5?�}:=S<#̈́<�ǯ=G��=��S=��P=���=ؙ�=�i=��>�U>039��2�G�5�g;��R9`�1⟎0�C�8Q)�@L?@_�@V�H@WT�@W��@XO!@X�9@YD@Y/4@YGZ@YR�@YWaA�G�;ŏ�{@��                    E���5��G�!�G�KFa��Eê�E��D9p�CJ�:BH�                                                A�dC"�IB��A�^�AL�0@�C&?�:x>�]=�F                                                                    EA�N@���F��A���2��dF\�B,��<�D//��                                                    {@��A���A���@Z�[7�RK%��$�暦�M���:��ݩ�>��l���@�}H/��`41''    >�y�AP�T1�1�@k�A��m<�C�    B�F1B�F1C��fC��fC~S�?w��@Zk�3��C��5�[>6�lhA)�1C~p>�y�C$�B��MCJd�A��&BL��C&�A�H�B�<XA|����:� �Z A|��    A|��BK�B(�?yU�?yU�A�]3A3U�B�@b��A�MAդ$@��&F��U2��d0rNr-��2�,�8+�K5Ɇ�4��rGK<G1�G���G���>��A��gCl�@�.�?�M?O��?*��>��`>�D�>��}>�P&>��_8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�=�    ?��?��{@��6��7,F�5��6{��7�ݜ7YJ�6wU&5pr:4�aG                                                A��[�Au�Au{@��C�^�B$�B��$@DG2>�O�{@��{@��    C�E	BQ65��*    6Du�7��-@���3��7�j    =��C��J7�-�Bv�%D�>�C���C
�B���A���@�:�@;!m?�)                                                ?���A���A��@7+?�KU>��g>0G�=��~<ps                                                                    Drz�F���F* D��gD
�CS_�B���A�=z@�(@                                                @暱C�'B��AO�z@���@*W??��>��=�TY                                                                    C��IE�EyEY��D	!:C^B�B�[�B,�A/;@%Ub                                                @WBynpA��@�d.@?�`?�;�?
�.>ME=>�8                                                                    6��6�dA�w>Bu�A���?"�?���+�v)���,q��+�o:+�v)ߠ�,��j,���,dq�)��x7X�54+��2�-��            7��K4P3(2 r:2)j2�-�����K����7��K3~�    '�ʅ4!�U4!�,/�R�    =���7���7e�2 �*0���.@�                        +v�21x1��Z/TҸ    2 �2�ƃ            4>m�:?�->L��>L��>L��?	+?y;�?<{?�f?:
�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CB��G��>*VKA�P                                                B�                                          B��                    A�                @�p�    $�(A �ZD4>�N4L��2���2�,�?6��7�GR?���7]*�;b�I;a_�                7K    5{�D3�m�    8�KS{@��2!n�6ڤ8�㉴�L�{@��{@��{@��4 ߮7!��            8Ơ�    8Ơ�{@��    8Ơ�    5��7���{@��8ڣ7{@��6w3�6��5�5    5�6l{@��6-�06���6�dg6�dg<��    B���6���F�(C���A~�@<a`B��            ?̚�4+��:h�A4+��?=Q�?8Mm@�=�@��U@h�y?ֶ�? b>�CQ�"� �����������C����􊀾��^�V�D.��/O��tV��?�ۡ?=N�?8LP@�E�@���@h�?ֶ�? a>,��"� �����������C����􊀾��^�V�D.��/O��tV��=�    .Y� 8@<)7�j>�%9�&�7��3��3�-=��                                                �!���DU4��z�Ċ��ĹG.ħ�Đ )�hd�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7L�5K��>,<9�?'7��B`rA#��B3��A���>E-�8��<6��.G{��6��@?-Cb]]A��?�r6��@K8kT8к�7*G�CC�;D�C�D���D7�BC��C�cBRֲA�h�@�^�                                                Aê;CQC�C>��B��BB)�A��c@�ֲ@h�?^�                                                                    E��gG
�LG��F��*F>��E�I�D�OqC��C1                                                 C�VgEI��E=��E��D��C�5C��B-��A;�]                                                                    E�ѐF��F�UFzosF\��F.��E�>�E��EG�4                                                D)R�D�A�D� |D�"�D��OD~\�D0��Cߘ�C�k                                                                    2 �65�GT*�GJ��G�F�WRF|��F�E�+EPS�                                                A��A��5A��@+�                                                                ??:@�3�A@��A�VB)�VBYq�B���B���<�<�<�<�<�<�<�<�<�<�<�<�E�aE���Ee�DE=�D�i�DW�-C�Z	C��                                                {@��{@��{@��{@��{@��A�-B���$m<�f��7wKB�n�@+Z�    @⽎@⽎{@�ξW�ѾW��C���{@��C�87C��bC���{@��@Z�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��>��D��6G[qD��D���D���A;�0A;�0F\v�B7a�B�#�C�kFK�$FK�$D��D��F]��B7\�                @G�uC��JC���C�w�? C�_OC��2C��2C��rC���C��rC��C��fC��/C��SC��C�6$C�T�C�vAC���C�ÃC��C�%�C�Z�C���C���C��C��C�1�C��C��C��C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���D�v>�j�>�m�>�t?>��>���>�Ϝ>�!B>���>�)�>���>��>���>���>��>�i�>�3r>�"�>��>�͝>�*�@!p>?蹣        A����!�]�!K[�!<Ł�Ǒ�~Ǒh]Ǒ������#���#���#�ęA        >_��B$��    �j�F�j�F{@��@q&�FP�6M�'7A�        0�/6=���<d~C>2p�?4{E?6��?:fC??��?G�_                                                ���7�X�A�  ?k�    <#�
&e      nA     �06E�`4U�~5�.K6��6\!/5`2�41`�3{r�                                                                    5z 3��4��5��<5�|4��\3`�2��                                                                    G�!�G�KFa��Eê�E��D9p�CJ�:BH�                                                1�i
/���1�2.�1�n0�Қ/��s/�;                                                                    1��/R�0>�f1\ٹ1�a0O.��.&K<                                                                    4T� 2�%�4O65�U�5O��4j�3^u2��>                                                                    5�j3�'�4�z�5�
o5���4�Ư3���2�u                                                                    30y,1���2���3��w3Y�S2��1�ݪ0�v�                                                                    5��4&w4��N615�.�4�H3�/�3�n                                                                    6��'5 ��5�4�6�!m6_��5~�.4n�D3�i�                                                                    3���2 M2��:3��3�[�2潻1�j1#ag                                                                    6��'5 ��5�4�6�!m6_��5~�.4n�D3�i�                                                                    6�_4q4쪵5�3�5�~c4���3��I3O                                                                    3q�1�u�2*�(3\{:3.�]2aC]1qb0�Q�                                                                    6�_4q4쪵5�3�5�~c4���3��I3O                                                                    �,�$���"2L�54�=4��d3��u2�,2J,B                                                                    9�4 9ѤQ;M�;��;D=:tE�9�+p9�                                                                    3��>26036�4��3���2�yv1҇�0�sv                                                                    5 ]�3���4�^5��4��3�(�2�e�1�"L                                                                    � ��µ��T����|�X�5�`K�=�ұ��-                                                                    ��l���,�`���[�A��� ��<b��/bG                                                                    4O�2���3m�4{J3��V2�.71���0�6                                                                    �&մ��Q5����5�j���%E�������                                                                    ����Ez$��T
��B8����� M���b��                                                                    ��z� �~���ܲ�[�����y��Eٯ�0                                                                                                                                                                        0� �/��                                                                                            5�k�4�85n)�6�l�6V��5p��4VAF3���                                                                    2�0�K1��-3<2�B~2��0�&0GRv                                                                    4��G3O;4iO�5�dH5^y4|�a3j��2��k                                                                    1 �/��0��2�`1�;1	�0�/Z�                                                                    5���4=5duF6�[�6e~P5�j4�83��
                                                                    2V�0���1��3�2�3�2��1�0mɭ                                                                    5h3� t5,�6o��6/��56�4��3ngZ                                                                    3	B1���3D�74��44H�3PH2/Δ1�;                                                                    3��2[=4
5A�534l2���2@�D                                                                    1���0���2��3];{3"M�2(W1�0\+�                                                                    5b�3� 5RL6�d�6V��5^uJ4<	3���                                                                    36'�1ܳm3pV�4�N�4u|J3~<�2V�
1��                                                                    1���0v��25ӕ3���4�&3�*o3޽3Q�                                                                    0�d�/GC`1�52��2�Ei2�]1�1��                                                                    1��0��Z2^;}4µ4?�J3�PO3:�Y3-�B                                                                                                                                                                        &e      nA     �068C�'�ʅ                        6�E�3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ۈ�    @ې@    16:49:49        F�� @ې@    @ۘ     &�      n`     � 5���6�w�B'�fB'�f7�*�D"g�    (A�+xN6��?�  ?{��4��3�k�2&�-1�G2`~10NN    .;�I*9�r>n
>(_=�d�=���=�Iu=E<j�5�:@V�<#�B<���=G�b=��4=��0=��8=ؚ�=�N=��>Ʊ>0�@�	I@(ʵ?�ՙ?}?�?D'? [?��>��@UB@\�@S�B@TP�@T�@U�@V8$@V��@WhR@W�"@X6�@XgA�;<��b{@��                    E���7)�G�v�G& FbC>Eô�E�?D9�HCK&BI(w                                                A �C"�B�_/A�$AM�@�>�?�h�>��Q=Ҥ�                                                                    EA^�@��AFЫA���3oVF[��B,�@��h/��W                                                    {@��BA'�BA'�@WQ8(�ޤD#�Mç��/Y��&=6�>�7(A�A4�N/kr�4*��    >��BA�@1Νl@�AA�r�A7�>h�    Bd�Bd�C��8C��8C�'>��@�̐3�LC�M�6d�O6���A(LAC9�l>��BCN�gB��yC^=�A�׹B�'�C:رB�B���BYc�����    BYc�    BYc�ByW�AXg�>9D,>9v�A��A4K�Aېj@d�H@m��A4�@X�F�l3oV1{p.6e�3�o�8�5L6�q�5�xG%<Gk^G��PHd�>YuA�YAM!,@e,>ӌ#?%$?y�>��	>���>���>�d >�P8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�    @�N�@�N�{@��7���96��97�8��*8)7�7��6�s�5�5
�g                                                A���@z;@z;{@��A�֤Bf��B�q@9�d>�g�{@��{@��    C�v/BP8,6g�    6A��7���@�i�3ִ7#��    >	BC���7��Bv��D��hC�aSC&QB_��A��U@��@:$?��                                                ?�A�F�A�@7��?�S;>ʢ�>.p�=�%p<t�                                                                    DsSF��F�D��D
-?CR�^B��A�r@Ë�                                                @�.JC�B��AVW@�d\@)�R?���>Ԇ�=�=d                                                                    C���E��1EZ��DG�C]JMB��OB=A1n�@'�m                                                @X��By��A�g�@ˡL@>.�?��~?5j>P�`=B�k                                                                    70��7���A
�L>?gA�"_?.�?��h+��*A9L,�M�,=�>+��*���-g-'�,��@)��7l�58�2�h�0 ��            8��5Hr2*��2�n2�h�0 ������ގv8��4�L    ,��u5
Ul5)[2K    =�\<8l�58A��3ׁ,2�|1�`                        /1B�3�3�.2'88    3���2�2            5��;T-?�D?P`�?B3�?8��?Sr?xX?dd�?R�*?d4l?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CXA"G��m>+ϽA�                                                B�                                          B��                    A�                @�p�    -�cv)�,5!��61�?3�z�3��u?�=8A��@E&�8 ��;�f�;�$G    +��        7�ۛ    5z��6(o�    8.y{@��    4�^w7�#ϴ�M{@��{@��{@��7T�7��d            81U    81U{@��    81U    6K�5|g�{@��8* {@��4S��4)�&3(��    66�
{@��6��W7J/7���7���<���    B��6��GEꊶD��Aĵ<?�BQC!�            ?C�x59M9ޏ59M<�D�=�/�@�@3�?�U?fm>�X�>c��"� �����������C����򊀾��^�T�D.��/O��tT��?:��<�Bh=�%�@�T@�?���?^3>�<�=��V�"� �����������C����򊀾��^�T�D.��/O��tT��=X�    -��~5��7�&;�^�;>��;Y=M =�=�9�                                                �,..�0��(�Lę��Ď݋�}���S���Do̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6A�5\$P;���8���85�QB��A/�pBC��A�̤=��7|��6���EX�%4H�<��AH4O>;��<�]4H�J.چ8Q5��=CE`�D�EcD�F�D:;eC�~#C�BT�A���@�G�                                                A�`�CREcC?F�B�;eB,~#A��@��@��?G�                                                                    E���G
�YG��F�ӊF>��E�A�D�M{C�rC2?                                                C�L�EI�;E=��E��D��aC��C�qB-�SA;��                                                                    E��(F���F��Fzl_F\�=F.�E�?|E��EG�Q                                                D)RLD�>�D��D� ED���D~]ED0�ZCߚ7C�k�                                                                    3���7�HGT(GJ��G��F�b�F|��F��E�-�EPU�                                                @OJ,@��@Պ�                                                                    @�h�Al4A�L+A��B3<Be�B�^9B�@;<�<�<�<�<�<�<�<�<�<�<�<�E�jE���Ee��EI�D�k4DW��C�^�C��b                                                {@��{@��{@��{@��{@��C�$MBB9-��q�f(7��C5E�@X�    A<)A<){@�ξf>��f>�C�K;{@��C��C�q�C�K;{@��@WQ,{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��GT=�@4D��tG�D��KD�-�D�-�A<��A<��F]�B9�gAy��C���FK�eFK�eD��D��F^��B9�G                @��3C���C��C�tm?   C���C��C��C���C�=C��C��C���C��C���C��lC���C��C��>C���C���C��zC��C�7�C�iLC��C���C���C�&�C�WC�C��C�z�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��6C��C>�=�>�'�>��~>�4�>��L>��>��F>�zb>��>�t�>���>�^>���>�U�>�:<>���>�1[>��>��e>��@+��?��x        A�q�m���m!W�l�nŪĻ��H8��v�����BZƿ|qƿ|�ƿ|pĎЪ        >��BY��    @K @K {@��@q��FP��6��;7%�        ;H�?<,[?.s]?2Q6?@�?Bƌ?G�?O??],5                                                ��g7w]/A�  ?k�@��@�Ȏ&�      n`     � 9�#8T��7�
7=�V6��5�t�4�)�3��                                                                    87զ7�J�6�56o��5��4�4�3�j�2�|"                                                                    G�v�G& FbC>Eô�E�?D9�HCK&BI(w                                                4�~3�ˈ336S2�Z2	U�1*(400�/@��                                                                    3��K3�q2b_�1�b�1-z0V�k/^u$.sr�                                                                    7�6��B6x9�6��5y4���3�(�3J8                                                                    7S?@7�_6�P�6u��5��;4��4-��3!��                                                                    4�J�4�D"4�X4?��3o��2�PU2��1[                                                                    7�n73o<7��6�!�5���5
$�4TGC3E�"                                                                    8}{y8VT�7���74o6��05�|d4�=v3�E                                                                    5P6�5Ay�4���4���3���3(��2l�1s��                                                                    8}{y8VT�7���74o6��05�|d4�=v3�E                                                                    7�?�7���7AZ6�&�5عA53X4E�3NJ                                                                    4� 4�n�4NC>3��q3Op�2��1��V0�X                                                                    7�?�7���7AZ6�&�5عA53X4E�3NJ                                                                    �p|��/�4�2Z5*�|4�
�4$3#x�2���                                                                    1���8q��:�+G:�/#:e��:!�9`�49�<                                                                    5�l�5D��5�/4�Ų3��3�|2C�c16�7                                                                    6�t�6�%�6�5���4�p#4Q'3);�23`d                                                                    ��
���Q�z�+�(����H��_>��iױ֪�                                                                    ��+��*\���v����fײ#Y��Fn��&                                                                    5��5��`5Y�4�Ô3Ӄ�3d2*�p1454                                                                    ��֮��@��0�?����,�}�p����6R����                                                                    �x���Ei~���#�X������u�*&�V�,                                                                    �B����m���"�����ԧI��5�扰0�K                                                                                                                                                                        +�P/��                                                                                            8G�c8|�7�T�7'��6���5���4���4g0                                                                    4��M4�=�4�D3��J3	��2>M�1h-?0���                                                                    7C{{7B"6�h�6&\5�H44�HE3�3�                                                                    3Ť�3���3�Y2�}2g�1G��0~�k/��e                                                                    8?j�8j7�}37%)a6�¤5���4�4`'                                                                    4��x4���4m3�B�30�2Q(O1�}�0��1                                                                    7��G7�/�7K�}7X�6QT�5���4���3�2r                                                                    5��-5��:5h�F5��4o<23��{2�cZ1�                                                                    6�ǩ6~	�6$��5�[^5)'�4T�k3kK�2���                                                                    4��.4�*44<�3��"3ARB2s^11�t}0���                                                                    7��7��7x�}7&�6�D5�
�4��]3��R                                                                    6�)5ۏ�5�@G5>s�4�33�2�\�1�R�                                                                    4�%4uU�4W8�4s��4:��3�43��v3O��                                                                    3e��3F?�3-�3D۸3��2×F2iͭ2(b                                                                    4���4��(4��24���4d14�g3��X3~2a                                                                                                                                                                        &�      n`     � 65� ,��u                        6��o3�s�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ې@    @ۘ     16:49:51        F�� @ۘ     @۟�    &-      n~     ��5��E7�?B'�fB'�f8#&D"g�    (Ħ�+V�7cR?�  ?s��0��3�r�)M߻/h|���_�        /�+�>��>M�%=�!w=�ǒ=(��:��z3|�}2^��:<<#�W<��0=G��=���=��1=��=؛=�=�	5>��>0�@�.@yժ@:�?��u?��!?|H�?T�?=^@_>�@_>�@U{]@Uw@Uv8@UzQ@U�<@U�6@U�3@U��@U�@U��>!J=K6?
p{                    E��p7~"�G��tG�CFa�SE�=	EQRD9D}CJ��BH��                                                A�C"Z+B��A�|{AL�M@��?��>�b�=�g�                                                                    EAkE@���FۈA���1���F[�GB-A6��/M��                                                    {@��B�}�B�}�@\��8���%�*    �R
�/H'����~?   (Ħ�ATp_��֛1�%Y    >���B:��1�k�A3�*AhDA�Y        B��yB��yC�&�C�&�C�m�>��P@��3��C�Ch7/O�6���A*G�CK�D>���C_zB�;aCnf�A���B�-�CJؓB�B�StBbzs?�ֵ    Bbzs    BbzsB]�;qV�    6M�A���A;��A�R�@u�u@ �(@�]t?�کF�6�1���//-y,mh1q��8���6��5�\)G��GR�`G��-H ��=��W            >;i>��~>�T>��^>�OI>��W>�4@>�
$8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���MΰM    A#�hA#�h{@��7�a-9���9i �8�a�8TE�7���6��=6%��5=8g                                                A���A��A��{@��9��^B���B�H@:
x?
��{@��{@��    C��BM�_6��    6B{�7��E@���2ǘ`7B�a    =�<�C�m7�>�B^|�D��JC���B�SB,�TAiI�@�d�@�-?2
                                                ?�ϑA��A�@*�?q'�>�o)>g�=pJ<V�                                                                    Dqf�F��KF�ED��hD"�CO%�B��.A�W?@�K�                                                @㺖C��B�ĻAL<�@��$@$@)?�&>���=�f�                                                                    C���E���EX�:D;CWAaB�3�B"A,Uy@$~                                                @T%@BwxgA��@���@5��?�&�?ǁ>I^==�>                                                                    7}�7�k�A��>A�A�S?!{�?��,/@*�l�-X,��,/@+�.`Z-�s�->�*?��7�W4���2B�f/��.��.��(�i�7�w4Ր\2Q�2HF�2A	G/穊��w1��W7�wU4��    3I��4�j�4�3i2&��    =�ߨ8Z�A8@"M6�J75Ȑ2��                        5U�6�L�6��3��}3�?�7���2�}�            4̲>Z@-MsS?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Ck
�G�>+#rA��t                                                B�                                          B��                    A�                @�p�    .�=�*�z�5TBs6�z�2��J1�j�?��C8W�|@|�D8&x�;�ŷ;��5�W�            7��"    2AOC6�@�    2=��{@��            �jFU{@��{@��{@��7�hX8O�A            7Uݯ    7rH�{@��    7rH�    6��L    {@��2=��{@��2=��2%0Čh    6��L{@��6�"�7��d8Y�f8Y�f<���2	Q�B��z7;yE��PC�v�A�9�    C0�            <�wr4��/7�P4��/<�+<��=Z�f=4�<��<u�<Jކ<% ��݊�4����W��*������;��4��\�ˊBq&�-�����E�<��<��<���=ZKC=3k�<�!�<`��<�8;����݊�4����W��*������;��4��\�ˊBq&�-�����E�;ܼ-v:N)��6 t6qW�9#"C9 ��8˲�:���;���;��j                                                ď�č#Ć���{]c�`���=��
�ò��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         9��^        {@��{@��{@��{@��                                                    CG��DұND��D?�kC��XC�BU@�A�ZF@�A                                                AǺ�CR�NC?�B��kB/�XA��@�@�@
ZF?A                                                                    E���G
�%G��F���F>��E�,bD�9
C�pC+�                                                C�I�EJME=��E�@D��C��wC��B-�A;�8                                                                    E��4F���F��Fzo�F\�>F.��E�>�E��/EG�Y                                                D)SD�EPD�"?D�"�D��ED~\�D0�Cߘ�C�k*                                                                    7���7k�jGTF�GJ�0G��F�c�F|��F��E�. EPT�                                                                                                                                @jCHArq�A͞�B#{B8
UBmC,B�x�B�y�<�<�<�<�<�<�<�<�<�<�<�<�E�-�E�Ee��END�fDW��C�_�C���                                                {@��{@��{@��{@��{@��D]B��.S?#�fz�8uC~��@|�    AOy\AOy\{@�ξ}6�}6C�x�{@��C��kC��tC�x�{@��@\��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��;���D��lG�BD��CD�!kD�!kA7@�A7@�F_G#B<�(    C���FK�jFK�jD��D��F_�9B<�                @���C�v�C�D�C���?   C���C��kC��kC���C�;VC��EC��zC�^C��C���C�e(C�&WC��C���C��'C�a/C�H�C�@�C�J,C�b|C��WC���C���C��C�C��C��C�}�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C��>���>�!�>��>�ݱ>��>ό�>�y>�f�>���>���>�l�>�A�>���>��.>�%1>��;>�6�>���>�x>��
@5�@ �        B�ǚ��Ǚ��Ǚ���⍽�%G����vݧƿ݋ƿݦƿ݋ēU�        ?	�BC49    @ػ2@ػ2{@��@�{�FP�7a^�70�        ;�!6?A9?A��?CQ@?E��?I�?P
�?ZCr?l6z                                                �*M73&A�  ?a[A��A z�&-      n~     ��9Sz8��Q8'�7u$M6���5㕹4��3��u                                                                    8��t7�e�766�ӏ5�
w5��4�a3 ��                                                                    G��tG�CFa�SE�=	EQRD9D}CJ��BH��                                                4ݖ�41ٜ3�	�3 i�24ҫ1n`�0�Х/�\q                                                                    4�w3`�2��2"4�1dhm0��/�=r.�t�                                                                    7d��7I:N6՘�6If/5�ѽ4�d4�3D�A                                                                    7���7T��7�6v#�5���5��4Vty3Dő                                                                    5.�i5,�4��%4>�3�-�2�B�2;<�1+��                                                                    7���7��q73�B6�k85�/�5'�4�J3p�                                                                    8���8�t�8��7c@6�U	5�p~5%� 4#��                                                                    5���5�Z�5*Lw4��|4��3e�2�T�1���                                                                    8���8�t�8��7c@6�U	5�p~5%� 4#��                                                                    8Mk7���7c;�6��i6��5T�Z4��23��                                                                    5	[�5	xI4���4�}3�_�2���2&L1!��                                                                    8Mk7���7c;�6��i6��5T�Z4��23��                                                                    �J�ް��5��5���5��4B~t3��b2�
D                                                                    5�5�E�7�7��{7_7!6倈6�~�                                                                    5��5��c5;4��3�	�3A�2q�[1]�P                                                                    7%47�6z��5��75
-�4Fw3v��2u�                                                                    �,)��'����ϵY:��� ���D�zP���                                                                    ��dl�������+-���#�QOɲ�[���z                                                                    6#�6��5y�~4��q4	W'3E��2xQ�1u�O                                                                    �3�o��.��w�
��e-���������$�                                                                    ��~���v�����PY���|�8K��x���J                                                                    ���m��l�.�w��ڳ5H�;�0�\ �u#�                                                                    (&+[(�eE+�y\+p�+��,+e�-�F.xL                                                                    +�e,�C                                                                                            8��p8~p\7�w�7\��6���5��_5+I44P:                                                                    5Q�5 �U4{6�3�-C36mR2�4q1��N0�?�                                                                    7��`7yA!6�g}6[5��85��4)�3E�I                                                                    4Qw3��3v�2��Z2<�
1��f0�|/���                                                                    8��Q8t�7�Wi7Y�26���6�57��4WW                                                                    5Q 4���4p�%3��r3C�2��Z1�U�0�L                                                                    8��7��27�m70F�6�s5��4ؠ�3�OP                                                                    6��6�5Ě�5Iu54���3�Cv2�� 2	�.                                                                    6�ɖ6�{y6�j6q�5_�4�`�3��2��r                                                                    5z4�4��U4"�R3~�2���1�t0���                                                                    8*�s8�-7�A�7Wr�6��5��B5b:4w�                                                                    6C%6/��5�Kd5v9�4��D43K�2(��                                                                    4�>�4�o74�ע4��'4vS�4)��3�H�3�a                                                                    3�xh3��"3��{3~��3G�3	@�2���2h�                                                                    4�L�4�&4�@q4���4���4O��4�3��                                                                                                                                                                        &-      n~     ��66T�3I��            H́    =�?�6��H3٧@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ۘ     @۟�    16:49:54        F�: @۟�    @ۧ@    &�      n�     �p5� 7��B'�fB'�f86ND"g�    (��+"�17S��?�  ?w�1��647o)�`�/��Y��{        .�}-*脐>�^�>!�=�#t="�;�;g7F�q/�82e��:.�=<#	3<�D=G��=��?=��d=�ǹ=؛|=�5=�	�>�/>1@˶p@meU@7�?��?��/?Up�?%
�>�Z�@_ov@_��@U�@U��@U�q@U�x@U��@U�/@U�2@U��@U��@U�=���=�5?��R                    E���7��xG��GSF`�`E�l�E
�cD8��CJWBHW�                                                AΣC!v�B�E�A�i1AK�.@�by?�\>Ӥ=��&                                                                    EA��@� F 	WA��Q1��F\��B0VA4A.p.                                                    {@��B��B��@vY�8�)b%�i�    ��o�/Jd��~F?   (��Ay�Z���?2gA    >��BAE2M�@��AGLAF��        Bt2�Bt2�C��2C��2C�x�>�@@�]3��C��P7Rd6�?�A0�)C7j >��CIkB�ČCLϧA��EB�]0C+N�B�B��BE�����    BE�    BE�B)!P            A��A0O�A�ބ@j��?��@�V�?��F�4�1��.��+� /1!<18�z�6���6V<G��G8H��G��=>*N�            >C�L>�o`>�7B>�%N>ɑ�>��>t>ؐB8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�    A�^A�^{@��8�9�A9�La9t�8���7��7�6E�Q5fB�                                                A���A�� A�� {@�ΠC	B�!>B�bZ@R%y?�.{@��{@��    C�*�BLxK6��-    6E��7���@�ʒ2�m�7�Q7    >*N�C��7�p�BDA�D�f�C��kB�&�A��}A>�m@�D?�
v>�uh                                                ?b�FA���@�?�K�?08�>�ף=�=<Y<2j�                                                                    Dn�{F���F��D�6�D��CI"3B��A�y�@�h�                                                @��CW�B}	�A=��@�q�@y{?��>���=�sR                                                                    C�JbE�2EUL$D;�CN�B�C�A�|WA$}@��                                                @MWlBt^A�"�@��A@)@�?�b�? t�>=
=4к                                                                    7���8��A��>H	�A���?&��?�C�,�"<+���-�k-0+�,�"<+�y.̺.>=-��y*���7׋�4��2Q�,*��0*s�0*iL*(�7
�s5}1�K�2F~2Ft, ��
�s1ѹ�7
��4�A�    3��4�xo4�n�.��v    :��b8:��85��6��%5ZS�2P�<                        5mx6��26���3��W3�U7ZCL2Չ�            4�Ԝ>C�-h?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C[ղG���>+�A��*                                                B�                                          B��                    A�                @�p�    .�k�*�,�5X@6��2��1;��?�I�8L�q@{�I8$;N<	�<j?7#            7��v    3,��6W^�        {@��            ���{@��{@��{@��7�%8��            7w1    7�M&{@��    7�M&    6T��    {@��    {@��                6T��{@��6�@�7��8 N8 N<�@�1�A�B��;7�E��D"'�AR�    C*;?            <4��4�|6��H4�|=0�<��{=9�'= �4<�;�Ud;H�:�|�ֻ���Ǌ߿͊�_Ҋ�K��*��{�ǊY���?���+a��?�P�<4x=-�<��1=9?�= Z#<~�3;�7�;G�8:�a�ֻ���Ǌ߿͊�_Ҋ�K��*��{�ǊY���?���+a��?�P�8&T�+Z��)���6N6��90�9�8��h8�87�m�7�-                                                �tԣ�p���d\K�Se��:sP�jv��m�S��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �C	        {@��{@��{@��{@��                                                    CH�2D��xD�I�DB�dC�3`C{�BT��A��P@�H�                                                A��2CR�xC@I�B��dB13`A�{�@���@
�P?H�                                                                    E���G_G�KF���F>]uE�	�D��C��C�                                                C�F�EJF\E>*�E�D�r�C۱Cm�B-�aA;��                                                                    E�ӣF��F��Fzu�F\ýF.�6E�<pE��oEG�                                                D)TD�PD�+BD�'$D��[D~[�D0�#CߖsC�i�                                                                    7ZCL7�y�GTr�GJ�xG F�UgF|�F��E�+<EPQ�                                                                                                                                @t��Az�8A�#�BwQB<��Bt�_B��B�=<�<�<�<�<�<�<�<�<�<�<�<�E�NCE�"mEfEE=D�X.DW�^C�[�C���                                                {@��{@��{@��{@��{@��DɑBMM.=�u�f�r8>��C�@p     A�.xA�.x{@�ξm�m�C��o{@��C��!C��!C��o{@��@vY�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�;ŻmD���G�{D��D�VcD�VcA/��A/��F`�iB?[�    C��FK�QFK�QD��D��Fa NB?V�                @��C���C���C��?   C�1�C��!C��!C�:C���C�c�C���C��eC��C���C�*C�˟C�l�C��C��C�iC��C��C���C���C��fC��3C�ЯC�wC��C��C��C�(�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�6YC�j?�?Qb?��?��? c�>���>�c�>�(�>�v>��p>Ђ�>�>�	>���>�f�>���>�Y>�%�>��>��
@4C?�F�        A��hǒ�ǒǑ����\f�F��
���
@��Z"�ƿlƿzƿk�ze        ?�B��    Af7�Af7�{@��@Ch�FQ�7w� 7]"        ;K��?F<k?F�r?H�1?K�?P$�?W�?d��?|:d                                                ���7ъA�  ?<��@�^$A z�&�      n�     �p9�a�8�I\851�7��6���6Mg5��41_                                                                    8�j8�p7d�6� !6�519Y4B�L3Gғ                                                                    G��GSF`�`E�l�E
�cD8��CJWBHW�                                                5��4a��3���3 .�2`�1���0�D�/��Q                                                                    430�3�z2���2JU�1���0��//˵;.�A�                                                                    7�z�7�y�7[N6}�f5̽�5�T4E,3t��                                                                    7�6�7w�7�X6d�%5�UY5
{{4P:43L{�                                                                    5S��5 �4��c4/3�V.2�S�25�126P                                                                    7�B�7�Bz7!��6��5̄�5)A�4~�3y�                                                                    8��'8ԁ�857'7�Q6�6{5D�K4E��                                                                    5���5�\5KR	4� �4 ӫ3�~�2�"81Â�                                                                    8��'8ԁ�857'7�Q6�6{5D�K4E��                                                                    8/?D8��7�ua6��6'�O5}� 4�H$3�+�                                                                    5.u�5.I�4�bG42��3��P3�2Abd1?�{                                                                    8/?D8��7�ua6��6'�O5}� 4�H$3�+�                                                                    �0u2rp�5�ʥ5�s�5/.k4�Mw3�4�3	jY                                                                    5��5��7�`$7|�`7�6l�5��`5Wc                                                                    5���5�9�5*6�4�3�T�3�,2j�[1f�^                                                                    7U�7<�h6���5�~o5'�(4n�3���2���                                                                    �]�H�E�ݶ�굇�A��`��\�%2�9)�                                                                    �$m^�Uس~��8����&���P�$��/�                                                                    6SJ:6;d�5���4���4&j?3m�`2��1���                                                                    �gM��H+߶�F@�/S
��P����ƴ��$�x                                                                    ��܅��b��Fq岳��J�cp
��������                                                                    ���y��4��\����.?��h9ٱ������                                                                    )��S)��,be,0��+�==+=�:*�/�*:�                                                                    +�P�,��                                                                                            8��8�z�8Q7�Q�6�E6��5?+�4`�                                                                    5=ӑ5$F�4�%4�3b��2��l1،!0�|�                                                                    7��7�+Z76�R�5�l�5 ��4Q��3vW�                                                                    49�`4 ��3���3��2jq�1��0�q;0�{                                                                    8�T8���8�7�S�6��6( �5d74��                                                                    56/5�
4���4��3rR#2�8�2+*1̟                                                                    84�8 ��7��7\�z6���5��5�M4AN                                                                    6Mŉ67�P5�#5|!�4��4��332+�Y                                                                    7~�7Ҫ6�sM62F]5��J4���3�}�2��&                                                                    5&G�5^y4ȃ�4K�!3�hU2�Q1�kP1
Ã                                                                    8\�8D[n8�27��6Ѥz6��5#��47�&                                                                    6{�6`h~6��5��4�g4��3;[2Q�P                                                                    5m�4�z�4墴4�^.4�3�4R,�4*%3�|�                                                                    3פk3�h03��h3�}<3w�3)�r2�7/2�
e                                                                    5#W5�5U54�:84�>�4�p�4"��3�_�                                                                                                                                                                        &�      n�     �p68��3��            GB�    =�Y6���3�-d@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @۟�    @ۧ@    16:49:57        F�x @ۧ@    @ۯ     &�      n�     �@5��7���B'�fB'�f8��D"g�    (�w�+GdH7@�?�  ?{��2*c.4-�)�N/�s���         .�ލ*���>�K>(Q�=��=4�5;���7�/�#�2X��:>�<"�<�=G��=��f=���=��*=ؚ�=��=�	�>�E>1'@̖�@oN@@˺?�`;?�ӯ?_��?3��>���@^.�@^q)@T��@T�c@T�&@T�:@T�@U�@U"@U<Q@UPz@U[�>.A�=��n?���                    E���7�l�G�G K�F_e�E�c�E
�D7�:CI�BG~Z                                                A��C k�B�o�A��AJ��@��u?�o5>ң;=��;                                                                    EA�-@�T�F ;�A���1Z��F]=�B1�OA4�b-4��                                                    {@��Bf:~Bf:~@|�^8����K�    �n/�/8�`��3?   (�w�A�s�ǐt2��    >��B#T2r�@Я|@�b	@��w        Bq�BBq�BCĘ�CĘ�C�[w>�3�@�P43)�C��7��6�%�A3�gC}�>��C%�6B���CE�1A�=bBI>�C!�,B��B�%wAډ��B    Aډ�    Aډ�A�q�            Av�-A.wAϪ�@Q-)?��p@��?���F�_1Z��.��,ۛ1h�8�;:6�]6&G��G*�gHߡG��><A            >Bc�>�By>��>��>�"e>���>��k>�<8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    Ax0#Ax0#{@��8 �9�b9��9�8x�e7�[�7;6Fp�5jC{                                                A��	A��$A��${@�Ο��OBB�B~@@_�>�B�{@��{@��    C�GBMNl6�r�    6Io�7�k�@��2�Ц7��b    ><AC��_7��mB-�UD��C�=,BO4	A�d�A	@���?��>��E                                                ?E/�A���@�]+?�,�?�>f�=ȶJ=c�<�                                                                    Dk�TF��+F6�D���C�KCBH�B�N�A��+@�=w                                                @�SOC��Bu]A.w�@�U@�[?vָ>�f=�v�                                                                    C�3YE���EP��C��!CEG�B���A�S�A��@�`                                                @E�=Bpy�A⦻@�q,@)�?��>��>/v�=)��                                                                    7��{84�A��>M�A���?+g�?��[,�H�+F2�-�XU-	,�H�+�!�.z��.��-�� *��87��T4�;�2f�2,,�03�03V*1_5���5 �1���2[r2[g�,#wx����1ӭp5���4�    3@
�4��4���.�:�    :��8��8��6i��5�2`�                        4��6E�>6B�3E��3$�o6��2��            4�U5>X�,�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CW��G�.>+I�A�                                                B�                                          B��                    A�                @�p�    .��s*�x5'1�6��2��)1�w?���84n@B�Y7��<�[<Z�6�o�            7�l    2��6/��        {@��            �g�!{@��{@��{@��7���84fj            7.{y    7s�s{@��    7s�s    6.��    {@��    {@��                6.��{@��6֭�7k�86��86��<���1���B�*r7r�F O�C�$�A�z    C�,            <s�4�=6ݚe4�=<�q�<��=�	�=3��<���<$G;�"�;�
������f�፰�������3��^u�}�f�[u�AT��,í�'�to<r��<�l&<�x=��>=3Ix<�ɧ<#�`;�Ui;�̋�����f�፰�������3��^u�}�f�[u�AT��,í�'�to8aAi+!��*"p65��6�pF9}�9*Y8�Py8>�F7�j�7e=�                                                �}/��y'��l�Q�\���D�q�"����DÂ9Z̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         ���O        {@��{@��{@��{@��                                                    CH��D���D��iDC�8C��$C/pBR�A�uw@��                                                Aȯ�CR��C@�iBÈ8B0�$A�/p@��@	uw?�                                                                    E���G>�G�kF���F>'�E���D��C�o>Cr                                                C�AjEJ��E>`�EȷD�K�C�nCCJ�B-h-A;��                                                                    E��+F���F�%8Fz}PF\ŏF.�FE�9�E��EG��                                                D)U7D�\�D�6#D�,�D���D~Z�D0�Cߒ�C�g1                                                                    6��7�;�GT��GJ�BG�sF�9(F|��F�yE�%vEPL�                                                                                                                                @r�rAyA��B�PB;�>Brk�B�ЬB�C<�<�<�<�<�<�<�<�<�<�<�<�E�scE�CEe�DE0D�A�DW�fC�R�C���                                                {@��{@��{@��{@��{@��DV�B~z.HN}�f��87�EC��]@��    A��KA��K{@�ξ)+�)+C�6{@��C�$C�$C�6{@��@|�^{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��;�r�D��.G�D��D�D�D�D�A(m�A(m�Fa��BBk�    C�m*FK��FK��D��D��FbBBf�                @�,_C�B{C���C�(?   C�7�C�$C�$C�ӿC��_C�u]C�DC�
C��=C�wUC�)�C��0C���C�3�C��C�vC�sC��2C�X(C�
C��C�ٟC�سC�yC��C�C��C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�:/C���?Og?	�j?Uk?S�?<G? �N>�^>��b>���>�GF>�;�>��>׏�>���>�Ly>��>���>���>�~�>�[�@bb?ǋi        A��Q�n���n6!�m��ś/���O���b/�����<vVƿVƿV$ƿVĀP        >�U�B��    ANO�ANO�{@��@L��FQW$7�%�7��        ;O��?Eb�?E�B?G��?J6?Na�?Uk!?aY�?u`�                                                �6�6�@�A�  ?R�F@ĮiA z�&�      n�     �@9|M�8��t81$7��6׳26�5��4$s�                                                                    8�Y�8�7_�~6�͑6;H55h4I��3O�f                                                                    G�G K�F_e�E�c�E
�D7�:CI�BG~Z                                                50�4WҺ3�~g3�2a��1�n�0��/�8�                                                                    4&�3�O22�N�2H�[1��C0�/�h.ي�                                                                    7�_�7w��7A>6|�05�!i5�74K�53~@                                                                    7�7C7X�|6��T6-n�5��/4�d94,�V33P�                                                                    5:��5f�4��4!W3c�2��P2,�1��                                                                    7���7�x�6�d�6S��5�IV5v#4SI�3[)�                                                                    8��K8Ϙ8,!<7��6���6�5C��4Fn                                                                    5��z5���59�U4�Oh4��3��2���1�M                                                                    8��K8Ϙ8,!<7��6���6�5C��4Fn                                                                    8%$�8�a7�b6��6#+35y34��<3��Z                                                                    5"C�5'̙4�C4$��3�p�2�#E2:��1;��                                                                    8%$�8�a7�b6��6#+35y34��<3��Z                                                                    �J�n��@�6!O5�A�58��4���3�wJ3�                                                                    4��-5i87�;)7��57=��6�U�6's�5��                                                                    5�35�4�V$4W��3��3�2C��1J�R                                                                    7I��78�J6���5��B5%��4nB#3��\2�m.                                                                    �PZ�>�M�_j����������+�-�A0�                                                                    �>����
?���t��V���/ϲ*�
����                                                                    6G�67�5�R�4�Oh4#��3lg"2��*1��                                                                    �X���A���pƶ.w�����ղδ�d�*�                                                                    ��_��.:�C�f���J����g����Ŋ��?�                                                                    ���B��E˴X���߄Q�0d��n�O���3��B                                                                    (���)Cm,7�+��?+��+4m*��*:h                                                                    +�
A,?�?                                                                                            8��8��:8@7���6��76�x5EHe4h�                                                                    52H5m�4���4�3c�B2��u1�x�1��                                                                    7���7���7�56���5�c�5#4XQ�3Ai                                                                    4.pm43j3���3�2ky�1�3�0��0��                                                                    8��N8�Oz8�)7��46�ˢ6+6�5kZ�4���                                                                    5*ϓ5� 4�H�4�3sc02���2La1J�                                                                    8)98д7�6"7\?6��5���5d�4�                                                                    6A:�60�5�5{t�4�q<4:�3N�23*Y                                                                    7�X6�4�6�J�61�)5�P4���3�H]2�]�                                                                    5$�5��4��4K2x3�[�2��2 ��1Ƿ                                                                    8N��8=7�8K�7�u�6�96�
5*^�4?��                                                                    6l+6X?�6��5���4�f4$�3B�f2Z��                                                                    4�T�4�)*4᧶4�\4�D44X'�4
�X3�^8                                                                    3�I�3���3�X�3�=�3z�3.��2�\2�h�                                                                    5�5`D5	�}4���4��#4�V4)s�3�a                                                                                                                                                                        &�      n�     �@6<y�3@
�                        6�l"3�w?@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ۧ@    @ۯ     16:49:59        Fݴ @ۯ     @۶�    &Y      n�     ��5�"�7��=B'�fB'�f7�eD"g�    (�Vh+\�7=�v?�  ?~is3�4�o+`�p/[_گ�!        .��1*���>���>H�'=��=�Y$= �9�1��)2B�U:�<!>|<�j�=G�\=��.=���=��+=ؙ}=�P=�	\>�U>1K@��!@j��@ R\?�ҹ?���??�d?��>�<�@]u8@^�@Ti@Tv@T�@T��@T�R@T��@T�@T��@T�@T�p?�i=o�q?�W                    E�1M7�g�G�P�F�H*F^RE�}�E	u�D6��CH+�BF�]                                                A#\C�"B��TA��2AI��@�?��6>ѧ�=��                                                                    EB�@Ɛ@F q3A�&�2�� F]�B2�A0��-�bD                                                    {@��A�\+A�\+@�%@8d3�%��?    �>�x/`T3��Z?   (�VhA�����m3(�C    >�T^A��T2�n@�����R�j�        BE�GBE�GC��fC��fC��=>���@`s�3 �C��\7�f6�+EA6�B��>�T^B�·B9C*�ANVA�ܑB��7A��}Bke�Aqg����]    Aqg�    Aqg�Aq�;�?�    6���A&�@�1CA��6@�?���@bȔ?g�}F���2�� 0	�-+�2�;{80��6�)5���F���GZ,G�T�G�Q>j��            >NT�>��>���>��'>��>��$>�e�>�қ8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    A/<A/<{@��7�c9���9v�A8�F8Z�K7��E7�I6:L�5`F�                                                A���A[KA[K{@��:c�Bz�Br2q@nef>�!�{@��{@��    C�5+BM��6�S3    6Mw�7���@���3
9^7Z�    >i�QC�ź7�IDB ��D���C�+�B*ԯA�0AA�@p��?��U>�	�                                                ?4��A�o�@���?U�>���>R]�=�}�=�[<��                                                                    Dh��F��DF
�D�6�C�fC=/�B��kA�P@�#-                                                @ϔ�CT�Bo*�A#��@��@
_?i�6>��,=��                                                                    C���E�U\EM\�C��_C>q�B��kA䤺Ak�@v�                                                @@@Bm��Aܡ�@���@��?�Or>� w>$B= ~�                                                                    7y��7��A�2>S�TBV�?0�\?��+�$�*�H�-a�,�ɥ+�$�*��-��E-kH-v�*+��7�iE4�U 2#j3*�3�/�Q/�6(��r�N54��P1�4u2�2�*�ŋ7N51�K��N�4��    1��/4�9�4�7�-d~�    :���7�W07�bF5ZN4��0��n                        3��>58I56�e1�TX1Ƀ�5�u�2��I            4��=�x,��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C�G�c�>*�3A��                                                B�                                          B��                    A�                @�p�    .l��*i	�4�'�6�e@3l22�U�?<_#7�t?���7f�D;��;�+7A��            7S�<    0�o5�C    2Ȅ{@��            ��=R{@��{@��{@��8u�8h�x            7�	�    7��L{@��    7��L    5�	
    {@��2Ȅ{@��2���2���1�G�    5�	
{@��6�p�6�J8k��8k��<��E0M.B��X6�s�F.cYDY<,@��E    B��            =���4�ĥ8@5�4�ĥ<��}<�"w>��>��>~=�[�=a<zg��&��j���B<�ʻ���Mw��,(�}j��[+��A3�,���ǊD�=�_�<��<��>��J>���>M�=�,6=+<y�C��&��j���B<�ʻ���Mw��,(�}j��[+��A3�,���ǊD�8�d9,�Y*(@5��6V�:h=9��9A`�8��8W��8$.                                                �X$��Tb�HC��9Q�"����RøY��;8�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         :c�        {@��{@��{@��{@��                2zj�                2zj�            CGˊD��7D��YDBB�C���CvmBO�=A�^
@�@�                                                A�ˊCR�7C@�YB�B�B.��A�vm@��=@^
?@�                                                                    E��Gc�G�F�x�F=��E��pD˴mC�D�C �                                                C�7aEJ�GE>�DE�
D�&UC�,FC&!B-ILA;�|                                                                    E��@F���F�+|Fz��F\�.F.�nE�6�E���EG�                                                D)V D�f�D�?AD�1HD���D~Y[D0�Cߏ�C�d�                                                                    5�u�7��]GT�GK�G�rF��F|f�F��E�XEPF�                                                                                                                                @��A���A�{HB
�BA�zB}�mB��_B�q<�<�<�<�<�<�<�<�<�<�<�<�E���E�]�Ee��ED�)�DW�dC�GnC��r                                                {@��{@��{@��{@��{@��D�@B��.�ѧen�8�C�$�@��    AVL�AVL�{@�ξ\ﹾ\�C�Ds{@��C��!C��/C�Ds{@��@�%@{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��Gʩ;�yD���G�hD���D���D���A#PA#PFbIBD�,    C��GFK�XFK�XD�D�Fb�BBD�                @-aC���C�	�C�]m?   C�c�C��!C��!C���C��YC���C��C��-C�ˑC�C���C���C�vGC�KC��C���C���C�3�C��_C���C�L�C�DC���C���C�jC�:C��C��1{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�c�C�ɫ>�J�>�Y>�f>�+Y>���>�2�>��>�T�>��>��>�X>ߋ�>���>�r�>�>>̟�>��>�">�(v>�c�@ [U?�        A��r�&���&A4�%���0�(ǡ�.ǡ$Ǡ�����8ƾYSƾY`ƾYS�^p        >�GLA�Ĉ    @��@��{@��@j��FQ�&7��y7��        ;n�=?K��?Lv?Nc�?QM�?VE
?^�?k�I?}��                                                ���6��A�  ?&D�@���A z�&Y      n�     ��9R�B8���8iS7��6�c�6�p5\4�                                                                    8�77�;Q7DO&6�®5�U�5*`"4A��3Jt                                                                    G�P�F�H*F^RE�}�E	u�D6��CH+�BF�]                                                4��490�3��'3x2Nʹ1�I�0���/��f                                                                    4lv3i��2Ͳ�23��1���0�xc/��.ӜK                                                                    7hH7U�6�6a��5���5ξ4B��3v\@                                                                    7�d�7-�6�w]6߬5nB4ʖ�4�R3*                                                                    54�Ʌ4=��3�l�3?��2���1��1��                                                                    7�%�7T&�6���6$�|5�{�4��42�3AO�                                                                    8��8�v�8��7h�6�z�6
�@55t�4;UB                                                                    5�w�5�Tt5�%4�!Y4|�3g�`2�ѫ1��J                                                                    8��8�v�8��7h�6�z�6
�@55t�4;UB                                                                    8&�88V7^�A6���63�5b��4���3�e�                                                                    5j�57]4�e	4Y73���2��2)�1-�J                                                                    8&�88V7^�A6���63�5b��4���3�e�                                                                    �O�x��^ 6DT5��m5,	O4� �3�!�3��                                                                    3�NT3�>}9[��9��8���8 Zn7���7��                                                                    5���5lCI4��E4(��3�Jk2�o2%Y,13�G                                                                    7*k7 �	6�5��5t4[L�3�@N2��                                                                    �/	ֶ$�s���s���x)�k�%���<��                                                                    �y� ,Ͳ�RQ��{"��⥲{Kx�$���ያ                                                                    6(�d6,�5{�e4�<�4�A3X�2��w1��                                                                    �6�&�H����������$���H�	�$��                                                                    ��{���#�+j���`��
�ݱW�]���{����                                                                    ��KA���8�>����!��`����t��w�                                                                    )mݰ)��-��-`l�,�o�,L��+���+��G                                                                    +��,(:�                                                                                            8�ə8�J�8	�7w)L6�n6��5;�4`��                                                                    5k�5��4�#�3�K�3N�o2�d1ԅ�0��                                                                    7�Ʋ7��#7�86ud�5�'�5!4M��3vm`                                                                    4`�4 W3�M�2�|12U��1��0�_0��                                                                    8���8��i8]7s��6��\6O=5_��4�r                                                                    5U56�4�x3���3]
m2���1��1��                                                                    8
8�+7���7E��6�n5�fY5��4	�                                                                    6"T�6��5���5a�)4�4�3�A3�2.�                                                                    6��6��e6��6�F5��4��3�b�2�U�                                                                    5-$4��O4���46�3�<�2�w�1�p�1U�                                                                    8-��8#V�7��7q��6���6L�5$e�4;                                                                    6Fgz6:�m6Gk5��4޲+4��3;�2UĤ                                                                    4� �4��4��4�J�4��4K�h4��3�R                                                                    3���3�.3���3�D'3f�V3$Ȓ2�Y�2��                                                                    5 ig4�?�4�O�4ذ�4�q�4y<)4#�]3߲d                                                                                                                                                                        &Y      n�     ��6@A�1��/                        6�A�3�&)@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ۯ     @۶�    16:50:02        F�� @۶�    @۾@    &�      n�     ΰ5���7s�:B'�fB'�f7��BD"g�    '�w{+���7!�"?�  ?�3̩�3�q/��1|�R���R/IZw    -ν)�%�>���>yK>%}=�K�=�d�<���64~j210�:�p<!�<�UJ=G��=��=�ͺ=�¨=ؘ�=��=�	N>�s>1@��j@d^�?�U?�	?XO>�n>	K>I�I@Y`�@]�2@TR�@T_3@Tr3@T��@T�\@T�{@T�J@T�y@T�l@T�@��<Ї�{@��                    E��J7Y�G��;F�c1F]��E��lE�$D6H,CG��BE�                                                A�}C�B�JVA�)pAH��@��s?���>���=�Q�                                                                    EBQJ@���F �zA�ST3��F^%�B3��A&��-琇                                                    {@��A%܅A%܅@�8.�)�����S�X�^/73%��~?W*'�w{@,H,�$3�QI    >��4@�@#2V�@����f��͈        BqFhBqFhC�CfC�CfC��>���?�3L*�C뒶6�Ϩ6��tA8M!BMn�>��4Bj��A���B�C�AtAb��Bw-�Ar�:BM,��Zb�?�X�$A��Zb�    �Zb��k��?�ݍ=�^�=��@�v]@WcAU�)?΢~?���@���?�;F��3��197.;�3��T7��5E�5�9�F��UF��Gn�G{�>���    ?��?n�>q��>�o:>�p{>��>���>�g�>�X:>�*8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    @:��@:��{@��7��9U��9I*8Ō�87Δ7��{6�16�P5BZ|                                                A��?�@N?�@N{@��?�,@A�etB`HR@y��>(��{@��{@��    C���BO��6���    6�A�79��@��T3j�7'��    >1�}Cɑ�7��B\WD�;�C�89B)� A�A��@��?��B>��                                                ?1�RA��@��$?V�l>��>_j�=�;�=��<�x                                                                    Dg��F��F��D���C��C;C	B�t+A���@���                                                @��C�vBk��A�@�e�@�*?e�>�� =���                                                                    C�o�E�#dEKLC�O�C;�FB��A��`A�x@Rd                                                @=�Bl�vA�c�@�K@D�?���>�5�> W0=<�                                                                    7L67��A" �>X��B ?4�?��*��)��3,�QH+��*��)���,4:,���,[R�)�3_7�S�3�571l2'��/�/�%�녷���4%�2 ��1b��1b��'Ȱ67��̰�#����4���    +�9�3���3�ԭ*6��    ;��6Jz�6��2���1?),,Z�                        .�a2k;�2j�-��M    2���2�.            3��\;1m,�|�?�  ?�  ?�  ?�  ?�  ?np�?Q�_?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B�c�G��>++�A��%                                                B�                                          B��                    A�                @�p�    -���)��3��X6So.3�2�3��`>���7���?E�F6�Z�;��;���7���            6��+    4c�.5��    7 ��{@��0w��3�A    ��{@��{@��{@��7�KS8�I            7&��    7�ٳ{@��    7�ٳ    5k��5*�{@��6��{@��6���6*�5�[    5��m{@��5�Ci6>k7�[
7�[
<���    B�>�6�ƲFFyD�U@	9�>���BD�7            >��3�5�9�kY3�5�<�֥<���@�e?ͣ�?POQ>���>=|���U��^���7��ʲ��DɊ�%�}^��[!^�A�,��o�>>�W�<��$<���@��?�]c?P7>���>�=|Y��U��^���7��ʲ��DɊ�%�}^��[!^�A�,��o�>:�X-��{)�z�5�'�6�;w��;J:he�9��88M�8Ph�                                                �j�$��z�����Iç���^��<�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 3�>�3�ݑ>؟6F
5�/�@�~�@��@� A��=���5!&�4�:_C�+}6��;�
q?�ֹ=��/<���6��G���6���6�AQCFӫD��D���D@�sC�V�C
��BM[OA�m�@���                                                A�ӫCR�C@��B��sB-V�A���@�[O@m�?��                                                                    E���G|%G�F�`uF=ɟE��
Dˇ�C��C ��                                                C�*-EJ�E>�JE�>D��C��C�B--A;o�                                                                    E��F��NF�/�Fz��F\�{F.��E�4�E��EG��                                                D)V�D�m�D�E�D�4�D���D~XoD0�Cߍ"C�b�                                                                    2���7S��GT��GK7�G��F��hF|K�F��E�TEPA�                                                >���                                                                            @���A�ńA�9tB�1BMN�B� �B��B�Ph<�<�<�<�<�<�<�<�<�<�<�<�E��E�pNEe��E�D�LDW|�C�=�C���                                                {@��{@��{@��{@��{@��D8�B#�Q-'Y�e�o7���C���@��>    A��A��{@�ξZ�ӾZ��C��{@��C�)�C�<FC��{@��@�/{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��;�}�D���G��D���D���D���A!A!Fa��BDs|>���C�^6FK�OFK�OD� D� Fb��BDn\                ?��C�C�C��?g;C���C�:JC�:JC��C���C���C�4mC�s�C���C��TC�C�?oC�[�C�oC�w6C�s,C�_�C�;(C��C���C��MC�]�C�"GC�$C��C�JC��C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C�c>�ă>�)%>���>��>�8O>���>�Ë>�q�>�M�>͉Y>�Y�>В�>�t>���>�x�>�H>���>��>�gx>�1@!_�?�R�        A����|\��H����ġl��0���0/��/�m�5��ƽUƽbƽU�2��        >�~A�Q    ?��>?��>{@��@zxFQ�M7\v�7�        ;@�?S�V?U�?X �?\.�?c+?n��?|�?�                                                  �<�6���A�  >�f>��@�b�&�      n�     ΰ9'��8���8 ��7e8U6�@i5�)�4�B4
�                                                                    8S�T7�F�7"S>6��C5��5��4ǈ3/a�                                                                    G��;F�c1F]��E��lE�$D6H,CG��BE�                                                4���4c�3��[2�+�22]�1g��0l�/�jL                                                                    3��83=�P2�2��1aM�0�R�/��~.���                                                                    79c�7.h�6�dv6=��5�M4�7�4�^3T�                                                                    7K��7	�6_�5�35X�h4�-E3�P3��                                                                    4�j�4�M�4J�3��J3.�L2��i1�#�0�Y                                                                    7x��7'v�6���6ԭ5�Z�4،�4
(�30�}                                                                    8���8��|7�Ew7B36�Fy5�*5 �4!�                                                                    5m]E5~��4��x4qG%3���3;S�2yx31�3                                                                    8���8��|7�Ew7B36�Fy5�*5 �4!�                                                                    7��7�~c74�6��G5�0t58
�4]�A3��`                                                                    4�\�4���4h�r3��%3Zp�2��y1��]1w                                                                    7��7�~c74�6��G5�0t58
�4]�A3��`                                                                    �RNp��+;5�/�5���5`4U��3��k2��{                                                                    0��1k�:���:i@|9�`�9NVd8�-�8ly                                                                    5���5:a'4���4N3�B2��!2 �g1$_                                                                    7�7�O6PK(5�
;4��s42�G3I��2vG�                                                                    ���KƵ��M�Z��঳�:������$#�                                                                    �ί���l�ː=���ʲ�;��N�{����mH                                                                    64T6A�5M�4�w�3��U30��2HȻ1t��                                                                    �I��䧶���OԵ^L1���ȳ�yT���                                                                    ���M��e\�䲅�б�h�/ǰc寠�E                                                                    �fy��]jX�Z������
 �8�;�I�q��O                                                                    )���)���/8�{.�.+�h-@�_,�t,!�                                                                    +E�+�HT                                                                                            8k�8\�<7��7O)6��-5�f5Ʃ4A��                                                                    4�x�4��4e3�>30�'2I1��0�5�                                                                    7g�7Xw6��I6M�)5�_@4�5y4�E3T1�                                                                    3��3ځG3`Wq2ҹ�27 1��C0��/�\�                                                                    8b@k8S��7�E�7L$(6�!S6pF5#+�4f�;                                                                    4���4���4[��3�5?3=A2�I�1���1��                                                                    7�߹7���7��47'�6���5���4�ʧ41�                                                                    6�j5��5�%�5>�c4��d3��C2�0�29                                                                    6�U6� :6{�C6"5^�K4�	3���2�C�                                                                    4х�4�Ig4��K4P&3~�2�
e1���0�9                                                                    8
�*81�7��7L9F6�|/5�d'4���4"�2                                                                    6s�695ټ�5if4���3�)R3
֮2:�                                                                    4���4���4�/�4��4w�4'��3���3�eE                                                                    3��R3�3�{�3rW3G��3�2���2���                                                                    4��04�84��4�E4���4M �3��3��                                                                                                                                                                        &�      n�     ΰ6u&V+�9�F�P�    =%�T            6���4:�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @۶�    @۾@    16:50:04        F�. @۾@    @���    &!      o     �P5��7;pB'�fB'�f7P��D"g�    &=��+���6�qC?�  ?�i4�Dv2�1)�0i�Ů��        ,{�(x.�=��=�d�=��=�d=�Hf=h8�;���3��:4i<!N�<�S�=G��=���=��c=��H=ؘU=��=�	�>Ǫ>1�?��l?��F?6��>�e>��>�0=��2=U�@T�@]�O@TR@T^�@Tq�@T�G@T��@T� @T��@T�&@T�@T�fA0x*<���{@��                    E��Q6��G�ĥF�!�F]Q_E��
E��D5�CG)UBE��                                                A~%C�B�)A���AH��@�@?��>М�=��x                                                                    EBm9@��F �RA�m�3p��F^N�B3�@Fu�-ɧ�                                                    {@���*���*��@�S7�˘�������~�/QU�%V�?��&=����0�-ʜ 4S�    >�D�>��2:@Q�����9= ]    B�OMB�OMC�>�C�>�ChU�?<�>��c3 u�C�5�6��6��A8�gA��>�D�A���@�S�B��@�L@�}nA��@�"�B����G@���#/b����G    ���G���@��>��%>���@���?�<@���?`1�?g�9@��?�g�F�xo3p��0��.��3v�6�#2�.�5.��F�x�F�<�G��G,w�?e��?b~�AB�@5��?B�U?+vp?w�>��9>�V�>��A>�d�>��K8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���Mұ
    ��o���o�{@��76�T8x��8�S�8$:]7�T:7CbV6���5��_5"�                                                A��I������{@��B~XjA%��B<�<@~f=�{@��{@��    C�iBQ�-6fI�    6d�06��@�D?35L6��    >���C���7�1B!��D��kC��B>�ZA�ΐA,�@�1?�#�>Ơ�                                                ?8A'A��5@��9?t�?�>|,�=��9=�V<�O                                                                    DgW�F�զF?YD���C�ӷC<$B��8A���@�n�                                                @͝�C##Bk�A �@��@	��?g0�>�1�=���                                                                    C�\E��EJ|RC�C<y�B���A��A�"@ް                                                @>*�BmSWA�3�@�ð@�~?��A>��^> �7=�                                                                    6���7,�hA#{D>Zp�B;n?6�?�+�)���(�dG,��+W��)���(�0_*�,�+���+��a) ��7P�815��/!=$�!�+ґ�+ґ�!��l��	�3+��1ޢ�/| /|$��w7�	ѳ����	�4/\L    'w(1-{�1-{�'��    =1�ڶ�ޔ4m�M/��N.>��(�MN                        *�
/\$�/[�x*!JB    /���2���            2�E�:���,��!? ��?ټ?v]?#�s?@%0?g?c�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  BB�SG��>+եA�ӵ                                                B�                                          B��                    A�                @�p�    ,!3�(�1I��5���3��3w�=�ʾ6���>[��5;Ċ;<P);A��7)�k            ��Ã    4���3�    7T{@��3"�25D
5d#��r�{@��{@��{@��6�7r3            6��N    7��{@��    7��    4�hH5��{@��6��({@��7�6?T�6�Z�    5�T{@�ζ'݋47=�7=�<�e    B�D�6�JMFZ�uD�#?W�?mA�W�            ?d��15��:�15��=j�=�A@���@LE�?ɸ?)zu>���=�{��U��^���7��ʲ��DɊ�%�}^��[!^�A�,��o�>?dZ�=j�=��@��@K�?ɥ�?)rH>���=��!��U��^���7��ʲ��DɊ�%�}^��[!^�A�,��o�>:��-���)�*)5��53��;�^�;�y�:�[9ɾ8;:���                                                ����q?���5��T��|��"���c�És̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 5h_o4���>���7��7!��@yT�@_<�@ev�@���>}6�\#6>�%E"��6�=1d�A�e>Yp�=��6�I9��6�ԝ7_
CFD�D���D���D?��C�lC
�BK��A�V"@���                                                A�D�CR��C@��B���B,lA��@���@V"?��                                                                    E��1G�gG�F�QRF=��E�j�D�dC���C �"                                                C��EJ��E>�DE�;D���C�ɒC�B-�A;T�                                                                    E��F��iF�1�Fz�xF\�<F.ܑE�3�E��EG�9                                                D)V�D�o$D�H/D�6(D��ZD~XD0��Cߋ�C�a�                                                                    /���6���GT��GK?KG�7F��F|:.F�?E�eEP>�                                                AH��A��>AG��A1L?��C                                                            ?�%�A jA���A�O�B?ԁB~��B��BǕ<�<�<�<�<�<�<�<�<�<�<�<�E��RE�u�Ee�/E�WD��DWrxC�7MC��                                                {@��{@��{@��{@��{@��C@��Bqǂ+���g��7�d�Ce�@j��    @[}�@[}�{@�ν������C���{@��C�k�C�WOC���{@��@��-{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G^�=?��D��YG��D��0D��dD��dA!ݷA!ݷF`�BB�<BV�C� �FK�6FK�6D�[D�[Fa�BB�                ?H��C�@C�?�C��t?��C�eC�;�C�;�C�|�C�ɃC��C�SwC���C���C�D�C���C��C�C�=VC�o+C���C��C���C�ЦC���C���C�|�C�E�C��C��C�KC��C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�[�Cί:>��j>��k>���>��>�Z`>�U�>��I>��>�	x>���>���>�OO>�M�>��M>��>�?�>�Ik>�;�>�[�>��@'|?���        @�H3�����J����H�n���1�����̾3�{�	ƽa�ƽa�ƽa��>Bn        =T@�A�a�    ���R���R{@��@V��FR�7(g�7�        :8^>���>�
�>��j? �?S0?]+?g�j?s��                                                �t|�6DGA�  ?.�:>=�?���&!      o     �P8>��7�<&7Up6��96U9*5�0�4Ȁf3��e                                                                    7pݹ7
6�̈́6�u5���4��3�C�3�h                                                                    G�ĥF�!�F]Q_E��
E��D5�CG)UBE��                                                3��3^i�2ߨz2n�01�d�16z�0R�/q��                                                                    2�j2�x�2A�1��G10f�/��.��                                                                    6SAa6�5O6!%�5�q5I�94���3��y30��                                                                    6q
86K�5�5�5xXP5�4��3�_3#1                                                                    4I24\a3�dq3=#2�ע2�B�1���0�(                                                                    6�M[6x*�5��55��157*�4��4��3!�X                                                                    7��$7ۼ7FD6�3�6D�i5�b�4��W4��                                                                    4���4���4M�i3���3�h13��2]�1}�q                                                                    7��$7ۼ7FD6�3�6D�i5�b�4��W4��                                                                    7�7#836���6�f5�K�5�*4EI�3a�                                                                    4  R4.Ʒ3��H3e,#3
!>2�:�1���0�[                                                                    7�7#836���6�f5�K�5�*4EI�3a�                                                                    �@�t�R�56�5��4�u�4$�^3xg�2� �                                                                    7�x7�k�;٩:�4r:eH�9��9+h8n�                                                                    4��4�݈3��3��>32Nj2��#1�l1�M                                                                    6!�-6D	o5��B5#"e4��\4A�33g72M�p                                                                    �Tr�G ��:��嚴Uӳ����ٍ���K                                                                    �넪��A�)U{�7�Q�?	�#��ؽձ���                                                                    5 ]5A�V4��4 �/3�Fg3��22zp1Lk�                                                                    �%���IXR����=�
�b�}5Q���0��х                                                                    ��a��˒R�i�2�Kx��#�	7��H'ܯ���                                                                    ��X3�������(�p���вƗ�3WI�a\6                                                                    (�*;)"tJ/,3�.�E�-�B�,�#r+�Br+C��                                                                    ,Rx�*�#�                                                                                            7�`�7�n�7;��6�)_6N�X5���4��4 5@                                                                    4�F4%=r3���3R@�2ܻY2H�1��0�y|                                                                    6���6�C67�Q5˱�5U�%4��n40@3/��                                                                    3�3!߁2��n2P�,1�h1Q��0��R/��:                                                                    7���7�ſ73�6�:x6\��5��=5�24?                                                                     4S�4��3��R3O?g2��2[߿1��0�~�                                                                    7J�7!{�7�	6�D�6,�l5�}4���3���                                                                    5�958��5#y4�~4E��3��j2��21���                                                                    5��6}�5�}'5�[�5��4g�|3��W2�C�                                                                    3���4"83�ju3��}3�12���1�	�0��w                                                                    77E^}7m6�7�6SC�5�b�4��4�y                                                                    54�F5a��55�4�?�4qr3�q2��2�e                                                                    3�3�*z4	i	4JQ4�4\�3��d3��/                                                                    2�p-2��]2��2�F�2���2��@2� '2V�.                                                                    3�w4}�4'�(46wF4=t54!ƶ3��3�R�                                                                                                                                                                        &!      o     �P6W='w(E��f>�9;��            6�}m3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @۾@    @���    16:50:07        F�l @���    @�̀    &0�      o6     � 2�
�7��B'�fB'�f7�D"g�        +�ڑ6���?�  ?��3nJ0���1�I/�-�ā                7��0���+=�.�IR5��a<�^�<�ل7;`�:#��<!�'<�dU=G�=���=��-=��z=ؙC=뽶=�
e>� >2G7��00��+<��+!��+S+��>2Ӥ�:��i@P��@]��@TQ�@T^�@Tq�@T�@T��@T��@T��@T��@T��@T�+{@��{@��{@��                    E�)[4ZG�r|F�.�F^@�E�X/E	-YD6.3CGn�BE��                                                A�C�_B��2A��AI�1@���?��<>��=���                                                                    EBo�@��F ��A�s14��F^�B2�=��-!��                                                    {@���T��T�@p��7Y�%����ƥ!��/E3~$(p?t�    ��#-��2?�    >�y�;�oP1�@a��ћB����;�y�    Bj72Bj72C��C��CRtm?�  >��3�.�C�{6���7$MA5�A�G>�y�A&��@�Ar��@>��?���Ac�@p�2A�&g���A�    ���    �����    ?ZE�?ZE�@�=>�]�@LXz>�{�>��,@*&(?A�F���14��.�i�,�1zz�5�F�    3�.F"��F"��FrQ#Fn�L?fA�=G�9A�L@y�f?v��?M�?*5?>���>�/,>��'>�2�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�\    �''��''�{@��5�$1��Z    31D�3�+�5ƞ5���5!�4ϕ�                                                A�aL�hٗ�hٗ{@��CH�j@�B3��@o<��){@��{@��    C��=BQ `6,�    6ڥ�7��@�`3���6��s    >�9�C^�J7�+}BE�=D��ZC��B��B#2A�)�@���@�g?>�                                                ?h+A�;v@ӎ�?�C�?ar>ӆ%>#�=_�<^��                                                                    Dk�}F�'�F
�@D�gAD��CJ��B�e�A�+�@��:                                                @�"5C6�Bu
5A6U@�x2@I?���>�.�=��                                                                    C���E�h�EN�C�fCL��B�@A��Ad�@e                                                @I��Bs�A�@�g�@'W�?�?O>�Yi>2��=/jh                                                                    4i�4���A"!8>X��B?4y^?���)Ē'��R+|R�*ƽ)Ē'���(�\+D�+n(lV*7�                        ���3�.�1��4            6�ኳJ�X���2u�                        >J�����E                                                                3O]�            3O]�{@��,�*v>L��>L��>L��>L��>L��>�!m>ٳ�?Y��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A���G�kj>-H�A���                                                B�                                          B��                    A�                @�p�                2�
�1��i1zz�=	q�5�ը=��4�=:��*:�i�                ��!                    {@��2��5_��    ��{@��{@��{@�εpk6
��            08    08{@��    08    5[��2!��{@��    {@��6��63y �6���    5��{@�ζ7��1+l}        <��~    B��X6zͯF�n�D�):>�TS>1:A 8�            ?� }    :t�    >3�y>L�\@���@a-�?���?c'�>�<�>�%��U��^���7��ʲ��DɊ�%�}^��[!^�A�,��o�>?�b�>3�y>L�\@���@^W�?�+�?bʊ>���>=ˋ�U��^���7��ʲ��DɊ�%�}^��[!^�A�,��o�>;���    ��0�;.���<��.=5�/;��A:�n:�: Yd                                                �!���O51��E���)��AhA�:y�d�ģ~�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6^P�5a��>�Ih8s7D��?�1?��?�`�@L7�=�d�7�N6��EE�I�6%/�=��bA�C7�'�=�(6%/�JJ��5;��6��CF78D��1D���D?�NC�x�C
 �BK�oA��@�O                                                A�78CR�1C@��B��NB,x�A� �@ˈo@�?O                                                                    E�~�Gf<GF�Z�F=�2E�wD�]�C���C �                                                C�EJ�@E>��E��D�$C���C�/B-�A;F�                                                                    E���F���F�-rFz�F\ȖF.݃E�4�E���EG�$                                                D)U�D�h
D�BD�2�D���D~YyD0��Cߌ�C�a�                                                                        4�I�GT�GK"�G�zF���F|AFѧE�nEP=�                                                A��A��fA�yPA�Y	Bg'A��&A6Z�?|�                                                =؈@JMA>�AHjA���B{�B_��B�d�<�<�<�<�<�<�<�<�<�<�<�<�E���E�`�Ee��E��D��DWqlC�79C���                                                {@��{@��{@��{@��{@��A���B�W-    �fn%6��C:�m@Q��    �uEB�uEB{@�ν�(��(C���{@��C�U�C�t�C���{@��@q��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��GB�>8k�D���G�xD���D�!�D�!�A-~	A-~	F_��B?+�C+qC]�FK�FK�D�	D�	F`��B?&�                >�мC���C��C�Y�?6GC�/MC�Y�C�Y�C��YC���C�4�C���C�*C���C�LC�[/C��fC��cC�-�C�sJC���C��tC�5C�_�C�w�C�}C�r�C�T�C��C��C�JC��C��-{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�"nCӔ@>�>��?>��>��0>��k>���>�|�>�k�>�y�>��>�˫>��>���>�y�>�[�>�à>�PY>��>��>�{�?���?�@        @0�.��W���Wz��WkĤ�J�����W���īN���>b��>b��>bĤ�        =L��A���    �,W&�,W&{@��@B�rFRy6�p7�        9,97��    ;PM</Y">,,�>���?��?>
(                                                ¡��4��_A�  ?k�    <#�
&0�      o6     � 0NN    2a��2��4%�4�V�3��j3�̢                                                                    /�LE    1���2H3Q2�3�7{3 V�2��                                                                    G�r|F�.�F^@�E�X/E	-YD6.3CGn�BE��                                                +�3�    -�w�.qT/���06�z/��)/y�                                                                    +�|    -Y	-�;(.�.�/f��.��.:H�                                                                    .c��    1)�)1�]�3+�3�]3�:2�;�                                                                    /�"r    1wh1�\3�c3�ܠ3D!2��J                                                                    -t�g    .�t/�^�0�AY1�e[1%�'0��                                                                    0 jb    1 �c1���37M@3�n3o�(3��                                                                    1�t    2V_�2ȇH4
�4��#4u�3�s                                                                    -��    /c��0B�1k��2$�1��(1)��                                                                    1�t    2V_�2ȇH4
�4��#4u�3�s                                                                    0^T|    1���2��3{��4��3��'3t                                                                    -Y4�    .��/w�0�}�1�ڝ1��0��                                                                    0^T|    1���2��3{��4��3��'3t                                                                    ����    0-�Y1�	2��3!�2�vl2g��                                                                    8��8��n;,k; �:���:3�9C��8�d                                                                    .�    /)k�/�YB12�1�S�1b�1ZV                                                                    /�
�    0���1(6�2*`3��2m.2O                                                                    �+��    �(A��s��%9I���ڲ	����                                                                    ��a    �2��9Rq�*�#:]�	!��G��                                                                    .��6    /��&0&@_1|2\2�T1mfw1Yq                                                                    �2�f    ��Z����ˎ�|�q���r��M                                                                    ���@    �vZ#�,9�eO�	 Y�|�:�"3                                                                    ����    ���ί*0���4��ɰb�а	\�                                                                                                                                                                         �]��ɤ                                                                                            /�
#    2E�02΅�4 H4�J[4�3�9�                                                                    ,��    .��/S�%0��Q1G�00�O�0]#�                                                                    .��    1A��1��3%��3�k|3&�$2��                                                                    +��    -��.R"�/��"0Q�</��M/rzJ                                                                    /� �    2=��2˒	4+&�4ʌ�45y�3��[                                                                    ,��    .���/P��0���1[�G0͐�0��j                                                                    /f8    2�22���454��j3�Z93�70                                                                    -P@    0C�0�v�2:�2��1�B�1��7                                                                    -�J�    0���1�`�2د�3h�2�|�2Z��                                                                    , ��    .��_/���0��g1���0�EI0y��                                                                    /*`}    2'�2��[4#ޖ4�w�4~?3�Ct                                                                    -B�j    0>�[0��2;G�2Ȉ�2��1��`                                                                    +�#]    /�)0b`1�v�3t�2ގ�3!��                                                                    *���    -�<�.��0�Pf1��1��!2�
                                                                    +��V    /1$O07͑2�43!��3�3E��                                                                                                                                                                        &0�      o6     � 6��    F��    =G�    >�P�    7�4jā@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @���    @�̀    16:50:09        