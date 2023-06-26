CDF      
      time       bnds      lndgrid       levsoi        levdcmp       cft       glc_nec    
   ltype      	   natpft        levlak     
   nvegwcs       string_length         levgrnd       hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      Sun Jan  9 16:23:28 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1919.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1919.nc
created on 12/10/21 16:45:43    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:28 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1919.nc had following "history" attribute:
created on 12/10/21 16:45:43
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�fF�  @ؘ@    @ؠ     $ѹ      b�     x 3�26���B'�fB'�f7<��D"g�    >��+�q6}c�?�  ?��4�-�1@��2�J\0yf�/Q}        }=1o .8��81�U)0���7w=���=?�l9�ա2�R�:*�,<"HK<�w=F�&=��=��V=���=�y�=듀=��X>��>�8��1�w�+�4?0��
7��0A:87�>���@cPI@iK�@_�^@`��@a6�@a�:@bY�@b��@b�^@c�@c*�@c+�{@��{@��{@��                    E�}*5rA�G�YG}�Fgo�Eǒ�Ep�D?F�CNBMc�                                                A��C$��B��A��AQ+�@�Q�?�l�>��b=�&�                                                                    E@��@ěF%�A�p�2�\4F[K�B)�%=��-�C                                                    {@����1J��1J@^@H7k�)�����!�G��/L���>�a~>�����(/tU�2���    >�y�=L��1�ξ�~���
����A<Y�=    B;B;C�	aC�	aC�� ?�  >���3zC�6-��6��A$�A*��>�y�Ab�a@V�A�d@oA�@E�AA��@��oA����x>�=!/4'�x    �x����@m	�?|��?|��@_޹?D��@�m?;�+?5u@s,U?��	F��2�\40[��-ck2�Ɩ6b6u    44E&Fm_;F�>�F���G�Q?�ehA#=B�~@�S�?�1 ?O�#?)Pj?�{>���>��U>��h>��]8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�+�    ���F���F{@��69'4��    6l�6�h�7C�16f�O5b�~5-�                                                A��+���^���^{@��C'�
@��*B.�6@[x�=2({@��{@��    C�O�BY#b5�g�    6GB�7w��@��x3	{�6�0    >��QC�>X7���BV\D�Y�C�ۮB�i�B!�Ay�@�@d�?&}                                                ?o�A���@�&?�R;?R��>���>i�=H��<I��                                                                    D|�F���F��D�ˋDD�C]�MB��A�<�@Ì�                                                @�xC	)�B�9�AF)@��!@0�?���>���=�                                                                    C��E��_EhyD
�@C^]B���B	Q.A%��@'�                                                @ZZ�B{[,B  @�T�@5k�?�Y:?&�>7�%=:X�                                                                    5�ܿ6%+�As�>9��A�i�?�,?���)��L(��+�7+G�|)��L(�Y�)���+�(d+�u(�z7<��                        �2�F3z1��{            72�F��a�2�F3�X9                        >g���G<                                                                2���            2���{@��,���>L��>L��>L��>L��?_F�?-�(?��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A�CfG��>,z�A��                                                B�                                          B��                    A�                @�p�    ���~    3�22�@2�Ɩ=�}�6M�>G�4���;>-D;C�x                ��
    5$�0،2    7��S{@��3T��4�SX7�'����{@��{@��{@��3|��8Vd            7ҩ�    7ҩ�{@��    7ҩ�    �uo46�w�{@��7:5�{@��7�de6��7Af    �:ʥ{@�ζr�2��7���7���<��K    BÞ�6~]�FtS�D�e�>�+>&�vA#�Z            ?�kT    :7�W    >�6t>C�@�ΐ@��)@)�??f&<>��>Gj.�f������~��)>��kP���T�s�R�(�9�#�&5����?�S�>�6l>C�@��/@���@)�*?f&<>�02=�ċf������~��)>��kP���T�s�R�(�9�#�&5����=B��    '���4yÒ,��=�XK<��6�|1�^(=C=�                                                �!���D�	���g���a����ķ�ğ�Ă_�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 61mx4��4?Jק7�Ǖ6��@
�V@@!v@���=�A�7~��5��+F��27���>gݝB��?��>j�7���J@�07���85+`CN++Dє�D�l>DN��C���C��B_cA�ِ@�zh                                                A�++CQ��C@l>BΙ�B8��A���@�c@
ِ? zh                                                                    E��hG	��GofF�TlF?m�E�P�D�[C��$C�_                                                C�,:EH&�E<D�ET�D�8�C݌bC�BB.p�A<�E                                                                    E�rF�]KF��<Fz |F\��F.��E�2�E��;EH/                                                D)9�D��D�ԴD��D�lD~;TD0�4CߙC�wh                                                                        5�;�GR��GIc�Gl�F��F}SoF;bE�A�EPq                                                A���A�@A˼jA���<��\                                                            >�L4@|Q�A>�A��|B'bBVd*B�ѩB�J�<�<�<�<�<�<�<�<�<�<�<�<�E�E�Ef�E�1D��&DXC�|�C���                                                {@��{@��{@��{@��{@��A�+uB��z�0�e�p6۴�Cȅ@X�    @A�~@A�~{@�ν�v��v�C�:�{@��C�.PC��]C�:�{@��@^@}{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��GVn>��D���G�D���D�:�D�:�A:�A:�F[�,B2��B��"C��FKךFKךD� �D� �F]R�B2ծ                ?�>C��C�}�C�,�?�*C��C�Z�C�Z�C�elC�wC��AC��C��8C��C�A�C�vC���C�בC�C�C�C�}�C���C��C�$C�HEC�]=C�b�C�XYC�*sC�C��C��C�jO{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�zC�>�͇>�L�>�>�.>���>��Q>�B�>�!_>��@>�F>>��J>���>���>��>���>��'>��2>���>�*E>���@��?�r�        @f���k���j���jrĜi���Ƌ���j�¿�Ĥ�qƿ��ƿ��ƿ��ę��        =M��A�X    ���W���W{@��@U�-FO��6���7�c        0�"l:Bm    >+�>�V�?2�=?6��?<J?C�                                                ��\5��A�  ?k�    <#�
$ѹ      b�     x 2�1�    5���5�Ȳ6FW�5R%�41��3�S                                                                    2	,�    4���4�3h5z��4��h3`L�2߬�                                                                    G�YG}�Fgo�Eǒ�Ep�D?F�CNBMc�                                                .c�=    1r�1�w1��0�2�/�	7/9}�                                                                    -��    0;�*0FLa1F�0�.��+.jM�                                                                    0�[�    4X�4y��5@O^4_�3_p�3��                                                                    2w�3    4\|�4�u5E=�4���3�3�3�                                                                    0�Y    2$j2I+�3*�2X��1[da0ܙ�                                                                    2��    4��4��5qq4�93���3"}�                                                                    3��w    5��5�u�6M��5qǬ4`�g3�s                                                                    0~۵    2���2��.3���2��1�2�1\�                                                                    3��w    5��5�u�6M��5qǬ4`�g3�s                                                                    2�Ñ    4��l4��a5���4ʦ�3��}3A��                                                                    /��    2՜27L3"�2Qʚ1T��0��                                                                    2�Ñ    4��l4��a5���4ʦ�3��}3A��                                                                    �*=�    3��3��4�K3��@2�xs2�k�                                                                    9T��8���;A-�; 8�:��@:ߺ9C��8��3                                                                    0��    2�j�2���3oPy2��71��@1��                                                                    2��    3���3�I�4���3�2���2+l�                                                                    ����    �K�C��s�?=��M��;���6�                                                                    �� !    �Vs���
�)K��)A�9�o�q��                                                                    1й    2��2�?-3���2�޾1���1-�q                                                                    ���/    �׆�.ga��P�!������k�                                                                    �6�
    ���^������ư�x믳�w�OĿ                                                                    �[�    ����������˱�yް����'�7                                                                                                                                                                        *ķ�#��                                                                                            2¸    5��5���6JV
5hI4X�.3�HK                                                                    .�b2    2��242�/{1�t0�>�0���                                                                    1�,    4{T4��5Q^�4s�3me:3	7[                                                                    -�\.    0���1J1߳�1�x0t</�nG                                                                    2Ǡ    5u̲5��d6Xg�5~��4�$�4J�                                                                    .�V+    1���2�`2�7�2
.�1I70��                                                                    1�"!    5%V�5S�v6005&��4V�3��"                                                                    /�'    3<�t3q�41[�3>�2.�1�u�                                                                    0bzd    4�s4+�4��,4�|2�42�A=                                                                    .�j�    2��2C�i3Q�2�1�0�J�                                                                    1�F)    5J�5�aH6=��5K�4:13ɌF                                                                    /ý�    3f�3���4X�W3i�2T�L1�W+                                                                    .Q�y    2.92�0�4	`�3��m3�z3D
o                                                                    -)XM    1�_1��2�*2uu1�s22j�                                                                    .�.    2T�2�h4'��3���38�3o��                                                                                                                                                                        $ѹ      b�     x 6:�I                    =�    6�?�3ӮA@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ؘ@    @ؠ     16:45:43        F�8 @ؠ     @ا     $�      b�     }@3>��6�}�B'�fB'�f7��D"g�        +���6>s�?�  ?�a4��0�΀1զg0�7�/sJX                9�M22ԙ,��137s�9�sw=zJw<3�]4d#Z:2�W<"�H<�w=F�=�׭=���=���=�{E=��=�׉>��>93G26uO+yL�+!��+S/M��6�e�>w��@_w�@i1b@_�?@`��@a6�@a��@bY�@b��@b�F@c�@c*�@c+�{@��{@��{@��                    E��4�}G��HGmFhpsE�mzE��D?��CN��BM��                                                Ac�C%UAB�zYA��ARy@��6?ȷ�>�X�=�Qr                                                                    E@�F@�x�F
�A�T�26[FZ��B(�b<�D/-��                                                    {@�ο�ٔ��ٔ@Z�E7,�'%զ$����/B3�%��)>츞    �$.e.��2]��    >�y�>ey�1�X�?����F���n;���    B@�B@�C���C���Ccȱ?�  >��39�C�R�6r66���A"y�B	�b>�y�B5GhA]�Bz��A��A,7B3]ANb�BJ+_�z��?�Y    �z��    �z�����>i_�?}��?}��A-�@?�cAP��?�}�@3�AIJ�@Q�zF��X26[/�ȳ-p@�2Ä�6�(3    3��F���Fo�RF��IF�p�??�C?��VB�+�@�3S?�"�?O��?+	?>��>��>��9>��8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�(�    ��+��+{@��5���        5�6�5���5��_6d>m5^!�4��                                                A���Nl1�Nl1{@��Ci��Aq�+BA#�@T)U=�S�{@��{@��    C�w|BY)5��    6D��7w]�@��3��6��    >��Cx�7�ҨBfn�D��C��B���B<�0A��]@�a�@�7?	zn                                                ?�	�A�U�A b�?��?|��>�=>#��=`?{<M�<                                                                    D~�F�U_F��D��Dd#Ca5"B��BAƍ�@Ĝ&                                                @�>�C	��B��qAQo�@�]@5b�?�67>�
=��                                                                    C�	 E��nEk-DRQCf.�B�6�B �A)�c@(/                                                @^�B|�B^%@��@?�?�%k?�3>>L,=;�H                                                                    5��5�dA)d>4Z�A�C?Jx?���)uIi(,+���*��|)uIi(D:�*R6S+���+L&(��7��                        �� 39�1�g�            6� ��J��� 2�c�                        >�P���G%                                                                2���            2���{@��,�@�>L��>L��>L��>L��>L��?>(�?�W?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B�_8G��>>+�rA�J                                                B�                                          B��                    A�                @�p�                3>��2[�2Ä�=�qm6> >+�+5<I�:�n&;�                ��d    2���/=mM        {@��3��4�os0w%���j{@��{@��{@�β�	 6��            0w'�    0w'�{@��    0w'�    4��E5>�Y{@��57T�{@��6��4gj6щ�    4ʖr{@�ε���3�k�.��e.��e<���    B�C26v�EF`,D��K?Ɇ-?�YB�1            ?���    :F�f    ?(
0>��s@�g@�cI@B�u?�s>�B�>q��f������~��)>��kP���T�s�R�(�9�#�&5����?��G?(
+>��s@�SK@�0�@Bkf?�s>��d=�� �f������~��)>��kP���T�s�R�(�9�#�&5����=s��    !�6~4�X'�*�>9C�=̖G;�01Q�};'��>'��                                                �!���C�������)���!6��+�����̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6�J5k��>��m8e�7e�xA*~@��`A��A7o�>�A8�b6�ƸF�1�6���>�]QB�"�<�yV>���6���J�4�Ƿ71�CN1�D�]�D�-�DN�C�� C�B_��A�@��/                                                A�1�CQ]�C@-�B��B8� A��@߶�@? �/                                                                    E���G	~�GVuF�]F?��E�\!D��C��C�                                                C�',EG��E< �EZ�D�O�CݝHC�B.x�A<�C                                                                    E�JF�X1F�ݻFz�F\��F.ɉE�4&E���EH�                                                D)8�D���D��(D��zD�kkD~<jD0�Cߚ�C�w�                                                                        5,�GR��GIF�Gp�F���F}Z&F=FE�C�EPq�                                                A��$AСpA��A��>A���                                                            ?��@�)A#��Az�PB֚BHNBo��B�2�<�<�<�<�<�<�<�<�<�<�<�<�E�E�NEf��E��D��DX�C��C��p                                                {@��{@��{@��{@��{@��A�;�B�wt    �f�x6�7B��@��    <L�h<L�h{@�ξ4��4�C�R{@��C�0�C�p�C�R{@��@Z�~{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��GAy>�	@D��lG��D��DD�: D�: A?�XA?�XFZ�dB1>�B�ڢCv�FK�bFK�bD�	D�	F\��B19�                ?�|�C���C��C�R?-�C��8C�WC�WC�&C�>�C�\�C���C���C�؝C��C�9�C�b|C���C���C��C�#�C�]C���C�΋C���C�"&C�8�C�F�C�/�C��C��C��C�,�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C�]L>��L>�K�>�p�>�ܕ>�>��u>��>�eE>��2>��~>�I>��%>�3">�>�D�>��_>���>�{;>���>���@�m?��        @�J���ix��f��c���ʶ���������ڢ�����z���z���z��ħ�        =^��A��    ��7��7{@��@N��FO�P6��7~        7�O<        =[�> J>�{$?%�y?*��?0�                                                ���5�sA�  ?k�    <#�
$�      b�     }@        4���4�]#4髳5N|4+�3���                                                                            3���3␘4��4�,�3X�2���                                                                    G��HGmFhpsE�mzE��D?��CN��BM��                                                        0:��0;��0t�}0��/�*�/$2�                                                                            /k��/muP/���0f�.�Q4.OhD                                                                            3�^03��?3�5�4Z�3V�-2�                                                                            3��,3���3���4���3���2�(}                                                                            1f�`1��1�1+2_=:1k;�0ƪ�                                                                            3ʱ�3�qu4�A4�I�3�_3&�                                                                            4�8�4�Q,4�^�5o`�4\�3�^�                                                                            1�!-1��2EG�2�-�1���1D9!                                                                            4�8�4�Q,4�^�5o`�4\�3�^�                                                                            4NG4
LX4I��4��3�
�3,@�                                                                            1K�01f�1���2Q(�1S�0�i"                                                                            4NG4
LX4I��4��3�
�3,@�                                                                            2�_2�ڶ3(|D3���2�A2~`�                                                                    9� �9H�;@g�; �2:ݭ:!�`9B�>9	o~                                                                            1�[(1�?�2��2���1���1�                                                                            3$T\3`�3D3��=2��2P{                                                                            ����������I_�4f����[                                                                            ��,�
G��m��N�2\�U��                                                                            2$�!2�w2E�P2��1�)?1X7                                                                            �E̺�PH������S����fv                                                                            ������s�*����(������7�2                                                                            ��:߱ RG�9r±����1�-�                                                                                                                                                                        $�b����                                                                                                    4��4���4� 5c�?4PP�3ݙv                                                                            1"T�1)�L1~J`1��`0���0{�                                                                            3�H�3�S�3�F^4n��3djN2��[                                                                            0�0(g�0���1k~0^1/��8                                                                            4�4�$�4���5z�4x�4.�                                                                            1�1'1|1���2�M1��0��{                                                                            4N��4|��4���5#i4dT3��$                                                                            2l�;2�\�2С3:�=2'N1���                                                                            3':x3L%�3��4�2��2k�R                                                                            1??1iO�1���2�12-0��S                                                                            4|��4�b�4�5G�S42�g3�G,                                                                            2��z2�q!2���3dA�2L{�1˿                                                                            1Y�L1�j�2���3��U3�3-g�                                                                            0/��0�X�1���2pf[1��2                                                                             1��2	$�2�v�3��h30��3S��                                                                                                                                                                        $�      b�     }@68��                            6�w�3���@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ؠ     @ا     16:45:45        F�v @ا     @خ�    $ҁ      b�     �3~6XZ�B'�fB'�f7:�D"g�    
[�|+��U6��?�  ?�3��0�I�1�x0a��/e�        ��n�� 8ԣ�2�}+R�B.�#*5�D�=�<�5�4�#�:9�8<#+h<�=F�==���=��C=���=�|�=�=�ؔ>�>�8��2%�+d�#+!��+S.��N6a��>�@[�@i	�@_��@`�C@a6f@aӾ@bYg@b�w@b�1@c�@c*�@c+�Ap)�;�^){@��                    E��64��G�k G�&Fiv�E�TE]mD?�CO{BM�8                                                A�C&CB��A��7AS�@�I�?��>���=ך�                                                                    E@['@�L}F�7A�0f1Ӏ2FZ��B'��<�D/.a^}                                                    {@���,u@�,u@@V�.7G�%{�U�������/.���B>�n
[�|�#�.cѴ2E�y    >�y�?�O�1��V���@�>���
o<-��    B�^�B�^�C���C���Ca�8?�  ?���2���C��5į�6�J�A�YB�u6>�y�B�[�B(qoC�A[�*B �B��5A��TB���A$3���w��i3A$3    A$3Ax��AU� ?}A�?}A�A��@�/�A���@,��@�JaA�@��WF���1Ӏ2/��-c'p2���7��#3\Q3R^�F�~qG�GW�AG�?�A@�KB�>�@���?�/y?O��?+�?m>��	>�t�>���>�7�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�7�    �P���P��{@��5_4�5 .    5FA5xx�5b�6��5S��4�ǃ                                                A������y���y{@��Cf?�A�L3Bt(6@I�>[�{@��{@��    C�1�BY'35��    6A�Q7v��@�N<3�36�Q
    >5��CtT7��Bw�D�G�C��B��B^89A��}@��@$&?��                                                ?�GXA��}A	�@?�2�>��>.I=u�d<Ux�                                                                    D��F��dFQD�h(D��Ch��B���A���@Ɓ�                                                @�I!C
��B��A\�P@�ct@?�8?���>ȌG=�a�                                                                    C�(!E��BEm�*D�Cn��B���Bq'A.Vj@*�                                                @d	�B~�B�@���@J�&?���?�>D�G=>��                                                                    4�{�5M��A>.[A�͙??�r *T��(�-,��+z��*T��)*s+��y,$�%+��)�|79)01��E/��e,a��            6��3Z�1��I/��%/��e,a�ڶ����^�6��2���    %_;g1�W�1�tK.q�    >�e6���4��/�6.�R�,�V7                        )>$/���/�Y�-� �    /�:�2��M            2��^:885-�*>L��>L��>L��>L��>L��>���?Aa?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CEfG�Y�>+p|A�G~                                                B�                                          B��                    A�                @�p�    3�.�1�E�3~2�2���>�X�7E�<?$�6�;I�-;Pe@                ��T�    5��1��I    8O��{@��4�:58��7ެȴw�n{@��{@��{@��3Qb7��l            8( _    8( _{@��    8( _    �o�7L��{@��8'��{@��6�~�6O��6y\    �c�{@�ε��l5A��7�q�7�q�<�B    B�lH6u�F0��D@�@�,?鋈B�D3            ?��_1��E:O�G1��E?nc6? ��@�Gq@���@BAo?�$�>�r�>��Ƌf������~��)>��kP���T�s�R�(�9�#�&5����?��?nc? ��@�1b@��@;��?�#�>�r�=�#ڋf������~��)>��kP���T�s�R�(�9�#�&5����=�m�    *d�w6)���N/>�`�>���=�/B8(�5��>O�                                                �!���C�������:������/o���̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7�5QC1>�i�9�7�aLA��%A�2A�#2A�w�>@��8E|6�sF�Մ6��>y UB��w?ܠ�>{��6��J��8Ro8|CN?~D�//D���DO,�C�C�C�*B`.+A�q�@��                                                A�?~CQ//C?��B�,�B9C�A��*@�.+@q�? �                                                                    E���G	b�G=�F�e�F?��E�uD�&lC��-C�}                                                C�%+EG�&E;��Ea3D�hfC��C�fB.��A<��                                                                    E�XF�SF��<FzvF\��F.��E�5�E���EH�                                                D)8@D��<D�ǞD���D�j�D~>:D0�8Cߜ�C�x|                                                                    /�:�4���GR��GI)�Gu#F�>F}hpF?�E�FvEPr�                                                A��A�δA�FA��HB  u@��                                                        ?ͫ@�@JAW�Ar�*Aȫ�B>V�Be�B��f<�<�<�<�<�<�<�<�<�<�<�<�E��\E��Ef�}EΦD��DX�C��C��\                                                {@��{@��{@��{@��{@��A�7wB�r&�@�gw6�t�B�vs?��+    @�>@�>{@�ξ�,��,C�j{@��C���C�غC�j{@��@V�f{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G<l?FlD���G�?D��~D�[�D�[�AD��AD��FZ1%B/��C��Ce�[FK�qFK�qD�iD�iF\pB/�l                ?�i�C���C�o�C���?�C��7C�9�C�9�C�CC�QRC�b�C�x�C��!C���C��4C�0C�6 C�\�C���C���C���C�UC�S�C���C���C���C�LC�,�C�.�C��C��C�	C�F�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C�!�>�CU>���>�W�>�(�>�4V>��>���>�ƞ>��4>���>���>���>��>���>���>��c>��>�޴>�l�>�b=@]?�A�        A<����n���(���=�f$�����3��=���ź��ź��ź����I        =��A��    �m���m��{@��@V��FO��6D�7Z�        6*��;?.    =
A�=·�>[��?�?/[?#��                                                ���7(9�A�  ?k�    <#�
$ҁ      b�     �3�    4b߅4�f�4]�5��4 s=3�D�                                                                    3�    3�I�3��h3� I4%:3J��2�;�                                                                    G�k G�&Fiv�E�TE]mD?�CO{BM�8                                                /`8    /��0�/�L,0�L/�/~                                                                    .���    /0�/5�l/��/� �.�U�.>�i                                                                    1�~p    3,c3c��3Uݒ4
��3I�j2��                                                                    3��i    3qC3��G3�u%47�3�O2�j                                                                    1�    1'A�1|�i1dj�2m�1qA�0��J                                                                    3�f�    3�p3Б3��4_��3�V�3Y�                                                                    4��    4y�4��4oeH5]4S4�3���                                                                    1{$6    1���1��"1� b2��1лc16�S                                                                    4��    4y�4��4oeH5]4S4�3���                                                                    3��    3¦�3�(3��4�t3�u�3 &                                                                    0�P    1��19)�1?��2�1L�0�{S                                                                    3��    3¦�3�(3��4�t3�u�3 &                                                                    �&��    1���2v4�2�Q�3J�.2�Ek2g�,                                                                    :�(9��;B�,;"'�:�2::��9_v�9�u                                                                    1��    1��)1��[1���2Y��1���15E                                                                    3	c}    2�$�2�+d2�%/3m�=2�2!�                                                                    ��TN    �!��m�7�T���?�(�Y��;                                                                    �|Y    �*^ί��n�<U�hOg�&�ʱDk                                                                    2�    1�C�1༇1�[2p��1���1(_                                                                    ��0�    ��y��߳7��(γD����                                                                    �2#�    �~A鯣诡]��V�C��IW�)X}                                                                    �
�f    �����뤰�J��Ri�2��5�                                                                                                                                                                        -2����                                                                                            3�    4L�4|�/4a�5s�4C~f3��T                                                                    /���    0�a�1� 0�p'1���0�q�0g�                                                                    2�    3G�&3{e3h��4��3V[p2߮k                                                                    .���    /�->0 ��/��60�Q/���/}_�                                                                    3	/�    4C�`4yE�4p��5�s4i8y3�]�                                                                    /��C    0���0�us1 �#1�z1�0���                                                                    2��X    4j�4@��4,�h4� I4	�3��                                                                    0��    21
2\i2E3�2�/2�1�<4                                                                    1\�    2�da3؉3o�3�_�2�s�2X�:                                                                    /|K�    0��12
1[1�I	0�h0w�                                                                    2��    4 �4k��4R�4�'v4'y3��b                                                                    0�̬    27�E2��2q�3��2?e�1�I�                                                                    /L��    1
@X1�5�2�03<��3r�3cq                                                                    .%Dm    /�o�0�Y�0��[2U�1��a2 �z                                                                    /y�~    1(�O1�A�2:��3fg�3%��3B��                                                                                                                                                                        $ҁ      b�     �65��%_;g                        6�F3ͩ�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ا     @خ�    16:45:48        FŲ @خ�    @ض@    $��      b�     ��4���6Y yB'�fB'�f7w�D"g�    #�`<+q3#6�?�  ?B5mH3.Wz37��1Y1��1���1\��*�&��=�?=�@<ԫ�;'9I5D��<T�N<L*�4�$�:>�U<#��<��Q=F��=��
=��t=��C=�}�=�=�ى>�x>�>�)<�L�7<i�0�8�+�*.���5��<=���@Wӷ@h�t@_��@`��@a6@aӀ@bY8@b�V@b�@c�@c*�@c+�AS6N<��{@��                    E�Yx6��G���G�Fjt�E�C�E*D@��CO��BN-+                                                A�C&��B��
A��IAS�2@��,?��K>ق�=��I                                                                    E@(�@�#VF��A��2��xFZ.pB&�Z=A��/��                                                    {@��A~i�A~i�@RD7�4��1Q�$��$��j�/W�x&�8�>�^�#�`<@�R(0%�B4ܻH    >�LA��1���?٨{A�
NACM>��C    B��B��C�C�C�C�C��W?&�@`ɳ2�
|C��5�F@6���A��C�>�LC��Bi}�C0�0A���B7;�C�A�p�B���A��K�;s���UA��K�;>A��KBW�Ae�>�;�>ن�A��-A	��A�^@9u~@�#�A�}@`�IF��2��x0��*-1��2��T80�C5�T4�9G$Gim�G�^�G���>|��A\�%A���@d"�?p��?N�	?*��?�?0�>� T>�?7>�|d8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��b    @$xV@$xV{@��7LO8�B8���8G�6���5T)Z5�@�5A�D4�f                                                A�pv?E�?E�{@��B�LHB!�rB�ڹ@< �>�#{@��{@��    C���BXp]5��*    6>�7u�@��3*\7�    >Q|pC���7��B���D��ZDC�B��B��A�A	@1Ջ?�                                                ?���A�;�A��@�h?���?,�>HA�=�L+<^�(                                                                    D�3F��FC�D�h�D(Cp�}B��-A�w�@�Ћ                                                @��C
��B�ɣAg��@ݨ�@K_�?��>�\i=��                                                                    C��E�5�Eo��D_#Cwd�B��-B7A2�.@,c�                                                @h��B�CB~�@�_�@V66?� F?	�>K�$=B=                                                                    6�FT7޵A�P>(�CAޔt?{D?���+�.)��,�D�+��+�.)�+�,�f,�O-,sO�)�VI7[v�4:��2/;8�            7���4ZA91�802�2/;8������SH7���1�)    )@@�41?'4.B�1?�    =���7�T7y<�2o�1#�3/��                        ,2r�2A�"2/��0���    2o�2��v            4J��:S�v-�?	�*>���>�o>p��>L��>�8�?ч?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C(p�G�8>>+R�A�$5                                                B�                                          B��                    A�                @�p�    )�l%���4O�5�˧3\2��?A<s7�"W?�6:7m�^;��;� ;    /���        6�G    6|B5Ю�    8=F{@��2j*58ͽ���{@��{@��{@��9��7��            8Z�uJ�8Z��{@��    8Z��    �[Kb6�J�{@��84�j{@��5���5`KO4���    5.<{@��5��6�ȅ7���7���<��t    B�n6���F��C�	�Ap46?�twB�K�            ?��[4:��:9[4:��>���>�o^@���@x�U@*��?��->�">�l܋f������~��)>��kP���T�s�R�'�9�#�&5����?��$>��>�n|@��@c��@�`?���>�
�=��f������~��)>��kP���T�s�R�'�9�#�&5����=�w    -1��7R�6��>E[@>���>[��8��?6%|t>oV�                                                ��®�A�� ����j��G����K��@!�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6`�q5��=,�R8�f8=R�A�e}A��A��A� =��H7��e6�!�E�g�5(��=k�A�:�>��!='o�5(��Jo�8��6��"CN�VD��zD� {DO��C��_C��B`��A���@�A                                                AΥVCQ�zC@ {Bϰ�B9�_A���@��@��?!A                                                                    E���G	NG(�F�nQF?ғE��*D�:�C� �C�D                                                C�(1EG�YE;�|Eg�D���C���C�*B.�?A<�L                                                                    E��F�N�F��NFz�F\��F.�E�8KE��eEH�                                                D)7�D��;D���D��9D�jD~@,D0�CߟC�y`                                                                    2o�6�RKGR�uGI�GzNF��F}xFD!E�IpEPt�                                                Ad�AT�5A�ǥA���B�X@��                                                        A�A���A���A���A��jB=��B_��B�t�<�<�<�<�<�<�<�<�<�<�<�<�E��gE��Ef�"E�JD��#DX#rC���C���                                                {@��{@��{@��{@��{@��A��B�_�)m�`�f��7ZMB���?�g>    A �)A �){@�ξmb��mb�C���{@��C�3.C��1C���{@��@RD{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G\�?`D���G��D���D�MLD�MLAIH AIH FZ%-B.��Bǃ�C�{OFK��FK��D� �D� �F[��B.��                @C�rC��C�D�C��$?  C��!C�%C�%C��C��dC���C���C��yC���C�ݙC�� C��C�=�C�b9C���C���C��C�UC�U�C���C��	C��C�nC�(�C�C��C�C��k{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���Cƽ�>��I>�P>�+�>���>�6>>�ݷ>�S�>��>�Ć>�?c>���>�̑>���>�*x>��m>�˼>���>��U>���>���@/J?��'        A��������zGţ�ǃx�ǃ7�ǂ���N$��ۆ��ۉ��ۆļr        >�B'�m    ?͔�?͔�{@��@���FOX�6D��72        9��J?�/?ʞ>�qm>���>D?A�?r?��                                                ¾V|7�FkA�  ?k�;�vT<��L$��      b�     ��8�f�8�<7E�5��L4H,.4��4��3��x                                                                    7���7& K6<��4�S3|�}4N36��2��+                                                                    G���G�Fjt�E�C�E*D@��CO��BN-+                                                4G�d3��m2�u�10ʚ/�Ĥ0W�]/�sj/b�                                                                    3|2�#�1Ţ40_P�/|2/�W�.�Nk.5�                                                                    6��{6��E5�	�4�"|3@=[3��35��2ћx                                                                    7@�6��66�4��#3��4)�m3�*W2��                                                                    4�Б4^�q3��2�G�1~�2	�1k!�0���                                                                    73��6�A6_�5�3� ~4O1L3��3	��                                                                    8+H8��7*V�5��4`L{4��4A��3��                                                                    54텸4GM�2��t1�A2d��1���10h�                                                                    8+H8��7*V�5��4`L{4��4A��3��                                                                    7E7D�6�x�5Ց3�/�4R�-3�A�3��                                                                    4��E4\�$3��m2sh.19$�1�ݐ1>]�0�                                                                    7E7D�6�x�5Ց3�/�4R�-3�A�3��                                                                    ����3��3{`B2c&3U:2�o�2Y%�                                                                    93��9�;(6�;g:�>�:0i9q<w9/�l                                                                    5S�4��x4oj�3kX1�؝2I`�1���1=�                                                                    6�n�6g��5�m4/H2��+3?�Y2��r2��                                                                    ��I��e�1���ٳ��?���7x��V����                                                                    �_�s�5X���x� ���(�%�6%��u�:`                                                                    5�85g�4���3�$1��D2B��1��1	                                                                    ���kn��k�B�.�ǒ��̓�9Ų��                                                                    ��Ҳ���'@@�Ǆ/��Ű(����a'� ��                                                                    ��!������.�ֱ�2ܰ��v�%��z0��                                                                                                                                                                        /���.���                                                                                            7��7��7F	5��?4JPv4�%40T�3���                                                                    4���4A6�3��)2�80�)�1v/�0ǽ0[a                                                                    6�ѿ6�4�6��4��"3QY13�e�3AX�2�[�                                                                    3|��3=F�2�� 1i/߭�0�0�/�/p�(                                                                    7���7�O�7 �65�~4Xa�4��N4R\^3��                                                                    4wT49V�3�9�2F�0�1~1�I�0�H�0�ۛ                                                                    7rI*7:��6��
5kv4p4�x�3�GQ3~�                                                                    5�r�5U:�4ĉ03���211�2��X2�
1�(�                                                                    6C�X6��5���4>E�2���3�J%2�T2MFG                                                                    4_��4,N�3��2Yt1/�1��0�q�0j�v                                                                    7�D7d	�6�/)5��4=�4ƓW4��3�=&                                                                    5�7)5�N�4�5�3�r�2X�2��2,1�jt                                                                    45v�4�3���2��32	�3�2�f�3�B                                                                    3�P2�og2�"�1���0�sV1���1İ1��r                                                                    4]�543�r3��2�7Z2'x�34��3��38��                                                                                                                                                                        $��      b�     ��62M�)@@�                        6�͖3ʫ�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @خ�    @ض@    16:45:50        F�� @ض@    @ؾ     $�I      b�     ��5���7$�B'�fB'�f7�3�D"g�    (�ž+F�6�C?�  ?tj�3��4�61��/�O1g=-cM    /g(+
�!>��\>7`H=�lo=Ȧ:��];\�R7�\�2\�V:@R<#�Q<���=F�`=��N=��"=���=�=��=��N>��>J@�6�@t��@{ ?�ѐ?���?l?E%L?�6@eC@f��@] �@]�@^�@^��@_@@_�3@`s�@`�x@aS@a��@�m=h�?^�K                    E�@7���G���G�rFj6�E�\�E@hD@�2CO�BN@�                                                AmSC&.�B�D�A�}AT@�6�?��i>ٰc=��                                                                    E@^@��F��A��"3�4FY��B&H�A3w�/EV�                                                    {@��Bj�(Bj�(@O:C8cV�$�B�    �mJ/.y�&x+�>���(�ž@m��.�#3��y    >�_�B0��1��A-\^A{1Ay1[        B��B��C�HVC�HVC��[>�G@�P|2��;C�sx6nqu6��AgkCG��>�_�C[�B�hLCs>�A���B��CO�}B�B�hB.���@�A    B.��    B.��B:��=cB19 l�9,�A���A<h�A�j@`ĸ@53J@蔭?�f�F��\3�41&�.5^P3��8��#6�E5�#�G�G.�%G�zG��<=�E}    94�    >EH�>��>��^>ӕi>ƚ�>�<�>�Z�>�;�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�v�    A=t�A=t�{@��7�y�9���9g�*8��i8ZH�7���6�>6\�5,<�                                                A�e
A1�HA1�H{@��?s�cBf=�B�5�@1ہ>���{@��{@��    C��^BV��6�N    6;��7u>>@���2��7L\F    =�C�%*7��Bu�RD�?C���B�|�Bc�^A��A ��@)-�?��                                                ?�b�A�*�A��@ V�?���?~�><��=~)Y<U8@                                                                    D���F���FR�D��DHCq�B�_�Aϋ�@�+�                                                @�j7C
#rB�$IA`�<@�)q@Lus?���>�sc=��B                                                                    C�E�EnIfD&Cu��BɮsBփA2�)@+��                                                @f(6B}�4Bߚ@�[�@S��?�5?q>K��=As                                                                    7�v&8}�A�>&�A؟2?
c�?�i-,���+��->-�,��e,���+N .Ty]-�7e-t��*o�Q7�=L5 4>2�i�0��X            7���5
':2��2��"2�i�0��X��������7��H4B�    1�P�4�04�w|2���    =���8[��82�65G8�3�?�2q��                        3���5%��5~3�U�1��5�ON2�8�            5b�<�-!�p?�  ?�  ?~b�?{�?p��?wa�?u�p?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CiPG¶�>-��A�D�                                                B�                                          B��                    A�                @�p�    .���*��5�i6��X3��d3�/�?�
�8>:S@V��8��;�K�;�+`5�n�-�^        7� g    3�|�6���    42��{@��    +�	    ���d{@��{@��{@�ε��6���            5��    6F�{@��    6F�    6�^    {@��42��{@��4Ҩ4�"2�+    6�^{@��5�G�7�4h6�K6�K<���0a�YB�t6�)kE��9C�}ZA��    C*�U            >�-�5 �9|��5 �<��C<��^?���?�/?>Iz>��D>t��>P*��gs���
���Њ�*n��ld���5�s�
�R�p�9�D�&8� ���[>ƥ�<���<���?�"�?��o?8%)>�>S��=K���gs���
���Њ�*n��ld���5�s�
�R�p�9�D�&8� ���[=D@�/�d�()L�6	�6:i';i�Z;5>�<Ċ$<��=��>+�                                                ĝ��ę�Ē��Ć�y�mMb�a!��=�����̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         ;gLt        {@��{@��{@��{@��                3�%�6��        823�%�            CPa)D�H�D��YDS&�C�
CWjBb3.A�
%@��"                                                A�a)CRH�C@�YB�&�B<
A�Wj@�3.@
%?!�"                                                                    E���G	d�G91F�tF?�E��D�@�C��C�\                                                C�;^EG�E;�Ek�D��	C���C �B.�1A<��                                                                    E�8F�R�F��Fz8F\�F.̾E�9	E���EH*                                                D)8�D��D���D��oD�jDD~AD0�Cߟ�C�y�                                                                    5�ON7h5GR�0GI*�G�tF�&�F}��FF�E�KpEPu-                                                8o�w            ?o��                                                            @v��Ay�gA��xB=aB6 Bh��B��B�(<�<�<�<�<�<�<�<�<�<�<�<�E��|E���Ef��E�D�DX(<C��C��                                                {@��{@��{@��{@��{@��D�B1E.h�k�f�8C�1@(�    Aa�qAa�q{@�ξl8�l8C�W|{@��C�7�C�8C�W|{@��@O:C{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�i<ҧ^D��0G�D��D��D��AF�AF�F[;>B0ϸ?o�wC��FK�	FK�	D�,�D�,�F\=�B0ʠ                @�/�C�g	C��_C�Va?   C�n�C�7�C�7�C�~C��C�q	C��SC���C��C���C�n�C�9�C�vC���C��_C���C��C�C�4�C�c�C���C���C��>C�C��C��C�.C�ig{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�}RC��Q>���>�r>� �>�fQ>�G>�^>�JP>�,�>���>���>�@>��E>��>�{�>�7'>��>�<>���>�J:>��@f<?���        A��Ǘp ǖ��ǖt��-(�U����>�Ɓ��ƿ��ƿ��ƿ��ė�k        ?R�B*�k    A�A�{@��@?�xFO/6���7�        ;t|�??�s?@�t?B��?E�.?G�b?K�?Si�?b��                                                2�7w3�A�  ?b�@�-@��$�I      b�     ��9\�8�p�8
}7ne:6�]_5���4�z3�!�                                                                    8�P�7�o75�L6���5�M�5@�4��3��                                                                    G���G�rFj6�E�\�E@hD@�2CO�BN@�                                                4�7�45Ί3���2�٪2' �1]%0j��/s'�                                                                    403e��2���2��1R�H0��
/�J5.��:                                                                    7hv#7G�P6��>6F�z5��4�Z�4so32�                                                                    7��7V�(7"�6���5�z5��4N��35�                                                                    5-gu5��4�R4`��3���3Ͱ2,jh1��                                                                    7�M+7�[70*�6��~6Ē5C�4|v�3]�                                                                    8�k
8��8 ٭7o9x6��5�8�52�4?�                                                                    5��5���59C/4�,�4R�3gע2�1���                                                                    8�k
8��8 ٭7o9x6��5�8�52�4?�                                                                    8
�7���7{�6�6�6b*5V3�4�|93��c                                                                    5u[5
s4��4'G?3�tN2�H2V�1�U                                                                    8
�7���7{�6�6�6b*5V3�4�|93��c                                                                    �G�糪5F��5Oi~4���4��3^��2���                                                                    3|˔3�+�:L%�:%�9�L/9�d%9'�8�(                                                                    5�mu5�~�5<��4�]E4�3=j2p�1SO�                                                                    7&n�7�:6��G5��B5X]4CC3a^�2d&s                                                                    �-��4��糵M=봘L��ɲ�Բ��                                                                    ��²�jβ�y����Z��:�@��Rڱ��                                                                    6%��6��5��"4�ڃ4�s3F3c2eX�1g��                                                                    �6�P��F��+��u��Z'���ʳ��;��U�                                                                    ���^�����#����뫱���.��d�����W                                                                    ����z�(�+��1>�����)q�A��Z�p                                                                    (��(.�.	�V-�^�.�f�.=r�.��U0�zk                                                                    *q*BϞ                                                                                            8�@S8}8�j7]��6���5�-�5	�l4$�O                                                                    5�4��4�3�$3-��2}�h1��[0���                                                                    7�;7w��7 ��6\:5�-?4���4��34��                                                                    4�%3���3�`�2�o23�!1�10�� /́�                                                                    8�5�8r��7��?7Zx�6���6 ��5$954Dm�                                                                    5�E4�t4U.3���39��2�zx1��0ހ�                                                                    8�-7�׋7�F�7&�y6v�85�a\4�Ų3�Sv                                                                    6 ��6t5�q5>Wf4�i3�& 2��1�>                                                                    6㨰6�w�6�].6��5GhB4�r�3�(�2� 6                                                                    5@4���4��4ϭ3c��2���1�.�0�۬                                                                    8,*�8�7�9�7K�?6��?5�Z�4�*�4��                                                                    6D�*6,e�5�A�5h��4�XH3�g|3<�2b�                                                                    4��4�G�4���4��	4Z�4[X3��b3�U                                                                    3�k:3���3�K13n��303�2���2�W�2Nѹ                                                                    5 �4�tl4մ�4��(4�@�48��3�k[3�h4                                                                                                                                                                        $�I      b�     ��60�1�P�            G�^    <��6��3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ض@    @ؾ     16:45:53        F�, @ؾ     @�ŀ    $ӭ      c     � 5��7ek�B'�fB'�f8!:GD"g�    )��s+6�7 ��?�  ?sYx    4�T))�m/�����        /ƹ"+�1+>��:>5T=���=c�:�<
5ni�.	؍2d��:5&�<#2�<��=F��=���=��=���=��=똁=���>�>�@ޙ@�9@3�r@�R?�N�?�� ?��?���@in�@igP@_*~@_@_�@_@_�@_�@_6@_�@_+�@_5c=���=��w>��                    E��d7j��G��uGP�FiuBE���E�WD@F?COY>BM��                                                A�nC%��B���A���ASyg@��e?�v[>�9�=��=                                                                    E@�@�&&F�YA�:0�Y6FZ�B'�AF�n//�                                                    {@��B��B��@RR�8��`%�\.    ���v/*�צ
?   )��sA(d_���1�,f    >�y�B?mB1؄zAH�2@�T @ԭ        B��^B��^C�R�C�R�C��G>��@�*c2�6�Cζ6�D�6�IA�CO�:>�y�Cd?�B�tfCr�%A�,wB�b|CN�B�B�.fB[�?�    B[�    B[�BSh.            A�zLADi�A�Z�@v�@(�@�ެ?��F��U0�Y6.b�Y+��=0�,8��e6�� 5˃�G1aG2WnG쮮G�7�=��,            >�q>��>��>�n6>�#�>�-�>�>�K�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    AM�AM�{@��7���9n�I9Pu8��8H�o7��b6�	6�5 z[                                                A�q�AGYMAGYM{@�Π'<?B�b�B�B(@00?	 �{@��{@��    C���BS��6��k    6;��7uK�@�z
2�S37B�    =��,C�l�7�v�B\�~D���C��B�N�B*��A�:P@ч�@}>�U~                                                ?t��A�Z�@�J?�n�?bt�>�ߵ>bw=Q<�<>�m                                                                    D~�.F��$F��D�X0D��Cj<�B�m�Aɫ�@�}�                                                @�SC�B���ASg�@�cZ@A��?��|>ƙa=�Xz                                                                    C���E��Ek�tD+YClN�B�%�B�%A-�@('                                                @`gB{OB^<@�Z<@F�?�P�?�>B��=;��                                                                    7h��7��A��>)��A��?�t?�uT,���+�Y-GO�,��,���+Q8w.Y��-��-��*y��7���4�{�2D-�W�*�5�*��L&�<7�/N4�(�1�{2�2B�-�U���/N1�,f7�04s�~    3��4w�4v��/��    =^8Df�8*]�6���5j��3#O�                        5��6��96��4L%)3�J~7�R�2̪�            4�pq>k��-$I�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Cm��G��>+��A�H                                                B�                                          B��                    A�                @�p�    /u��+[�5?��6z߉1�X`1��?�.�8I~@:c8�L;��;�~�                7��        6�N        {@��            ���M{@��{@��{@��6��H7�'�            6�]c    6�]c{@��    6�]c    6�N    {@��    {@��                6�N{@��6�3w7��'7�#7�#<�w{2��B�	�6���E�1�C��NA܉�    C4            ;���4���6Fx�4���<��?<��1<�V�<�O<-�;�ǵ;�!:�Խ�������ڨ���x��曊��b�u��T�5�;k~�'{K�`��
P;��l<���<���<��8<��<+s�;��;	�:����������ڨ���x��曊��b�u��T�5�;k~�'{K�`��
P8�,    )��6�\6uU�9K��9\�?9��9u�8��	8��
                                                �4_��4�(�1��(˻�S��
4*��If��<S̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �'<?        {@��{@��{@��{@��                                                    CQ�6D�Y�D��aDU{C�f>C�3Bc)�A�Ͻ@��G                                                Aщ6CRY�C@�aB�{B>f>A��3@�)�@Ͻ?!�G                                                                    E���G	�!GRF�wFF?ͼE��D�7�C� �C�A                                                C�I�EH�E<TEnD�~ZC��ZC�B.�A<�^                                                                    E诇F�W�F��yFz�F\�OF.�zE�7�E���EHF                                                D)9�D��DD��SD��D�k"D~@�D0��Cߞ�C�x�                                                                    7�R�7X�1GR��GIGyG��F�*F}��FFE�J�EPt                                                                                                                                @*��A@�5A�|�A�	�B��BJ�OBs�B�F�<�<�<�<�<�<�<�<�<�<�<�<�E�&E�bEg�E�WD�
�DX'�C���C��T                                                {@��{@��{@��{@��{@��D%�A���/��e�"8sCV�@F��    An�lAn�l{@�ξ�L���L�C���{@��C�C�C���{@��@RR�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��;'�@D���G�cD���D�xPD�xPA?|�A?|�F\��B3��    C�K'FK�;FK�;D�5{D�5{F]_B3�|                @��|C�|eC��OC�ޮ?   C��C�C�C���C�^[C��C��OC��C�6�C���C�zC�.�C��C��C�JC��C��C���C���C��C��2C���C��C��C�C��C�@C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�C�j6?�>��7>�>�>�F�>�8�>�w�>�&n>�=1>�E�>�V>�W�>|>��>���>�SC>���>�WZ>�:5>�$:>���@A��@�Y        BrǡCǠ��Ǡ&,�J���/�����Ə�|��E��r��E�/,_        ?��B0�    A�A�{@��@���FOA87F�7"�        ;�\Z?wK?��?��?!I�?$W�?(ș?.�?5
                                                �.��70KsA�  ?k�A	��A z�$ӭ      c     � 9Cic8�#�88�7b�"6��C5ҵ�4�`3�f                                                                    8v��7�:�7(G�6� �5���5R4��3
b�                                                                    G��uGP�FiuBE���E�WD@F?COY>BM��                                                4���4#�O3��e2�2(g�1\�"0k�/e�f                                                                    4^3N�q2�aX2 �1T��0�pK/�u�.���                                                                    7Nچ74z�6�2=6>��5���4�2|4j^3(H�                                                                    7���75>�6��W6Z�g5�1�52�4+�3�s                                                                    5�/4�-�4���4!P3��2�m�2�a1 ��                                                                    7��7]�j6���6��j5�<�5 Ze4R#�3<7                                                                    8���8�Ĕ8*�7\�6��5��55��4cD                                                                    5�N�5�T�5w�4��4B3]��2�b�1��                                                                    8���8�Ĕ8*�7\�6��5��55��4cD                                                                    7��7߄#7[6��6�v5Or�4{�3p+�                                                                    4�K�4�CU4�/c4��3�X}2آh2�A1�                                                                    7��7߄#7[6��6�v5Or�4{�3p+�                                                                    �,У�]5�t/5nIt4�^�4,�3r�_2�                                                                    5n
5o�~7N@�7&w6�O6$>�5�za5�                                                                    5���5z�@5+z4�7`3���3�j2H12��                                                                    7��7�6r��5�
-5
y�4?93]�2Tj                                                                    ���	�D�� ϵC������3����ʱ���                                                                    ��2���[���-��� ���(�:�;��뱕_                                                                    6)\6�;5s��4��4B�3Am32`]1WFl                                                                    �"�������N���^�ﴡhٳ�o8��$z                                                                    ��M�����w��؊���y�.�fpd��6�                                                                    �~ȋ�b�������>��\7�)(��Bf2�N��                                                                                                                                                                        +]�:+���                                                                                            8��8d��7�iS7U�(6��5�S5
˾4�;                                                                    5]4�v	4y"�3�n311	21�8X0�P�                                                                    7�7�7`C�6�c�6T=�5��/4��m40>3*�$                                                                    4�t3⾪3tV2ق�27Z(1���0�c�/�S�                                                                    8}-8[��7�^�7R��6�_E6Q�5%��49�                                                                    4�۔4�K4n��3��3=�H2�)1���0�WG                                                                    7��h7�3�7�H�7�6yf�5��4�\�3�z                                                                    6P�5�x5�S055e�4���3�M�2ڲ�1苇                                                                    6ʫ(6�]w6yY16 B�5I�[4���3���2�m                                                                     4�
4�!�4�|4�w3fS�2���1���0��I                                                                    8D�8f�7���7A��6�iy5˄�4���3���                                                                    6/)�6�5ׂ5]�D4�/f3��3�B2a                                                                    4��l4��44�=�4��v4\nB4��3��3q�l                                                                    3���3��[3��3c��32 52��E2���2Cs�                                                                    4�Mh4֕4�K�4�4��(49+w3�&�3�Ϧ                                                                                                                                                                        $ӭ      c     � 60;�3��            G��    =���6�?3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ؾ     @�ŀ    16:45:55        F�j @�ŀ    @��@    $�      c5     ��5��`7��B'�fB'�f80iD"g�    *!9+7׵7GKI?�  ?v=    4<z)(j�/�_��'^7        06N�,.�]>���>%LS=�z7={�:��5��.9�2_5�:%�}<"eY<��=F��=���=��V=���=�=�r=��>�>>�@�q�@�x�@D�$@�?�H?ѳ�?¬l?�ي@jm.@ji@`$W@`C@`n@_��@_�@_ѫ@_�d@_��@_�\@_��=�[�=���? �S                    E��7gaG�M�G��Fh�E�H�E�
D?�[CO�BM�3                                                A\�C$��B�H�A�CAR�@�z�?��>���=פ�                                                                    E@5�@�F-F�|A�-�0M��FZ|?B)�AOf�.]�                                                    {@��B���B���@`�8��~%�;v    ��,f/>���"`�?   *!9A28� �	2$    >�y�BAp�1���A&aA�[AJ        B�)=B�)=C��C��C���>���@�X�3K�CԄ*7��6�ںA!��CB%>�y�CUEuB���C`��A��%B�~C<��BM�B�7fB"GT�$��    B"GT    B"GTB,�}            A���A;��Aޢ@u3�@rW@���?��MF��0M��.#//+���0��8�Z6��5�:�G 
�G7"�H�TG�i>��            =�$K>��w>�o>��3>�ƙ>��>��g>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�    A�(A�({@��7�C&9o�9M�x8��8>�O7���6ްd6
��5��                                                A��A�\kA�\k{@�Π>!�B~gB�8�@=��?w[{@��{@��    C��+BS�T6�U    6>`b7u�@��2��7�$4    >��C�)G7�vBL@7D�#PC�SB~�BB
�)Az�@��?�5>�C�                                                ?_�TA�Zi@�po?��?6��>���>|%=9�<7FX                                                                    D|�F���F��D�S3D޸Cd_CB��4A�l>@�x                                                @�.�C��B�I�AH��@�4�@9u?�>�=x=�b�                                                                    C��lE��Eif�D�eCd��B�JLB	,pA(֋@&�                                                @Z�Bx��B �@�	�@<�?��Z?��><+=8�L                                                                    7_�\7���A��>.�3A��?~�?�ܢ-	bY+��f-���-<(�-	bY+߄�.�.G-*���7�=(4��K2|I-��/�Z�/��E+hӥ7��4�nh1���2��2e-�r���2$7��p4�1�    3��4��}4�*�0U�    :k{�8.K8$��6�&P5l�3n�                        5�L6���6���4�_�3���7i�2���            4́(>f��-�p?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CjH|G���>+�sA�%q                                                B�                                          B��                    A�                @�p�    /�!+˴�5C��6vc51h�1��?���8D8 @wCS8ҳ<�"<�)                7�]{    2��:6�ǫ        {@��            ��
�{@��{@��{@��6(l7]�*            6'�
    6'�
{@��    6'�
    6�3�    {@��    {@��                6�3�{@��6�p�7�4C7b�7b�<�4)2x�B��o7/�E�C%Dt�A�l�    C1$�            ;���4��*6AdE4��*<��k<�&#<�C<��<&}=;�� ;
�:�5?����q^�߁7��'~����� ϊ{q^�YvΊ?�T�+1���)0;�#<�К<��<���<�UK<!�n;�;��:t�ˋ���q^�߁7��'~����� ϊ{q^�YvΊ?�T�+1���)096��    )v��6 �6���9�Ʉ9ӳ�9���9Y_�9�8��                                                ���Ÿ�cŵh�Ţ�lŎ���u-�Q���8p�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �>!�        {@��{@��{@��{@��                                                    CQ��D�G�D��DT�C��jC Bb��A���@�}�                                                Aч�CRG�C@�B��B>�jA� @⮱@��?!}�                                                                    E��G	�7GhPF�x:F?��E��ED�/C���C��                                                C�U�EH-gE<:�En�D�u�C��|C�B.�
A<�                                                                    E��F�\�F��jFz�F\�fF.�NE�7E��CEH�                                                D):�D��PD��D��D�k�D~@qD0�;CߝxC�x�                                                                    7i�7U�GR��GI_�G�}F�&7F}��FDE�I�EPsf                                                                                                                                @	��A&UA��XA�/ B*/B:(OB_h�B���<�<�<�<�<�<�<�<�<�<�<�<�E�'�E�EglE��D��DX$�C���C���                                                {@��{@��{@��{@��{@��D(��A�&�/�;��gm8�Cyl�@G�@    A��[A��[{@�ξ�����C�G�{@��C�6�C�6�C�G�{@��@`�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��:���D��GvD���D�� D�� A:L�A:L�F]��B5��    C� XFK�FK�D�8�D�8�F^a�B5��                @��C���C���C�a?   C��C�6�C�6�C�l�C��C���C�M�C���C�ysC� �C��C�FC��#C���C�B�C���C���C�J�C��C��
C��C���C�հC�\C��C�C�SC�j+{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C�%_?��?��?�?��?��>���>�2>��?>�V>�K�>�6�>�F>ͣ;>�Y>���>�2>�Z�>���>�t�>�8@�?ğ�        A� �Ǟ��Ǟ)ǝ�'�3���������ƚ�Q���#���C���Ů6        ?*�uB�    Ah6�Ah6�{@��@LݰFOkD7^G7@>        ;kB�?��?��?`<?
�?�
?�]?�?�p                                                ��7�:A�  ?k�A	��A z�$�      c5     ��9C��8���8��7\sH6�l�5��4�e�3��:                                                                    8wEn7�4=7#٪6�;V5Ŗ�5|�4	M�3��                                                                    G�M�G��Fh�E�H�E�
D?�[CO�BM�3                                                4�14!(k3��2�
2#�1Vѓ0c�+/[�                                                                    4�83K�_2��X2�c1O�0���/�ׯ.��                                                                    7O�)72P�6�3m6:b 5���4�4�4
I3!�                                                                    7�27)��6��26-��5�� 4��K4�3v�                                                                    5N4ψx4W��3��3q��2�;�1��_0���                                                                    7�F�7Ol)6�M�6Th�5���5	3.44�83,�                                                                    8���8��q8E7Q��6��5���5>44�                                                                    5�#�5�6�5!�4�d4vO3P+D2���1��A                                                                    8���8��q8E7Q��6��5���5>44�                                                                    7���7�mz7P�6�n6�X5Dҽ4n*�3cW'                                                                    4��V4��4��4
�Y3��2�/L2��0�ҏ                                                                    7���7�mz7P�6�n6�X5Dҽ4n*�3cW'                                                                    �+�\��q�5�
�594�u�40|3uU�2��                                                                    5��5p1o7P�Y7!�6��6{5�n�5
�                                                                    5��%5k��4ӵ�4\��3�Nn3d\2,RD1$��                                                                    76:72�6i)5�� 5K+46��3R��2I�4                                                                    �}��9���ӵ>�w���b��u���ٱ�B�                                                                    ��y�����8���7V��~�6��.P��Ϋ                                                                    6j�6��5ibs4�4�4��38�(2Ux�1LD                                                                    �#S��Nඔ��l=�Y\Ҵ�޳�@5��6�                                                                    ��"[���@������I��=��*1��_o>�w.5                                                                    ��-��`"ִ���2���Ĳ$ⴱ<�L�Fb                                                                                                                                                                        +A��+�{�                                                                                            8���8b'u7�c}7P�,6���5���5��4�"                                                                    5
L4䧱4s3�K3-�2x!�1�p�0��f                                                                    7��p7]��6�}�6O_�5��M4�H&4� 3#>�                                                                    4O3���3n72ԅ�23	�1�6E0�&�/��(                                                                    8~W�8X�7旜7M�6�2�5���5 �I41�M                                                                    5 ��4�V94i$h3���39�2�[�1��0�/�                                                                    7�Y�7ܢK7��j7�,6s<a5�O)4���3��$                                                                    6305�'15�)�50�W4���3�
2�8�1��N                                                                    6��K6�JB6svS5��5D��4�(�3�T2��<                                                                    4��4�4�4�3`�&2��]1�~0��                                                                    86�8��7�I7=!S6���5�`�4���3�J�                                                                    60>�6�5�k�5X&4��3ⷷ3��2*�                                                                    4��4�t4�k4�P�4Wg4��3��M3g��                                                                    3�q�3�*
3��3]�>3-��2�2�P�2;I�                                                                    4�4���4��N4��B4�h44��3�S�3��                                                                                                                                                                        $�      c5     ��62R3��            G���    =��86��53�4]@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�ŀ    @��@    16:45:58        Fƨ @��@    @��     $�u      cT     ��56\7��B'�fB'�f8��D"g�    *���+H�?70�?~�?v�    3��))�/'Z���        4��0���>�_�>Yށ>�=�6�=\�< [{5t��2L��:<!�J<�m=F�,=��A=��h=��=�~=��=���>�N>�@�R�@��@ke�@7e_@I:@�?�G?���@nlc@nI�@c�&@cFh@b�@b��@bT@a��@aN�@`�Z@`�s@`��=>��=�G/>��                    E���6�UeG�+G�[Fh��E�+CEx�D?��CN�DBM��                                                A9�C$�&B�!�A���AR԰@�d�?��e>س�=ר�                                                                    E@e�@�k�F�A�L�/�FZ�B*E�AZ�-=܃                                                    {@��Bd� Bd� @g��8P��%tO�    ����/3ɋ&���?   .M��A/s)��01�<�    >�y�B�f2��Aw�?�~�?�^p        B��B��C��C��C�4�>�3@��3lHC��6�>Q6�w�A$��C�j>�y�C&2�B�i(C>%A�&�BKmC}�BF_B�o�A��l����    A��l    A��lA���:�
L    5�Ax׹A�A�ژ@O��?�Uu@�#�?���F�0/�,ݚD+�0c�)8���6[�]5^ �GgXG4�G�f�G���=��            =���>E=�>��>���>uH9>rP�>r\B>}�C8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��U    AG�"AG�"{@��7o��8���8�I8m��7��7.��6��j5��4�.                                                A��.ABp�ABp�{@��8_��BD��B���@K >�!m{@��{@��    C���BT��6��6    6B�#7v�@��2���7V\    =���CZR�7�_YBHQD��QCͺ&Bv�BȧAq[o@���?�| >�$�                                                ?Y}�A���@���?��v?,��>���>��=3�l<9��                                                                    D|$�F��F�D���D�Cbd�B�P�A�f4@��                                                @�vC�B�D�AEu�@�vl@5�f?���>�*�=���                                                                    C��E�^Eh�vDL`Cb��B�H�B�LA&Ξ@%��                                                @X�rBw94A�/V@��<@96?��[?s0>9l^=8                                                                    6���7p�iAK�>3�1A�z�?�3?���,AI�*�]j-=�,�g�,AI�+��.��-��E-V�;*^I7�h�4-�1���-_�/TK/S��*���7�.4�}�1�h�1�Κ1�b�-W�\��.1�<�7���4��/?�3��Z4'4&��/���    :Jlk7��y7���6�r5B3��                        5�l�6���6�m
4QS�3�u�7L/�2ν�            4�V�>~�-	�D?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CM��G���>*�}A�oq                                                B�                                          B��                    A�                @�p�    4,1�0�g5'�6!�/� �0�g<?�VO8'�}@4G7Ծ�;��;�[                7���    1d9M6z�N    1_*{@��            ��Kh{@��{@��{@��7�q7���            6�QT    6�QT{@��    6�QT    6z�@    {@��1_*{@��1_*1M��/���    6z�@{@��6��7w��8=�8=�<��l2�B�}6��]F-�C��&A*3�    C2            ;�/14�j�6\4�j�<@�;<B-�<���<j�%;�h;^O�:�WI:^��m����	��qϊę����5���5�u�	�T�i�;<��'Q��;�	�
;���<@�N<B'�<�<f� ;↭;U�t:� /:MV��m����	��qϊę����5���5�u�	�T�i�;<��'Q��;�	�
8��p    )ΛY5;1�5�&19Ym~9�d�9E�c9TJ8�q�8��                                                �-dƽΡƾ.:ƼGƫ*dƃ!�'������̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         8_��        {@��{@��{@��{@��                                                    CQ#ZD�G:D��JDT�C�SwC��Ba��A�)�@��                                                A�#ZCRG:C@�JB��B>SwA���@��@)�?!�                                                                    E��G	��GqF�y�F?�6E���D�*�C���Cˬ                                                C�[�EH:E<G[Eo�D�u?C��/C�B.��A<�(                                                                    E�iF�^#F���FzF\��F.�}E�6�E���EH�                                                D);5D��JD��?D��D�lVD~@�D0�-CߝC�x�                                                                    7L/�6�GR�8GIh�G�bF�$dF}~IFB�E�HpEPs(                                                                                                                                ?��K@���Aq1�A���A��BrmB=XBm��<�<�<�<�<�<�<�<�<�<�<�<�E�.�E�pEg�E�D��DX!�C���C���                                                {@��{@��{@��{@��{@��D7+$A���0=�E�e�T7�g�C��@OP    Ak� Ak� {@�ξ�)���)�C���{@��C���C���C���{@��@g��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G)::U2�D��dGF�D��<D��QD��QA8W�A8W�F^��B7?�    CZ)FK��FK��D�9D�9F_01B7:�                @{�C��C��$C���?   C���C���C���C���C�t5C�d*C�QdC�8>C�>C��lC���C��WC�b�C�(?C��C��LC�FJC��0C���C�T!C�'C��JC���C�	C��C�C�fC��_{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��iCZSp>���>���>�ς>�z�>�>��>��>�:�>��>���>�<�>��>؆$>�m�>��1>�(>���>��>�r!>��l@@:�@��        A��3ǧ��ǧ��ǧ6���ʘ����e��
�
���3�[�3�s�3�Zư�<        ?bPB\�    A�6A�6{@��@�҃FO�7ppt7��        ;��^>�m�>��>�4t>�-<>�-t>��>�G>�c                                                ���6�~PA�  ?k�A	��A z�$�u      cT     ��8�A8*z�7�ۿ6��A61��5|��4�`i3�*                                                                    7�
m7WW�6��]6�5`��4���3� �2�                                                                    G�+G�[Fh��E�+CEx�D?��CN�DBM��                                                4G0�3���3��2�f71�l�1��0!�i/*�                                                                    3{��2Ṥ2Bn�1�s�0�{�0'i/LI�.W�v                                                                    6��6��-6cO5А�5-��4��3�с2��                                                                    6��q6��G6*��5���5%�4��t3�b2���                                                                    4��[4a�3��3�S3 ��2W�1���0�e�                                                                    7�a6��6P�5�٦5J�c4�.�3��3�}                                                                    8$�\8$��7���6��267~�5�p'4�v�3̓�                                                                    5��5�g4�n4�3�r�2�^2;IV1E��                                                                    8$�\8$��7���6��267~�5�p'4�v�3̓�                                                                    7t �7v��6��Q668\5���4�R4%'3.�&                                                                    4q	�4���4�3���30G2r�\1��90�L                                                                    7t �7v��6��Q668\5���4�R4%'3.�&                                                                    �3\	��~5U�5?�4�43ߠ32��2���                                                                    4�G�5-�,7$L6��,6nx�5��5mP4�<��d��.�˥��H%������ê�e.�F${�.�ڂ�
��k� �k                    54�C5ZO4d�3��3J�72��1�A�1 �                                                                    6�7a6�z55��15=;�4��N3ܜ�3�[2?                                                                    ��µ��=�RX��ջ��*���v�����k���t                                                                    �a�f�mM��]3��=����q������p�^m                                                                    5�j�5�ג4���4>�@3��2ެ�2n�1�@                                                                    ���Y���ƶ'����ƴ����A���j��G�                                                                    � s��Oq�����������Bv�<��?�d                                                                    ��(���䄳�f��01����˖2�u���                                                                                                                                                                        +��+�	5                                                                                            8 Л7���7��$6��67��5���4��h3���                                                                    4�=F4}�t4	'+3o^�2�Z2�1W��0��5                                                                    6�aU6���6��5��S5>)�4���3��2�F�                                                                    3+�3x��3[�2m��1�-�1 �0l�|/�r�                                                                    7�!v7���7��6�;�6D�>5���4�fs4	�D                                                                    4y�4sj4�N3k��2�S2'�(1��0�l                                                                    7uDS7u �7*��6�pT6
e�5Hg�4��?3���                                                                    5�'5� �5CX4�7<4*�3e�2��$1�?~                                                                    6F26E��6	�5�'4߫�4!�3UC�2t�I                                                                    4b�U4bD4�3�,�2���29.1s�a0���                                                                    7��7��w7P�[6��6)&�5t�4�H13�Gu                                                                    5�L*5�5nl�4�C�4AP�3��B2�R�1ӿa                                                                    47O�4@�43� 4��3�� 3�a3�o�347(                                                                    3!�3�
3�2��2ž92�`�2R�y2��                                                                    4`a4k\4[c�4<�4�3��z3�l%3\Cj                                                                                                                                                                        $�u      cT     ��65i�3��Z                        6�"N3��0@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��@    @��     16:46:00        F�� @��     @�܀    $��      cr     �`5��b7� �B'�fB'�f7݌�D"g�    )���+g��7)�9?�  ?~i&    3��-(���/�t���        0'��,!F?>��]>N!J=���=�y=�7:�2�w�2?�j:�"< ��<�=F�l=��7=���=��0=�|�=�`=�ں>�[>�@�&@���@J,B@C�@��?�3?��?��X@p��@p��@f/�@e�@e��@eIX@d�@d�l@d(k@cΙ@c��@c\�>�Y=n��?�4�                    E��a79TG�MG�YFh~�E�#
E�D?��CN�BN.�                                                A DC$��B�A��AR��@�k�?�H>��/=��                                                                    E@��@ďqFA�jx0���F[\�B+�AN�-���                                                    {@��B �B �@k�8F�V��n�    ���p/*�&���?   )���A����1���    >�y�A�12��@�q
=W3�=��Q        B~��B~��C��C��C�:>�?@f�3
�C܊%6��6��$A&}B�H�>�y�B�	�B i`C�hAb��A�QB���A��-B|�6?ܒ���    ?ܒ    ?ܒ@j�u:��u4��5b��A.@�{A��[@��?�o�@]�?q��F�>-0���.l�|+�	�1*�85�6�5�E�F��*G'��G��MG��>%��    6)�    =�y�>���>���>��s>�Z�>��a>y�>tl8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M̵�    A,'�A,'�{@��7��9:�9&V]8�=D8�77z�6��5��?4��                                                A���A,|jA,|j{@��:���B	hBi:@Xqc>��{@��{@��    C���BU�o6�'�    6E�7w]�@��32섖7Z�    >%��Co�!7��BFpID�OPC�]�Bqx�B��As�@�@�?�ߖ?�                                                ?V�A�c�@ۄ�?���?(�>��P>	,p=<��<Io�                                                                    D{�F�ѪF��D���D�Cb�B���A�6�@�f�                                                @�5C�B�q)AB�h@�E�@5$?�e->�S�=�@�                                                                    C���E�&�Eg�CD4fCa.�B��$B*�A'��@'�                                                @WoZBv�A���@���@6�l?��L?�>:�=;��                                                                    73�7��EArR>8GaA�?��?�u�+�J�*��-f�,���+�J�*�k�-���-f`�-wb*05N7���4W�1�#.��/}��/|��+�u��Wc4���1�h�1�\1�=.��6�Wc1��䶈R%4M�U    2�24O{T4N{�0u�    :v��7���7���6*��4�}�35�                        4���6�J6��45y�336��12�~            4�e�>ȯ,���?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C�G�c9>+yA�{K                                                B�                                          B��                    A�                @�p�    /�[�+��4��6BԴ0�{�17_�?Q�a7�ɘ?�|�7q;�U;���                7W�     2��6�{    1a�:{@��    0�    ��{@��{@��{@��7Jq7�Ϗ            6���    6���{@��    6���    6w_    {@��1a�:{@��1a�:1@�?0��    6w_{@��6[�6��7�'`7�'`<��#1�ئB��6�w�F#�>C�B�@���    B�Mg            <;`�4q��6��4q��<��J<�3�=J��=H-<���<k�;wf9:���Bf��2�ؠ-������=<��jȊs�2�R�\�9�݊%�������<4[P<��<�.�=F�
=\<��G;�U�;c^�:�7�Bf��2�ؠ-������=<��jȊs�2�R�\�9�݊%�������9�;    )� |5�L&6�:��F:{3:T3:9�<�97?r                                                ũ�)ů֬ů�:Ŷ�j��h}��Đ��q�#�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         8��P        {@��{@��{@��{@��                    3r��        :��q                CP��D�E�D���DS�C�CU�BaiA��w@��                                                A���CRE�C@��BӚB>A�U�@�i@�w?!�                                                                    E���G	��Gx�F�{xF?��E���D�,�C���C��                                                C�b�EHF�E<RSEqD�v�C��zC��B.��A<�t                                                                    E�QF�_�F��0Fz2F\�hF.��E�7�E��\EH�                                                D);�D��GD��D��kD�l�D~A5D0�Cߝ�C�y`                                                                    6��17+��GR�GIp�G��F�$7F}~�FB�E�H�EPtD                                                :��                                                                            @�0A%�2A�`>A�uPB	5-B*��BB��Be%N<�<�<�<�<�<�<�<�<�<�<�<�E�5_E�%UEg�E�D��DX![C���C��U                                                {@��{@��{@��{@��{@��D.p�A�`^/�+�f��7�o�C�9h@O��    A`.A`.{@�ξ�/���/�C���{@��C��;C��AC���{@��@k�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G5�:�D��
GXQD���D�knD�knA7A7F_�B8:��Co��FK��FK��D�:ID�:IF_�_B8�                @)q�C��C��AC��?   C��LC��;C��;C�ƿC��mC���C�� C���C��nC��C�lC�S C�4C�zC��}C���C�x6C�5qC��MC���C�n5C�>C� C�uC�C��C�yC���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�ԹCo�>�ms>�k>��N>�m|>���>�`>���>��>���>��>��H>�Z_>�t$>�+>��>̶�>�7�>��c>�2>�I�@N	 @��        A}�y�@/h�?��?l���ǰُǰKjǯ�P�I+D��d���d���d�Ż,        >��A�    A�UA�U{@��@��4FO�7g�)7�        ;��Y?�?��?�?��?�C>��>�Q~>�'�                                                �	��6���A�  ?k�A	��@�w$��      cr     �`9��8wxg7�T�75ŉ6|��5�"�4���3�k�                                                                    8?�7�L7��6e�35���4�r3�p2��                                                                    G�MG�YFh~�E�#
E�D?��CN�BN.�                                                4�Mv4�H3a�f2��Z2b1&��0%\/_                                                                    3�9m3#��2���1�!1'8�0R��/P{g.A�i                                                                    7!�*7�6�|36^5vc�4�U(3�cu2��<                                                                    7D�'7��6B"6�?5sa�4��3��!2��                                                                    4ן�4�4'�3�ou3<){2�@�1�l0�r�                                                                    7pB7#H�6���6$��5���4��4�i3G�                                                                    8���8pSM7���7+3�6��o5�HS4Ȝ�3�[A                                                                    5Pt5T�4��4go3Կ�3�)2B�d16t                                                                    8���8pSM7���7+3�6��o5�HS4Ȝ�3�[A                                                                    7��<7��7-#�6��U5�R]5�4+^�3 <~                                                                    4���4�:4d��3��z3N�[2��1��0��                                                                    7��<7��7-#�6��U5�R]5�4+^�3 <~                                                                    �G�3��tl5��!5U�$4��Q4�N30e�2np�                                                                    4d4�,x7�4�7�n�7�6�qP6�r5��.                                                                    5��,5:�44�J�4,��3�M 2���1���0��3                                                                    6�36�Һ6B��5�S4�B�4W�3��2�;                                                                    ���N��`ɵ�N�[�r/.��H������g�                                                                    ���ᲬdĲ�D����V;�Ou���G �                                                                    5�8�5��A5B��4�S�3�Ee3��2��1��                                                                    �������t�ҵ�	l�.�
�r$X��%ײ���                                                                    ��n��cq��_w�]f�����8T�!Z�+�                                                                    �G}��4�z��g����C�ǥe� б�S�
 8                                                                                                                                                                        +T��+�П                                                                                            8N7�86Q7Ƙ�7+��6�:}5��z4���3�׏                                                                    4�v4�U�4H�63��3$k2?Ok1[�E0jL�                                                                    7J�72�x66*q�5���4�A�3Ԯ,2���                                                                    3�?�3��3D��2���2��1Hɒ0p�/�t                                                                    8E�H8.�q7��7)7�6�H�5� }4�eo3��H                                                                    4� I4�ҋ4@�h3�j�3�V2RC�1�h0���                                                                    7�_�7�7zm�6�b�6D��5|�4��m3��7                                                                    5�m�5�s�5�I5�K4`�3�
2��}1�@?                                                                    6���6�ھ6J^5�_S5�34K��3Y��2[��                                                                    4�[54�g�4gF�3�ڨ352h�z1x��0z�                                                                     7��7ٔs7�
U7�6p2�5�.4�~�3�C                                                                    6	&�5���5��<52]\4�A�3�~2���1��M                                                                    4���4��4���4b�x4-�X3�V�3��3!|�                                                                    3m'�3b
3Tι37!/3d{2�S92W �2~�                                                                    4�Ym4��4��4�}�4TW�4&�3���3E_e                                                                                                                                                                        $��      cr     �`68�k2�2                        6�x53�W4@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��     @�܀    16:46:03        F�" @�܀    @��@    $�=      c�     �05dlz7^i*B'�fB'�f7�� D"g�    )/oW+��I7�t?�  ?�I    3x'�	,0&\��v        /k]�+b�>zx�>l�!>"	=��=�s�<��.7XB$26�`:v{< �2<�=F��=��W=��=���=�|=��=�ګ>�v>/@��?@��@I�C@i?���?�?���?ܠ�@uz�@uG�@j6�@i��@i/	@h�B@h@g~�@f�r@f��@f)g@e��?�D?=�A�J                    E��70�G���Gs�Fh@�E��9Ej�D?�0CO+BNo�                                                A��C$k�B��$A�m�AR�@�U�?���>�ܞ=�;�                                                                    E@�>@ĭ�F>RA��32ȰF[��B+�VA?�~.(PD                                                    {@��@��@��@my�83���!�Z�'�k�/E����?o)/oW?�ܬ�@c1�(    >�y�@�]�2�@����R-��4        B�F�B�F�C�NC�NC�{�>��[?ώ�3'�ZCްB6�D�6�G�A'��BE��>�y�Bb�A��B���A��AW�Bv��An��BT`U�yV:@=�p#(�V�yV:    �yV:��dK?-!=Z`�=[�@ꂢ@Lk3AV�A?�)?�w�@��C?��JF�$2Ȱ/�j2-q��2ĭ7�a95E��5��F��<G��G��>G��=>���    >���>���>*��>�צ>�"�>�)g>��F>��>���>�uw8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    @K��@K��{@��7���9
�9��8��8
L7f�z6���5���4�>s                                                A��?̕�?̕�{@��?سuA��|BV�7@ce>&�{@��{@��    C�aBW��6�؂    6n�|78Up@�/?3Q4�7��    >iN�C�]�7�ψBD�HD��KCǁ`BnG�B �~As�@�mG@u?�?                                                ?T��A�ϵ@��?���?'+>��">�e=DeV<Y��                                                                    D{XF��iF4bD�O6D5DCa<�B�֧A�/�@�t�                                                @��lCi�B��VA@O�@�JZ@4`?��z>���=þ                                                                    C�!E�$&Eg 9D	��C_��B�BLjA(�{@+�                                                @V�=Bu�A�c@�u�@4�S?���?�>;��=@�                                                                    7�7�K�A15>;��A�Mq?k�?�_+
#n)���,���+�#�+
#n)ื,q�,��8,r�)���7~�Y3��.1PR+-��.®6.��+n�r4	�1�%1Kq�1JE{-�%q7r1�(�r14<
    /�T�3��3�_0J�    ;��$6d9�6��4/�2Ҝc0�	�                        1��4�v4
n2��+�O�45.�2�׫            3�m<�|,��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B���G�>,��A�                                                B�                                          B��                    A�                @�p�    /�+b�3��A6�2�2��>�G�7y�f?S�M6��;�C;��2                6I�    4I)�5>r�    5U�{@��/K�k4��    ����{@��{@��{@��6�g47DE�            6&I(    6&I({@��    6&I(    5-Oo3@b{@��6�{@��6A�5��V5�    5\��{@��4�o5��/7�67�6<�f�*�
eB���6�u/F@�gC��N@%�?0�fB;4v            =�zb3�!�8��3�!�<Z��<WI?�>���>c�}=�Ƀ=.�<�te����D���=�S�����(��sD��Rd�9X�%� ����U=�BB<Z�C<WBg? �*>���>Gc(=���=O$<m����D���=�S�����(��sD��Rd�9X�%� ����U<I�    (�͉5R�.5ԑ�<ۛF<��E<�B�<���;�_�;g��                                                �)���y�>Ł/q�~*bŀT�Ń��Ōp3Ŗm̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 2�k	2�k	=ĊJ516�516�A'c�@�^,A
o�Ar�=�+�4k4kB��g6?w�;��>�/�:!�=q��6?w�F�a�295�5q[CP�MD�?tD��DS'�C��4C+B`��A��w@�=B                                                AбMCR?tC@�B�'�B=�4A�+@��@�w?!=B                                                                    E��DG	��G��F�|�F?�E��D�*hC��]C�H                                                C�iEHU�E<`�Eq�D�tnC�ۧC�LB.��A<��                                                                    E�<F�a%F���Fz�F\�F.�E�7�E���EH�                                                D)<�D��D�؅D��D�mD~AwD0��Cߝ�C�z                                                                    45.�7GR�dGI{�G��F�!xF}|�FA�E�H�EPug                                                ?��&=��#                                                                        @
��A/&A���Aϳ�BL�B7A�BUc`B}�@<�<�<�<�<�<�<�<�<�<�<�<�E�=1E�-�Ef�rE�.D��DX 5C���C��0                                                {@��{@��{@��{@��{@��D){A��C.�t��f�I7�+Cp��@F��    A�*A�*{@�ξ�������C�5h{@��C�\YC�\lC�5h{@��@m�6{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�9�;�D���GHD��xD��oD��oA6S4A6S4F^�pB7��?���C�-�FK��FK��D�;QD�;QF_mB7��                ?ǩ�C�LGC�o@C�*�?'�C�	-C�d~C�d~C��LC� C�:jC�sC���C��C�^C�F�C�`�C�t+C�qC��HC�x�C�bC�>zC�C���C��3C�pmC�4�C�
�C��C��C��C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�C�88>��G>�CV>�� >�>�8">ħ�>�!�>ˑ�>�>>ϰ�>���>ѣx>ѻy>�"�>ϩ�>�M#>�<�>��>�b�>�T2@6w�@�A        Ar�������}��Y�ŕG\�8��7�5�7g��ɨ���ő��Ş��Őń�B        >�A޳I    @�?@�?{@��@��FP7I{�7�        ;}�Q?��?�+?O9?b?��?��?�W?��                                                �(��6���A�  ?k�A	��@�&$�=      c�     �08߮�8O�7��s7!�o6iK�5�<94��3��                                                                    8E�7�96�a�6L��5�XQ4���3���2��                                                                    G���Gs�Fh@�E��9Ej�D?�0CO+BNo�                                                4jv3��$3B{2���1�x�1!�s0+/'�}                                                                    3��3	��2u��1�n1g(0L!h/X��.S��                                                                    6�Q�6��X6�do6��5cP�4��^3�2��{                                                                    7�6ޤR6X��5��5`|�4�!_3��3 ��                                                                    4��A4�Y�4��3�~r3-�12�Z1ķ)0ոf                                                                    70��76��6�5�/�4ٷ4�Z3i�                                                                    8A��8J��7�kP7�t6tU�5�ٛ4�k�3ύV                                                                    5	�52{�4�TG4K�23�E�3��2K��1K�R                                                                    8A��8J��7�kP7�t6tU�5�ٛ4�k�3ύV                                                                    7�\{7�Ŕ7�6r�5�͐5/#43b31�&                                                                    4��4�'^4B��3ú�3=�$2�>�1ǘ>0��u                                                                    7�\{7�Ŕ7�6r�5�͐5/#43b31�&                                                                    �1*Գ��x5q��5@S4�d�4�d34�m2�                                                                    4w�2H�:9�}�9qi9 ��8k;�7�+�7'J�                                                                    5Q�h5wZ4���4E�3���2��2
�1��                                                                    6�Oc6�R56&��5|�\4Ä	4	�3�K2��                                                                    ��U���w��\�S��_�?������OǱ��                                                                    ��:A��������y�E̼��Q��R��YC                                                                    5�N5�xt5&�4~)�3�J�3
v�2 ��1�+                                                                    ��$����>�R����B��!.��jN峓�Q����                                                                    �=6j�>������E[��7
����'!5�;�l                                                                    �ʹ���tt�g[���p��0U�ϰ�U                                                                                                                                                                        *��:*�(�                                                                                            8�8J�7��D7�;6p*F5�͑4�S�3��=                                                                    4���4���4,�>3���3 MJ29�1d�0��                                                                    7�F7+B6�y6��5x��4�C�3��2��`                                                                    3�t�3��b3)S�2�sp2�C1BJ�0z�/�V�                                                                    8��8�7���7�6�n�5���4�.G4��                                                                    4�Sh4���4%��3�T�3	9;2KvI1�0��n                                                                    7���7��b7X�6�I65�D5tO�4�Ӡ3�7                                                                    5�c'5��5v�5#�4O~3��02��I1�c                                                                    6i�6q��6.��5��5�24ElD3a��2o�                                                                    4��u4�=�4G�#3�Sa3'��2a�M1�"0��6                                                                    7��_7��7�7--6]��5�M4��53�&�                                                                    5�#�5��5��54}�3��:2�P<1��                                                                    4X8�4k3�4c,4J�4 �w3�W3�6�30.�                                                                    3.��3>�37��3#Q�3�k2��+2__�2^�                                                                    4�"�4���4���4wG4D-�4��3��X3WU�                                                                                                                                                                        $�=      c�     �06bK�/�T�F�P�    =%�T            6��e3��F@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�܀    @��@    16:46:05        F�^ @��@    @���    $ա      c�     ��4T��7*�LB'�fB'�f79Z�D"g�    $N`G+�m6��?�  ?�n49d�1�5//L(�1/�0�ͳ        *�C&�3.9���7��<��-=@=:�=�)<�)40k :��< �5<�O=F�f=�х=��i=���=�|=�6=��">��>�@@o�@\x@A0@�?��<?��?�b�?��B@{h@{%�@p�@o�:@o2@n��@m�#@l��@k��@k+�@j��@j4�@`�<EC{@��                    E��+6�
G�:rG��Fh�SE�@�E�&D?�CO �BN��                                                A?iC$�B�@A���AR��@�sn?��>���=؆u                                                                    E@Ӈ@Ļ�FN�A���/GRF[��B+n�@���-��                                                    {@��������@i�b7��!{��"�l&��/)���L�?�$N`G���j-��3�_    >�y�>�Y1�c?|�>���D��<N.Y    BT�BT�C���C���Cp��?8_�>��3E\C�ʥ6�v>6�	�A'b�A�}�>�y�A��g@�}B@��@��xA�_O@�VB����4$@��#g�r��4$�����4$����@<�:?G�?G�@�C?�o�@�]R?\Ms?��@��(?��JF�C/GR-;+|-�QO2�I6�K�,.{)4��nF���F�ЧF�~CG}J?�L@�2A��a@|*,?��?z}>��J>�:�>�/�>�J>���>�ټ8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��s    ������{@��6�m7��]7�!67G7	+�7m6z�X5���4�o�                                                A��!���/���/{@��B��XA"#�B7	`@e�C=�J�{@��{@��    C�s>BY��6P��    6��6��{@�"q3��A6�`�    >��Cf�,7�S3BO��D���C��B��nB�A�>@��(@wZ?b                                                ?bReA���@�'�?�J�?>�n>�(�>�R=R��<k|g                                                                    D|$�F��sF��D�:�D/Cc�)B���AǼF@�l�                                                @�C7�B�5�AFw�@�[}@7�?�A�>�X=�}�                                                                    C�/$E�DdEhV�D��CcHB���B
0�A+#�@.�                                                @Yu"Bw��A��@��t@9�?�=�?��>?�T=EĀ                                                                    6Ժ6��_An7><ÜA�^?J�?��)�s�(a?�+�!�+6��)�s�(�*�*,a�+�`k+���(�+h79Z�*���(?dP%�@&            �4�;3E\\1�k (C�Q(?dP%�@&74�;��jǷ4�:3�{�    !?b%*��@*�
�'���    =Y�;��iZ-�F#)���(G@�'�                        $x�#)���)pC�(%�    )���3f            3�;E�,�-�?C��?S:�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B5�LG���>+�;A��o                                                B�                                          B��                    A�                @�p�    *, &�\*���4�[�/Y��2�W�=���6�\�>5��5xq;1�6;6��                ���
    4�8�0��    7͵�{@��3�ɗ4���6�KF��8�{@��{@��{@��6�|�8            7lV�O J7lV�{@��    7lV�    ����6\�{@��7Zu{@��7�%�6ȕm7R �    4�{@�ζ
 �3wu7�Q�7�Q�<���    B�l�6��F\VuD�?H��?;�A�5�            >��*���9h�*���=���<���?ǈ�?�<?O��>�	7>=r=[ި���E[��=��㊰젊�)�sE[�Re�9Yj�%��������>���=��u<���?�ߡ?��,?��>�A�=�K=�����E[��=��㊰젊�)�sE[�Re�9Yj�%��������=�3�    *���4�605�A>5H�>Zx�>HP9=��=.�4<�aW                                                �^��Ù͇Ŋ� œ�Ř��ž{�ţ� Œ?%̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 5�o 4�w*?n�8)I�7>�@��@aܤ@s��@�I�>%D+7��6�E�h�7��1=�$
A�Q�?	��=�t�7��1I�M�7���80?�CP��D�CHD���DS�C�}�C�B`��A��@���                                                AТ�CRCHC@��B��B=}�A��@���@�?!��                                                                    E��2G	��Gy�F��IF?�E��aD�,�C���Cӕ                                                C�g�EHB{E<S�EuLD�}&C��uC�)B.��A<֪                                                                    E��F�^�F��)FzhF\�F.�zE�8yE���EH�                                                D)<[D���D��D��D�m D~B&D0�ACߞkC�z�                                                                    )���6|GR�GIr	G��F�'F}eFBuE�IFEPv�                                                AZ#Ato�AB@���@t��                                                            @S�@�EA"�
A��jA��B-ƜBP�B~�A<�<�<�<�<�<�<�<�<�<�<�<�E�2�E�&=Eg�E�D��DX �C���C��I                                                {@��{@��{@��{@��{@��C�3B4�3)ݱl�gK�7
��CB9�@,+�    @>�I@>�I{@�ξ8��8�C�Q�{@��C���C���C�Q�{@��@k� {@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��=4K
D��VG0�D��-D��D��A9[pA9[pF]��B66nB6|�Ca��FK��FK��D�:�D�:�F^ʾB61V                ?D{C�VAC��vC�?�?U C��C��UC��UC��xC�,tC�{�C��BC�){C��C��C�a"C��gC��#C�1�C�kaC��fC��5C��C��uC�̩C��AC���C�S�C��C��C��C��C��o{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�[C��>���>��e>���>�ic>�CW>��D>�c�>�u�>���>���>�ʾ>���>�
>�%�>ŉ�>��G>��>�*>�1�>��@ܳ?�L        @����z��L���ū���������9Ų5��� t�� t�� tũ��        =O"�A�#�    ��������{@��@oFP'�7�7�        8��%>��>-�K=���>.��>٫�?�
?�?F�                                                �WÏ5�
UA�  ?k�@�UU<#�
$ա      c�     ��7�-�6�ϔ6.�6�S6��5dR4�i�3�                                                                    6�b!6	�5[��5I�\52ۘ4�3�3���2���                                                                    G�:rG��Fh�SE�@�E�&D?�CO �BN��                                                3@p2t�1�ly1�7�1�`�0�:B0��/��                                                                    2r��1� k0�n0�900�l�0R/,��.4�                                                                    5�5�#:5i�5c5	�54r{�3���2�I�                                                                    6�~5��R4Ӽ�4�GU5 �4��^3���2�ͽ                                                                    3��3(L�2���2��;2�Tl2e��1��-0�                                                                    6!+)5���5e<5�56d4��3��3��                                                                    7*C6���64�_6k�6��5��B4��83�Ȕ                                                                    4�3Ԓ3B�13Kz�3q��2韮2%�12I                                                                    7*C6���64�_6k�6��5��B4��83�Ȕ                                                                    6|_$63�z5��5pB�5v�g4���4��3��                                                                    3z�3D�E2��o2Ë,2��2d�1���0��9                                                                    6|_$63�z5��5pB�5v�g4���4��3��                                                                    �"����,J3޷�48˶4M�3���352T5W                                                                    8507���:b��:H�9��9M�Q8��7���                                                                    4>�)3�@<3¢3n�36WV2���1���1
˵                                                                    5�5U/`4�sY4y�U4oF3�ę2�$2GB                                                                    ���=�N)V�����
,L��c�^�Ʋ�2���<M                                                                    �Y��"=�K�u6���ʲ&�����8��                                                                    4�I4TRZ3�h�3{��3q�M2�62s1	�B                                                                    ���Q�S���aԴ�#}��a��-B��k���*w                                                                    ����uE�G���A�A�P�����ð!%��                                                                    ���5���P�M�|�cκ�_�1������3����N                                                                                                                                                                        ,��,7E�                                                                                            6��6�^�6 7N6L6��5|�'4�^�3���                                                                    3z�3-C�2���2�2���2�x15��0Z1�                                                                    5�`5���5�595�~4���3��N2�6                                                                    2u�}2)�1���1��31���1�v0G/�/o?�                                                                    6�6�b�6��6&6�5�І4�Q�3��;                                                                    3p��3&4:2�b�2��S2�W2st1X�50�&�                                                                    6l*�6'K�5ʢ5�?5�A054�%4^��3{��                                                                    4��4?1�3甬4 $I3��73N�*2~z1��,                                                                    5>�d504��t4�5�4���433��2K�m                                                                    3Z�3�2�"�2�2�h�2&�r1M�O0h��                                                                    6�R�6Lx�5���6	
\6��5\�f4�3���                                                                    4��4i��4��4� 4�"3|wP2���1��6                                                                    30��3x#2�3F��3²�3�{�3\R3�u                                                                    2��1�y�1�02 ˞2�U#2��T21��1��7                                                                    3W��3 �G37�3s3�3��+3��3��37 V                                                                                                                                                                        $ա      c�     ��6�Q!?b%F���>�9=�            6ێ�4">@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��@    @���    16:46:08        Fǜ @���    @��    $�e      c�     ��4�T6��B'�fB'�f7)�D"g�    �a�+x�6���?�  ?��4���1�S1��60�i?09�k01U]/��y%��! �9�2�7�6�E}<�<=�Z3=�W<�A�4��7:#��<!�m<�˟=F�r=�ӗ=���=��v=�|�=뗼=�۲>�>�9'�8ש�>�:>!ځ=�??/�z?���?���@vg@@y$�@n�q@o�o@o�N@o}m@n��@n:\@mbW@l��@k�|@k�{@��{@��{@��                    E��R5��IG���G3%Fi��E���E�TD@*'COi BO$�                                                A�ZC%o�B���A�� AS�T@��Y?�X/>�I^=��                                                                    E@�w@Ĵ FL0A���3 &�F[�B)�<�D/-2�                                                    {@������@^�N7c�եR���'�X~/Cb�����>�H�a��3�.b&4�+V    >�y�=�91�#?X4� fP��N�<	�    B�PB�PC�Y�C�Y�Cw=?Xפ>JB_3�sCڟ�6g+6��6A$=�@ۉ)>�y�ASV?��;AL:�@+1�?��A2#@V0cA�v���n?���!(2���n    ���n��I�>�ŷ?}{�?}{�@;>�@A�?�4>�6@,$�?d�BF�3�3 &�1U�-�&`2�! 65    4X��FJ^>Ft�F�&�F�@N?�v{@j��B���@��?sԎ?5�8>�ص>���>�R�>���>�T�>�B�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�7�    ����{@��6jS�6�6}6��6�{�7C�7JU6|7�5���4�R                                                A���0��0�{@��CH�@�<B,�@]
L<�s{@��{@��    C���BZ�6
�    6ka
7m�z@���3?�6�.�    >��(CjLf7�D�Bi�D��yC�G�B��B> A�EI@�E@)K=?>�                                                ?��A��wA Q�?�Z�?{|�>��>7 =z�<�+&                                                                    D`�F��aFӭD��D�Ck!TB���Aϯ\@�[S                                                @��C	a�B��?AY.�@� 8@B�?�� >��=���                                                                    C�k�E�#rEk��D��Cn�/B���B�A3�@8�i                                                @a��B|B�E@�8@I�=?���?8E>KF=TS                                                                    5��6_��A|>:�A��?��?��)���(+A�+�H�+�`)���(dW�)HZt+���+x$(��7)�                        �;�u3�s1�a            7;�u��(�;�u3v��                        >9vطI}                                                                2���            2���{@��,���>L��>L��>���?	1?j��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A��G���>*�A���                                                B�                                          B��                    A�                @�p�    $�>9 �    4?)�3(w�2�"'=<�6��=�:�4:�M;~�;;�                �u�]    3�la0Q��    57F�{@��3���4��25�1e���c{@��{@��{@��1�pH7���            5�V�    5�V�{@��    5�V�    3���6��
{@��5�${@��7~�6�@�7-��    4N� {@�ε���2�M�6W�6W�<���    B�H*6}0�F~!�D�Ts>���=���@�	0            ?��    9��|    <��;c�?�Q/?�P�?�|�?9�>a3�=�R����E[��=��㊰젊�)�sE[�Re�9Yj�%��������>�,<��;c�?�x
?��g?u��>��{> �X==�����E[��=��㊰젊�)�sE[�Re�9Yj�%��������=�'q    %���1�e�2�9>�%=���=�.3>m�=��=q                                                �!��©����Ŕ�g��N�����ŵ]UŬ�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6l��4�&N?<{b8��U6,�?��?�'?���@E�=���7��5��F�܋7��!>���B�'?:W�>�FT7��!J���5��o7�'BCP�#D�a{D��&DS]NC��]CUBa>�A�P�@��Z                                                A��#CRa{C@�&B�]NB=�]A�U@�>�@P�?"�Z                                                                    E���G	�oGc�F��)F?�E���D�4�C�C��                                                C�c�EH �E<48E{�D��pC���C��B.�lA<ߵ                                                                    E��F�Z?F��FzeF\�}F.� E�9�E���EH
:                                                D);�D���D��*D���D�l�D~CD0�Cߟ�C�{�                                                                        5�	GR֚GIY�G��F�5�F}��FD�E�KEPy.                                                A��A��3AB��@V�d>;�R                                                            ?3�-@�v�A1A�mJB��B-7JBM4oBw�{<�<�<�<�<�<�<�<�<�<�<�<�E�!EE��Eg�E�aD��DX$�C���C��"                                                {@��{@��{@��{@��{@��A�7�Bm�O$XTA�e�6��`CF�@UX    ?�<*?�<*{@�ξb�"�b�"C��O{@��C��C�q~C��O{@��@^��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�!>���D��UG&,D��,D�X�D�X�AA��AA��F\��B3pBXgCe]dFK��FK��D�:KD�:KF]� B3j�                >�PC�u.C�)C��7?$c�C���C�6�C�6�C�E�C�`�C���C��HC��YC�+
C�o�C��lC���C��C�N�C��8C���C�EC�:�C�e�C��C���C�~ C�`C�%C��C��C��C�Mv{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C�O�>� �>��k>��>��>�]�>���>���>�L>���>�\~>�R3>�|Z>��#>�:g>���>��>���>�A�>��V>�*8@-f@�f        @+�)�ߗ#�ߖ�ߖ#ţ���������ťW��F��F��Fţ��        =L�IA��,    ��u���u�{@��@�-<FP"�6�p7��        7>�=�zB=ǌ�=�L>�#�>�`:?��?��?�]                                                n516�A�  ?k�    <#�
$�e      c�     ��6!�,5�L�5��"6Z�69Q5^�s4t�>3w�0                                                                    5L4�H4ۤa5���5?4���3���2�[�                                                                    G���G3%Fi��E���E�TD@*'COi BO$�                                                1�S�1Fl�16?�1�J1�w�0�h�0 -�/��                                                                    0��]0z�0f5�1]0�+]0j�/!��.#ǟ                                                                    4+G�4[,4��057P 5644l��3��2�τ                                                                    4�l4�n4�	5gp�5>��4�[l3���3��                                                                    2A0�2-��2]`I3* 
3��2�|1���0��j                                                                    4�u�4���4�ԙ5�o�5iC4̌/4D%3�%                                                                    5��5�p�5���6V $6$�H5���4��B3�>�                                                                    2��2��2�V�3��3��2��M2#с1-3s                                                                    5��5�p�5���6V $6$�H5���4��B3�>�                                                                    5M5*��5S5��]5��?4�\Y4��3iF                                                                    2z�2=6�2MD�3��3^(2kr�1���0�6
                                                                    5M5*��5S5��]5��?4�\Y4��3iF                                                                    ��E�fd3i4Zr4Jk�3��2�Rj24�`                                                                    6���5��\:q|�:hCs:"e9�]19 �8C�                                                                    2�nr2�	�2�e�3�`�3i�2�:�2�1�+                                                                    4@XM4I�4#�4�Z}4��3�x�2�3�1��D                                                                    � ,�'%U�x!��<E<����Y²�������                                                                    ��h��e���kI�������R�ł�~	�'w�                                                                    3?�A3IP�3$ml3��J3��2Оo1���1��                                                                    ����+6�Di{����ئ�)Mz�]����}                                                                    �������ƕ���B��_$*��~���^���                                                                    �S���ɤ�̌���3��n�,�������U(                                                                                                                                                                        '%�'�[h                                                                                            5Zj�5���5�m�6M�6��5v�
4���3��!                                                                    1���2�+2!12�2�2�>2� 1*�~0G3                                                                    4U�44�$�4�.5K�Z5! �4���3�Y�2��                                                                    0�T�1	�61�<1вE1��1cW0;M/Zk�                                                                    5Q��5�O�5��w6J,�6&iu5���4��3��                                                                    1��2�B2�i2�1�2��j2�1Kȅ0m�I                                                                    4��?5�~5IZ�6��5��50#�4P��3d��                                                                    2�s�3�3fK4.�#4?)3IM�2nsV1���                                                                    3��63�5(4"��4��}4��4U�3(��28�`                                                                    1��b1���29�G3�2���2"�I1@��0S6%                                                                    4��5%ƙ5vI6:�-6��5WHI4u3��:                                                                    33=uB3���4Ug4$N3v	x2��C1���                                                                    1�S1�@Y2S��3�{@3ϟ3�E�3N4�3�                                                                    0{�0�S1+
�2Z��2��_2��n2&��1ۡ�                                                                    1��N2Q�2�Y�3���3��p3��3|�3&�                                                                                                                                                                        $�e      c�     ��6^��    C���    ;��    >��{    6�O73��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @���    @��    16:46:10        