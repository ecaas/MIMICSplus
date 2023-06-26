CDF      
      time       bnds      lndgrid       levsoi        levdcmp       cft       glc_nec    
   ltype      	   natpft        levlak     
   nvegwcs       string_length         levgrnd       hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      Sun Jan  9 16:23:30 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1948.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1948.nc
created on 12/10/21 17:00:27    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:30 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1948.nc had following "history" attribute:
created on 12/10/21 17:00:27
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�fG� @�w@    @�{     )>�      ��     8�4���6�B'�fB'�f7'�D"g�    'K�r+��p6��T?lS+?lKl    1�~�$z��1���͔J        -�@�)u�[:f�W9.*>S��>0�>R=���=��6,��:/��<%z<�X�=G��=�<=�	�=��=��-=��=�t>�>v�@�P@���@��:@M��@-�@@�@�@E�@��@@��-@x%�@w»@wS�@v��@vG�@u�;@u'�@t��@tQO@t�{@��{@��{@��                    E�b�6%�G��IF��HFZ�HE�nEfD8��CM
6BP3?                                                A�vC�B�OA�iqAH	E@��%?�K>�ҩ=�=�                                                                    EH.�@͆�F%�GA�T�%��Fd*�B0B�A<��,�M`                                                    {@�ο�������@c�7�$��̤;��"/G�q�#1v>�]d'K�r�U�?���1��E    >�y�=Bo�1���T�Y�����e        A�ǹA�ǹC�8qC�8qC{��?-�s>G`�3T4C��6q�16���A'��Aަ>�y�AQGO@=�UA�x�@d�J?�A1�@��NA�*�� >��,    ��     �� �9R<���?��?��@a��?5�;@y�\?9vw?5)�@|�?� �F�kc%��#�6�-�g�3d�6�    4���F6�FA�F���F�ݍ@N>    CrI@n�W>I�!>T3>b��>t)A>k�X>cZP>U�l>TB�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�-    �Mg�Mg{@��6��F7�<�8%C�7��7�76熫64r�5m�f4��                                                A����V���V��{@��C$��@��{B-z�@`��=2��{@��{@��    C��ABX�6!�    6KKu7Y|�@�kl3>�6�cH    >���CBd�7��DBM`rD�{*C�ZB�$pB^��A���A8,@z��?�&L                                                ?idA���@���?��^?��}?�$>�4�=�T�<�3                                                                    D:_YF~�E��D���D��C�&�B��+A�A��                                                @�0�B�;�B2*AY��@��b@q��?��?$�>FK                                                                    C��?E�{�E
�JD��C}�B�8�B.�A`m;@o�T                                                @=/�BL�	A��@�ܣ@f��?��?Gů>�S�=�p�                                                                    6`D�6�m�A��>8�JA�?Ǚ?�r�)s�"(PG+���+��)s�"(C�)X'C+xN�+T2�(���7'�                        �i��3T41��z            7i��1��E�i��3_�C                        >w�����                                                                2��            2��{@��,���?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A�o�G�S">)�lA�                                                B�                                          B��                    A�                @�p�    --�X)�n    4���%�13j�=3��6��=��$4q�;��;UE                ��c        (�-        {@��4?��3�'�    �:�{@��{@��{@�έ���8b+            �Pm    �Pm{@��    �Pm    ��&Y0���{@��3I�:{@��8��7��7���    ��3{@�δ�;z2�
D0�x50�x5<��?    B�$�6�'�Fh��D�I>��e=��AO            ?��    9���    ;	� ;|�@I@
�4?v�n>�fL>=K)T�&+I�e��>�����|芬RЊ�e��{Ȋ]��E��2�#@�>���;	�;|�T@
�<?�+�?'p�>ts`=��=
��&+I�e��>�����|芬RЊ�e��{Ȋ]��E��2�#@�>:iE    S>2�4�C?zWN?0��>��=�o=.�<�=T                                                �q�U�O�7�	�����ƽ�ƶ	ƴըƼ9R̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6�H�4cՃ@�j8�ʭ5�m�@��@؄?���@���=מ�7��5��G%�8L�L?��Cr<9S�a?�8L�LJ�w�3�'R8M��CT��D�#vD�sDH-�C��JC �Bq��A��@�~x                                                Aԁ�C^#vCHsB�-�B>�JA��@��@�?6~x                                                                    E���G>lGz9F��SF8��E��D�?�C�Q�C�*                                                C��TEPZ�EB&$E&D�4?C�f�C��B-R�A<zk                                                                    E�;DF���F���F{i�F]5�F.�xE�J�E���EG�X                                                D)�wDé�D�_�D��^D��!D~��D0�Cߋ	C�I�                                                                        6GI�GYq�GN6�GuBF��Fz��F�E�0�EP9�                                                @<�@��                                                                        ?=`B@�	AT��A���A��BrB&�$BF�u<�<�<�<�<�<�<�<�<�<�<�<�E��E���EcӭE�D���DW<�C�mHC���                                                {@��{@��{@��{@��{@��D++A���,���f�y7+��C W�@+    ��FԾ�F�{@�ξBѾB�C��]{@��C���C�C�C��]{@��@eq�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G?�:~3�D�Y�G�D�Y�D���D���A!��A!��Fd��B8Z�@�~�CA�FK�%FK�%D��D��Ff,B8U~                ?�)C���C�JqC���?
�C���C�GC�GC�u�C���C��zC�	C�<�C�v�C��#C��C�ZC�@{C�kDC���C���C��fC��C�6gC�N9C�ZjC�Z�C�L�C�"
C�jC�$�C� �C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C���>�ٱ>�z>��>���>�D�>��>���>�M>�/>�^�>��z>�~>�t>��>�u>���>��L>�B>��&>�!@�p?�`J        @crJ�:@��:@��:@\ƺ҆�D�{�D�l�D�Gƾx�������ơ�v        =L��A�n!    � �� �{@��@f��FX�/6���7v�        8H=��g>B�>��>�4>�{>��i>�X�>�_                                                ����5P�A�  ?k�A	��<#�
)>�      ��     8�7kl�7Y�7 �J6�m}5��25�4�I3��                                                                    6��R6���6"�]5��4��u4+��3?^�2>x                                                                    G��IF��HFZ�HE�nEfD8��CM
6BP3?                                                2��?2�`2��y2
[~1OY�0�/�w�.�ȴ                                                                    2��2�_1��;1.�i0��Z/���.�o-ȑ�                                                                    5�|�6
��5�$5h��4�]*46�3P%#2yh                                                                    5��N6+��6"p5���554���3�kN2��                                                                    3|#V3��3�KP3���3Ԫ2x�^1���0�)8                                                                    6	q�6Q�6%*5��5]]�4�E43Z3�e                                                                    6��}7j�x7Ѓ6��5�m�5D�44z'}3��                                                                    3��4f�W4,}C3���3jE�2��2
�1.�                                                                    6��}7j�x7Ѓ6��5�m�5D�44z'}3��                                                                    68s-6��o6X÷5�
�5Xwp4�13޸�2�k                                                                    3D�3�qM3�E�3V�d2�P2E��1�G�0�e�                                                                    68s-6��o6X÷5�
�5Xwp4�13޸�2�k                                                                    ��I����4��\4u��3���3'TT2|]z1�i�                                                                    c,y��:�c:���:Zy9N5H8�M{7��\                                                                    4!�x4m��4293�R�3]gz2��1�N1P                                                                    5Z`H5��5d��4�4A�F3���2��1�s                                                                    �K<m�Ҥ�����g���z������*n�                                                                    �9H���/�ϔ���Z۱�α{N�/�ʲ^                                                                    4[H�4�g74h�\3�k3H2��91���0£�                                                                    �X����r"��Fs�&b������c����/8�                                                                    ��	ڲ[ٷ��b��1��Ȓ�zX�����z�                                                                    ����-����^����ص�\�@���%��;                                                                                                                                                                        W�Ie�o                                                                                            6�ٗ70�x6��6���5��5(6�4Pp�3o�                                                                    31�j3�t$3q�_3>+2n�1�Pb0��0jj                                                                    5�DN6,��5�7^5�)4�p40��3d��2��                                                                    2.+�2��F2l�m2F�1v[�0�Xs0r"/{l                                                                    6��7)O6�X76��5�R�58�4x�Y3���                                                                    3*��3�.g3g�{3O-2~��1�`�1�I0!�m                                                                    6$��6��6�l�6<55�aS4�b�4��3
M>                                                                    4<{04�Y@4���4V�3���2�q2J*1"                                                                    5E 5� 5{1!5��4{�3��m2ѳ�1߄�                                                                    3N�3�ۄ3���3-��2�11���0��/�r�                                                                    6I��6��6���6e�5��5��4�m3)	                                                                    4f]�4��j4�4�^3�	�3�L25>31A.�                                                                    2���3�Tc3�ri3��3��y3K��32�{�                                                                    1��2T?�2�e2���2h'�2$�e1Ф�1��y                                                                    3[3���3��+3ԙ3��[3ys3�Z2�	
                                                                                                                                                                        )>�      ��     8�6>�[        <�g�    G�O    =D$6�3�i�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�w@    @�{     17:00:27        G� @�{     @�~�    )>�      ��     =�4�i�6��cB'�fB'�f7-�[D"g�    (ֲ�+�Ԇ6Ab7?l]?lL�    3�%�ׄ1����ށ1`_G    /+ />l��>�2>a�C>)��=�$�=�Am=��5��x:7��<%�<�s�=H=�>�=�{=��s=��H=��=�t�>�>w(@�g @�E-@���@K%�@,�R@�4@�@�l@���@��_@w�c@w{�@w:V@v�n@v�h@v� @vk?@vI�@v0�@v"�{@��{@��{@��                    E��q6�p�G�BmF��hF[��E���EΖD9!CM��BP��                                                AhPC86B�jpA�4AH�$@�?���>�yN=�]                                                                    EHs@�]	F%��A�2
)?�Fc��B/�YAO�-6�c                                                    {@����z���z�@_�7�l%Ӝn%�'��'k/�$U�p>�]d(ֲ��u���+�1��f    >�y�>�a�1���r��`��Vy        BJ��BJ��C�O�C�O�Ck��?�>�!�3,C�;b6J�6���A&k�B��>�y�B?��Ax��B�"A��A2��BF��AL��B@T������,�#2u���    �����Q�?\�m?�  ?�  A=S@W\4Abq�?���@S1'Ae��@\F�e�)?�'օh-�l3LQ6��    5r�F�iqF�� FրGwS?�58    C��K@��m=�n.>>*�>u�>y�>m�~>c�0>U)n>T��8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M� �    ������{@��7*�8���8��	8�(7�z6��B67�=5n_4�t�                                                A��������{@��C��*AdD�B<X�@YL=�m{@��{@��    C�\�BY5���    6H�7X�@��3
W(6�'�    >�\�CH4�7��7BS�D�ЌC��IB�;zBdΕA�bA;a@~�L?�iM                                                ?pl�A��@��@ ;r?�dK?>��_=���<�:�                                                                    D;d�F�E���D��2D`�C���B�JxA�6}Aٌ                                                @ʝ�B�6B4�A^M�@���@uR�?�J�?j>Uj                                                                    C��E�g�EFD�C�1B��B0Y�Ac��@s��                                                @?�@BM��A���@�$1@j�O?�pQ?K�>��=�                                                                    6���7"NA[s>2�.A�|�?�7?���)���([z+���+(OC)���(��*v�+��+�j�(��7-�[                        �k��3,1͟�            7k��1��f�k�f3s��                        >�>J��k1                                                                2���            2���{@��,��4?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B�uYG��>,�QA�bQ                                                B�                                          B��                    A�                @�p�    .��h*��f    5jT	)B�43Ll�=Ͱ�6��>DD�5a�`;М;%�                �ے'    3���/��$        {@��4r	'3��s    ���{@��{@��{@��$��7�H|            $)6|    $)6|{@��    $)6|    �Zn7��{@��6-K�{@��7��,7h�7-sl    �Zn{@�ε���3��C6� �6� �<��U    B�Y)6�^F@�!D^J�?�4�>��{B
��            ?=    9���    ;�CY;ժ>@~k�@3��?��m>�s>.K�=��*�&+I�e��>�����|芬RЊ�e��{Ȋ]��E��2�#@�>�J0;�@�;դ�@�?�?CR�>�]�=��=e��&+I�e��>�����|芬RЊ�e��{Ȋ]��E��2�#@�>���    K4�Ȥ50�o?��?�Ҝ>�4J><��=���<��g                                                �I1�ƞ|�ƵǉƲ�mƳ^ƵƷ �Ƹ�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7
GT4�C�?q�&8�N�6MS�A>�@��A(�'A=V>=�86��5��XG��7�-�?�wYC��g?��!?�wY7�-�J�@H3ǎ�7�ECU �D�6�DȆxDI�C�ʎC!�`Bs/?A�i@�ކ                                                A� �C^6�CH�xB��B?ʎA��`@�/?@i?7ކ                                                                    E��NG'qGi;F��F8��E��D�XMC�m�C�                                                C�ǸEP9_EBmE0�D�F�CԅaC��B-gA<�g                                                                    E�;NF���F���F{f�F]5LF/ VE�L�E��EG��                                                D)�Dä,D�[�D�։D���D~� D0�Cߍ2C�K�                                                                        6��GYY�GN$�G�F��JFz��F��E�5EP>B                                                ?|΀>�'�                                                                        ?���@��AfH�A�8�A�fB6�B&�^BG�<�<�<�<�<�<�<�<�<�<�<�<�E��pE���Ec�E��D���DWE�C�t>C��+                                                {@��{@��{@��{@��{@��D:=FA�c$.f�-�g��7\7}BЃ�?ᾈ    ?���?���{@�νٌ��ٌ�C�|.{@��C�D�C�ˣC�|.{@��@eǵ{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�:spvD�Z�G�pD�Z�D���D���A$=-A$=-Fc�LB6��?��%CE�uFK��FK��D�0D�0Fe��B6�N                ?�N{C��?C�f�C�[�?��C�=�C���C���C���C��C�FC�3�C�OBC�o�C��kC���C��UC���C�VC�ANC�kC���C��sC���C��C�-KC�;�C�AC�(�C�IC�$6C� �C��i{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�BC���>���>��o>�7�>�L>�q}>���>��X>��=>�m�>��>���>���>���>�>���>�?�>��>�Ƃ>�HA>�#�?�w7?��        @���4 �3�*�3�Ƭ���A�	�A�!�A�ZưJ��I��I��Iƚ�X        =\2A�Z�    ��������{@��@;�FX_56we�7W�        7�NT>�?�>��O>�!�>��9>�7>�Or>�w7>�*                                                ��W�6{�A�  ?k�A	��@�7�)>�      ��     =�8a!7�V�7%l�6���5͍5	�N4�E3c�                                                                    7�#�6�r6P�"5��F5�b4-��3>{J2?:�                                                                    G�BmF��hF[��E���EΖD9!CM��BP��                                                3��3E.�2�#2��1Xcu0��o/��Z.�j�                                                                    3��2y�1��Z1?n�0���/���.��-�^�                                                                    6���6n$�6�35k4ר14]I3O��2z�d                                                                    6��6�|�63��5�?N5@[~4��93�U`2��                                                                    4_��49Kt3�23��3�c2�	1�@�0��.                                                                    6�:6���6[�5�|5k~4�')4^�3�                                                                    7ӫ�7Ź�7/�w6���6Ő5H��4{]�3�/�                                                                    4�Y�4�u�4ar3�T�3w�.2���2H�13,                                                                    7ӫ�7Ź�7/�w6���6Ő5H��4{]�3�/�                                                                    7�I7�6��>6��5dV4���3��2�5
                                                                    4*��47S�3��Q3nO�2�T2J��1�}�0�i�                                                                    7�I7�6��>6��5dV4���3��2�5
                                                                    �"����R�4�3�4���3�e3&Z�2wݑ1�H4                                                                    �mg]3H$;�j:̕�:7�9��!8�.�8#g                                                                    5r4�S�4mY�3�q3k�2��n1�i�1*x                                                                    6<��6-J�5���4��g4K�T3�h�2��51���                                                                    �Acf�4a�����}HV��Tͳ&��ױ*��                                                                    �+�
����q��q��-T�}8ɱ�t�ʿ[                                                                    5=��5.�54���4y�3Rm�2���1�}{0ŝ�                                                                    �N^$�:c��1�6H۴��J��KL��[�0��                                                                    �Ц]��rò?H�����'7��}�3���W���                                                                    ��k����
�I��˅�$Q�^�R���߯���                                                                                                                                                                        Z��E��                                                                                            7��{7�J�7�S6���5�25*��4P[�3p�B                                                                    4)\�4��3�D3r*2y(D1�0�W0w�                                                                    6��6�5g6pE5�и4�M�43@�3dv�2��                                                                    3%�A3��2�(2b�1���0�G�0e//��                                                                    7���7� 97O76�ǣ5�iG5;��4x�{3���                                                                    4"u�4�3��>3R�2�=�1�s*1�40"͞                                                                    7�
7`=6���6M�V5���4�W	4 �l3
gq                                                                    53Y05'Ik4�.�4j��3��#2���28{1-                                                                    5���5�25��U5&A4���3�>�2�0�1߮�                                                                    4�4.f3���3=�J2�a�1ʐ�0��/��[                                                                    7?͚72�h6���6{4�5Ų�5
�4q�3))                                                                    5[45Lv-5
�r4���3��30�23�1AS�                                                                    3�}�3��3�\[3��E3�t�3M{�3 �2��L                                                                    2���2�׸2�f2���2q�/2&1�/1���                                                                    4�4	��4��3�.73���3{%p3�,2��                                                                                                                                                                        )>�      ��     =�6<�                FQ�\    <��6�*]3�V�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�{     @�~�    17:00:30        G @�~�    @ႀ    )?Q      �     C�55�6Td�B'�fB'�f7_�KD"g�    (�r_+�Dj6��?l]d?kČ    312�'�=�/.���~        /��*�r>��>>�!a>?��>�k=�G=r@*<�[�54jz:>8�<&|�<��=H"=�A�=��=��F=�ػ=��=�v	>x>w�@�e�@��@Z��@*�@�W@
��@Y@@��>@��@x�@x1�@xF�@xY�@xg@xk$@xdA@xS�@x@ @x1�@� <	J�{@��                    E�&�6ܼG��}F��F\>E�"E	,�D9��CNA�BQ��                                                A��C�5B��QA�BAI[@�g�?�l
>��=��                                                                    EG��@�,F%��A�	�3��oFc{RB.o4AMV�.GS�                                                    {@��@\p�@\p�@[�7䛕��%y%%�ا<`/PK%��>�]d(�r_>>*���Kp1��/    >�y�@δ81ȣ��L�RA���C�        B�^�B�^�C�1C�1Ct�
>�|J@[ 2��)C٣�5��6�?A#��B���>�y�B���B ��C-WAYLPA�7BㅧA���B��3�����Q	%#�4�����    ����@��A�.�?��?��A���@���A�@0w@���A���@���F�Q@3��o1�o.3lV3��7˥#5 ��5L�F�PGώGf�G��>e�O    C�L�@�R+=�h�>�7�>���>�u�>���>��>elE>[8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�s�    �8Ui�8Ui{@��7i��8�@H8�78bC�7��7,�x6u��5��4���                                                A����+4 �+4 {@��C���A�zBsȅ@M��>\P?{@��{@��    C��RBY �5��F    6EM-7X�@�V�3�7��    =�q�Ch�Q7��BT�{D�j�C���B���Bc�A�ҽA;��@�<x?�t'                                                ?p�A�I�@�;�?���?���?��>��=�۶<�tw                                                                    D;�/F~��E�
�D�9DS�C�&B��B ��A	��                                                @�ƷB�-�B67A`7�@���@w�~?У�?G�>m                                                                    C�wE�=�E-LD�C��B��_B2�Af!@w"                                                @@ƼBM�A��l@�O@l��?��7?MZ�>�ua=�4b                                                                    6�`87w@A!E>,yvA��?�?�5L*��x)�T�,T�+�m%*��x)��,��;,���,?)f�~7W`�3[@0�I�0M��            �G��3��1�L�1�0�I�0M��6G��1��/�G��3��N    (a��3Q��3
A�2���    >g��77��6�5S1��0K�e0"�                        +Թ�1�r[1mg�1+}    1�#2���            3�U�:��Y-
L�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C�G��>,�-A�o�                                                B�                                          B��                    A�                @�p�    .�
O*��3s�&5���3��V3�%�>�ki7��?XB�6ؐB;b�V;f�=                5��3    5���1A�l    8��b{@��3��;4;n�    ��D^{@��{@��{@��8NW�7(            7���    7���{@��    7���    ����7�	�{@��8�e�{@��6j�6w�5�]    ���X{@��3�u�6-S�6���6���<�:�    B��[6��F.C�p1@�e�?�Y�B�#�            ?a'�3[@
9��=3[@
<@�N<&�E@�1@W�h?���?��>Wv'=�4�&+I�e��>�� ��|芬RЊ�e��{Ȋ]��E��2�#@�>�s�<@�8<&��@*�K?��?U�->�Np=�G!=.m֋&+I�e��>�� ��|芬RЊ�e��{Ȋ]��E��2�#@�>��    *0I85E~ 5�TO?�h.?��?3��>���=ȥ-=��                                                �CQU�Mf�N�!�[�0�j����ƐZ�ƣp�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7;��5C&/=��923	7K�OA���A�A�;FA��^>'�e8}6�)�G�3�6?0+?V� C��DA2�?V�c6?0+K&d8�+�7C��CU�D�V�DȶDDJ<1C��:C"�WBt�JA�F&@�P                                                A��C^V�CH�DB�<1B@�:A��W@��J@F&?9P                                                                    E��0GyG_F�lF8��E�0PD�qFCC��                                                C��FEP'�EA��E:}D�W�Cԣ^CƿB-|&A<�X                                                                    E�;�F�F��UF{e>F]5%F/"E�N�E���EGʄ                                                D)��DávD�Y6D��DD���D~�ID0�wCߏJC�Mw                                                                    1�#6���GYM�GNkG�vF��xFz��F�}E�9tEPB�                                                                                                                                @�|A!��A���A�SmB�BɔB3<�BMV�<�<�<�<�<�<�<�<�<�<�<�<�E��E��PEc�E��D��DWOC�{[C���                                                {@��{@��{@��{@��{@��D8Z�Aœ�.Yxݧgh�7��^B�)�?���    @�@�{@�ξR�T�R�TC��{@��C���C���C��{@��@a�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�:Z�OD�Z�Gx�D�Z�D�$�D�$�A%#�A%#�Fc{B5+�    C_s�FK�
FK�
D��D��Fd�B5&k                ?�=�C�NtC��C�|)?   C�T�C���C���C��C�3;C�E4C�X-C�n�C���C���C���C�ؒC���C�cC�(C�H�C�m�C��`C��~C���C�qC��C�.yC�*�C� \C�#�C� �C�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�v�C�6t>�B�>�Mp>�e�>�O�>�H�>�o�>��g>�k�>�܁>��>�Q�>�µ>�T�>��>�/�>�o�>��>��&>��T>�K@#^�?�2        A.��3kk�3P��35��n��Ǉ�,Ǉ^WǇ,hƈ����_��_��_Ɖ�!        >YA��    ��zݿ�z�{@��@x�>FX&q6@q%71	        7�8_?�1>���>�(}>��>�<3>԰�>���>��I                                                ��D7Z<�A�  ?k�A	��A z�)?Q      �     C�8���8?87�� 6�f06�S57�w47c3% A                                                                    7���79�6�Y�6h�58[�4g�3gE�2Pl                                                                    G��}F��F\>E�"E	,�D9��CNA�BQ��                                                4AG53���3K�2a��1���0�B/���.��g                                                                    3t$2�Ė2*�f1��Y0��/��.�n�-�x�                                                                    6�Nx6��6I�5���5�t4UaM3|�c2��:                                                                    7!��6��6�x�6��5�"�4���4�P3�R                                                                    4��4�}=4;�3Ә`3RJ2��V1�{�0��o                                                                    7E�47��6��61�5�*h5 ��4�)3��                                                                    8*Vb8,07�D�6�T6=L5��84��,3���                                                                    5��5��4�ȧ47��3�ݧ3
��2*��1$�2                                                                    8*Vb8,07�D�6�T6=L5��84��,3���                                                                    7���7l�-6۲�6BZ5���4�54	93�9                                                                    4���4�,i4(�3�l�3-^�2�_M1���0���                                                                    7���7l�-6۲�6BZ5���4�54	93�9                                                                    �7�γ�/A5
.74�n�4&Z�3^�I2�*X1ϝ                                                                     /��00i�;$چ:�L�:_�t9��W8��C8:�_                                                                    5i	G5"�I4�h�47�c3� W2�
I2"�1(m                                                                    6��6�s5�5<"�4��3�c�2�	�1Ф�                                                                    ��f���q��=񃴼6�7w�3��@@!�9Kj                                                                    �d=�Y���M���-�n�g����/�?���&O                                                                    5��5�9c4�F�4A{+3�|2ϟ�1�|�0�z|                                                                    �� ���Z3��嵈A��8T�'�:V��AF�                                                                    �)�������5���ñn���Tv��ͮ���                                                                    ��c��37�����'߲h�[��Z���{�����                                                                                                                                                                        +�_a,Y�                                                                                            8xB7��7o�@6؋�6&M5e�4~Cx3�ݭ                                                                    4���4p7�3�:}3]�2��1�Ts1�0^�                                                                    7�d6迧6j��5��{5,K4py53�f2���                                                                    3�*�3kR�2�J�2\U�1�ܠ1P�0�5/#�~                                                                    8�7��P7e��6�s261�|5{��4��}3�PQ                                                                    4�[4fm�3�[j3Z�o2�
"2w}1+�s0223                                                                    7t�7e�j7"i6��@5�05	�4�3\�                                                                    5���5�,�50'S4���4��3&�Y22K1+נ                                                                    6Nm�69�5�*5vִ4�3�3�2�!�1�a                                                                    4k�?4T 4X�3��2Ө�2�X1H0
��                                                                    7��7�H�7<b�6���6%52}u4>��37ƙ                                                                    5�i�5�S65WL�4�V�4 N3K�2Y��1R�                                                                    495�40�4&��4O3Ӥ,3��c3�2��b                                                                    3�"3��3�L2�2��2\��1���1�\E                                                                    4b^!4W�N4L4,p�4V3�!�3=��2�X[                                                                                                                                                                        )?Q      �     C�68�(a��                        6���3���@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�~�    @ႀ    17:00:32        G2 @ႀ    @�@    )?�      �2     I`5���6��~B'�fB'�f7�?=D"g�    '��+��R6O�o?�C?|�M4p3��v1��03��.$</        -�U�)�n�>��>B�U=�h�=b�:���3�X,�2S�G:B��<')1<��=H�=�C�=�7=��Z=�ٍ=��=�v�>�>w�@�1�@S�K?���?n��?M�>�>���?P��@�g�@�@W&@}s�@{��@ze@x�@x#�@w�@w��@w��@w��@��<��{@��                    E��7b�UG�b�F�ӇF[�rE��
E	+�D9�MCNG�BRN                                                A�|CY�B���A�e=AI_@�f�?Ì�>��=�"a                                                                    EG�$@��F%ptA�� 4��FcGB-�A$i�/Il^                                                    {@��A���A���@\B\8C}֤�����(ڧt��/7�!�08W>�1�'��@��G.F�34?�"    >���A���1�j�?���A�fE@��        B�k�B�k�C���C���C�Е>�!�@�2��zC��K6Q-6��\A!��C�e>���C$�B���CD+�A��BK� C �A��=B�<PA����9�!Ȼ|A���    A���B+�A��$>���>���A��NA<8A�y�@@��@H�&A2
@zF�"�4��1�P-� h3@�8b��6.R5ӊqG��GJ�G��2G��>6�    A��L@L��>u�=>� �? 
>��z>�P�>�b>���>��8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�uW    @�r�@�r�{@��7��n9|�l9bB�8��z8V��7�oo6� W64�R5` ,                                                A���@r��@r��{@��A�� B"��B�c�@D֡>�Y�{@��{@��    C��(BX:5���    6A�e7WG�@�S3�~7-@�    =��C��7��~BA�	D�@(C��~B�ɢB>�A�JaA%+�@Zl?|�O                                                ?Z	�A��U@��?�q�?y��?i,>o��=�͵<�`�                                                                    D:l�F~��E�nbD�b�D�C���B��A�EA]j                                                @Ǽ�Bܟ�B3!GAW��@�ٍ@o��?��<?
��>H�                                                                    C���E�ߓE��D�WC{�B�UAB.�sA`O�@r��                                                @<�%BL'}A���@Θ�@c�_?�?H�{>��I=�po                                                                    7{��7��dA�L>'�xA�;�?�g?��_+��*H�,�u,>G�+��*�j�-e��- �,�V�)�b�7��i4���2I�+/4
�/G�/F�b+�ۇ6�4��v2!jN2IR�2F��/3���A��6�64g[�    ,���4�l.4�k61�>    =��8A^7�g3n�S2"�0(�                        .�V`3E��3:��1(gu    3o%~2�Ig            4���;m�-��?�  ?�  ?�  ?z?j`I??rp?n��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CGP�G�
�>+*A���                                                B�                                          B��                    A�                @�p�    -n��)Jw�4��N6:W�4�d3@��?wV�8�o@,�7�m;��9;�H7��            7O�X    5�ł6[�    8��{@��0.43_9e    ����{@��{@��{@��8n܍7�o            7җ�    8I»{@��    8I»    5Q�G5�r.{@��8xc{@��5���5G��4��    5U�S{@��64Dz7<�7�-H7�-H<��;    B�6�Y�F�LC�qoA][�?'zCܫ            >�ho4��N9��4��N<� j<���@ �#?�j�?M�5>�,@><�Z=�I��&+K�e�>����|늬Rӊ�e�{̊]��E��2�#@�>�3�<���<�z@ �J?��?6��>���>�f=I%�&+K�e�>����|늬Rӊ�e�{̊]��E��2�#@�=A�11��*6s<5��6R�;d-=L^�=�r�=��=g�=1�M                                                �(���I��Ea�֡�ó��Ë�L�]0^�;C�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6W��5w��<�E�8�â8��A��$A	�vA�lfA�/�=L��7�;�6�8 E��5���=��A��?(#�=*X�5���Ju�F8�a7�ACV�'Dއ�D�(�DK��C�]�C$gqBw��A�=@���                                                A��'C^��CI(�B˞�BB]�A�gq@���@ =?;��                                                                    E���G(gGj�F��F8��E�?aDȋzCC�B                                                C��EP:�EB�EAVD�_�CԹGC��B-��A<Խ                                                                    E�<�F��[F��aF{gDF]5�F/PE�OME���EG�H                                                D)��Dä�D�\0D�ֽD��GD~��D0��Cߏ`C�N                                                                    3o%~7_��GY]PGN+�G�8F�GFz��F��E�=�EPFS                                                =v�z                                                                            @�Z�A�@|A�'sB�BZ�|B��]B�^�B�Lt<�<�<�<�<�<�<�<�<�<�<�<�E���E��Ed	�E�D��DWYnC���C���                                                {@��{@��{@��{@��{@��D�B,Va-�_�f��8�mB�6?�P�    AW�AW�{@�ξ�������C�RE{@��C��#C��qC�RE{@��@]�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�<H��D�MNG`�D�M$D�h�D�h�A ��A ��Fc_B5�=v�zCΚdFKХFKХD�%/D�%/FdóB5L                @EC�7�C�wvC�Z)?   C�baC��&C��&C���C���C���C�gyC�H|C�-:C��C�dC��C�	{C��C� �C�6�C�S�C�u�C��mC���C�߅C���C�C�(C�!GC�#�C�!	C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�r1C��>�	>��0>��>��H>���>���>��`>��.>��>���>���>�"�>��V>� c>��&>���>���>���>��>�Mw@=&�@	��        A�u��D�w�D�7�D3��^t'��e����b��V����ƻ�ƻ��ƻ��ƻk        >�j�B�U    @K �@K �{@��@��tFW�[6�׻7�        ;I\)?\��?]݁?a��?f??m&?u�o?{�b?}                                                ²J7��A�  >�)?^��@�2g)?�      �2     I`9:`�8���8�h7`v�6���5��!4�;�3��                                                                    8kl�7���7&�36��5��u4Ӿ4��3x�                                                                    G�b�F�ӇF[�rE��
E	+�D9�MCNG�BRN                                                4�.�4�83��2�K2�	10t%0v�/��                                                                    3��{3I�}2�p�2<�1Jj0^�/���.�ވ                                                                    7X��7AR�6�y�6H�5��N4�.�4#=�3Q�a                                                                    7��a7_�{6��E6}X[5��L5*94�Π3�a                                                                    5025	��4�?C49*�3�S3	%�2j�1���                                                                    7���7��=7�U6��p6��5O�4�Q�3�v�                                                                    8�KA8���8
�7p6��5�U�5@Q+4Y�                                                                    5�P�5���5-�z4�+\41�3w
�2�Ư1��                                                                    8�KA8���8
�7p6��5�U�5@Q+4Y�                                                                    8�|7���7\	�6�b
6%;�5T�4�$�3L                                                                    5
)�5w4�[�43��3�M�2�2P��1q%�                                                                    8�|7���7\	�6�b
6%;�5T�4�$�3L                                                                    �jH���G5���5k�V4�:3��3Z�2���                                                                    1�^�20��:��v:YD9�Z9UӇ8�{�8W)A                                                                    5�;)5�(5]�4���4٢3K�2�,�1��A                                                                    7�'76�6h��5��-5w)44��3�_�2��                                                                    � 9�m��=
�E<W��!��QL��ݲ٤                                                                    ��3���4��i���*O�����/�����7�                                                                    6a�6<l5l�y4ƶ 40#3:��2�p�1�j�                                                                    �+4�h����ж�>�j}贐ֆ��4��                                                                    ���a��=�Jղ�S����̱�����]����                                                                    �����pԣ� �+�����򒿲tǱJ�#�h8�                                                                    )��*@�}.��010�n14��1P�K1$9O                                                                    +�� ,_-�                                                                                            8��
8u�X7��q7dy!6���5�*5$�s4K,                                                                    5\4��14zz�3�$�3:�2dݚ1�h�0��"                                                                    7���7p�06�(6b�5���4ݠ 44qz3^�                                                                    4	�3�s�3u`52�xU2A%k1p40�eZ/�'�                                                                    8�*�8k�7���7a4�6���5��5DR�4r1�                                                                    5��4�c�4pE�3���3G��2{�s1�a�1	,                                                                    8��7�Q7�nb7 K6n�t5�V�4ǳ�3��                                                                    6n�6w�5�p56�4�r�3�c`2�;2��                                                                    6��`6���6�t6U�5@��4W�3�`2���                                                                    4���4���4�P�4ϱ3\�h2vH�1�m�0�=                                                                    8�a8��7��\7C�x6��G5��4�v4�                                                                    65j�6%�05�L�5_��4���3�@�3yh2��                                                                    4�D�4�7�4���4�Ҙ4\lg3���3�g�3��f                                                                    3�#3�?73���3n�32�2ɱ2�S�2[MA                                                                    4�]4޶4ѮZ4��4��4�z3�~�3���                                                                                                                                                                        )?�      �2     I`65Щ,���            Eٿ]    ;3	B6�z3�+�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ႀ    @�@    17:00:35        GQ @�@    @�     )@      �Q     O05�Pv7PMB'�fB'�f8&�D"g�    (��++w_�7
'#?�  ?t~2T�{4ԫ)߁3/�^����        .�8�*��>�E�>#H=���<`��8�6�1�Qk+!�2f.�:@X�<'��<�$=HH=�DM=��=���=��=�=�w,>>x6@ҏ=@w�@�r?��-?�A?d+�?*-�?%5�@���@��@{�@z�M@z6@y}�@x�@x|7@x7+@x�@w�?@w�?�=P�??��                    E��C7x�}G��$F���F[%OE�h�E�D9yoCM�BQ��                                                ASC�DB��A�tAHq�@��?�<l>؝	=ܩ*                                                                    EG��@��F%m�A��;2��Fc
�B.�A5�/W_i                                                    {@��B�mmB�mm@_�w8��%{Y�    �H�/2�,'3��>���(��+@�&����2�sH    >��B@|�1�S�A,�y@��(@��        B�jB�jC�s~C�s~C���>�X6@��22�(%C�޷6�K6�cA"C�CO(V>��Cd7�B�C�Cw�A��ZB���CS^B�JB��RB?zr��u    B?zr    B?zrBE�~=�5D��8MM�A�yXAM� A��n@v�@%��@Ś ?��F�2��/�~I,�K!2 ��8��06�/G5酋G�~GV��G�#`G� �=�͏    6���    >>��>�l>�5<>��5>���>�N�>��>�f�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�F�    A(��A(��{@��8  �9�(9s�]8�z[8^�#7�Y37�I6CS�5ma                                                A��DA�CA�C{@��;���By��B���@?	]?�b{@��{@��    C�v�BTF�6|d(    6AS�7n9@��?2���7R�    =���C���7�ϡB'K�D�~�C��	B�OBxWA�ZA��@&��?;I�                                                ?:'�A���@���?�P?> �>�">>X�=x��<��                                                                    D7r�F~sWE�D�FDcPCx��B���A�2OA [�                                                @�!�Bٿ�B,��AI�N@�2J@_��?��{?�>}                                                                    C��!E�L�EhD�2Cp
B�m�B'm�AT\@d��                                                @5BI$�A���@�~v@T@�?��X?=M�>��=���                                                                    7�w;8�'A3�>)�Aߔ�?�'?��W,r�g+�=-7�B,��,r�g+C�Q.G:V-�t -k�*b}7���4Ϫ#2v��/��/��b/��,`B�7�ل4��J2A1i2s~�2qH�/j���ل1��*7�ٽ4��R    3�.4��4�3�1W��    ;[P}8X��8/5�6���5-Z�2qm                        5��6d�w6`&N3�%/38��7$32��            4�>Y>�9-1�Z?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C�?G�T�>,PLA�y�                                                B�                                          B��                    A�                @�p�    .�0�*�m�5;�6H�2�P2$uW?ͺ�8O��@h;�80�;�V�;�5!6�w            7�i    1�|�6�!�    4���{@��    /�8�    �ci{@��{@��{@��7�7��@            6���    7C;u{@��    7C;u    6���    {@��4���{@��4l8�4*�3��    6���{@��6C��7�g�7�7�<��1�G�B���6�R-E���C�[3A��i    C<i�            =Wn�4�a�7�@4�a�<�n8<�+�>U@M>[0=��Y=2�6<�(<�b�&����܋�6��-��Y�����܊|���^v��F��3���#�=N��<�g�<�Q>T�s>F=�a�=2A�<���<MG�&����܋�6��-��Y�����܊|���^v��F��3���#�;�x.���)�5�6M��6��:��9��#96�R8��8D}D<~�                                                Ć$�ă"N�y�m�f��J<��$�����dÇV�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         ;�;�        {@��{@��{@��{@��                4���4[6r        :ʱ�4���            CW3�D޲yD�{�DKp�C�d&C$��Bx�A���@�=@                                                A�3�C^�yCI{�B�p�BBd&A���@��@ ��?==@                                                                    E��EGH�G��F��F8�TE�G6DȤPC�֔C�                                                C��{EPi�EB8�ED7D�_�C�ĬC��B-�TA=                                                                    E�>!F���F��)F{k�F]6�F/ �E�NlE��hEGɰ                                                D)��Dì�D�c%D��D��"D~��D0�7Cߍ�C�L�                                                                    7$37w�GY��GNL�G��F��Fz�kF��E�?�EPGp                                                :{�                                                                            @n5�AvSGA�Q�B!B;G�Br�B��?B� ]<�<�<�<�<�<�<�<�<�<�<�<�E��E���Ed�E�cD�
�DW`C��C���                                                {@��{@��{@��{@��{@��D�B��.:ҍ�f�8kC �@{�    AP��AP��{@�ξc6�c6C���{@��C���C��C���{@��@_�w{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G3;9�X;D�LQGnKD�L'D�YD�YA�#A�#Fd��B7��:{�C���FK�FK�D�,�D�,�Fe��B7�h                @�j0C���C�]�C��R?   C���C���C���C�e�C��C��lC���C�`�C��C���C��lC�g[C�5�C��C��^C���C���C���C���C���C�˦C��8C�3C�".C�!�C�#�C�!C�[/{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��hC��9>�8�>�]X>�W�>��a>�l�>���>�W|>��z>��>��>��T>�!V>�ͪ>��3>��>� >�M�>�:�>���>��@$��?��A        A�:ǚH(Ǚ�ǙL��ޟ��&�������z�1��~���~���~�Ą�B        ?J�B3�`    @�c@�c{@��@kKgFW�=7<��7�        ;�Y�?CCZ?D3?E��?H��?M�L?U��?b�&?tq�                                                �l?K7|x)A�  ?C�w@��	@�_w)@      �Q     O09NN�8��78�f7qZ�6�V�5��5;z4
�*                                                                    8�L�7��71��6�o25ٱ?5�j4'�3/|k                                                                    G��$F���F[%OE�h�E�D9yoCM�BQ��                                                4�(�4-\�3��2��25k�1v\�0�.�/�Ii                                                                    4	'[3Z�2��2 w�1e)�0��	/��.�Ȅ                                                                    7q6�7R��6���6Xk�5�,�5	��49nq3j!9                                                                    7�r�7T96��q6P��5м~5<~,4t�3��*                                                                    58]�5�4���4d�3��S3>�2J�1kS�                                                                    7�6�7��6�ы6~�E5�E5faR4�kJ3�@3                                                                    8�Ҫ8��8r�7x�,6�V&6#%�5OV54d��                                                                    5��5��5,��4�Pj4=-�3�:+2���1��                                                                    8�Ҫ8��8r�7x�,6�V&6#%�5OV54d��                                                                    8��8��7_��6��62ƌ5��<4�O�3ː�                                                                    5nD5"�4��44i�3��3 n?2\��1x��                                                                    8��8��7_��6��62ƌ5��<4�O�3ː�                                                                    ��L�3�C�5�v�5��35��4;�3��2�?>                                                                    5^.5�Q�8��j8��t89Y�7��7m�7E3q                                                                    5��5�6�4��4�_�3���3`�^2�1��{                                                                    7,2/7�H6o��5�H5#
�4t��3�a�2���                                                                    �1�¶(���ř�T�U��I�����X��rg                                                                    � ����O��ϲ�jN��d��W�
Y����                                                                    6,�6��5r�4·�4'�B3|X@2�N]1��7                                                                    �>"��%����� y�����ˮ;�	%�&RW                                                                    ��=���[�$����V����\�,��Y���fU                                                                    ������4�,2��N7�	矲=���e������                                                                    )a��)���-���-,�4,��,A�+�b�/�
P                                                                    +Q�B+��                                                                                            8�U�8��8P�7vWr6���6y�5;#64b�G                                                                    5
p5~4�ǋ3�t�3T��2��A1���1 ��                                                                    7�0y7�G�7��6t��5�ێ5��4M1�3x�Y                                                                    4�|4�J3��2���2[�1��0�n�0��                                                                    8�k8���7��7r��6��B6#/�5_@v4�_�                                                                    5��5��4�S�3�ج3cV�2���1��1X                                                                    8h08'K7�i�7,��6��5���4�'3���                                                                    6%	\6��5��y5E7�4�$�3�5@3:�2�                                                                    6�b�6л�6��K6r�5[d�4���3���2�-                                                                    5\�4�I4�t�4^`3z�h2�J1��D0�|�                                                                    80W8ڱ7��7R�6��5��5
4K4*                                                                    6I�64g�5�5q�4��{4K53�21T�                                                                    4�.�4�h�4�r\4�b�4z�4.��3��.3���                                                                    3�	k3�T�3�=b3��3J{�392�p�2sY�                                                                    4��4�@4�6q4��(4� �4UL�4	6�3��                                                                                                                                                                        )@      �Q     O0653 3�.            G� �    =�y>6�,3�s�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�@    @�     17:00:38        Go @�     @��    )@}      �o     T�5Ò7Ue�B'�fB'�f83j�D"g�    )� +N%p7�Z?�  ?q��    47�)0�u/�3߰``        /</�+.�>��o>��=oM;�� 7�.O0�G�+��2i�X:6�c<'�B<�=H�=�B�=�!=���=��=�>=�wq>3>xf@��@���@��?��>?���?�K/?�e?b�@}j�@}��@r��@s{{@t@t�f@u1�@u��@vQ�@v�h@w�@wM-=�b=�Ȇ?�0                    E�j7�^aG�9|F� �FZi�E��.ErYD8�!CM7fBQq                                                A�HC 7B���A���AGԭ@���?°�>� )=�&�                                                                    EG�0@�!wF%|A��0��VFc1�B.�eA;�D.�ѱ                                                    {@��B�[B�[@b�G8����rM    �zr/B�{&ݡ�?   )� ALMv���o1�si    >�y�Ba >1�X�AB	�A�eA��        B��&B��&C� C� C��v>�q@���3 ��C�<u6���6�=uA#�CR��>�y�Cgx�B��.CuKA�2+B���CO��B#�B�-�BQ����'    BQ�    BQ�BXh@;���    6fqA��GAL�A��(@�|'@��@�(�?ԥ�F��0��V.�+���1��8�xw6�4�5��G!a�GWNG��fG�"=��            >+��>��>М?>ȦO>�h>���>���>�:�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���MױM    AX7�AX7�{@��8�|9�_x9zu�8���8_��7�J�7�%6;,55_\�                                                A���AT=LAT=L{@��9z�!B��B�8�@=:�?�{@��{@��    C�`�BS�6��    6E7	(v@��32���7cX    =֗gC���7�yBr�Dp�!C],�BQ�A�Ag[�@̗I@�?@q                                                ?!��Ao��@q�x?}��?n>�@>�=E�<\p�                                                                    D4��F~X�E�6D��DCna�B��AA�-@�r                                                @��^Bօ�B&SA<��@̤�@P�?�c>�>�                                                                    C���E���E�C��5CeqRB� _B��AH��@W��                                                @-�QBE�6A�J�@���@E�;?ˬ?0�>mڢ=��u                                                                    7��i8D�A��>-�)A�l�?�;?���,�N�+B4s-g��,�D�,�N�+�xM.���-�ZS-���*�V7�pb4��2b�,�]t/���/���)�Э7��5 (2%�2^�2]�1,��1���1�si7�K4��    3�� 4��4���/7E    :���8_�!8H<6�F5��g2�b�                        5�ŷ6�@J6�4G�3��7�X2�^~            4�"u>a�-8��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C}XiG�O>+,�A��                                                B�                                          B��                    A�                @�p�    .�>J*���5k�"6Lr�2-��1#�?�_8^c�@�^�84�f;�V�;�gU                8#�    .1��6��@    2K�g{@��            ��X{@��{@��{@��7RD7�%^            6�]�    6�]�{@��    6�]�    6���    {@��2K�g{@��2K�g29|0�V    6���{@��6�P7��7�7�<�F�2��B��K7	|Eū�C�e�A�3    C?��            <s4�|�6jl�4�|�<�p�<ܗ�=!��<��<>Լ;�
;xe:��Z�(�M�r�ʊ����x������r��`��H��5y��%�x<;r<�h�<܃X=!"<�
�<=��;�̖;{~:��(�M�r�ʊ����x������r��`��H��5y��%�x8b �    *F�]6}�e7%1�9=�,9ag8���8�h7�s�8P�                                                ķ#ĵ<�Į��Ĥ>rĔ�)��t�M�t��>̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         9z�!        {@��{@��{@��{@��                                                    CU��Dޜ�D���DIVC���C#d
Bw�A��u@�B�                                                A���C^��CH��B�VB@��A�d
@��@�u?<B�                                                                    E�͒GlG�jF�9F8��E�J�Dȸ^C���C�                                                C��EP�EBbECD�[QC�ɞC�sB-�_A=,f                                                                    E�?^F���F��UF{p�F]8-F/ qE�L�E���EG��                                                D)�sDõ�D�j�D�ݵD��
D~�HD0�Cߋ�C�K�                                                                    7�X7�4�GY�eGNftG��F���Fz�PF�E�?QEPG                                                                                                                                 @V��Ac�=AÒ{A���B0PQBb�tB�
�B��><�<�<�<�<�<�<�<�<�<�<�<�E�0�E��Ec��E�D�DW_�C���C��                                                 {@��{@��{@��{@��{@��D�>B�,.�x�f�f8$B�Cc{�@<�|    AvNUAvNU{@�ξ��]���]C���{@��C�+tC�+|C���{@��@b�G{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��Gn�9���D�R�G��D�R�D���D���A��A��Ff&.B:�E    C�p-FK�CFK�CD�-�D�-�Ff��B:}                @�H#C��?C��C�2{?   C�Y�C�+tC�+tC���C�wkC�4�C��.C���C�R�C��uC���C�uNC�5�C���C��`C��EC�J�C��C���C���C��C���C��7C��C�!�C�#~C�!#C��a{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�[�C���?3�>�W�>�}�>�>��N>�%�>�E�>�=�>�1>��|>���>��f>�[>�q�>�>�A�>�9|>�>��c>��@>h�@�        B�ǫi�Ǫ¿Ǫ9��|��,�+t��*� ƍkg�Ģ��Ģ��Ģ�ĵ��        ?$oRB5�W    A�sA�s{@��@�KEFW�7@�]7�        ;�c?7�U?8'c?9~Z?;��??]#?D��?MHf?Z+�                                                ���7/z�A�  ?k�A	��A z�)@}      �o     T�9W��8���8D�7z �6�h�5�985�E4�                                                                    8�=c7��78��6��5��5t�4&�W3+�M                                                                    G�9|F� �FZi�E��.ErYD8�!CM7fBQq                                                4�44��3���3�p2:��1z�}0���/��                                                                    4g�3dB�2�|s2&7J1k�0�_�/�@E.���                                                                    7}5�7[�Y6�]6_r55���5�48f�3d��                                                                    7��p7?�6�h
6-V�5�k�5NV4@�3Z�                                                                    54+ 4��4S��3��!3�m�2��2 tn15��                                                                    7�m7j��6�b�6S��5у�57�4k�3���                                                                    8ʾ�8��48�7x�K6���6�5E8V4U �                                                                    5��5���5&�V4�N�45�=3�*2�c�1��                                                                    8ʾ�8��48�7x�K6���6�5E8V4U �                                                                    8�8%7];�6�+�6/^�5�)�4�'�3�^�                                                                    5��5��4��-4./Q3���3M�2M9@1b΍                                                                    8�8%7];�6�+�6/^�5�)�4�'�3�^�                                                                    �Awb4�"�6F5���5�24U��3�0�2�                                                                    5��6�h7��7<�,6� �6*r5�=15 �@                                                                    5�`5���4�4\��3�i�33-�2ak�1F�                                                                    75^�7��6pC5˟�5"}�4p:*3�g2�}z                                                                    �;��&�7��՝�\����`�꿾�_�o�                                                                    �V���ߏ���Ĳ�޲�c�[��
Kα�ҝ                                                                    65��65�5r��4��V4&<�3v�k2�@1���                                                                    �G�'�,'��S��������Է�x�"%
                                                                    �ɾ��l�*/��������`+C��0���                                                                    ��(���e̴3�����a5�A�I�e��|�                                                                                                                                                                        +�W�,��
                                                                                            8���8���8�A7}��6˶�6��59��4]@1                                                                    5#�U5HJ4�$4�3Y�a2�io1�d�0���                                                                    7���7���7؍6{�5���5=�4K��3r�y                                                                    4 k#4
gQ3�Sw3�2a9�1���0��0	f�                                                                    8�]8�
�8�7z76��,6%�5]�4��`                                                                    5�5�X4���4 '�3h��2��$1���1~)                                                                    8�U8@57��l73F�6�'�5�~84�2�3�*                                                                    6-n6��5�5L�4�-\3ٴ�3A�2L                                                                    6�A/6ږP6�ם6ޒ5b�`4��3��\2�:z                                                                    5%@4��\4�?�4%��3�o�2��1�� 0���                                                                    89yL8%Nz7��7[p6�L�5��(5
;z4-�                                                                    6S�W6<��5��5zj�4�ũ4
�3��2-�:                                                                    4ۥ�4Ϣ�4�nn4��m4�V�42 3�r�3��n                                                                    3�~U3��'3��(3���3Q)3��2�_82nt�                                                                    5:�4�ƥ4���4�h�4��4Y�m4	)�3�U                                                                                                                                                                        )@}      �o     T�69K
3��                         6�Z3ӝb@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�     @��    17:00:41        G� @��    @��    )@�      ��     Z�5�oJ7��~B'�fB'�f8#vgD"g�    )3�U+9��71v�?�  ?u��    4b8))��0�j�*C        /j4�+Zs1>�A=�5�=.�<:��)4Yv�-c��*��#2fsF:))L<'�<��=H=�?V=��=��b=��J=��=�wV><>x|@�L^@�5@D�?��?�2q?��A?���?���@w38@woq@m�@m��@nN�@o@oƼ@p�@qWD@r0@r��@r�=�h]=��?a�e                    E�o7��G��1F���FY��E�*rE?D8d�CL��BP��                                                A��CbUB���A�YAG!�@�,�?��>�^I=ۭ3                                                                    EG��@�A�F%�+A�1"Fc�DB/� A>�S-��1                                                    {@��B�.YB�.Y@e�8��$�9��    ��yg/�z&v�?   )3�UA{3��=1��    >�y�BI!�2�SAA�dA
�	A
T�        B��B��Cȁ?Cȁ?C���>�Bg@�3J�C�
6�:-6���A%r�C@��>�y�CS��B���CdA�A�	�B��TC>elBQCB�3�B���N�"    B��    B��B%�U            A���A<6�A�ӝ@vR�@J\@��J?ĿpF�%1".�{j+��j1C��8�\�6��e60�G1G6πG��ZG�� =�E�            >!�I>�Q\>�yk>� �>�� >��U>�@�>�"r8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���Mܘ�    A�ΣA�Σ{@��8
�9�},9��#8��88iU}7��76/6</5\�V                                                A���A�  A�  {@�ΠH�B�$B�'6@B%�?��{@��{@��    C�?�BS�6��s    6Kcd7
�:@���2�$�7y�!    =�E�C��57��iB��D`��C>YCB0TPA��AB'6@��I?�e0>���                                                ?��A^/Z@NQh?SE�>�]�>���=�\�=&v�<<��                                                                    D1�0F~<-E��D�pD�Cd��B���A�bE@�!                                                @��B�/2B^�A0�@��7@B�y?�j�>���=��                                                                    C���E�E`�C��C[^ B�J|B��A>s�@MS�                                                @&�BBK6A�y�@�$�@8b?��S?$�>^O�=s��                                                                    7���8	�FA�T>2m/A�k?�c?��D,��\+���-���-,�~,��\+���.��
.,2�-�L�*�(E7��4�K�2rA�-A��/��%/��s*u�n7L��5��2�'2m��2m[j-=ʌ�L��1��7L��4���    3�a�4Ǹ+4ǐ/�1�    :�Q<8=C86�~6ĳb5v��2���                        5nX6��m6��4$)�3��79��2��            5z�>4}Q-��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Cn�G�\�>+s�A�>�                                                B�                                          B��                    A�                @�p�    /	�+NE5Hf6X�25=O1Y@?�!�8F6�@p�	8��<,'<�                8�7    2a.6��        {@��            �=��{@��{@��{@��7p-�8Ȁ            6��x    6��x{@��    6��x    6�p{    {@��    {@��                6�p{{@��6ң�7��[84W84W<��1�B��7֖EٷxD"GA���    C0)d            <.�o4��6�G�4��<��E<稾=FU�<�<i��;ϹU;;��:�:�(gM�Zc�����Ҋ�T튮���Zc�~ቊ`���H���5[��%r�<,�<��%<�=E�<<g��;��E;7��:�p�(gM�Zc�����Ҋ�T튮���Zc�~ቊ`���H���5[��%r�8�3�    *� m6�7=Q�9�O�9@�}8Ѓ;8u�8��9�                                                ��J������Ӗ\�ȡĹ �Ĥs/ċ��^�)̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �H�        {@��{@��{@��{@��                                                    CT�D�G�D�
�DFr6C�KC �ZBtV`A���@���                                                A��C^G�CH
�B�r6B>KA��Z@�V`@��?9��                                                                    E��G�G�2F��F8��E�HD�ŀC��C&%                                                C�'eEP�eEB�IE=�D�RC�ƊC B-�A=N�                                                                    E�@�F���F��/F{vfF]9�F.��E�J�E���EG�U                                                D)��Dÿ�D�s-D��D��D~�pD0�C߉zC�Jl                                                                    79��7�	�GY�LGN� G��F��<Fz�F��E�=EPE�                                                                                                                                @Jx�AY�A���A�B*T.BZCkB�r�B�@K<�<�<�<�<�<�<�<�<�<�<�<�E�M�E��FEc�"E�D��PDWZ`C��#C���                                                {@��{@��{@��{@��{@��D�B�.�߿�fV�8.��C��_@b[�    A�nA�n{@�ξygx�ygxC��{@��C��C��C��{@��@e�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�\9XYD�[�G��D�[oD}�LD}�LA A Fg}yB=�    C�ʑFK�FK�D�+�D�+�Fg�'B=��                @���C�B�C��C��@?   C��/C��C��C�)�C���C���C�=�C���C���C�*qC��C���C�8�C��C��5C�H&C��XC��ZC�i�C�8<C�^C�RC� �C��C� �C�#YC�!/C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C���?�%?�+?	w�?-T?��>��0>���>��@>�p>��>�K�>ضN>�Tl>��6>Ȱ|>��>��>��>��/>���@*�7?��        A���ǟ0�Ǟ�3Ǟ"����2�������zN$��&���'��&���{3        ?9�7B�    AL��AL��{@��@w.pFX)�7o��73�        ;�&�?1��?1�?2��?4MI?7�?;J�?AC�?I�&                                                ��\�7��A�  ?k�A	��A z�)@�      ��     Z�9iC�8�hW8�B7�Q�6�!�5�~�5�74
y                                                                    8�SU7��7E�6�gV5���54�4,�83.q                                                                    G��1F���FY��E�*rE?D8d�CL��BP��                                                4���4Bq3���3W�2F
�1��p0�	/�i�                                                                    4�3u+2�L'21F�1z(�0��!/��t.��$                                                                    7�k<7l��6�Z
6l��5��'5G4=�3f��                                                                    7�"n7:�6���6�^5��5R14)��3>�                                                                    56��4�/D4@�Y3�B3m�"2��2�31                                                                    7��7d��6���6?��5�t5 ��4O2�3hB�                                                                    8�a�8Ů[8�7��6�k�6 A�5D�4O��                                                                    5��5�m�5*4�4��44�93��_2��g1�
T                                                                    8�a�8Ů[8�7��6�k�6 A�5D�4O��                                                                    8%s�8�7g�6Ϣ�61��5�J�4�3�}                                                                    5)-C5+��4��G4.'3� �3��2G�O1X��                                                                    8%s�8�7g�6Ϣ�61��5�J�4�3�}                                                                    �Q]�3�(c67�5�?�5"_�4oUg3�p�2�?�                                                                    5?ZO5�
7��#7m�6�We6VoD5ȋ�5L�D                                                                    5���5��J4��4H]�3�b�3��2F'1^�                                                                    7Fe87.��66�5�k�5'? 4t,{3�cN2���                                                                    �Ke��3�o��*�lc����|�����.0�u�                                                                    ��
	P���
��no���U�g�\��k��D                                                                     6F]6.�N5�CC4��4*\E3y��2���1��                                                                    �Xv�9��`��({���M״�3���i�#4H                                                                    ��۱��ѳ4Y������i=Ⱎ-Ư��V                                                                    ����[C�A	���ް�vβL)��o\V���s                                                                                                                                                                        ,�,��@                                                                                            8���8�-�8�87��6Ն�6�5>��4^�N                                                                    51��5��4�d4	4�3d$(2��1���0�B                                                                    7� 7�?7��6��5��*5$��4P�w3t/M                                                                    4.H4��3�h~39�2l�1�r�0�p0
L�                                                                    8��y8��8
�7��6�_�6,j�5cQF4��&                                                                    5*h�5��4�l�4>�3t_2��2 �1x=                                                                    8%B8	�7�v7?�56��5���4��3��                                                                    6<�p6&�5�*�5[:�4�e�3家3��2�                                                                    7`6��6��w6�5p�14�f�3�f�2�1Z                                                                    5m�5ޚ4�}41'�3���2��1ٙ�0ꁌ                                                                    8I��82}�7�?�7js�6�:�5���5��4-S                                                                    6f�m6K�|6��5��4�C\4\�3$�{21X^                                                                    4��4��4�� 4�@o4���4;��3��3�V�                                                                    3���3�3�U3�;�3^`�3�2��(2r��                                                                    5�P5�4���4ؤ4�,�4e��4�E3��?                                                                                                                                                                        )@�      ��     Z�6?3�a�                        6�_�3�"�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��    @��    17:00:43        G� @��    @ᕠ    )AE      ��     `p5��K7�1B'�fB'�f8��D"g�    )h)+:67Ff�?�  ?yк    4K@2(޶�/�ƪ�e�        /8��+-�>��>O=T�s;s�6��1/��T+t2W�~:�7<'�<�=G�3=�:o=�	�=��W=�� =�#=�w>5>x�@�8[@� �@Yj?���?�*�?��Y?��?t@t��@tӪ@j_�@j��@kA@k��@l�@l�L@m!�@m�@n�@nB�>4�6=���?z                    E��k7��G��MF���FX�|E�jVE��D7�nCK�BP�                                                A>wC�mB�`�A�eAFW�@��o?�go>֟�=��                                                                    EH�@�t0F%��A�G�1i�FddB0�A:��-��                                                    {@��Bo�NBo�N@j��8�`�$�1    ��c�/'���?   )h)AY��ѰN1͓�    >�y�BM�2��A��?���?���        B�:)B�:)C�a�C�a�C�S`>��@��23��C�B�7�6���A(O*C��>�y�C)�mB��CK��A�ֽBQ�SC'��BR�B��A͕B��3�    A͕B    A͕BA��z13�p    +Ɠ�A~7TA��A�Æ@Okj?�Ɣ@�1?���F��1i�/�,!�B1��\8��?6�:R5��G �PG7!G�[�G�'�>�:            >/�>��>�,3>��>�kd>���>�k->��x8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�'    A^'A^'{@��8&�9�:^9�@�8���8e>�7�U�7z26?z65f$\                                                A��~Ac��Ac��{@��/H��BBQB���@M��>��{@��{@��    C��BT(-6���    6RR7C�@�IL3�7`�h    >�:C���7�#SA�DR\�C$�rBgA�+BA*c~@���?�J_>�7�                                                ? ��AO�@1��?9ą>�1>m�X=�=�<&�~                                                                    D.j:F~:E���D��.DG�C[{�B��AЌN@߅.                                                @�ܠB��TBA$z�@�=+@5`�?�5>�q=�1�                                                                    C��UE�OAD�t�C�-CQ��B�SB�A4��@C,�                                                @��B>�6A�e`@���@+MP?�t-?��>O_�=dv�                                                                    7���8_�A�>9�]A��?�"?��,���+;-f�,���,���+R�l.St�-���-�)�*���7���4�92u
-N�/��/���*1�6��{4�٫2
�@2p�2ps�-�^���{1͓�6���4�9_    3;�R4ȏ�4�q /t�@    :�8k8�8��6u<f5��2�1?                        5�6Mz�6G�A3���3�"6�)�2�j0            4�r=�a�-�}?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CY��G���>*�VA�jH                                                B�                                          B��                    A�                @�p�    .���*⪲5�B6L:T2aϚ1�T~?�o8(�G@;�S7�\p;�9�;��\                7���    0�f	6��<    (5�{@��            �d�z{@��{@��{@��8�i8`-�            7]mQ    7]mQ{@��    7]mQ    6���    {@��(5�{@��(5�'�C'Q    6���{@��6��7u[�8b�%8b�%<���1�ءB���7
��F��C���A9�@    Ct            <po�4�u6ޙ�4�u<�U<�#=���=.E<��<V4;��;�q�)-��;���������Qʊ�r��;���a���I�a�61��&6'<ny
<�~�<��=���=-z<�"`<R�;�ϖ;$n�)-��;���������Qʊ�r��;���a���I�a�61��&6'8�l�    *а�6m�R7�9���9K �8�3Q8��8� .9Y�+                                                İFDİa�Ī��ģlyĘ��Ĉ߫�h��4-̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         /H��        {@��{@��{@��{@��                                                    CQ�D��rDǟDB��C���C�Bo��A�Y!@�ۚ                                                A��C]�rCG�B���B:��A��@��@Y!?6ۚ                                                                    E���G�:G�RF��F8��E�?xD���C�C9i                                                C�6�EQ�EB�`E3�D�B�CԹiC1B-�UA=j�                                                                    E�BFF��TF��F{|xF]:�F.�(E�H�E���EG�N                                                D)��D�ʩD�|cD��)D��D~�hD0�C߆�C�H�                                                                    6�)�7��GGY��GN�Gn�F���Fz�DF�:E�8�EPCS                                                                                                                                @[_�Ae�A�	oA�]�B/��B`��B�{�B�l<�<�<�<�<�<�<�<�<�<�<�<�E�k�E���Ec��E~dD���DWPC�x�C���                                                {@��{@��{@��{@��{@��D~4B�.���g=�8(ɃC�'?@��X    A�N�A�N�{@�ξf��f�C�DA{@��C�[�C�[�C�DA{@��@j��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��9��D�b�G�tD�b�Dx5]Dx5]A�6A�6Fh�B@�    C�bFK�(FK�(D�''D�''Fh��B@�H                @{��C�AAC�YEC�aJ?   C�o:C�[�C�[�C�C�C�/]C�:C��C��C���C�@C�FkC��C��aC��[C�B�C��hC��^C�=�C��C���C�d]C�;�C�yC��C� 
C�#/C�!:C�@F{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�o�C��d?��?��?�?��? �r>�a�>��s>�|�>��>���>�?�>�5>��p>�H�>�>v>�E(>���>�+�>��|>�$�@*Y�?�<p        A�@�ǀ�"ǀ�jǀ(�ŵ�o�G9����[��N�j���F���y���Eĸ��        ?�B��    A&o�A&o�{@��@y��FXn�7�.*7bk        ;�U�?:�?9�&?:�?<0h?>��?Bz8?H��?R��                                                ���6���A�  ?k�A	��A z�)AE      ��     `p9[(e8���8�7��V6�l�6995U�4��                                                                    8�j[7ᄉ7C�6�)5�̡5#:�45�3:��                                                                    G��MF���FX�|E�jVE��D7�nCK�BP�                                                4�+4;�3�>a3�!2I�O1��0���/�|�                                                                    4� 3m\[2�3�22*1~��0�ѕ/��i.�g�                                                                    7��7e�6���6k�5�'
5Kb4D��3t�S                                                                    7��97$�6~
�6�u5�b�4���4��34�                                                                    5!��4�V�4(��3�0N3T��2�
�2"`1�                                                                    7���7I&�6�?h6*��5��?5�4?f�3\Õ                                                                    8�"�8�(�8�E7z��6�E~6�:5D��4U                                                                     5�5���5;4�2B4,"g3�g2�>�1�>�                                                                    8�"�8�(�8�E7z��6�E~6�:5D��4U                                                                     8��8B�7\��6�V�6,�5��m4�~�3���                                                                    5�;5&I�4���4#�F3�v�3Ez2C��1Y�f                                                                    8��8B�7\��6�V�6,�5��m4�~�3���                                                                    �Z}R�R�6}�5�^45*"�4�q3��Q3gE                                                                    4��5�6c8�7�4�75(6��6 .�5��                                                                    5��5f|�4��.42��3�Q�327U1S%^                                                                    7<��7,�t6vO�5�P5%<�4s��3�:2��j                                                                    �@A��/:{��m�n:����� :[��Q�'�3                                                                    �
�I�Jq��2`��;���*��o�a�V\�Ɵ�                                                                    6<A�6,�g5v��4ѝK4'��3xD~2�c1��                                                                    �K�_�3�^��Hɶ&�ٵ��A�چ���,�                                                                    ��9γ��Ƴ06-���W�l�}��죯�܀                                                                    ��}&��s[�?�ճ�bƳ���Si4�{����]                                                                                                                                                                        ,]�k-L�                                                                                            8��H8�Ş8w�7���6�׸6L�5DC4j��                                                                    5'd�5bK4�&4�3dz�2��1�P�1�^                                                                    7�0�7��:7
��6���5�G5'0�4W3!3��#                                                                    4#��4a~3�3�_2llg1�4�0��0�b                                                                    8��8���8��7��>6�f6/�5j#%4�O                                                                    5 ��5`�4�3�4�3t^2��02��1�g                                                                    8281�7��m7AQ;6�;'5��4��r43�                                                                    62L�6"��5��5\�4��v3��Z3��2��                                                                    6�$16���6���67L5v�4�*�3�*a2��                                                                    5�5Q�4�h�42�W3���2�11�¸0���                                                                    8>��8-��7��7lF�6�i5�Z�5`4&x                                                                     6Y�6F��6��5��4ԣ�4XT3- 2>@                                                                     4�hz4��4̦�4���4��&4B�3���3�-�                                                                    3�%�3�{3�_�3�jQ3cm3/�2�y�2�>�                                                                    5	��50�4�!-4�m�4��/4m�b4�3��                                                                                                                                                                        )AE      ��     `p6E�v3;�R                        6�d43��f@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��    @ᕠ    17:00:46        G� @ᕠ    @�`    )A�      ��     f5Æ�7��B'�fB'�f7��SD"g�    (�_>+I��7A,?�  ? 2|�4-�)��,/�;���a�        .��*� f>��2>�/=��<C�8�4�1�f+ �2Fi:�t<'c3<�|=G�F=�5�=��=��"=�֋=�F=�v�>)>x�@�^H@o�'@E�?�:�?��(?]�?1��>���@u@u�@j��@j��@jڛ@j�@k�@kR�@k��@k�@k�@l�@"ʼ=|_�{@��                    E�1<7�PG� �F��iFW��E���E�<D7	dCKBORq                                                A�uC�6B���A�AE��@��?��,>��_=�MS                                                                    EHO�@ͮ�F%��A�x�30�Fd�
B2E~A36�-��-                                                    {@��A��A��@oے8a�"�!a�    �x�Y/=� ��{?   (�_>A',j��Չ2�z�    >�8�Aj�2^	@���B~��7/F        Bbn�Bbn�C�D�C�D�C���?'�@K�L3
��C⺉7�r6�t�A+�B�V1>�8�B�#LBm@C c�AJ�-AٓTB�"A��CBi��?ǽ����    ?ǽ�    ?ǽ�@)�G:
    4���Ah�@��mA�M3@��?~�@X8T?aK�F��L30�0��s-�Aa35��8 1S6�5�Q9F�9�G�G���G�i�>2s�            >Fy&>�f^>��>�LM>�&	>�ۃ>��#>�a|8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    A,KA,K{@��8m�9�ST9���8�4�8_u7���7��6F��5}0�                                                A��}A2��A2��{@��8�B	��Bb*�@\+�>�x{@��{@��    C�1BT/�6��    6Xv�7�w@�b�3,ۙ7Y�_    >2c�C�<�7�"iA�D�DG��C-�B{A��A
�@��x?�a>�w�                                                >�AD�	@�?1g>��>]�=��=`�</B                                                                    D+��F}��E�RuD�6KC���CS��B�9�A�RK@��                                                @��B͈�B�mA~�@��@*Z>?�.>�7=�g�                                                                    C���E��3D�$C�y�CI�MB�=�BqA+\u@8��                                                @�B<d�A�5�@���@!%a?�n<?�f>Ao�=T��                                                                    7��$7�tAp>A:fA���? ��?�k�+�a*��]- ,�;�+�a*�o'-w�h-[%�-�*0��7���4�M42)�*�C�/���/��l(2ϯ��)�4��2	�u2$t2$r�*˭7�)�1�л��)�4�*Y    0�4~	d4~�-,��    :�Ě7��T7�&�4M!3��/yA�                        1�8A4+34)��0��    4T>�2�F0            4�'�<w�,�Bu?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CL�G��9>*��A�`)                                                B�                                          B��                    A�                @�p�    .�(*}b�4��
6BjJ3p�36��?9P7Ռw?� (7U�^;�Z;���7��            7J�.    3�|�6�    0��#{@��            ��|{@��{@��{@��7�m�8�            7�    7�	S{@��    7�	S    65    {@��0��#{@��0�P0�n�/ڶl    65{@��6�-�6�Ia8{k8{k<�SU    B�W�6��4F6p�DT�@��=    B���            ><]m4��68��k4��6<��L<Ʊ?d�?�R>�*>[F=��/<�͝�)�鶋zi��#��%f��Mና鶊�r�ad!�IiҊ6�&�>;�D<��c<Ƥ$?c�^?2�>��k>=��<�	4�)�鶋zi��#��%f��Mና鶊�r�ad!�IiҊ6�&�:F)-���*EO]6=&�6�;��:�D�:OM9�}o9�DY:�F�                                                �kn��g�C�[�@�L��7I?��[����]n�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         8�        {@��{@��{@��{@��                /"2                /"2            CO�cD�\zD�ekD?zC�"C��Bj?A���@�                                                AπcC]\zCGekB�zB7"A���@�?@��?3                                                                    E��G�G?F���F8�cE�0�D�ĩC��CFp                                                C�?UEQJBEB�EE%�D�/�Cԣ�CdB-�vA=}�                                                                    E�CbF���F�mF{�CF]<cF.�iE�FsE���EG��                                                D)�_D��;D���D��_D��D~�SD0�kC߃�C�F�                                                                    4T>�7���GZ�GN�!GXcFͻlFz��F�E�3oEP?�                                                                                                                                @xpA|��A�ƪB/�B=��Bt�dB�kdB��d<�<�<�<�<�<�<�<�<�<�<�<�E��E�Ec��Eg\D�׫DWB�C�nzC��j                                                {@��{@��{@��{@��{@��D&YB�X.'�P�f��8�?C�4
@��0    AdO#AdO#{@�ξ�$A��$AC��={@��C��LC��MC��={@��@oے{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�A8�g�D�d1GųD�dDs��Ds��@�V�@�V�Fh��BB��    C�FK�xFK�xD� �D� �Fip�BB�d                @*|�C�l+C�VC��?   C��MC��LC��LC���C���C��vC��C���C��C��C��'C��C��C�|vC�QXC�?C��(C��#C�B�C���C���C�y�C�=/C��C��C�#C�!CC���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C�>x>��8>�B�>��>���>���>�_�>�g>ꪁ>�K�>�>�O>�3>��i>�=�>׿0>ҿ>͘J>ȭQ>�Z�>��@T��@�h        A8�� yv� )Q��)�(	�ǖÿǖS�Ǖ���:���9��\��8�s�        >� �A�4�    As�As�{@��@��FX��7��=7��        ;�J?G�?H8�?J:�?L��?P��?W�-?dI!?z��                                                ���6���A�  ?;b�@�oA z�)A�      ��     f9P�n8��L8�*7��6���6�.5� 4&��                                                                    8���7�y87?��6��g5�Y5'�4@�63R�C                                                                    G� �F��iFW��E���E�<D7	dCKBORq                                                4ۈ�46
�3��3�v2J��1���0���/�{q                                                                    4
�C3e�2��l20�1��0���/��.ݩg                                                                    7wR�7_{6�\6gZ�5��5��4N��3���                                                                    7���7�r6kt�6n�5��4ݙu4��36��                                                                    5�F4�^4
�3�"3HM2��F1��Y1V                                                                    7��s75�o6���6j�5�]�5l49��3_�o                                                                    8ɫ�8��x8p*7q��6���6Rj5Id>4f$U                                                                    5�T%5�15,�4��4#e�3���2�2b1�C                                                                    8ɫ�8��x8p*7q��6���6Rj5Id>4f$U                                                                    8�8�7Q�M6��86'/5��*4�&�3ǀ�                                                                    5�v5 ��4��'46�3��L3L�2C{91e��                                                                    8�8�7Q�M6��86'/5��*4�&�3ǀ�                                                                    �f,����u6��5�5e5-�64�/"3�-�3��                                                                    2�?�3Se�: �(9���92��8�0�81�7�X                                                                    5��5O��4�o|4&U	3�v�31R21�D1U�                                                                    75A<7)�R6l/5�5"�4s`m3���2���                                                                    �7�C�*�T����l�?��n�ak�#F\�>-�                                                                    ������c��0�����w|��!�3���}                                                                    64�B6)a�5kd]4ʦm4#��3v�I2�W1�Y�                                                                    �B�N�.����XD�"�ѵ����ܤ�Gͳ>��                                                                    �Ĳ��v��*4ز�	r��&�o#��:���                                                                    ���}����<�j��X�)�ZA㱆�����N                                                                    )��G*l".kP. D-~M�,���-4f�-�@A                                                                    +���,x�	                                                                                            8��8��d8�r7��
6�~�6 �5M�4�)	                                                                    5��5<s4�(J4��3a��2�R�1�HE1p/                                                                    7���7��f7ܩ6��5��H5(��4`�/3��9                                                                    4g�4QJ3�W�3�2i�1��c0��0!�*                                                                    8�z%8��h8�7U�6�3�60ƌ5t��4�Gj                                                                    5&�5	f 4���4�-3q�]2��2
��1/�&                                                                    8B�8
za7���7@J�6��h5���5~W4S4                                                                    6*��6B�5�.�5[� 4��.3��3k�20_                                                                    6�;6�͋6�C6c.5xo�4��#3�!�2�j                                                                    5	ء4��V4�v41��3���2�jq1���1��                                                                    86nK8)@>7���7k6��c6K�5!��4<��                                                                    6P~6Am�6��5�L�4ָq4L39 2W��                                                                    4�{G4�K4�J�4��N4��4H�S4]T3���                                                                    3� c3�b-3���3��f3ei3"IP2�Nm2���                                                                    5��5��4��4�e
4�}�4uuI4 ��3�t                                                                                                                                                                        )A�      ��     f6K�0�                        6��23��<@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ᕠ    @�`    17:00:49        G� @�`    @�@    )B      ��     k�5���7u��B'�fB'�f7��}D"g�    (�+���7%��?�  ?��2�A(3���*�)1����        .=��*1
>m*<>IT�=�O�=a�U;�A�5�ַ.w�26��:[�<'h�<�=G�=�2D=�=��e=��Y=��=�vK>*>x�@��w@rX�@�?��<?�m?Rqp?k�>�1@tk2@u1m@j��@j�>@j�J@j�@j��@j�@k@k$J@k=�@kL�@�6�="�{@��                    E��=7LgG��F���FW KE�#�E�ND6s�CJ^FBN�H                                                A..CyWB�jiA�s�AE 3@���?�7>�/=كi                                                                    EH��@���F&'�A���3ҾpFd��B2�nA ��.3�/                                                    {@��Aq�)Aq�)@oÛ8*'�%�uH!�ḧ��/7ᕤ��?
v�(�@�V4���3�y    >��5@��2T�@j�$� �����        BoeaBoeaC���C���C� ?Z�?ˁ�3PTC�� 6�h6��ZA-1�B9��>��5BS�!A�N'B��*Ae�ACƺB\l�Ac;�BG���g�@M_��ޗA��g�    ��g����>��=��F=��@��@?"[ADƩ?�"2?d�o@|�?�dhF��3Ҿp1S�.��{3��7�0�59�85�O�F�KG�Ge��G�|�>��T    ?y��=���>���>�,�>�U�>ۗ>ͯ^>�L�>�G�>�(�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    @E�n@E�n{@��7�P�9Q-�9Js28��86�7�-�6��/6/�5_U                                                A��3?�
�?�
�{@��@"�kA�ԘBS��@e\�>&}�{@��{@��    C�l�BU�6�E    6�ש6���@���3�-�7�9    >K��C��7�GyA��oDE2�C�B#'�A���A,G�@��?�r�>��%                                                >���AE7~@#?L��>�jY>suy=���=O�<ɉ                                                                    D*�F}�E�5�D�ajC��cCPI'B��NA�R�@�٧                                                @���B�NUB��A�@��e@%��?�C�>��a=�̘                                                                    C���E���D�C�6�CG �B��,B��A&_�@2u�                                                @�zB<�A�7�@�[]@��?��?	�O>9�=K                                                                    7Lg87�s�Am>GDlB��?&�?�]*��)���,���+��r*��)�gh,9b�,��<,j�)��>7~%�3�z1j��+�?%/�/�)����4�2o1b�1a��+�?�7��1�g(��i4�1}    +�\\3��3��. �    <$P�6Y��6�=�2c�=1w&.��                        -���2;N�28�R/2��    2d+�3��            3� ;�O,�`p?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B�E�G���>+�A��_                                                B�                                          B��                    A�                @�p�    -��`)���3���6|3�<�3��>��z7tUh?P�Y6��;�==;�K7:Rp            6��    3z��5��    5Z{@��1���3dA�    ��_k{@��{@��{@��7�r8@            7tv    7��s{@��    7��s    5���    {@��5�N{@��6)��5���5�IG    5�h�{@��6+Wo6�o8	��8	��<��f    B��&6�ʒFB��C��@ݶ>�MXB0�            ?3�9��'3�<�w�<��@+�H?�L�?_k>ŋ�>/Z=�O�)7��Ћz���N��%���N���Њ٠�adJ�Ii��6�&?�<�th<��i@+X�?���?^�@>�/@>$��=�Cm�)7��Ћz���N��%���N���Њ٠�adJ�Ii��6�&<ŀ/���*J&5�L�6�a�<]\;w��:ק�:8��<#g<�#�                                                �,�,�W�,�K�g�;��%��!�õ~�4�x̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 1��}1�(s=��t4���4��+AE8@��b@�ʆA'�=0�3F�3C��B�660w.<�J�>���    =��60w.E�f�294�(CM�@D�oD�5�D<�zC�>�C�OBefKA�KP@�@                                                A͖@C]oCG5�B��zB4>�A��O@�fK@KP?/@                                                                    E���G�dGaF��}F8i�E�!yDȺ9C��CM                                                C�A+EQoEC^E[D�0Cԍ�C��B-�A=�f                                                                    E�DF���F��F{��F]=�F.��E�DxE���EG��                                                D)��D��D��pD��D���D~��D0��C߁EC�D�                                                                    2d+�7R:TGZ7�GN��GD�Fͤ�Fzm'F��E�-�EP;9                                                ?�;�y                                                                        @o\|A���A��(B	>gB@�hBz`B���B�Fo<�<�<�<�<�<�<�<�<�<�<�<�E��kE�'�Ect%ES�D��DDW5�C�d1C���                                                {@��{@��{@��{@��{@��D�BFd-���f��7��jC�s�@q�x    Ak�Ak�{@�ξ�A��AC��{@��C�wBC�{�C��{@��@p/}{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��GrG8�� D�_�G�D�_�Dq^�Dq^�@���@���Fh�QBB??��0C��FK��FK��D��D��Fi8�BB                ?�)�C��eC��C��?Q�C��C���C���C��C��C�`4C���C��C�$�C�]8C���C���C��$C��
C��rC��[C��rC��bC�\�C�#C��\C��0C�eOC�&RC��C�"�C�!KC��;{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C�`u>��>�Y>�.�>��>��$>��v>��%>ψ�>��>���>�tv>�^V>֥�>�7d>�ۇ>��>�G�>˅�>Ǧ?>�@OI@�        AQ���m���67��Aļ���5��5���50��G�u������������dR�        >2*A��    ?ǌ?ǌ{@��@�/FX�C7\Et7�P        ;���??]=?J�q?M?P�?T�9?]l�?k��?}ֳ                                                �7B�6ǶA�  ?(�F@�w?@�5/)B      ��     k�9p�8�n7�7Y�6��U5�'5��4��                                                                    8@�r7��7�66��p5�~�5j4+��3;�\                                                                    G��F���FW KE�#�E�ND6s�CJ^FBN�H                                                4�u�4�R3�@?2�C!2*l<1q»0��/��&                                                                    3ʯ�34M2� O2�1WEY0���/���.�ͩ                                                                    75M7/n�6��`6<-5���5iQ45��3p��                                                                    7DU�6��Z6Rx5샧5k�04�g�4�	3(�                                                                    4��#4�K41o3�G�39.G2���1�*71m                                                                    7o�@776�f6�X5�*�4��P4.��3N9�                                                                    8��t8��	7؝�7D�6�@�6��5/T�4H	                                                                    5p��5�ӡ4��4��c4
P3n�2�i1��I                                                                    8��t8��	7؝�7D�6�@�6��5/T�4H	                                                                    7��7�e�7%�R6���6
15`pQ4��3��?                                                                    4ܐ�4��4[��3��3�!2�FD2'�'1DbM                                                                    7��7�e�7%�R6���6
15`pQ4��3��?                                                                    �s��伽5�G-5�7�5�|4cx|3��3Vy                                                                    0�G_1;��:�ΰ:~u&9��o9`1b8ư�8I"                                                                    5�ڦ5"��4��4��3�R2�;�2'�1E`�                                                                    7۪76;��5�f>5�%4OP�3��)2�	;                                                                    �엶"d���9�B�G�����n��Џ�*Z{                                                                    �¬ܲ�P��Ǐ���XѲ�nݲV̲����tK                                                                    6<�6*B5:��4���4�B3Qѱ2��1��{                                                                    �a��Ͷ�c;��*�jι��D�+\�'(�                                                                    ���,��y���ز�?����%�J�����̯�Y                                                                    �^s{�]&����������=")�ph
��n�                                                                    )Wr�*��/�T".��.N�o-��/�[B02Q�                                                                    +�h,'��                                                                                            8g#�8^W�7؊�7Q�#6��6��52��4d�                                                                    4�,4��4Z��3���3:�b2��/1ʾ91/e                                                                    7bn!7Y��6�!d6P>N5� 5E�4DA'3z�                                                                    3��%3�8z3Vy�2�iu2Ag�1�H�0�NU0�                                                                    8]�o8UH7Ϸ�7N�x6�'6	5U�`4�                                                                     4�,4ף�4R�3��3G�2���1��q1�                                                                    7���7ٱ�7�K7!}6��d5�-�4�j4
<�                                                                    5�A�5���5���54��4��3��o3:�2�a                                                                    6��"6���6qH)5���5QN�4�rr3�;|2�j                                                                    4�:p4�a4��4	�3o5�2���1ڍ 0�T�                                                                    8�x8�7�x,7AE`6�J5�5��4(�                                                                    6�e6
45ЉW5\�4���4|3%G�2A                                                                    4��b4���4��p4�X�4oF�4.�3�3�y�                                                                    3O�3��3��'3l�D3AZ�3��2��R2��                                                                    4�\4ˡ�4���4��b4�94T��4f!3�͑                                                                                                                                                                        )B      ��     k�6�Z�+�\\F�P�    =%�T            6�!�4$%�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�`    @�@    17:00:52        G @�@    @�     )Bq      �     q�5�
7;�VB'�fB'�f7Z3�D"g�    %�+o+�i�6��#?�  ?�!55�2�J�2�­1p��/��I0|�    +��'��j<��=5O�=l�=��u=U��;&��3��H2+�F:	f<'�<�	=G�5=�1�=�q=���=���=�Z=�vV>L>x�?�vv?g�e?+u�>�&�>��>�m�>�R>n��@r95@t{0@i�@j8 @jj%@j�@j��@jł@j��@k�@k"�@k2�ACw�<;��{@��                    E���6�?,G�� F���FW�E��Ek]D6%�CI��BNx                                                A&�C�B�b�A�X�AD� @��?���>Ԃ=�'                                                                    EH�]@���F&=A���3��|Fd�GB3�?���.v�                                                    {@��@9'd@9'd@nqJ7��`%Q��y�ͧ�K/)~��|Q�?f%�+o��|*/S04v�"    >���>�sy2��@���c�S�`�]=�au    BG<tBG<tC�q�C�q�C��!?N�?�S3b	C��6�n*6��A-�	A�)�>���A��;@���A�2@�e@]��A�;@���B%�����@�J��Fx$����    �����ԳL@6� >y�|>y�V@F�8?��@��?LP�?��@C6�?I�yF���3��|1�.@��3��I6���0���59 F�'�F�l�G'�BG9�e>���?���@�}?��?[�X?5r?6�>�U�>ƗZ>�3>�xv>��8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M���    �P#��P#�{@��7?�8Vd�8~��8I!8�[7d��6�w�6�53��                                                A���V���V��{@��B.0A%aB6t�@jJ�=���{@��{@��    C��FBV��6k=�    6^-7
�@��f30vz7y�    >a�8C��k7��wA�vDJ��CUB@p�A�qA?g�@�b�?Ƹr>�λ                                                ?3SAM{�@4B;?t�{?4�>�6=�c�=Ć<#�                                                                    D)ѧF}��E��	D��C�KCO�B��A��d@��{                                                @�\7BΊxB��A�@�9�@%��?��>���=�^�                                                                    C�zE���D��C۔lCF�B�nB��A$9@/��                                                @|UB=d^A���@���@P�?��?��>6O�=F�4                                                                    6�Ց73��A��>JaB�?(O�?���*�V(��5,&!+�|n*�V) d#*�.�,��+�]R)#q�7Z+�/Y2�-&��'�q�            ��|�3;<2 MT-'-&��'�q�7�|������|�4L��    $T�l/N��/N�n*-��    =��ބ2���-�#�,F�)Y�                        '���-aΩ-]�n*�'e    -�- 2�H;            2���:c��,�_�?�b>�?^^?�?d��?\�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B3��G�>,.A�6�                                                B�                                          B��                    A�                @�p�    +�'tZ�/qJ�5N�P3�S�3���=�� 6�Y>��5N�*;`�z;cI6�;�            5���    5$��4��    7R�{@��2?��5�I�5n�S{@��{@��{@��6�78}�            6�v    7'X�{@��    7'X�    5��c5�{@��7�7{@��6z�6"p5��O    5�8B{@�εB'4B�7g�7g�<��0    B�1�6��	FXh�C��?|'p>]�bA�H�            ?n�/Y2�:B�/Y2�<+��<�R@�b@YYp?�x?�9>�
>
S&�)8��ыz���O��%���N���ъ١�adK�Ii��6�&?iw�<+�
<�O�@��@Y*%?��+?0>k��=�8ڋ)8��ыz���O��%���N���ъ١�adK�Ii��6�&<��//�1(��4� �5���;�
Y;=)�:�5*: $�<�Ɏ=B��                                                ��]+��;w�o��Q�*�~��UwÑY̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 3�3�\�=���6*��6��@2I@@F�@A=@��=:\�4�\�4�*�D$g�6J�<<(-@ �=Z)�=P6�6J�<G�<�6|[�6f��CL�RD��qD�!cD;6jC��2C�BbK�A���@�w�                                                ÄRC\�qCG!cB�6jB2�2A��@�K�@��?,w�                                                                    E��G�7G$XF��bF8[tE��DȲ�C��CQc                                                C�;EQmgEC�EvD��Cԁ6C�iB-�LA=��                                                                    E�C�F���F��F{�NF]>SF.��E�CoE���EG�                                                 D)��D��D��D���D��D~�=D0�9C�vC�C�                                                                    -�- 6ɘ�GZ6%GN��G;�F͗Fza�F�]E�)�EP8r                                                A^x<A�:�@�B�                                                                    ?��gA0A���B�B:-�Bq`�B�^B���<�<�<�<�<�<�<�<�<�<�<�<�E��E�+Ecc�EG�D��DW-kC�\�C��c                                                {@��{@��{@��{@��{@��B�]fB}�<+&���e�	7��Cb,'@NN�    @��W@��W{@�ξ��R���RC�/r{@��C�_QC���C�/r{@��@n��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G7�=�*(D�X�Gq�D�XoDq�jDq�j@� �@� �Fg�WB?��B��C�*FK�]FK�]D��D��Fh��B?�                ?Hg�C�ƌC���C�+!?�9C�	YC�]C�]C��+C��C��C�K�C���C��C�S�C���C���C�3DC�udC���C��C�
C��C�!qC��C���C��bC���C�5�C�!C�"�C�!SC���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�C��r>�>�^�>��>��>���>���>���>��0>�t�>��q>��L>��T>�M�>�y�>��>�X,>�m�>�G�>�/�>Œ�@;@Â        @� ��������3����Z�r���Y�պ�դ�č�|ƿ�kƿ�qƿ�k�LՄ        =��A��    ��u���u�{@��@�"6FYi7'��7�        :�'.>�I�>��?.?Gef?Ls�?T>I?`��?rCj                                                �h��6EXkA�  ?M��>�q�?
G)Bq      �     q�8��7��;7FmW7+�K6f�'5��_4خ�3�e                                                                    7?�s6��f6z�6Y �5��4ȕ4�3�                                                                    G�� F���FW�E��Ek]D6%�CI��BNx                                                3��X31\�2��2��1��21''n0d�/}!�                                                                    2���2`	B2��1�q1e-0S$U/��.��C                                                                    64N_6ZzX6��6��5i[4���4hg3@�T                                                                    6P��6�J5��k5У5:�t4���3�),3�                                                                    3�13�n3�a3�>�3a62y�[1�> 0��                                                                     6~�6>kw5를5� 5d!�4��)4�)31	W                                                                    7��7�L70�W7$6o�J5��5
d�4�C                                                                    4|��4�4?�4LA3���3#�v2�i1���                                                                    7��7�L70�W7$6o�J5��5
d�4�C                                                                    6�E\7�i6�|O6uq�5�tA5>b4l�-3���                                                                    3��I4�`3���3�{E39)�2��#2uE1��                                                                    6�E\7�i6�|O6uq�5�tA5>b4l�-3���                                                                    �RJ$���?5>��5i̼4��4j;3���2�(H                                                                    6E6���;.H:���:cml9�)�9/8�>                                                                    4��&4U]3�9�4��3a��2�K�2
1)6!                                                                    6�6(�5��s5���4��4��3Q�Q2pB                                                                    �g��' :��Q��ʴc�C��߹�����
 �                                                                    ���������"����tȲRcĲ�y��.��+�                                                                    5:5'��4��4���3�+a3e,2T��1tt�                                                                    ��}�*�a��昵΃��%����[��N��]i                                                                    ��Oò�n/�ZIͲS�,�� o�
�R�i����u                                                                    �]�e���o�w3�}�]��������@��c��                                                                    (#��)��/m��.�F�.��-�Zv/ �|/�                                                                    *U��+$\(                                                                                            7f�7�n�7/?h7$6v�+5���5��45��                                                                    3�E4�m3�/�3�-�3�82JDy1�p�0�"�                                                                    6ak�6��n6+�C6"��5l64��K4V�3G�                                                                    2��3	V2��F2��2s�1TI�0��e/��                                                                    7\�57��S7(7!�n6� �5�5'�4YP                                                                    3�,4B@3��3��;3	�2^N�1�5�0���                                                                    6�"7�Y6�V;6���68�5��4��3���                                                                    4�Kq5�5�5��4S*j3�V�2�41���                                                                    5�D�5�.v5Ģ�5�ά5O4P[�3�Ä2���                                                                    3�sD3�~>3�3�3*��2n�1��N0���                                                                    7M�7%��7��7��6a�y5��Y4�4�/                                                                    5X�5=oy5)�5.kP4��3��32l�                                                                    3�X�3ϊ<44g3
4*��3��3�xm3�h�                                                                    2~LM2���2�W�3:��3	��2�Ҹ2�}�2U��                                                                    3�P3���4!W)4�I�4P�K4U�3�!h3��w                                                                                                                                                                        )Bq      �     q�6Q($T�lF2	>�9<%��            6���3� @� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�@    @�     17:00:54        G' @�     @��    )e5      �'     wP5C:�7�B'�fB'�f7]��D"g�    $+�p+�6�G[?�  ?��5	4�3g�2F��2p�2Z%c2��M    *]�m&J�J<�h�;z)�=���=Θ�=Qq;U53�}�2/�3:$�<(&-<�A=G�I=�3:=�T=��R=��-=��=�v�>�>y,@%�?�ρ?�H�?;�?�>���>kԹ=֊J@pPw@r�@i1l@i��@jL�@j�A@j�@@j��@j��@k'@k �@k0�{@��{@��{@��                    E�x6��G�8�F���FW��E��@Et�D6B�CI�aBM�g                                                AxCCB��!A��AD��@���?���>�]m=�ށ                                                                    EH��@� rF&<�A��+4Fd�B2�@b�/-K�	                                                    {@�ο�ᗿ��@jb�7�1=�at���6�ٳ�/[xv"��F?�t$+�p�6pH/�&4���    >�`=�$�2��T�|@W =c��=�    B�B�C�:C�:C�9�?3>���3CdC�Q6���6�̻A,��@�n�>�`@���?� �A&�\@�5?�=;@��@:�A�����m�\Z5�
i3��m�/k���m��M�@r�?��?��?�7>��S@Ac>�"�>l�6?�){? �F�x@41��F.�f�4�z6{11    5z.�F�?�F�T�G:�G*�?'n�?�5Aba@l'Q?0ϫ?'?�(>�Ƚ>���>�h�>ȶ�>٭08��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�;�    ��F��F{@��7�T�8�f8��K8uL�8Ծ7p��6�S�5�JF5-,�                                                A��:�-��-�{@��B?�@�i>B+�\@h�r<�:{@��{@��    C��vBW��6+o\    6��7H��@��]3W�?79�    >w��C��S7��=B�~DUC8�Bq-�A��HAY��@��`?ݱi>���                                                ?QEAY#�@Sa�?�{�?�]>�P�>�S='�P<4`�                                                                    D*�5F}�!E���D�jdC�V�CQ��B���A��X@˯�                                                @���B�eeB՛A#��@�P�@(ݑ?�N�>�r�=ʔ                                                                    C��6E���D��C��CG�DB�c,B�-A$�i@/L�                                                @�BB?Z5A�z�@�7�@ h�?��5?�g>7�=Fא                                                                    7�t7�'�A��>IH�B�?'�?��,*!Q8(���,,��+��*!Q8):�)� ",��+�S+))�w7]��                        ���3Cd1�p5            7����s����4[��                        =�����?                                                                2��            2��{@��,�Zz?3�;?.�8?I-?~o�?m�	?_��?k��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A�8%G�M�>,JA�|7                                                B�                                          B��                    A�                @�p�    *5�&��    5�Ü4��4��=��Z6g �>%��4�� ;n�m;q@*7!��            �\Q    5�ls3�#�    7<p{@��1���4���6
䲵;io{@��{@��{@��7�V]7���            6�[�	�B7�.@{@��    7�.@    �%�h6_
�{@��7>�1{@��7r��6�^7d�    ����{@�δ���3Q�k7><�7><�<��A    B�R�6�0Frc�D$�>�2'>�{@��            ?��U    :4+}    <*�<K :@�τ@�J@�?L�W>�Oa>53��)8��ыz���O��%���N���ъ١�adK�Ii��6�&?�:�<*�<J��@�l�@���@�<?K�>��=����)8��ыz���O��%���N���ъ١�adK�Ii��6�&<כ03s�(�Wi5FX5�25<E�f<&�e;YF;n��<�IG=�k�                                                bä7>�W��J
��2�@�v���h.�n�I̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 5M�4�Y(>���6�l�6b��?�w?��?xZ)@!{-=�}6,j�5���Em��7�1�=b
KAZ̉>��=�|7�1�I�*6���7�"�CL#`D���D��D;/%C���C\B`�WA��@��[                                                A�#`C\��CG�B�/%B1��A�\@��W@�?*�[                                                                    E��aG��G�F���F8SZE�DȵC�
�CUN                                                C�0�EQT.ECCE<D�CԀC�B-�A=�[                                                                    E�CF���F��F{�tF]>�F.��E�C�E��5EG�]                                                D)�'D��(D���D��lD���D~��D0�C�~�C�C,                                                                        7�{GZ#�GN��G>SF͏ Fz]�F�kE�'`EP7
                                                A��@�y4@���                                                                    @��7AupA��BvB?:DBx:B���B�]<�<�<�<�<�<�<�<�<�<�<�<�E���E�hEcgkEAD��]DW+BC�YC��                                                {@��{@��{@��{@��{@��CZ��B`F�)��#�fa7�WC7@2j;    @�9�@�9�{@�ξ�����C�wJ{@��C�#�C��C�wJ{@��@jnK{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��GB=�ƐD�S�GS�D�S�Dt�Dt�AoAoFf�
B=��A�P�C���FK��FK��D�DD�DFg�'B=��                >�C��C��8C�f"?T�C�5�C��MC��MC��C��)C���C���C�&C�]ZC���C�ՓC��C�G�C���C�ÿC��C�C.C�{�C��5C��~C��!C��TC���C�D�C�#ZC�"�C�!ZC���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�1C��{>��B>�u�>�b�>�:>�K�>��q>��r>��u>��>��u>��>���>��>��>���>� >åF>�&�>�o�>Ģz@J�&@��        @,�`��F`��D���C��;i���[���S���M�P�3ƿ)"ƿ)&ƿ)"�8$�        =Mt
A��M    �V2��V2�{@��@�mFY�6�I7��        9��?'Tv?@<�?H�9?N�?R��?Z�?h�\?wc�                                                ¡W5�M�A�  ?1�x>���<:�)e5      �'     wP8��8��7��7-�6o�5��4�w3��X                                                                    7���7A+U6�*P6Z�G5�y�4�
�3�n�3��                                                                    G�8�F���FW��E��@Et�D6B�CI�aBM�g                                                4'~�3���3%$�2�2�1�{1(_b0O��/sp                                                                    3S��2�R]2P�i1�%1v0T�`/�F6.�}T                                                                    6�VA6ń86m�g6(;5q&4��4Ȳ38@                                                                    6�_�6�B�6JMO5�h5^I�4�T�3�z%3L+                                                                    4z��4Tq�4-s3���31jE2���1�>.1 M�                                                                    7H�6�56wB6$�5�מ4��'4��3:$4                                                                    8�#8(��7�L7�6{X�5�~�4�=4N_                                                                    4��35?�4���4PA-3�o�3)%2v�E1�{�                                                                    8�#8(��7�L7�6{X�5�~�4�=4N_                                                                    7f�7~1%6�\|6x��5Θ5��4X��3���                                                                    4j4�Q44=�3ǩ^3D�c2�=D1��1��                                                                    7f�7~1%6�\|6x��5Θ5��4X��3���                                                                    �&�ȳ���5�ȹ5d��4ś�4�y3s�2ɦ�                                                                    ��B\5v��;Ms�;8'Q:�j9�ا92�8�J                                                                    5��4�Lb4�6�4�w3��l2ָd2�10Jl                                                                    6��6�@�5�cn5��l4ȼ�453?��2f^�                                                                    ��L��;�lhl�Ĵl�]���g��rH���                                                                    �J�T�gDE��8ֲ���Z����Ұı��G                                                                    5���5�ɟ4��4��3�9�3�2B�1j(I                                                                    ��K��1R�*���n�+X��5��*����
                                                                    �����ղ��T�T ���
���R�_����                                                                    ��OP��
k���_�����}��ñ/!��Z�X                                                                    (���)��:/�^/�%/ӳ.�}/�ib1<                                                                    *	�<*���                                                                                            7�fP7�Q^7�	}7$g�6~�5�|�4��4,�l                                                                    4s�4}4�D3�|O3(2J�1��=0��'                                                                    6�Q6�7�6�9r6#:�5���4ù�4Ew3=��                                                                    3n3w��3�C2�G�2�e1T"10��:/�ݦ                                                                    7�Q7��7�ig7"�6�L5���5��4N`�                                                                    4i'%4rŪ4�B3��3�K2^%q1��70��%                                                                    7c�7u&B7?��6���6@;5��`4�e!3�<�                                                                    5��5��5[>�5��4[r�3�S�2��1��C                                                                    67�6F�6�5�]V5*54Q��3�O"2���                                                                    4RH4bf�41*�3�j�31T�2o��1�5�0�ƭ                                                                    7�%�7��E7jxZ7�b6j��5���4Ҵ~4��                                                                    5��5�7*5���5/��4�V3�J2��G2S                                                                    4$T�4;��4P4i�41mo3�3��\3~t�                                                                    3��3�y3((�3<M�3`02�N2�(_2M�                                                                    4H�d4e^�4~W�4�g�4X�4QF3���3��I                                                                                                                                                                        )e5      �'     wP6r�M    FK�l    <�̎            6��40@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�     @��    17:00:57        