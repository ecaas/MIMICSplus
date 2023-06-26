CDF      
      time       bnds      lndgrid       levsoi        levdcmp       cft       glc_nec    
   ltype      	   natpft        levlak     
   nvegwcs       string_length         levgrnd       hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      Sun Jan  9 16:23:28 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1923.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1923.nc
created on 12/10/21 16:47:41    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:28 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1923.nc had following "history" attribute:
created on 12/10/21 16:47:41
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�fF�h @�@    @�     %m�      h4     ��3V�16�߇B'�fB'�f7@]�D"g�    u�+��_6��I?�  ?��4��1T�2.�0a�/Ń        ��h��8ʝA2�.��954K�;��=Q-N<
��45�:B��<)BZ<��b=G �=�ڥ=�=���=؆�=�v=���>�6> 8�E�2�+S\-�$4s�/���7D&?׶@g��@q��@g~�@g�@hR7@h�)@i+@iW\@i��@i��@i�D@i�[{@��{@��{@��                    E�/`4���G�#VG�5Fk�FE�I�E�D>�kCM��BP�	                                                Av1C%�zB�>-A��AT�@��?Ǻ�>׾K=ڪ�                                                                    E@�@Ě~FS�A�c�2?^�F[a�B(&#<�D/-d!                                                    {@����ӂ��ӂ@Zw{7Vm�&�[�$�O��r�//D%�>�mFu���Mo.�$\2t�    >�y�=���1ۈ�� ����F��`� <<{�    BW��BW��C� �C� �CxV�?�  >���3 u�C�E�60?�6���A �A2��>�y�Ah^�@`��A�R7@pe@�AS��@�U�AϪ��5�m�u�����5    �5� [�A
z?}?}@U��?N=�@� L?=�??$�U@g��?��F�:�2?^�0��-��2W�6r_y    3�t�FY�F��PGμG?�AW
�B�R@��"?y��?O��>���?	��>��?>�z=>�b�>�v�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�_O    �*R��*R�{@��5���3��6        5��6v�>6V�U5V4�B                                                A��y��lj��lj{@��CCu@���B,jn@W�r=1a^{@��{@��    C�78B[�<5��m    6Db�7�P�@�Bb3�6��    >��.C�1G7��Bj�cD��RC�KB�2B@1A��@�$�@q�?�                                                ?���A���A��?��?z/�>�l�>&Rk=^a�<a��                                                                    D&�F�F�D��9D��Cf��B���A� @��`                                                @�pCt#B��AV"@��@@7�m?��H>�1k=̱�                                                                    C�G�E�&@Em��D�ECi��B�`SB
l�A+��@5+|                                                @_�!Bz[�B�A@��@>X�?�`'?!�>;t�=H�	                                                                    5#Jl5�0A}�>5W�A��
?L?��Z)�m�(�ۑ,�f+_�{)�m�(�$�)�d@+�+�p�)�7@]�                        ���3 u�1�b�            7���	m���2�ȳ                        >Wt��š                                                                2{��            2{��{@��,Ӿ�>L��>L��>L��>L��>�y�?2�
?
�(?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A悩G��u>+:A��v                                                B�                                          B��                    A�                @�p�    I�-0 �    3V�12v�2W�=��6Zw�> �|4�=�;T2�;Y[                �:�K    5���0�v�    89/�{@��3�I4��7�_ѵ	_R{@��{@��{@��3�ǿ7��j            82y�    82y�{@��    82y�    �t
7P��{@��7ِ�{@��7�C6}<�6�.?    �es�{@�ζ��3dR7� N7� N<�I    B®W6vj4Fo�dDC��?�I><��A*�?            ?�m*    :�    >��W>�ch@t��@h@�]?iC�>���>YOR��}�͊�:P��⊽/���&���͊`�w�F1��1S� ��	�?}�p>��R>�ch@thQ@`7$@ِ?iC�>}Ǝ=�Sk��}�͊�:P��⊽/���&���͊`�w�F1��1S� ��	�=�    (�A�4 �+�T3;-O=��J8Lu71��=�mK>��                                                �!���D:��fL�@;��w�%�ѡ� ��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6�;4�<$>�O8Av6��u@�3@܊@ug@��=ʓ7��w6/09F�;73�>�Bр�@Z'�>���73�J�ߘ8;�8"̀CPy�D�.�D��(DT��C��C�B]u�A�~�@���                                                A�y�CR.�C@�(B���B<�A��@�u�@~�?&��                                                                    E�*�G	�`G��F�]�F@�"E��~D��AC�� CZQ                                                C��EH�]E<|�E|D�#C�c�C�B.` A=��                                                                    E躞F�b*F��!FzzF\�,F.��E�@�E��AEHd                                                D)A�D�&D��D��YD�u�D~TD0�XCߞ�C��I                                                                        5=�ZGS)�GI�-G
Fљ{F}�F:�E�E EP�U                                                A�@A�ZMA�B�A�@�@K��                                                            ?U�@`�eA��Az�4B��BGX�Bl�-B��@<�<�<�<�<�<�<�<�<�<�<�<�E�]mE�:�Eg�pE<�D�&�DX�C���C��j                                                {@��{@��{@��{@��{@��A�3�B�,~��g��6�U�C@U@�!    @��@��{@�ξ�mȾ�m�C�NP{@��C�C���C�NP{@��@Zw�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G͆>�[D�ռG-�D�ՓD�=D�=A?�A?�F\'B1�0B��Cy�FL4�FL4�D�f�D�f�F]��B1�                ?��C���C��C�>�?$o�C�C�Y�C�Y�C�`�C�neC�C��
C��=C��ZC��C�L�C�v�C��>C���C�
�C�CPC��%C���C���C�"KC�B`C�S;C�V�C�3uC��C�zC��C�d�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�SC��l>��Q>�C>��M>���>���>�iM>��D>�@$>���>��f>�K\>��>��>��>�=>��7>��>���>���>��@J�z@�^        @de���s[��r��qJ��uB����������֙��\��W���W���W���g�        =M��A��2    �@�!�@�!{@��@�{FP-�6��,7�        0�	9���        =���?"�]?$��?' �?*x                                                ��3�5���A�  ?k�;���<#�
%m�      h4     ��2g�x        4�5r��5@�4#(�3��7                                                                    1�W�        3�B04�r�4s��3Ng2���                                                                    G�#VG�5Fk�FE�I�E�D>�kCM��BP�	                                                -��        0>'�0��<0�P/��?/"Lc                                                                    -Y        /p1�0 ��/�;�.���.M,                                                                    0t[g        3��z4i��4L\3M�2��`                                                                    2��        3���4�uA4���3���3�M                                                                    /��z        1�r�2n�.2U[,1_XE0�g�                                                                    2&�        3�X4�r�4�O|3���3#T�                                                                    3$�=        4��5��B5`�4T��3�W�                                                                    0R�        1��.2΃�2�#N1�i�1L�                                                                    3$�=        4��5��B5`�4T��3�W�                                                                    2u��        4..4��4�13���34�r                                                                    /v�G        1d��2H�2>n�1Gp10�O�                                                                    2u��        4..4��4�13���34�r                                                                    �� �        2�8^3��3�s�2�.2pmp                                                                    9�fH9&}�;"p; �:��:��9P'�8�)                                                                    0D+o        1�~�2�q�2�L1�U�1!�                                                                    1�/�        3 w3��3�B�2�u2S(                                                                    �7W        ��0˳h���=#Ӳ,o�����                                                                    ���        �	��M2a��*�ұP^                                                                    0���        2��2�%�2�f1��c1"��                                                                    �?�N        �S=��$���E�xŲ��\                                                                    ���U        ��{װ�C_��r����8�9��                                                                    ���~        � ��������,��&�׃                                                                                                                                                                        +�:�#G�                                                                                            1��        4��M5uς5S�I4G�X3߭}                                                                    .k�        1,�2Q^1�W0�e40}^�                                                                    0���        3��4~[b4^I63[&Z2�B�                                                                    -6�        0*�1�.0��k/�=�/��                                                                    1�Z�        4�u�5�s�5h�#4no[4k�                                                                    .{        1)�>2r�1�J1
�0�!�                                                                    1��        4}}Z5<�*5{4�t3��%                                                                    /)�n        2���3W��3/h�2��1��t                                                                    /��        3L�4��3�i2�%B2fn#                                                                    .	6�        1j^2.N�2��19�0���                                                                    15�$        4��5f��5;��4+�3�CJ                                                                    /O�M        2�
W3��3Vc]2CtQ1�(U                                                                    -��!        1��s3&��3��#3
m�3(�l                                                                    ,�6�        0��2��2b�1߹42i�                                                                    .��        2ԍ3K�Y3���3)0�3NR�                                                                                                                                                                        %m�      h4     ��681�                    >��    6�oX3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�@    @�     16:47:41        FР @�     @�     %n]      hP     � 2�am6���B'�fB'�f7�1D"g�    
x�+��6AhL?�  ?�*3���0��&1~�0�'/���        @��6�U8�Ѻ2�5+H��/o�5��Q=�^<.x4]��:I�<*<�Y=G*T=��=���=���=؇�=릓=���>��> v8��2�+H �+!��+S/�*�7S�>׆�@d�@q[�@g~@g�q@hR@h�@i@iWL@i��@i��@i�=@i�U{@��{@��{@��                    E�Ժ4ij�G��zGf8Fl��E�1HE��D?�CN�BQ�                                                A #C&o�B���A�+�AU�Q@���?��>�N`=���                                                                    E@̦@�znF8�A�I�2+�)F[�B'Q@<�D/-�U=                                                    {@��>y�S>y�S@WL�7$`�%m�"{����/1(��"��>��

x��ʩ.QZw2\�i    >�y�>4]�1�y?�V��3���#e#;�]f    B-�rB-�rC��C��Cm�;?�  >�2�I~CԞ6x�6���ACA�Q>�y�B$�`AG�IBY��A�wA��BeYA<�%B4�j�m<�@&��!<(�m<�    �m<��b&?��?}x�?}x�A��@*;A2F�?Α�@�jA*:�@<FFF�U�2+�)0o�,���2O5�6�#    3C��FsOwF�X�F���F��?<�8@4N�B�Ƒ@���?z�O?O�>��8?4>ӵ�>�@�>���>��r8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    ������{@��5M��            5�5�ϊ6	\�5J�<4��                                                A����9D��9D�{@��Ci�Ap2B<?�@P�!=�Ω{@��{@��    C�tB[��5���    6A�-7�(@�{L3��6�z�    >�^lCd
67� �B|;\D��C�3xB��Ba�%A�z�@��@)e7?P�                                                ?���A�Ah�@&4?�u?��>0��=t�<di�                                                                    D��F���FkcD��,D��CmF.B�M�A�s?@�
�                                                @���C	*�B�c�Aa��@�ݯ@@��?�$`>Ż�=�a�                                                                    C�c�E��Ep�D�hCr�B��B�zA/��@6:O                                                @dw�B{�B!0@�S�@IM�?�gy?��>A�=J�3                                                                    4�"�5;WA`*>0�@A�%?�?�ki)�(4}+���+3�)�(Q�R*+�$+_B�(���7�1                        ����2�I~1�W            6�����2����2��O                        >�n����>                                                                2y5            2y5{@��-v9>L��>L��>L��>L��>L��?�"?�|?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  BqJ.G�C6>,@�A�i                                                B�                                          B��                    A�                @�p�    ����=    2�am2I5�2O5�=�lL6zeg>/�#56�|;
�i;'�                35��    3�b�/���    3<a-{@��3/.5.��4�Q��2:{@��{@��{@�β.y�6��(            4��    4��{@��    4��    5.�6�O{@��5���{@��6�n5�C�6��u    5"bN{@�ε)�3�=�5�5�<�
    B�]�6k�RF[�D�b ?��\?�A�<            ?�:    :"�1    ?Nv>�k�@zQ�@j��@	�=?��'>�ٲ>~��}�͊�:P��⊽/���&���͊`�w�F1��1S� ��	�?�(�?Nv>�k�@z&�@[��@	K/?��`>�F�=�����}�͊�:P��⊽/���&���͊`�w�F1��1S� ��	�=��    %�Ll4)i�;-O>mZe;��:�c=��>:�0                                                �!���D5l�3����X��,3�G��.R�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6�t�5|�>�
8�i�7��@�j�@���@���A#�>^�8"R6�IuF�%�6��C>�O�B�6>�y�>��p6��CJը58�F7%�2CP�aD���D���DUfC�[�CQ7B^)�A��z@��>                                                AЋaCQ��C@��B�fB<[�A�Q7@�)�@�z?&�>                                                                    E�))G	��G}|F�e�F@�E���D��C��C`�                                                C�SEHl�E<YoElD�;�Cބ$CŘB.iA=��                                                                    E蹤F�]F�ަFz]F\�EF.��E�BhE���EHX                                                D)A2D���D��}D��D�uID~U�D0�zCߠ�C���                                                                        4ܳ2GS�GIq�G�FѫF}��F=TE�G�EP��                                                A�A�lA�܉A��A��?@� �<���                                                    >��u@`O�A�Ae��A��5B.wtB`��B���<�<�<�<�<�<�<�<�<�<�<�<�E�D�E�%�Eg�hEI�D�1$DX4C���C���                                                {@��{@��{@��{@��{@��A�G�B�P8�1��g:�6���B�_�@�B    ?Ia�?Ia�{@�ξ+��+�C���{@��C���C�7�C���{@��@WM{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��>��?D�խG*�D�ՄD�[�D�[�AD�PAD�PF[!�B05�B���C`fFL3�FL3�D�fD�fF\�bB00t                ?���C�2C�Q�C��5?iuC��3C��EC��EC�ǷC��C�rC�X�C���C���C��C�#QC�IC�rAC��4C��KC���C�6�C�p�C���C��C�YC�%�C�>�C�5C�iC��C�	C�ӕ{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C���>�)�>��>��>�)>���>���>��7>�^(>���>���>��">��\>�]>��>�"�>�`%>���>���>��>��)@��?�S        @��N�˔��ˑ��ː'�������������q�F���Ȁ�Ȁ�Ȁ�ϻ        =Z�GA��    �=��=�{@��@V��FP6�}7�        7�-p            =m1`>�I�? d<?z?��                                                ���5���A�  ?k�:�<#�
%n]      hP     �             4(�E4�g
4��4x�3���                                                                                3U#�3���4��3@��2�R                                                                    G��zGf8Fl��E�1HE��D?�CN�BQ�                                                            /��~0$��0�P�/���/O=                                                                                .�N/PH8/��.ɴ0.=�[                                                                                3��3�&v4��3@l�2�.D                                                                                3O��3��44d�3�b�2��                                                                                1?�1��2B�1c�B0��Z                                                                                3}�j4
�:4\{3��3�                                                                                4%�4�d�5�94J��3?                                                                                1e��2
�G2|S1��k1>J�                                                                                4%�4�d�5�94J��3?                                                                                3���4
�4q�<3�Qh3(                                                                                0݅v1��>1��1@2�0�Ђ                                                                                3���4
�4q�<3�Qh3(                                                                                2�82�)�36J�2�#�2^o                                                                    9�ֱ9��;f�;~�:��M:�P9BW9u�                                                                                1��2��2[�v1��e1��                                                                                2���3��3`A�2���2�                                                                                ��䲖LV���C� �5����                                                                                ��Ͱ���Z���j�@m�                                                                                1��2;M2d�t1�tY1�1                                                                                ����U{۳��#������                                                                                �H�ү���H7m����, �                                                                                �n����e�E�����^� �                                                                                                                                                                        (�õ`�                                                                                                        4��4�/5��4:�73�Y1                                                                                0���1)�1��0Ӫ(0j�^                                                                                3��3���4Uu3L��2�[                                                                                /�j�0/ס0�.s/��/�đ                                                                                4w�4�q594^��3�\�                                                                                0�F�15�'1�i_0��"0�r                                                                                3��@4s�34ª�4u@3��                                                                                2D$2�b�2�y�2I1���                                                                                2��3E�3�NV2��[2U/�                                                                                0���1aEs1��>0��0s�
                                                                                4�4��4���4r�3�8�                                                                                2$e2�\�3�26:!1�@�                                                                                1Om*2W�31e3	�3&�                                                                                0'�1-�u2Yu1Ќ1�]�                                                                                1}�l2�s�3X��3��3>��                                                                                                                                                                        %n]      hP     � 65��                            6��U3ˌ	@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�     @�     16:47:44        F�� @�     @��    %n�      ho     ��25'6]FpB'�fB'�f7D�.D"g�    Cv+�]r6أ?s��?r4\2��0���0)6I/ԕ�/ZQ        E��= �8���1�TL+?Q+!��+���3/(�:�bX3&W{:N�t<*�.<�-/=G5.=���=�ǚ=���=؉�=�=��+>�A>!8�W�1�ތ+a7+!��+S/K~6��1>��o@`rq@q'�@g|�@g�D@hQ�@h��@i
�@iW<@i��@i��@i�7@i�OAB�;��X{@��                    E�~�3t��G��{G�Fm�9E�XE?KD?�%COEOBQx]                                                A ԣC')�B�mA�A�AV��@�3�?��M>��=�VW                                                                    E@��@�P�F�A�'�1uYFZ�dB&G<�D/.�                                                    {@��@@�@@�@S^�7H��$�~�$x7m&Я/M�<1>�XCv>���-,w2.8�    >�y�@M��1��M�#��A	����PR;�S    BS!�BS!�C�rC�rCtv?�  ?��%2���C�O5�(�6��CA��B��c>�y�B���B�C+APF�A��B�1�A��	B|ǨA`�x��,š(2A`�x    A`�xA�`�Aob?~�r?~�rA��@�&�A��@&u@��bA���@�,�F�t�1uY/0*�,��(2�z7��4EsD2kcF�wG%'GL��GnZ�?P�@�G�CM.@��n?z�t?O�W>���?��>�.D>��E>�[>�$78��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M���    ��1��1{@��4z�4g�;            3�zC4��F4�q'4��q                                                A��s��Q_��Q_{@��C� A�AB`��@E��>[�{@��{@��    C�AVB[��5�-    6>��7���@�A�3��6���    >�vCIE�7��B�@*D�e�D�C��B�KrAݻKAn�@>�=?#q                                                ?�n�A�R�A�?@ *?��?P�>P�=�-�<lן                                                                    D���F�Y�F �vD�e�D[CuU�B��mAҀ�@��                                                @��C	ݛB���Am$i@��@K�W?�Y�>ο�=э                                                                    CŚ�E��Es�D jCz�'B���B�YA5�2@8K/                                                @i�aB}��BP�@�J@T}�?Ƌ�?�>K�=M�>                                                                    3�ߞ4c�A	sK>*�A���?Px?���*���)-o{,%7�+�5*���)g?N,I�,CYO,
�;).��7A�c2�;�0�-~_q            6�Z3Om�1�}B0�{�0�-~_q��Z�{��6�Z1�[�    &��2��2���/lS    >r��6���5�`0߫r/��-�(&                        *Qd�0��
0�V�.�    0߱�2u��            3�h:d�-�>L��>L��>L��>L��>L��>L��>�*?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B��tG��>-�|A��a                                                B�                                          B��                    A�                @�p�    �[�Ϲ�2��625'1>��2�z>�4N7_��?%e6�\;T0 ;X��                5�*�    5"��1�s�    8�H{@��4wu�5@7h<5�K�{@��{@��{@��3y_7|B            7��    7��{@��    7��    ��X�7��r{@��8;��{@��7;5�6Fc7	�R    ��}M{@��4<�55�y46��76��7<��    B�d�6r�F#1�D�g@��?��B�:v            ?�N{2�;�:'�52�;�?�/?�@~��@lm@
v�?��g>�;R>��.��}�͊�:P��⊽/���&���͊`�w�F1��1S� ��	�?���?�?�@~�C@\.I@$�?�#�>�=����}�͊�:P��⊽/���&���͊`�w�F1��1S� ��	�=��    ,��6#f��;->���=�=-�;��]>Z��                                                �!���D������?��đ�Ń��Ŏ�@̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7$��5)�]?+�9��7L�UA��^@�N�AÐ�A�9�>9ί8v�<6�c�G �S7c%=>��&C7i@qQ>��7c%=K�L8�"7��CP�DѺ D�c�DUQIC��PCv�B^w�A�T�@�U�                                                AВCQ� C@c�B�QIB<�PA�v�@�w�@T�?'U�                                                                    E�(�G	�'GeDF�m�F@��E��D��C��UCh                                                C��EHC"E<64E!BD�T�CޫCC�>B.v�A=��                                                                    E�F�W�F��,Fz?F\�SF.�5E�EVE���EH�                                                D)@�D��+D���D���D�t�D~W�D0�Cߣ�C���                                                                    0߱�4BAGR��GIUxGsFѼ�F}��FA�E�K
EP��                                                A�cYA�ϓA���AꀀA�KLA��?��_                                                    >�>b@]4"A�vAO	kA�D�A���BQ�MBC�<�<�<�<�<�<�<�<�<�<�<�<�E�,mE�wEg�AEW+D�<uDXC��%C��                                                {@��{@��{@��{@��{@��A�4LB�F2�Dx�g-�6�9�B���@?    @�D�@�D�{@�ν�����C���{@��C�D�C��C���{@��@S_{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��>ӏ�D��(G*/D���D��gD��gAJ6AAJ6AFZzdB.��CT�CB"FL31FL31D�e�D�e�F\C>B.�q                ?�t�C�5�C�KqC�� ?{HC���C�T�C�T�C�WcC�\�C�c�C�lhC�z�C���C�ȡC���C��C�:>C�dfC��7C���C��$C�4_C�nVC���C��UC��|C�"kC�1�C��C��C�C�Y{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C���>���>���>��K>�6�>���>�J�>��>�j>��>��>��2>��>�=�>��>���>��>�;^>�O�>�>�F�?�L�?��R        A9�d���}��ob����au;�aI�a�ƣ����Ӝ��Ӝ��Ӝ�5(�        =�3eA�%k    ��z��z{@��@�fFO��6H}7�6        8Ay:���            <���>M��?	"h?�
                                                ��l�71�A�  ?k�<��<#�
%n�      ho     ��33�D            2��y3�8�3�W�3z"�                                                                    2cL�            1��c2��2�a2��                                                                    G��{G�Fm�9E�XE?KD?�%COEOBQx]                                                .��+            .MD�/��//@�/��                                                                    -�,�            -���.4:�.]_l.%l/                                                                    1;��            1�V�2��A2�۴2���                                                                    2߇�            2*-L2�(32g�2�'"                                                                    0x=�            0�&0�1��0�A                                                                    3�            2O�z3GC3Z93                                                                    3��            2��L3�{
3�>3��                                                                    0���            07�1��1a�D1(�                                                                    3��            2��L3�{
3�>3��                                                                    3<�w            28��3B
3F)�3�O                                                                    0@5�            /��x0��b0�t�0�G                                                                    3<�w            28��3B
3F)�3�O                                                                    ��Y�            0̕�1��2T�2?�B                                                                    :(c�9��8;�=;�:��Z:!g9=9"��                                                                    1!xn            0Rtr11U�1��                                                                    2b��            1/H2 @Q2,��2 �+                                                                    ���            ��:/���걯���f�                                                                    ��;�            ������4���b�'[*                                                                    1b�N            02�U14F11��1v�                                                                    �E�            ��PϲL�ܲ�1��t                                                                     ���            �_Q���(�* ��	4                                                                    �g�            ������9����!�                                                                                                                                                                        .�� ��:                                                                                            2oS            2�'+3�E3���3���                                                                    .�)            /R��0!�d0h%,0L�b                                                                    1j62            1��2��2�62�-                                                                    -��=            .Y��/)�/~��/`{�                                                                    2eW            2���3�3�}�3ם�                                                                    .��Q            /aK�01��0�x�0t<�                                                                    1�Cz            2�n3W}�3��n3e}                                                                    0o�            0��1vF`1�1�"�                                                                    0�t�            1t;�2."?2f�s29r                                                                    .��b            /���0G�0���0S�                                                                    2~�            2��j3��L3�W13�>;                                                                    0l            0�01��W1�?1�G                                                                    .��f            0�(@1�G�2��3�,                                                                    -�ۣ            /W40��1c�91�{�                                                                    .Ђ�            0��k1��2�d�3%��                                                                                                                                                                        %n�      ho     ��63(&��                        6���3Ȓ`@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�     @��    16:47:46        F� @��    @�#@    %o%      h�     �p4��#6�V�B'�fB'�f7p��D"g�    "���+{�^6/�?�  ?{)�5��3�2�4]1���1m{j0���    (��$��=��H=8��<���6��.�hT4�>E<e��4��1:S�B<+� <�QC=G?/=��=��=���=؋0=�C=��C>��>!x?�ё<���6�q 0Xx�+S/3�*6�o�>~�N@\�>@p�T@g{;@g�@hQ�@h��@i
�@iW+@i��@i��@i�1@i�JA�(;り{@��                    E�=6���G�G{�Fn�E�%E��D@��CO��BQ�d                                                A!v�C'�YB���A�O�AW�@��6?���>�˱=���                                                                    E@q�@�,`F�KA�	�2&LbFZN�B%e�@$޶.��J                                                    {@��AΝ%AΝ%@O#�7��A�5J�#����{���[��
��>�E"���@�g�/�?n4m �    >���A�Rw1���@��-AŇ�@ߕ�>T%T    B��2B��2C���C���C}�?4@@s�m2��Cͺ�5�Z6��A��C�>���C*� B�_�CS��A�x�BTm�C,�B�TB�&�A� ��fgi"4�A� ��+�A� �B�WA��?�V?��A�a!A(�eA��?@MF�@�sHA�<@�F��2&Lb0(�,�]�2	�8?�66��4���G�GA"�G���G���>�vA?ѬA���@���?B=�?N�?�_?؈>��>�g�>�/8>�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�
�    ?�KY?�KY{@��6�=-8��08���7�M�5�I�4/O
5'Z|51��4��                                                A��Ŀ.�*�.�*{@��C<�B ��B�@N@9�>���{@��{@��    C��[BZ�{5�h�    6;��7�K�@�_3�a7��    =e�Cs��7���B�$D���D\�CB�B��aA��:A+�*@Nh?){�                                                ?��<A���A�@1�?�Į?0�>t�=�T�<w�/                                                                    D���F��F"�D�o|D�C}�CB�;�A��v@׽�                                                @��CC
FhB��Awۂ@�|@W|�?�`	>ָZ=Ր0                                                                    Cǟ�E���EuCYD%�bC���B�E[B��A;D�@:�                                                @n��B~w\B
.�@�*S@_�e?�/�? �>R��=Q�                                                                    6���7čAK>&8*A�[t?
n`?�'+ 7\)�By,v�+�/0+ 7\* ,Q,҃;,�i�,r�)���7M�4o��2/~�/�XU            7̪+4�;1�(�25��2/~�/�XU�̪+��|�7̪+����    )���4dK�4^��1��    =a�N8��7��D2�#O1�C0%��                        ,ͭ�2�nG2���1+�=    2�0*2uK�            4~��:��a-	E?�{>�$>�V�>L��>L��>L��?-?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CK?PG��>+�
A�(                                                B�                                          B��                    A�                @�p�    (c�}$E;4�.�5���2l_,2\�?U��7��@?�*7�X;\�;W��    /�2�        7,(�    4�g�5B��    8`�8{@��1;��5���8 �� s{@��{@��{@��9i��7��            8CŞ�8C�
{@��    8C�
    5�o�6�R�{@��8x!�{@��77}26�5�6�Ĝ    6�{@��5���6ڎ�6�L?6�L?<��    B�C�6�HFb�C��Aeh�@5j�C�m            ?\ô4o��9�B4o��>��>� �@6��@'�:?�ƴ?`�q>���>�d��|�̊�:N�����/���&ߊ�̊`�v�F1��1R� ��	�?E�E>�>��z@6��@}?��??P�>��=�����|�̊�:N�����/���&ߊ�̊`�v�F1��1R� ��	�=�3x    -���8�|7y-:�p�>g+�=yΗ=�u%6��>g�J                                                {�±��ïaC�
bG���e��.}��V��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6�W�5܆�>��G8�ώ8D!�B �,A�gB$\�A�.�>D�7�P�7*�F ^�7e��=s��A�^?��=�ni7e��J���8<G|7��hCP��Dҏ:D���DU��C��6C�,B^ەA�@��#                                                A���CR�:C@��B���B<�6A��,@�ە@?'�#                                                                    E�+QG	�GPVF�v	FA�E�2 D� hC��Cp�                                                C�vEH%�E<�E'LD�oC�ԋC�B.�aA=�T                                                                    E踾F�S�F��FzMF\�XF.ݔE�H�E��IEHb                                                D)@�D��<D��D��OD�s�D~Y�D0��Cߦ<C��0                                                                    2�0*6HΆGR�~GI@(G\F��tF}�FG2E�N�EP�V                                                A�SA_"�A@I3A���A��-AL��                                                        @��kA���A���Ax�A���B6�BHI�Bs��<�<�<�<�<�<�<�<�<�<�<�<�E�#E�jEg�|Ed�D�HLDX&�C���C���                                                {@��{@��{@��{@��{@��C#��B��(��f�r7>�:B�?Q@(    @��@��{@�ξ���C��V{@��C��PC�n�C��V{@��@O#�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��>��D���G=�D�֧D���D���AO�AO�FZ�$B.;�Bʹ�Cg�FL5�FL5�D�hSD�hSF\ZB.6�                @B�C��KC�qrC��? �VC��]C���C���C���C��zC��KC��)C��?C��;C��C��[C��C�,�C�N�C�tC��.C���C�ZC�<	C�r�C���C��IC�C�)jC��C��C�.C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C�t>�a�>�w/>�@�>��#>� j>�kK>�ƶ>���>�!z>��2>��>��h>���>��>��]>�a>�m�>�|�>�]>�ߤ@^�?��4        A��j�Z
�Y�-�Yw��m�J��K�����?��*������������T}        >�{B/G    ��ٺ��ٺ{@��@E��FO�6o�7^}        9��?!+�?	�H>�z�>�N=*�>�K�>�ء?��                                                °h7�NA�  ?k�<��<�u�%o%      h�     �p8���7��R6���4���3E�4��4 x�3g<�                                                                    7�B6���6�Y3�l2I/q31�33"G�2��                                                                    G�G{�Fn�E�%E��D@��CO��BQ�d                                                4H!�3J��2�QS0O_J.�נ/�P�/���.� E                                                                    3|�2�.>1�YA/���-ҿ�.�*.��.��                                                                    6��6\2�5���3�v�2�J3��3!��2�߯                                                                    7��6�B�6+��4'�2���3�{�3�]=2��%                                                                    4�y�40�E3�Rq1�\�0t��1b��1kAK0��                                                                    7;��6���6Q�l4M1�2�k�3�^]3���3�                                                                    8,��7ç7�84Θ.39/F4/ے43
3���                                                                    5&4�`�4-2�p0��g1�c�1�1��                                                                    8,��7ç7�84Θ.39/F4/ے43
3���                                                                    7�6�7�6g4�4)5�2��
3��23�V-3
@.                                                                    4���4&��3��31���0�,1�@1.��0�)�                                                                    7�6�7�6g4�4)5�2��
3��23�V-3
@.                                                                    �,�����2��2m�L1ca2$,�2_
#2/*                                                                    9_��9G�:��w:��K:ckA: �95Ķ9(��                                                                    5]��4��4bԾ2V�0�^�1��1��0��                                                                    6��6-; 5y��3(��1��l2���2���1�,�                                                                    ��	�'�ȴ��ֲ�/��u����T��vr                                                                    �^�B��u���۰2կ���x�A���_�                                                                    5��5-�84}ѫ2,��0��>1��1�s�0�                                                                    ��ζ,��)d�d���0ò��5��I�uZ                                                                    ��v���8�����v�d��d�����@�
ǐ                                                                    ������J�����
���xy��a&߰]큯�k                                                                                                                                                                        0r G.��K                                                                                            7���7��6�K4�3�3��4��4��3�-A                                                                    4��l45t3b��19��/�0���0���0=^N                                                                    6���6�ъ5۹3��2%F3!V/3,�2�N�                                                                    3|�
3
T�2^'08_ .�l�/��o/���/O�	                                                                    7���7��6�'4��k3*��4(�4<Q3�p_                                                                    4wM=4tG3Y�$17K/�Z*0�P0���0a��                                                                    7q��7!r6�Qs4�z2���3ݢC3�u�3S�	                                                                    5�:�5��4�]2��;1�v1�K�1��21q�x                                                                    6Cz5�V5b��3\��1Ŧ�2�2���2+>                                                                    4_f�3�pc3��51|/��0̮�0���0C��                                                                    7��E7&a�6��4��3y84qb4��3�d�                                                                    5��5>&�4���2��d1*ӭ2ʹ2�I1���                                                                    44ϔ3�#(3��1���0�\�2Iʻ2��2��r                                                                    32�
N2m$0�[�/��1#n1��1�u                                                                    4\��4ܘ3�Vo2��1�'2v�s3��3�                                                                                                                                                                        %o%      h�     �p60
�)���                        6�Lk3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��    @�#@    16:47:49        F�X @�#@    @�+     %o�      h�     �@5���6�d,B'�fB'�f7��GD"g�    (&o,+F8�6��P?�  ?y��4�u�3�~.2�<�/�Ġ1�2U.���    ._S�*W��>���>G��=�.�=�L;4�7E��:���3�$:V6<,%(<�pa=GG:=��X=���=��	=،R=�5=��>� >!�@���@Dt�?��?�&.?a��?3�B?�>�F�@gY�@lwk@b�!@cm*@d�@dĸ@e�6@fCB@f��@g�@h^@hU�?��)<�>���                    E���7h�]G��vGMFnѲE�t�EVHDA:CPX�BR;?                                                A!UfC'YyB��A�8�AXJ�@�V?�a�>�09=� 1                                                                    E@SW@�NF�xA��U2r�cFZNB%�AV/�@r                                                    {@��Bg��Bg��@L��8d��%��    ��}/m@�&���>���(&o,Au��/��P4�&�    >� cB
G1Ĕ9A$�hAs�Aie=4=    B�d,B�d,C���C���C���>�a�@�j�2�Cˎ�6b�(6��gA+�CE'�>� cCY�B���Cm��A�) B���CH;�B6rB�VB[Zп�- ^�;B[ZЧ(I�B[Z�B`CA?3$�:���;�XA��A6�A�[e@`"�@Q��A	��@vDF���2r�c/�,�kZ1��8�e�6�E�5RG$ֱG_�ZG�qG���>,��?BG<W�\    >� �>���>�!>���>��;>�Q)>���>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���MŎ[    @�P;@�P;{@��7�09t�9W 88���8;��7�1�6ƺ56�5 ��                                                A�}�@��@��{@��@�UAB`��B���@/X�>���{@��{@��    C�oBY��6%p    69�m7�s@��b2�*j7&�    >$lC�+�7��3B���D�G�DB�M�B��A�nA,�(@I�7?&)g                                                ?�ʃA��A�`@��?��-?/�M>v�j=��<t>(                                                                    D��F�<F!x)DϓoD��C�c�B��jA�p�@�9�                                                @�^�C	~�B�S�Ar �@��s@[�<?��0>�=^=�u�                                                                    C�0E�y1Et�D#(�C�z�B�E{B&�A<��@;bF                                                @m"�B|��B�f@�]|@a>�?�\4?$mN>U�=R�]                                                                    7q�'7��AX->$|�A�sR?�B?�@�+�&M*���,�j�,Y��+�&M*�� -�p<-F��,��~*9�7}�	4�bM2I�P0�D            7��d4��2O�2\H�2I�P0�D���d�R��7���4��    2��4��4�k2��    =c�8BQX8�6J�O5(�2X�                        4��S6)6$�%3���3(O�7��2�y�            4��>�V-!��?}.�?~�?v&Q?^�c?RN-?R�?bګ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C`=FG�L4>+d�A�i�                                                B�                                          B��                    A�                @�p�    .�)�;_5�G6�t�2�@2ܓ?��p8;a�@A=�7�Ӣ;�);;�6s�4/VЍ        7�K�    3)Y�6�)�    6��{@��    2 �	5 ����ߺ{@��{@��{@��9r��8��            7D�3�7QA�{@��    7QA�    6�',    {@��6��{@��5�ܡ5H��4�'6    6�3g{@��6�*S7Xh)8	1v8	1v<�A/1�o�B�]6Ɠ�E�j�C��|A��    C'o�            =��4��e8��-4��e<��<��z>���>���>?��>�=���=�q��V������M���������t�����bd�G�1�>� �Q���=�,<��<��]>��>�k�>3VN=�Y=��<��T�V������M���������t�����bd�G�1�>� �Q���<���-��)Y܂5�Q36��9et�;�9<CS�<F	:�9\=~�W                                                �F��F\:�:O��.��� d�w�H��g��� �̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         <���        {@��{@��{@��{@��                5�l27��        ;��5�l2            CR��D��D�9�DY��C�9�C1>B`�A�_�@�o�                                                A���CS�CA9�Bٓ�B?9�A�1>@��@_�?(o�                                                                    E�:�G	�G[�F���FA+�E�DuD�1dC�Cx�                                                C�&�EHA�E<(E/D�|�C��5C�`B.��A=�                                                                    E�2F�V�F���Fz"F\�0F.�yE�J�E��PEHx                                                D)A�D��D�ŋD���D�s�D~Z�D0�aCߧ�C���                                                                    7��7N+9GR�-GIQLG*�F��gF}�xFL[E�Q�EP�>                                                =! ~:L͓    ?�Ǥ@�W�?n�
                                                        @��	A���A�~�Bm�B0�ABf-B��|B��@<�<�<�<�<�<�<�<�<�<�<�<�E�6"E�BEg�EuDD�VDX/�C��LC��k                                                {@��{@��{@��{@��{@��D. B1à-��l�gB�8�
C ��@/�Q    A)ݰA)ݰ{@�ξ�A��AC��Q{@��C�sC��C��Q{@��@L��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G �=�6D�ֹG[�D�֑D��=D��=AL�BAL�BF[��B0=@��LC��FLE�FLE�D�u�D�u�F\��B0!                @�F6C��eC��*C�A?   C�0xC��C��C�|C�C��KC�9^C��C���C��)C��&C�z�C�u�C�y�C���C��hC���C��IC�C�McC��C��XC���C��C��C�C�AC�f�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�>�C���>�3>��v>���>��.>�I�>���>���>�H6>��v>�1Z>��>��>��>��|>���>��P>���>�c�>�3B>�͖@DPA@A�        A� ��v8	�uv]�t�Ū�-��g��^��ے�E.KƽL�ƽL�ƽL��da�        >��QBEm�    @�) @�) {@��@�S�FOv^6���7FV        ;�<8?Rp\?SY?T��?T� ?<��?Fٱ?K�
?Y�k                                                �k��7aT�A�  ?U�@ƂB@�p%o�      h�     �@9I\�8�c8�7G�q6�M5��N4�$�3�+B                                                                    8~ZW7��7'd�6|�f5�y�4�y�3���3�                                                                    G��vGMFnѲE�t�EVHDA:CPX�BR;?                                                4���4*'�3�ض2рY2�1.�0E�/\a                                                                    4A�3V��2�b�2Q1)%�0\r�/x�l.��=                                                                    7S^�797�6��6%�u5u3m4�;�3�Q3"�r                                                                    7���7VK�7�6�W5�ю5&��4L�~3;�d                                                                    5 [�5X�4˩�4Y�3��3�2&r+1h}                                                                    7���7��X77u�6���6��5K��4z�3e��                                                                    8� �8�ާ8:7M�6���5��-5��4��                                                                    5�=�5��[5.��4�u�3�'�3A�_2�l1��#                                                                    8� �8�ާ8:7M�6���5��-5��4��                                                                    7�ڞ7��7k�36���5�~�53N@4d�d3{D�                                                                    4�+`4��r4�a4��3w� 2���2'141                                                                    7�ڞ7��7k�36���5�~�53N@4d�d3{D�                                                                    �HV0�u25&5�Q4q��3�y#3(%�2���                                                                    4���5��9J�9%$G8��8�';8)�#8Fʫ                                                                    5���5�`5F`H4��]4�3H�W2r��1^,K                                                                    7�07��6�7�5��74��4!�3E�W2Z!^                                                                    ���������)��q�䳡w.���3���                                                                    ���ݰ:��iF����T�k����ձ�$                                                                    6�6��5���4���3�\!3%�2Ke'1`GZ                                                                    �&!�÷������|Ե-�e�|���m��P�                                                                    ���ϳ�]��v�oHb���C���A\<�|�J                                                                    ���	�ho
�G����ײ�"@*�CG*                                                                    (�v�(魕+G%�+I6U*�l+*���,���-�{h                                                                    +6�+y_@                                                                                            8��f8k7��79��6���5�ĥ4��4\                                                                     5W'4�i4w�3�$�3
]�2G,D1��20���                                                                    7��7f946��686d5�'4��&3�f�3'�                                                                    4��3���3rx2���2-@1Q
"0���/�<�                                                                    8�Z�8aoO7��"76��6���5���5
�245�                                                                    5��4��4ml�3�l�3��2Z� 1�a�0���                                                                    7�Б7��N7���7	�i6C�o5�|4��*3�:                                                                    6.
6�v5�p�5�
4_�3��D2��/1۰                                                                    6θ6���6s�/5���5`14S�3�"2�U�                                                                    4�@4�M�4�X�3���35 82r_1��&0���                                                                    8T�8��7�j�7(��6o�5�$^4�533��                                                                    62�6�W5���5@��4���3��2�:2@�                                                                    4��4��
4�4r��4,�3��3���3cgz                                                                    3�q�3�^�3~�.3DV[3d�2���2!�27��                                                                    4��4��4��U4�z�4R�\4Ò3��3��                                                                                                                                                                         %o�      h�     �@6.0"2��                        6��3��y@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�#@    @�+     16:47:51        Fє @�+     @�2�    %o�      h�     ��5���7g�B'�fB'�f8$%�D"g�    )'x+7��7Xv?�  ?tK{    4	�)$�/�e�����        /`�B+X�>��r>4kA=�!�<ؓg9�i_2��p+��X2x�b:OL<,`'<�<=GI�=��<=���=���=؍=��=��>�n>"@ػC@��<@ n�?�VB?�%?���?���?kM_@lj&@l�E@b\�@b��@b�W@c2(@c�V@c�@d;@d��@dɱ@d�M=F�%=R��>�ֺ                    E�t�7u�G��GόFn�E��E�%D@�(CO�hBQ�                                                A �;C&��B�A�A�l�AW��@��X?��(>ټ�=��                                                                    E@]�@�)EFޮA��0�>�FZ(uB%��A?U�/�{�                                                    {@��B��TB��T@N��8�w	���|    ��T�/_���D�z?   )'x@�}Y��h�212    >�y�B; k1��kAM�K@�&�@�E�        B��B��C��bC��bC�N�>�P@�;@2�Q-C�6��6��/A�CM��>�y�Cb8B�>%Cl��A���B���CIB�B�'rBv�F@ߒ    Bv�F    Bv�FBl�M:��R    5:g/A�x�A?�KA���@t��@*z�@ԍm?�^�F���0�>�.ϸ+h#)0���8��.6���5�b�F�/G:!G�(zG�p�=lM�            > x�>�|�>��k>�P�>���>�ʒ>��>�e|8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�P2    A'GA'G{@��7��9to9W�48�l?8X�27���7 �6 +�55�E                                                A���A�jA�j{@��8���B���B��y@-F?�~{@��{@��    C��BV�:6�
�    69��7��@�:&2�$7689    =k��C�_7��Bq��D�-�C�3ZB�MEBO/A���A�\@ J�?h�                                                ?��sA�|A�H?�#?���?�`>?��=kws<O2�                                                                    D�^�F�UCF !uD�xD�yCxL�B���AҎx@�h�                                                @���CALB��PAd*�@�:�@O�?���>�<w=��                                                                    Cĳ�E�\�Eq��D�Cz�NBϾ�BnA5�@5�'                                                @f�Bz(=BtQ@إ:@S��?�Hx?5S>K~R=JN                                                                    7xp7��yA�4>'�A�2�?"�?�׿,2W�*�-4,��j,2W�+��.��-���->��*<�7��4|��2»,��%            7��!4��r1��72�R2»,��%���!2127���4oS    3�v�4s�4r��.�Ug    =B,�8C(�8,�o6�i�5�N�2�                        5��/6�Ӵ6��3�1X3�4�7��\2̐�            4�>��-64?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Ce�AG���>*�;A��                                                B�                                          B��                    A�                @�p�    /
V�*���5G�6��1�Q�0��[?�ov8J�;@zl#8� ;��J;��(                7�J>    0��16�w�    1���{@��            ��[+{@��{@��{@�ε�$�7 �=            5ƴ�    5ƴ�{@��    5ƴ�    6�m�    {@��1���{@��1�T`1f�i/��    6�m�{@��62-�7��7	�x7	�x<��c2 �B�6��E�zC��A�]    C0�            ;��4�H
6ic4�H
<Ňu<��<��H<��,< �;~I�:�L�:qx�*���U��	{���U���6�������U�a�	�F㎊1�� ����;Ԑ�<Ń	<��<�jQ<�t<=#;|9B:���:dK�*���U��	{���U���6�������U�a�	�F㎊1�� ����8&�g    (��+6�6b)8��A8���8C��8)�8>�8U�4                                                ��Z<��8���|���^HĮ�ĕ���q�O�5v̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         8���        {@��{@��{@��{@��                /��8                /��8            CTv�D�#�D�XVD[bEC�CwOBd�oA�
�@���                                                A�v�CS#�CAXVB�bEBB��A�wO@�o@
�?(��                                                                    E�F�G	��Gs�F���FA#KE�?vD�.C�C|�                                                C�8�EHkIE<K�E4�D�v�C���C��B.��A=�l                                                                    E�}F�[�F��Fz�F\�bF.�DE�I�E��gEHq                                                D)B�D���D���D��)D�t�D~Z�D0�CߦhC��;                                                                    7��\7b��GS�GIl�G77F��F}��FO.E�R�EP��                                                                                                                                @H�:AZ�A�&EA�$�B,KzB^=7B��B�$<�<�<�<�<�<�<�<�<�<�<�<�E�M7E�#vEh mE�D�a"DX5�C��QC���                                                {@��{@��{@��{@��{@��D!ݬB&.�z��gO�8;C["a@B��    ANm�ANm�{@�ξ�HȾ�H�C���{@��C�
YC�
\C���{@��@N��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��GQ><�1�D��uGy^D��LD�%D�%AE�XAE�XF\�B1�V    C�@FLS
FLS
D��CD��CF]��B1�:                @��XC��"C�c=C��L?   C�C�
YC�
YC��!C�fwC�(�C��C��0C�E�C��SC���C�V�C��C���C���C�n+C�J�C�8�C�:�C�N�C�n�C���C�� C�C��C�?C�UC��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�	�C�Z�>�aK>��>��z>�8p>ؖ�>�dk>���>��">���>�ڛ>�ײ>�'b>���>�E�>�E�>�E�>�`>�v�>�B9>�g�@N�@`        B��ǝ��ǜ�ǜ����:��C�����8Ƅv
�Ü��Ý�Ü��ԫ�        ?!EB>�    @��6@��6{@��@��`FO��7Q\�7O�        ;���?0��?1��?2�G?5�?9�;??��?G�X?S�_                                                � �7�!A�  ?k�A	��A z�%o�      h�     ��9H�z8�v"8k7p�<6��b5�;4��T3�YA                                                                    8}�/7�zF72D6�/�5�fa5Fv4μ3�                                                                    G��GόFn�E��E�%D@�(CO�hBQ�                                                4�N�4+GN3��J2�p22��1h��0z�/}��                                                                    4�J3XZ2���2o91a��0��y/���.�C]                                                                    7S�%7;9�6�A6Jb�5��4��<4��3<4v                                                                    7�ż7Ih�6�c6�Ƨ5�w750ry4R��39K�                                                                    5�4��	4��04D�3��3&62+��1�s                                                                    7��W7v*�7<�6���6�i5W�[4��,3bx�                                                                    8�tP8�9X8
�7mѴ6�p65"s4#3�                                                                    5�05���5*^W4��4 ��3w
2���1��\                                                                    8�tP8�9X8
�7mѴ6�p65"s4#3�                                                                    7�
7�W7js�6��6 ��5g&�4�f�3���                                                                    4���4���4��X4"�Q3���2�2�1�                                                                    7�
7�W7js�6��6 ��5g&�4�f�3���                                                                    �)&v3^��5�^5g�a4�f4>�3qM�2�?                                                                    5(x�5�\�7?�7 `�6�q�5�T�5tH�5��                                                                    5���5�v�5��4��j453T��2y:Y1Z��                                                                    7NN7h6���5â�5�r4R]�3t�'2vw�                                                                    ����2��m�M@��qq�ר߲�H$�˃                                                                    �믓�����υ����첎��DNͱ�`���A                                                                    66�5�F4�lx4|�3Wb2{iw1|�1                                                                    �&B��T���8ȶ9Ƶl�y���$��H��                                                                    �������$���б��^�9԰xd\��@�                                                                    ���j���#� ��3 �,�1�ñNR@�a�                                                                                                                                                                        *HX
*���                                                                                            8���8m�$8�7c�=6�E�5���5��40-�                                                                    5r�4�.�4�!�3�F�3<VX2���1�x0ǐ�                                                                    7�4�7h��7 �@6a��5�f�5�4$�3A-�                                                                    4��3�I�3�km2盥2B�1���0��//��N                                                                    8�t�8c��7���7`_6��}6��52{)4R-�                                                                    5�4�d�4i�3��3Io 2�I�1�,\0��                                                                    8 �7�M�7�x7&�46��5� �4��3��                                                                    6_�6,c5�@5>Ya4���3� �2��1��9                                                                    6��/6��p6�)�6�5S��4�j�3��2�g;                                                                    4�64՝4���4�F3q��2��T1���0�C                                                                    8��8Z#7���7K�]6��5��W4�:�4�x                                                                    62�6!��5���5h�!4��}3�r�3،2?                                                                    4�@�4���4�f4�ů4f�4L�3Ƚ23�M�                                                                    3��P3��3�f3m5%3:``3 ��2�6�2T5�                                                                    4��4�$�4��4�cd4��|4B��3�Y3�{�                                                                                                                                                                        %o�      h�     ��6.h�3�v�            G���    =�:}6�3Ŋ$@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�+     @�2�    16:47:54        F�� @�2�    @�:@    %pQ      h�     ��5��7�eB'�fB'�f8+�D"g�    )ٕ�+3T�7P�K?�  ?td�    4n)6��/�%f� �        0��,�>��!>'%q=��<�G�:b�4��,��62z�q:@�+<,?M<��=GF]=��=���=��q=؍=�"=��>��>"R@��@�k8@>�\@��?�͂?ʒ?�0�?�=�@l,H@l6@a�)@bU@b,@bM@bs@@b�;@bȒ@b�{@c�@c&�=F��=���>�FC                    E���7j�&G���G]�Fm`�E�YoE�vD@f�CO��BQ��                                                A E�C&�B��$A��9AW V@���?ɂ�>�S�=ۍ                                                                    E@��@�FDF��A�!Z0,��FZ�oB&��AL%*.���                                                    {@��B�CB�C@Vfn8�7|%�Fc    ����/L��#�&?   )ٕ�ARF���28֦    >�y�BJZ�2 9AVf�A,A��        B�1�B�1�C�߫C�߫C��A>���@��2���C�w�7^�6���A�CRü>�y�Cg�8B���Cow�A��/B���CJv�B"�B�5�BPZ��B�u    BPZ�    BPZ�BSd�            A���AL�cA�6�@|�@%��@��?ߥ0F��80,��.��+��0��f8��Z6�	�5ʨ�GNG(�UH��G���>��            =��>��>�ƅ>��>��>��>�3�>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�4    A�J7A�J7{@��7ު�9r��9M�8�y8D�7���6�� 6�5P                                                A��^A�ՙA�ՙ{@�Π2�oB���B�s�@4��?�{@��{@��    C��OBV�!6�}�    6<Z�7�[�@��82�/�7c��    >��C�q�7�
�B\�gD�.C�@�B�
B"�)A�E8@��@��? ZR                                                ?nّA�ҝ@�HV?��P?Q�&>�7G>&�=G�.<=�S                                                                    D�]	F��hF�'DÐ�DU&Cp��B� uA��@�Vy                                                @�U:C�B�agAX�c@ϩ�@E!�?���>��;=�A�                                                                    C�E�8�Eo^�D��Crw�B�`BA0>	@1{�                                                @a"Bw�FB�@�5@H
~?��?X>C2=D{                                                                    7f>�7�E5A��>+rA��?��?��s- �+�:-���-&�y- �+�\M.���./3�-���*�v*7�/4�S2��-�A1.��.3G)�~�7�V4ы31�a�2V�2��-�����V28֦7�{4u��    3��4zLn4y��0
s    <ꊞ80�08)\�6��_5�)�3ki�                        5�m-6��6���4�2�3��*7��2��$            4�Ƹ>ti-��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Cp�G���>+�#A�                                                B�                                          B��                    A�                @�p�    /��?+���5A��6|R�1]	1 ?Ƒ$8E��@z�%8�<4�;���                8Z8    11�6��        {@��            ��i{@��{@��{@��7n�7�y�            6��    6��{@��    6��    6���    {@��    {@��                6���{@��6�RE7���7΢}7΢}<���2B��b7ČE�~D�EAæU    C:N�            ;�1_4��6n�4��<�I�<��r<�w�<��~<g�;}�*:�0F:aTp�ޡ�B͊�3��jኽxo��b"��B͊aP�F~,�1`h� R7�B;�x�<�D�<�ϋ<Ը�<���<	t�;t��:���:O�f�ޡ�B͊�3��jኽxo��b"��B͊aP�F~,�1`h� R7�B8�(p    )эL6�6|�F9_[9��t9<�A9)s8�E�8�hR                                                Ŏ��ŊZ�ň���r��i;�N��16�cQ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �2�o        {@��{@��{@��{@��                                                    CT�wD�YD�LrD[9tCÇ�C��Be��A�W�@�=�                                                AԽwCSYCALrB�9tBC��A���@���@W�?(=�                                                                    E�Q�G	��G�TF��zFA7E�:GD�)�C��C�                                                C�H\EH��E<lE9D�p�C��hC��B.��A=��                                                                    E輼F�`mF��Fz�F\��F.�E�H�E���EH�                                                D)CqD� �D�хD��(D�u|D~ZLD0��Cߥ5C���                                                                    7��7W�GS-�GI��G:�F��F}��FOSE�REP��                                                                                                                                @f�A/mA�*A���B3�B?&�Bf yB�8�<�<�<�<�<�<�<�<�<�<�<�<�E�bE�4�Eh�E�OD�c�DX6�C���C��d                                                {@��{@��{@��{@��{@��D'��A�q�/i���f��8ƕC�$@QQ    A�WA�W{@�ξcFվcF�C�T^{@��C�r�C�r�C�T^{@��@Vfn{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G3<e"�D���G��D��D���D���A?��A?��F^(�B3Ќ    C�K�FLZ4FLZ4D���D���F^��B3�p                @���C���C�YC�t�?   C���C�r�C�r�C���C�MC���C��C�CtC���C�e�C�]C���C�Y�C�iC���C�f'C��C�ݕC���C��]C���C��LC��*C��C��C�VC�iC���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C�o?}?��?�{? \">���>��>�{>᝛>ڌ}>��>�7P>ɼ>ġ�>��J>�I_>���>��>���>�V~>� @'��?��         A��oǧ�9ǧA�Ǧ�q�*��!�� Va��Iƛ5��d���d���d�ń2;        ?>"WB$A�    AB��AB��{@��@pQ~FO��7b	�7qb        ;��?�(?6&?�[?W�?re?�?�?%�t                                                ��EY6�3A�  ?k�A	��A z�%pQ      h�     ��9G"�8�aD8`S7`!V6��O5���4ح3գ	                                                                    8{�&7�D�7'6N6��Q5���5 �W4� 3��                                                                    G���G]�Fm`�E�YoE�vD@f�CO��BQ��                                                4Ц�4"��3��02��(2$��1Us�0b�^/_�                                                                    4ǡ3M�2�0�2K�1PeV0�ϲ/�N�.�F)                                                                    7R��72�6�!�6=iB5��94�{�4.63%��                                                                    7�B�75��6���6H��5��5�4!�{3I�                                                                    5�4�f~4o��4�P3��2�Z�28�0���                                                                    7���7^%6���6us�5�0�5"�24E�z36v�                                                                    8�IY8�0\8B�7Y�6�7)5�N5�g44�                                                                    5��5��5e�4���4��3Y#d2���1�<                                                                    8�IY8�0\8B�7Y�6�7)5�N5�g44�                                                                    7�Uw7�>K7Y�6��65M�H4w:�3r�-                                                                    4��/4� �4�4�4
�3���2�3�2��1x�                                                                    7�Uw7�>K7Y�6��65M�H4w:�3r�-                                                                    �(*�2���5�E5q�4��4!� 3k��2�9�                                                                    5 ��5��a7 �+7F�6�m�5���5s��4�~u                                                                    5��_5~C4�a�4���3�33 ��2?@510�                                                                    7S�7��6q�5�M�5	��4=>3Y��2U�                                                                    �H>�ܵ���?`״�rH��J_�㉕��s�                                                                    ��y�բƲ�c1���������4�7������x                                                                    6��6�85sx�4�$�4�~3@�)2^��1Z�W                                                                    �%1��%���'�>o�\����\��ʃ��?I                                                                    ��S���"�O8���g���+D�bO߯���                                                                    ����`?���Ê��d�#vC�;�d�G
�                                                                                                                                                                        +�Ʋ+��m                                                                                            8�U�8b�.7��|7UY�6�d05�u5Od4
�                                                                    5��4�h4z��3ڥ�3/��2ycW1�g�0���                                                                    7���7]�6�6SӀ5�a4��4v}3* �                                                                    4y3�W�3ug/2�k25�1��0�M�/���                                                                    8��28YF&7���7RM6�ғ5��5"��48��                                                                    5H�4ۭ>4pLv3ׅJ3;�n2�[1�3�0у�                                                                    7�B�7ܾ�7���7L|6r�O5��j4��3��                                                                    6�5�G�5�Um51{�4���3��2��1��,                                                                    6�E�6�ad6uM5���5C��4��3�4r2�T>                                                                    4�84��4�P4k�3_��2��1��0��                                                                    8�68�t7�^�7=�B6�8�5ĪI4� 3�y                                                                    62G6+�5��05X��4�e\3�3 ��2�p                                                                    4�w�4�\84��X4��4U��4zt3�3g�{                                                                    3��3��3�y3]7y3,��2컏2�x�2;B�                                                                    4���4��}4��34�K�4�}�43q3�&3���                                                                                                                                                                        %pQ      h�     ��60�>3��            G�W�    =?�E6��j3�;@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�2�    @�:@    16:47:56        F� @�:@    @�B     %p�      i     ��5�s�7�d�B'�fB'�f8P>D"g�    )�G�+?�/7C��?�  ?{AB    3�l�(��M/�y��X        /���+�>�	�>3�=��j=��;I�5��.�x2h�:0w�<+�<�x�=G>�=��@=�ɪ=��=،e=��=��,>��>"~@���@���@7��@	.�?��?��L?��$?�#@m{V@mu�@c�@b��@b�b@b��@b�@b�-@b�E@b�@b��@b��=�.�=��?-}&                    E���7`aG� �G��Fl�xE���ET}D@9CO<:BQ�                                                A��C%��B�]�A�
[AV�6@�F�?�(�>��=�q�                                                                    E@�@�i�F+;A�?�0��kF[.B(��AG�t-E�N                                                    {@��BR?�BR?�@g 8{%n%;x�    ���/]5A&|Q?   )�G�A:�k���C2"h    >�y�B�*2�a@�I)@LF�@K��        B|7~B|7~C��C��C�M/>�s�@���2�:C���7
��6�U7A!�CV�>�y�C#�]BL�C:Z{A���BH~oCU�A�L�B���A�?����9    A�?�    A�?�A�<�            As�hA��A�L@J�N?�s@�@?��F��%0��k.NP�+��M1�`8�_>6i��5�K�G`"G?��G��GΘ>%��            >|�>�Xl>�4�>���>��0>��>��X>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    AG}AG}{@��7ժ_9`�w9F�_8�;�8>��7���6�\R6$�5�D                                                A��7AKJcAKJc{@�Ο��BFq�B�b�@JJ�>�>{@��{@��    C��\BW�E6��o    6@�7��z@��2�t 7^��    >%��C���7�`iBN"�D��C�{xB{oB+A�+�@��W?�W>�]�                                                ?\otA�է@���?���?.�>���>	� =4�<7��                                                                    D~�,F���F��D���D�+Cj�B�7ZAȯ�@��;                                                @�yqC��B�8�AN�.@��@<�*?� �>���=��D                                                                    C�ǆE�YEmB�Du�CkE�B�a�B��A,,@/�                                                @\1DBuU+B�?@�Q�@>X??�TY?�=>=,�=AK                                                                    7Z�7�0�A:�>0uUA�G�?�X?�:o,=�!*�=I-9��,�,=�!+ӆ.J�-��-Q�*Yg�7��}4���2TT-��/��/��|+c6�Į4���1��2:�2�-�%���Į2"h6��u4~�)    3���4�3k4��00v(    :z#w8n>8��6���5aW-3B�_                        5� �6�3+6�/�4�4t3}A�71(�2��j            4���>T�[-�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CJ�G��F>*�A�Q                                                B�                                          B��                    A�                @�p�    /�K�+��35r�6m2}1Q��18�?���8'(�@4��7�+�;���;��K                7�L�    1�з6N�        {@��            ��Δ{@��{@��{@��7���8$��            7!Lt    7!Lt{@��    7!Lt    6N�R    {@��    {@��                6N�R{@��6���7^��8'e�8'e�<�x1�ıB�Q�6���Fl�C�X*A'8�    C	�N            <�~4�!�6\	4�!�<�ш<��=q<Ø�<<��;�h�;:�V�ty��������C���T�������`�I�E�0����ٺ<�k<�͢<��R=�\<�{�<9]:;���;�d:��4�ty��������C���T�������`�I�E�0����ٺ9��    *��5�zD6W��9o��9�>Q9W��9M>�9;�8�*�                                                �_O=�a���_��ZL_�St��JVZ�@cA�4��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         ���        {@��{@��{@��{@��                                                    CTm#D��>D�=�DZVC�P7C��Be�DA��S@���                                                A�m#CR�>CA=�B�VBCP7A���@�D@�S?'��                                                                    E�[�G	�~G�6F��HFA�E�5�D�%�C�AC��                                                C�V�EH��E<�}E<�D�kkC�٘C��B.��A=�J                                                                    E� F�d�F��Fz�F\��F.��E�G�E���EH                                                D)D]D��D���D���D�vHD~Z"D0�rCߤSC���                                                                    71(�7P��GSGdGI��G:bF��lF}�"FN&E�QEP�r                                                                                                                                @&��A:��A�uA��=B�-B?aBb�B���<�<�<�<�<�<�<�<�<�<�<�<�E�t�E�D�Eh�E~D�a�DX4�C���C���                                                {@��{@��{@��{@��{@��D'��Aᓺ/@�(�g:48	�C�� @P�    Ar�Ar�{@�ξ|Y��|Y�C��E{@��C���C���C��E{@��@g {@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��;�	TD���G��D���D��D��A;A;F_8B5��    C�}FL^FL^D���D���F_��B5�                @x~�C��gC�5dC���?   C���C���C���C�y�C�W�C�9=C�kC��C��pC��5C�Q�C�!eC��lC���C�l�C�$�C��@C��*C�?\C��C��MC���C��C��
C��C�WC�}C�w�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C���>�JZ>���>�M>�)�>�`>��>��>���>��>ܙ�>��>ԫ�>�W}>˨E>Ƥ�>��/>�_�>��>�� >��@.��?�C        A����s`c�r� �r
���;���������S��\n~��4���4���4��dl#        >�6�B��    A.�A.�{@��@��FO�17w_m7��        ;���?e�?i ?�]?��?��?x�?�r?
�                                                ��w�6��A�  ?k�A	��A z�%p�      i     ��97�n8��8)�7]6���5��R4р33�i	                                                                    8h�7��07$j�6���5Ù�4�J�4P�3]I                                                                    G� �G��Fl�xE���ET}D@9CO<:BQ�                                                4��c4�f3�_!2��2"3�1P[@0[c/XO                                                                    3�*,3E��2�BE2D�1L��0���/��`.�~M                                                                    7B�7,N_6ʟs6;�5�7�4�kR4��3�-                                                                    7o��7'Ts6�o96&�O5���4��r4 3�                                                                    5=�4��k4H�I3�O3k��2�U�1��0�e                                                                    7�g,7L��6���6K�(5�"O5��4-��3+A�                                                                    8�va8� �8	�7S`�6�°5�l5E�4�
                                                                    5z=5|��5��4��@4�M3L��2��1�,,                                                                    8�va8� �8	�7S`�6�°5�l5E�4�
                                                                    7�_�7��:7S�6�|�6	��5C�b4i�3g��                                                                    4���4�E�4�b�4	%�3�"�2���2 w\0�Y7                                                                    7�_�7��:7S�6�|�6	��5C�b4i�3g��                                                                    �4V.��Il5��5~9�4★4&��3l�2�mp                                                                    4ڎ�5:9�7n��7?`�6���6-߈5�#35�                                                                    5�r�5j�04�iM4W;w3��*3
<�2'��1%y�                                                                    7�6���6l��5���5H�45t�3O
�2L��                                                                    ����Zg��ID�=-���7^��ᩲ�O���j�                                                                    ��L)��4����㲦�w���ٲ0�����F���                                                                    6b�5�k�5n4���4��38�2S�1Q-p                                                                    �2��3�����Q�Y`���" ���`���+                                                                    ��ͥ��M��^��޶��A��'ʰZ���w�                                                                    �p5 �X�A��ܳ��O�򷞲�ͱ5�k�@kP                                                                                                                                                                        +�K7,$�Y                                                                                            8x�{8Z�E7�K�7R��6��5�]5��4E                                                                     4�O&4���4u��3��3-}2s��1�0��                                                                    7s~�7VL6�V�6Q35��04��*4Uz3#�9                                                                    3�0>3�t_3p��2�dZ23�1�k0�~F/�f3                                                                    8nnY8Q�T7�a�7O��6�5�5���5	42q                                                                    4�W4���4k�`3��/39�2�׎1��s0ɶ�                                                                    7�t,7�. 7��17��6n�5�T~4���3�i�                                                                    6�5�I5��85/q�4��I3��G2�Q1�x�                                                                    6��6�D46q��5�Z5A4~ES3�wg2��                                                                    4�^�4��;4�T4��3\��2�L/1��0��1                                                                    8�8F�7�Ƒ7;��6�u5�J�4ځe3�                                                                    6%$p6�-5���5Vn�4��<3��;2��s2I�                                                                    4�T�4�T�4�Sq4�S�4R^"4;3���3` �                                                                    3�}{3�7u3|��3Z�3)�~2�{�2��G25`                                                                    4ׄ14̃�4��4�f�4���4/]3��X3��'                                                                                                                                                                        %p�      i     ��63�Z3���            F�J�    =G�6�%�3��+@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�:@    @�B     16:47:59        F�L @�B     @�I�    %q      i&     � 5�Dl7���B'�fB'�f7� �D"g�    )  ~+^ۚ7-�?�  ?(    3���(H�/�O����        /XA+N�>��><�=��^=9�:���3���,dI�2X�:'~<+�<�u|=G8l=��=��j=���=؋�=몒=��->��>"�@��[@��X@ G�?��?��!?��?��*?���@o@�@o47@d�D@d�5@d[�@d0�@d �@c�'@c�%@cz�@c\�@cL-?��=D��?�c                    E�f7l��G��*G�AFlTE�`NE�D?�{CN� BQr�                                                A_�C%WB��7A�^�AV@���?��E>ح=�K                                                                    E@��@đ�FX�A�a�2%=�F[��B)��A=��.w�                                                    {@��A�
+A�
+@lV�8R�A%r��    �{�/BC_�?��?   )  ~A?���ľ�1��8    >�y�Aau42!|@�jP����7z*        BU��BU��C��C��C�+�>�mK@M�W3 6�C��6�.6�A#f�B�o�>�y�B�uB�sC ��AJ$	A���B�JzA��&Bj�@�R���g!    @�R�    @�R�@�,�<��    7t�A�@��xA�F�@�N?�C�@U�?^*F���2%=�/��1,�2>h8""�6	V�5ђ F��G+g?G�?�G�>>q�            >%��>�HX>�g>�>>�8�>�[�>��>�!�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    A��A��{@��7�ټ9gE�9R�!8�Y8K�c7���6��q6X�5.\n                                                A�|A݋A݋{@��:���B
?�Bi I@X>��"{@��{@��    C�erBX�S6��c    6C�7�2�@�72��D7QW�    >p��C���7���BB�D��C�zB[�AA�Aa�M@��?�Ŝ>���                                                ?N��A�S@��?���?ȉ>�?�=��='�<2��                                                                    D}�F��MF�QD��[D\*Ce��B��A��@Ɋ                                                @��C�B�1
AER>@�C�@4�c?�&M>�T�=�P=                                                                    C���E��9Ek'�D��Cd��B���B	C0A(K�@,��                                                @W��Bsx2A���@�"�@5�w?��?
�>7��==�E                                                                    7f�17�/�A��>5�JA�J?Y�?��+���*L��- ��,j�y+���*�U�-8�-0nS,�0*��7��4i�M2��-Vz/��]/�s�*��^�C*�4�r�1�ε20n2M-�7C*�1��8�C)�4�6�    25�4aI�4a�/��    :yqe7�D�7�5ʔu4vK523��                        4F1 5�cg5��3m��2}��6H�2�b            4�<L=���,�d?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C��G��>+�%A���                                                B�                                          B��                    A�                @�p�    /{�*�	�4�t"6vx�2x��2BN�?:��7�R�?�ݚ7Y�;͉c;ˮT                7K˄    3,Ք5��~    3PJ�{@��            ��B�{@��{@��{@��7���8L            7�    7�{@��    7�    5�s�    {@��3PJ�{@��3N��3/��1���    5�s�{@��6��6�W8 ט8 ט<�5�0�VB�˂6�W�F(��C��@�h    B�i            =: �4t�7�F�4t�<�Z�<�Ӕ>I�>{�=�0�=Y�<���<	�ۋ�N����z��š������b�����b�ъG�i�2���!\�5[=6�<�WH<�͵>IX9>7=���=C�<��^;��ۋ�N����z��š������b�����b�ъG�i�2���!\�5[:O�4    )Жi5䍳6;�I:b�9���9�$�9�	�:��M:�Wn                                                ��2��ӵB��l��#�����ßmĻ6�īW�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         :���        {@��{@��{@��{@��                1ٞ                1ٞ            CS��D��!D�G�DYQ]C�d$C�Bdi�A�5@��                                                A���CR�!CAG�B�Q]BBd$A��@�i�@5?&�                                                                    E�dG
�G�F��$FA	�E�.�D��C���C�I                                                C�cEH�EE<��E?`D�dBC�ϲC��B.�A=ׯ                                                                    E�.F�h�F��@Fz!�F\��F.��E�GE��LEH                                                D)E9D�jD��D���D�wD~Y�D0��CߣWC��<                                                                    6H�7^�GGS^�GI��G99F��_F}��FK�E�OJEP��                                                                                                                                @O[�A\ZnA�i�A�͠B(�-BV2�B~M�B���<�<�<�<�<�<�<�<�<�<�<�<�E��^E�U�Eg��Ew>D�\mDX0�C���C��b                                                {@��{@��{@��{@��{@��D!��B��.��gb8�C�8S@T��    AR�AR�{@�ξ�䊾��C�7�{@��C�I^C�I}C�7�{@��@lV�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��;X��D��;G�-D��D�eiD�eiA6�A6�F_^�B7K�    C�a�FLa"FLa"D���D���F_�B7Ff                @'�ZC��C�|�C�J?   C�C�C�I^C�I^C�IMC�H�C�F�C�BZC�:dC�-bC��C�AC��GC�ϏC�� C���C�U6C�VC��HC���C�_�C�,�C�C��C���C��C�>C��C�I{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�DC���>�>���>ߙ+>�W�>��>�g�>�xy>�E>ڃ�>���>��C>ԝ�>���>�ȓ>��>��>�
A>�S�>�B�>��@A �@u%        A{w�����]�L��X��ǌݑǌu7ǌ��ۻ$�ę��ę��ę���x�        >��6A�~    @��@��{@��@�۾FP37mb�7�6        ;��.?4n?3�?4�?4J�?4�?6\�?8�W?=xM                                                ��Q6�A�  ?k�A	��A z�%q      i&     � 9<��8��B8S7o��6�ɓ5���4�-�3�w                                                                    8n=�7�H�7/�16�H�5�w�5�O4*f3�                                                                    G��*G�AFlTE�`NE�D?�{CN� BQr�                                                4ŝ�4%N33���2���21�n1gk�0wR-/y�d                                                                    3���3PΑ2�a�2{�1`��0�)*/�3�.��?                                                                    7H��76j�6ضR6J��5���4��4H�38&�                                                                    7pj*7*7�6��H6V�5���4�ԣ4ě3��                                                                    544�U�4@*�3���3c&H2�s1��1 ��                                                                    7��7P^6��6A��5���5+45��3@��                                                                    8��}8��S8#7a�n6�25�v�5�\4�                                                                    5�-5�#y5Y4��4J{3[��2��
1��c                                                                    8��}8��S8#7a�n6�25�v�5�\4�                                                                    7ﲭ7�C�7b��6���6��5TpH4��F3��                                                                    4�vX4�t�4�{�4��3�$ 2֚V2��1��                                                                    7ﲭ7�C�7b��6���6��5TpH4��F3��                                                                    �H�_��>e5���5�{�5��4?��3�<�2���                                                                    3�܆4 �8�O68��d8/ >7���7"56���                                                                    5�_�5oD�4�]4L�[3��b3�206R1:w�                                                                    7��7/�6/�5�k�5�p4Fx�3e�32j<F                                                                    ��<���p�Mq���/�׹�������g                                                                    ��Þ��DT�����$���o��Dj�������u                                                                    6�6�L5��4�A34��3I�>2j��1o�                                                                    ��?����'�P˵m�"�����ذ{���k                                                                    ��MK��p�!�]���h���8���ut���u                                                                    �w|��eZQ�#���[k�Up�1Վ�MV�^x                                                                                                                                                                        +���+��                                                                                            8ȡ8gW�8�97c��6�(�5���5��4+��                                                                    5N_4��4�b3�l3=H�2�ۘ1�v�0�                                                                    7z�b7b��6���6b$5�Qj5�_4"�3<W�                                                                    3�W�3�"v3���2���2Cݝ1���0��i/�X*                                                                    8u\#8]�"7�M7`�/6�z?6�U50^�4L��                                                                    4��4�^]4|�3��3Jr_2�81��0�c                                                                    7�7���7��y7&��6�@5�74��!3��?                                                                    666�5�s�5>��4� 3�?2�1�G#                                                                    6��I6�p�6���6�C5T24�u�3�2���                                                                    4���4Ѐ�4�>�4�3rl�2���1�i70�>�                                                                    8��8	�U7�1[7K��6�j�5���4��4�C                                                                    6*%�6�5�8h5hݳ4�Uy3�3�+2�                                                                    4��?4�D4��L4���4g*}4c3Ǆ�3�|U                                                                    3��3� 73��)3m��3:� 3 �02�:2QE                                                                    4��4�zT4�$4���4�D�4B�o3��3�B�                                                                                                                                                                        %q      i&     � 67C�25�                        6���3̻�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�B     @�I�    16:48:02        FҊ @�I�    @�Q@    %q}      iE     ��5���7\��B'�fB'�f7��5D"g�    ([7+���7��?�  ?�?3��I3���-��/�����        .V@�*L6|>���>S��=ۮr=b��</�5�-��72H��:%q�<+�r<�zT=G4�=���=�ũ=��q=؊�=�D=��1>�>"�@�m�@l��@��?���?��?S��?-a?CT@n�@orX@d�@d�@dβ@d�p@d��@d�@dkk@dR�@d>�@d2�@��7<�G5{@��                    E��77ctG�:G+�FkmE��lE��D?7QCN`�BQ (                                                A�C$��B���A���AU��@���?�C�>��=��r                                                                    EA#@Ĳ�F{�A�}�3}�mF[��B*|�A/��.��                                                    {@��AMʛAMʛ@m��82Rw�%
    &-�J/K�ǥ!��?�M([7@�U��۷3�=    >��@]��2
��@� �!l��=�        B$]�B$]�C��lC��lC�)�>�:#?�03%�C�),6���6���A$��Bs�>��B2��A^8�ByK�A3A>�B6{�AG,�B0"!�-�I>�z!�>��-�I    �-�I�4µ>��d8�|�9��@�@��Av@?�X�?X�@l�+?t�F�|3}�m0��R.�:3��7��^5�I5���F�@GN�G�c~G���>���    9��    >Q��>�6>��>���>�"e>�<1>�
k>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�\�    @���@���{@��7��9W3O9H��8���8F�H7�o�6���6%\�5J�                                                A��@9��@9��{@��=
�A�DjBI0W@c��>&u{@��{@��    C��2BZzR6��    6p�g7��@��3UQ�7+^�    >���C��7�r�B;v D�tC���BKfA���AO�@�+�?�]�>�                                                ?G<�A���@�m?q��?��>��w=��*=/<)��                                                                    D{��F��Fh9D�1_D}�C`��B��A��]@�[                                                @�vC�DB���A>?@�d�@.|�?���>�A�=�4                                                                    C�nNE��1Ei.�D�C_�B�$�B�A$0�@).�                                                @Th�Br�]A�� @���@.̓?�-q?(>1�D=8V�                                                                    7\��7�8�A�l>9�A��?�;?�5U*�)y��,�\�+�O�*�)���+��j,���,X��)���7�ԁ3��f1+^,VA�.���.���)�������3�n�1�w{1&��1&��,P[�7���1��.����4�Z�    ,�3���3�a�.��Q    ;��5�R�6��2�Z+1W�.�Ni                        .���2���2�2�/��    2���2���            3���;W-�,��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B�~�G���>)GbA��                                                B�                                          B��                    A�                @�p�    .|z)��y3�4�6h4�3�G�3�p�>�- 7b�?<*�6�d2;�|�;��7�\            6���    3P'5v_�    5i�0{@��    1[�    ����{@��{@��{@��8F
8��B            7���    7�k�{@��    7�k�    5i��    {@��5i�0{@��5UE�5�4�_�    5i��{@��6A��5���8�FS8�FS<��    B�i�6�"�FCeCɏ�@ �H    Bi�            >�?93��9H�3��<�Ί<�9�?�K?�\�?	@>�;=�+�=i9ϋ�ċ
k����ѪF���ˊ�N���
k�b�=�G�[�2p;�!G�"(>���<���<�5�?��?�(??�>��y=��=(�͋�ċ
k����ѪF���ˊ�N���
k�b�=�G�[�2p;�!G�"(;�/?�*��5��6u�;h�:�N:Hž9�<bNZ<��                                                �T�T�i�I|�=)T�*���&��B�Â��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         <�n�        {@��{@��{@��{@��                5�07"p%        ;�G�5�0            CS�Dҳ�D�TDW��C���C��Bb�A�x@���                                                A��CR��CATB���B@��A���@�@x?%��                                                                    E�jaG
)wG��F���F@�^E�#GD�eC��~C��                                                C�l0EH�E<�lE@lD�Y�C޾�C�B.�A=��                                                                    E��F�lAF���Fz%�F\��F.ݍE�E�E��EH�                                                D)E�D��D��YD���D�w�D~Y�D0��Cߡ�C��'                                                                    2���7Yh�GSu5GI��G4�F���F}�FG�E�LEP��                                                <E1�                                                                            @���A���Aڌ�B
�vBC ?B|�=B�$B�y|<�<�<�<�<�<�<�<�<�<�<�<�E��QE�g-Eg�MElrD�S	DX)�C���C�                                                {@��{@��{@��{@��{@��D�2B��-�.��f��8E0Cm��@FM�    AoAo{@�ξ��#���#C���{@��C���C��vC���{@��@m��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��Gf&;�D���G�D���D���D���A3�?A3�?F_ �B6�{<E1�C���FLaAFLaAD��D��F_��B6�^                ?�؉C���C�gC��R?�WC���C���C���C�C�/-C�QC�s�C���C���C��pC���C��C�"�C�*�C�*�C�!IC�HC��MC���C��KC�b~C�9�C��C��EC��C�C��C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��TC���>��Z>�h>���>���>��D>��>�\>ǧ(>�\�>ʟ,>ˏ>��>�	8>�l�>� >��T>�I2>�WN>�{T>�S@M�4@7$        A������߃�ڽ�ĥ���	L�̬��B�e�ƽ�ƽ�ƽ��cf*        =�4�A�c    @N��@N��{@��@��FPd�7H:�7�        ;��s?L�	?L�?O(�?Q�Z?U��?]�?g�F?ux�                                                �+��6�¢A�  ?2v�@���@�ѥ%q}      iE     ��9/p8���8��7m�6��5�5�u4��                                                                    8].�7�f-7*_6��/5�w�5�I4& �30��                                                                    G�:G+�FkmE��lE��D?7QCN`�BQ (                                                4�v�4�3�N�2�\?25) 1u#v0���/�N	                                                                    3羼3Fr�2�~�2��1d�y0��/���.��]                                                                    7:�M7-�6�O}6G�35�L�4��4'�N3WH`                                                                    7[��7�f6���6�)5��c4�A]4ی30Av                                                                    4�qG4��P4/z�3�EM3V��2�s�1��1��                                                                    7�AX7?H�6�{�60�5���5�4>~3Wlt                                                                    8��8�a�8i7[|�6�3;6��5$�V41׶                                                                    5nC�5{�52�4�vV4%p3b�2��1��                                                                    8��8�a�8i7[|�6�3;6��5$�V41׶                                                                    7�qJ7�#d7U��6�5�6ے5\Q�4���3���                                                                    4�Ƥ4�t4���4�H3��72�_�2�	1%�x                                                                    7�qJ7�#d7U��6�5�6ے5\Q�4���3���                                                                    �A����r5��n5�]�5�4O�H3�ek2�Mz                                                                    0���0��":g��:2F�9���9��8��8l�                                                                    5�m5[��4�$Q4:�3��3
�q28�1P�8                                                                    7��7i�6rE�5��V5i�4OX�3|ck2��                                                                    �ח��w����K�д�'c����
��v�                                                                    ��Z²�0���S��ʦ��{<�P��	D��g                                                                    6H6�5rVz4�!n4'�3RF�2�o�1�9F                                                                    �����*��WҶ%�p����j���J��                                                                     ��9���aֳò��߲ ���B煰��i��Z�                                                                    �f�7�Zǂ�Tn��
	�%D�<ъ�d飰���                                                                    )�)B .q�.�-�F�,�[2/���/�
                                                                    +Ҟq,iu                                                                                            8m�q8\��7���7`Q�6���6
5#��4Hd�                                                                    4�44�4}U53���3?̈́2�41�m�0���                                                                    7i%�7X�6�t�6^�5���5	�43~�3[��                                                                    3�3�{�3x+�2�>32Fx�1�0�R /��                                                                    8dL�8S�D7�Y�7]�6���6�5CJ&4o�                                                                    4���4��K4sI3♀3M#�2�'�1�6�1fs                                                                    7���7�\d7��7%i�6��5��C4�U�4 �                                                                    6��5� s5��5=K4�)�3�k�3 �l2��                                                                    6�fH6�m6{�i6��5X�`4�2P3��2�҉                                                                    4љw4���4��4�N3w��2��71�70��                                                                    8
�'8�7�}�7J,6���5�,5	��4*j                                                                    6�v6i*5ٴ5g�4�3=4�3`�23�1                                                                    4�4��B4��J4��~4l�4)E�3�Ɛ3��                                                                    3��<3�s�3��z3k��3>˖3�2�M2u�$                                                                    4Υ�4�b�4�"4�Qb4�I�4N�J4$3���                                                                                                                                                                        %q}      iE     ��6c�,�F�P�    =%�T            6�0k3���@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�I�    @�Q@    16:48:04        F�� @�Q@    @�X�    %q�      ic     5n��7)OB'�fB'�f7F�D"g�    &gZr+��Y6�	?�  ?�T4� �3Yƽ2S�0ڰ�/�1|/�    ,���(�<}>�v>��=�Gf=m9}<�8�;a�V3���2@d:)��<,U'<��=G5�=��F=���=���=؊�=�)=��Q>�->#@pm_@�&?�/�>��2>'��<���<�3 >`τ@j��@o�,@e�@e	�@d�c@d�o@d��@d��@d�A@di�@dV9@dJ�A"�$<�a�{@��                    E�O�7! G���G�Fj�E�@�EO�D>�YCNBPy	                                                A��C$:�B�)<A��ZAT�@�4�?��>׼B=�D�                                                                    EA)�@��"F��A��3<��F\B*��@���-��                                                    {@��@1:�@1:�@mQ8 ��&��w���*$���/#E�:gc?��&gZr��p�/&\�4K��    >��q>ې�2�,@z �ź���A�;�y    Br�Br�C�C�C�C�Co��?	RW? ��3�gCޤ�6��6�U%A$��A� �>��qA���@��B-5@��k@vV*A���@�F\B����F�A3}��Q��F�    ��F���C?��h>��g>�˳@e��?��@�m@?W}�?,5@}��?{��F�#03<��0�z)-�i3>�~6�z�2��#5��F��F��G&�"G-�?:C�=/��@�e ?Ǒ�?�q?�?��?��>�c>���>�8�>�r18��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M۵/    �V�}�V�}{@��7���9ߣ9�8��P8��7y�6�|5��5<��                                                A�2����{@��A�yPA ��B15
@h��=��{@��{@��    C�3"B\��6M"�    6Ivt7!sz@�`$3�:6�:�    >���CʅM7��B6>�D��C�X�BA
�A�'�AH�@�g�?�e>�j                                                ?B2IA��@��1?gr�?Ӄ>��M=ꐷ=�<��                                                                    DzF�(*Fb�D�j1D	��C]D0B��8A���@�!                                                @��Cr+B��A7�7@�:�@)�?�[v>�g=�V�                                                                    C�iE��Egj�D�|CY��B���B{A"@%K�                                                @Q��Bq�xA��@���@(�1?���?�A>.N�=2��                                                                    7�7��A�x>:�A�μ?�I?�*�(��J,�O+ms�*�(�$c*��, �\+Ȳ�)Z�7F�1*^�.���(�?�+��+Ǚ%�X�·%3YQ1�p�.���.��c(�4�7·%��D2�·#4b��    'r�&1"��1"��+��    <��Y��M94a�/���.G��*���                        *>S/f�t/e+��    /��]2|t�            2��M:�&c,�K?g��?y�?L3(?S�?B�f?�'?
�Q?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B:fWG���>*4gA��                                                B�                                          B��                    A�                @�p�    ,@��(+�"1=Q�6 3/3NVd3?$�=�O`6���>p�A5O�;O��;P��7��)            5�F    4i;l4o��    6�|S{@��0�B�5*9�    �O|j{@��{@��{@��7�=�8��            7-B�    83@{@��    83@    5+�"1��{@��6��{@��7�6��6��8    5�0{@�ε(�48�7��t7��t<�3    B��X6�A]F_�,D=�?~�>�P�A�	x            ?][l1*^�9�?�1*^�<���<���@�<@S�q?�f�?:C>_
'=��r��ċ
l����ѪF���ˊ�N���
l�b�>�G�[�2p;�!G�"(?\:<��l<���@�'�@S�G?�a�?9}>_&=�^W��ċ
l����ѪF���ˊ�N���
l�b�>�G�[�2p;�!G�"(;��3/��j)-�75&��5��,;#�:�Q�9%~7E=V7 \ <��o                                                �`_�ÜM7â�%ç�ð3'Î���R����̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 4^1R3�>���6[{(5� f@Z 6@_hq@B m@��=���5�t�5;�D�<7�B<�Xn@��:�ȝ=U�S7�BHA�p4�LR7�CR�DҌ�D�R�DV�C�7�CY�B`��A��@�^�                                                A��CR��CAR�B��B?7�A�Y�@���@�?$^�                                                                    E�n�G
:�G�ZF��zF@��E��D�[C���C|N                                                C�rmEI�E<��E?�D�M�Cޯ�C��B.w/A=�                                                                    E���F�oFF��=Fz)�F\�.F.�eE�EE��LEH�                                                D)F�D�8D��*D���D�x�D~YND0�pCߠoC��                                                                    /��]7ݍGS��GI�G.GF��ZF}�>FDfE�IdEP��                                                @�-�@�A$@ʓ�@hG�>��$                                                            @[g�Am��A�y�B]�B^.WB�E�B�,QB�<�<�<�<�<�<�<�<�<�<�<�<�E���E�v�Eg�E`!D�JDX#�C��	C��k                                                {@��{@��{@��{@��{@��C���BY��+��)�g�7�,C?�l@+�    @~��@~��{@�ξ�IR��IRC��{@��C���C�˃C��{@��@m�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�=�̀D���GS�D���D���D���A1��A1��F^*�B5FBA���C�/EFL_�FL_�D���D���F_/B5A&                ?B��C�]�C�^4C���?�;C���C���C���C���C�;�C���C��C�tC�d0C���C���C�
�C�4�C�Y�C�x�C���C���C��tC���C���C�s�C�UKC�*#C��C�2C��C��C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��mC��I>��0>���>�]�>���>�NL>��>�޿>���>���>�T�>��>�O>�ܞ>�b�>�mO>��>�X�>�2�>���>���@N�b@�<        @��oƾ��ƾ��ƾ�	���β��Ο��΋6�N?pƺ��ƺ��ƺ����0        =]�4A���    �v��v�{@��@��pFPzO77�<        ;E��?~?'D?7J?R��?mU�?u�\?{aa?�N                                                �d?V62ߟA�  >��>(��@R�5%q�      ic     8�A8WMm7��-7>�6�L�5���4���4��                                                                    8�i7��7 �6qB5�S�4�3��K3&\a                                                                    G���G�Fj�E�@�EO�D>�YCNBPy	                                                4n�3�R3T{J2���2oP1�W0b�/���                                                                    3�a�3y�2�2�1��/1+X0:~�/;oH.�-]                                                                    6��6�p�6��n6 j_5|<�4���3�?�3JN�                                                                    7�u6��/6=Y[5ڍ�5J�4��z3���3ւ                                                                    4�UR4���3�:�3�˫3JD2Z��1�گ0���                                                                    7,0 76gm66��5x�4���3�Kq3?��                                                                    8Ds�8Q�7�r7.�6�Cp5��g4��4$��                                                                    5��57�4հ�4^�a3�)�36�2$61��                                                                    8Ds�8Q�7�r7.�6�Cp5��g4��4$��                                                                    7�jP7�N�7��6���5�J5_x4�/3�M|                                                                    4�^74��4I�3��3Dҁ2���1��$1�9                                                                    7�jP7�N�7��6���5�J5_x4�/3�M|                                                                    � <	����5��(5p��4�k�3��:3&��2�>�                                                                    7&�-6o$�;��:���:_3|9���8��8`�                                                                    5M@5k04|�F4(3z��2��1�,h19�                                                                    6� 6��63QY5�C4ռu3�%a3�&2|%0                                                                    ���׵� ���Ļ�$x��wc���$g��ֽ��w                                                                    �������m��YD��3�Y��������)E��r                                                                    5�+�5�S53m4��3���2�f!2	�1�I                                                                     �����Û �f�G��V�3��Wo����	��                                                                    �@�s�E����-�g����`��~;��߯��i                                                                    �~�����w���3���ױ�����
��v�Z                                                                    )<��)�Q/J.�0�-& �+M*�+$�0��                                                                    +��+H�^                                                                                            8��8�7�}$73��6�q\5�S4�b|4;ܖ                                                                    4�d�4���4=��3�3��2*5u1F�v0�̏                                                                    7�7��6��%62W5�4���3�N�3M�                                                                    3�53�B<39�2��.2��12�0Y��/�U                                                                    8a=8M�7��%71P6���5���4�;4`�                                                                    4�L4���45�J3�sN3}�2;�1mA0�ݔ                                                                    7��t7�R7j�7v�6H�5`2�4s.Y3��                                                                    5���5�@]5�B#5��4em�3��2���2
�3                                                                    6n��6z��6=�<5׳.5"8�45+�3D��2�TH                                                                    4�_[4�;�4X��3���39e�2O21`�10�`R                                                                    7�{ 7���7���7#�6u\�5��4��S4yV                                                                    5�Cm5ؤ 5��5:m&4�4�3��.2��:2)�>                                                                    4[�%4s]4u��4k�S40��3�73puW3��                                                                    31��3Dre3F}w3>PC3�y2��2BON2h5�                                                                    4�fA4��4��4��4X0�3���3��`3���                                                                                                                                                                        %q�      ic     6<�'r�&    =��                6��Z3�C�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�Q@    @�X�    16:48:07        F� @�X�    @�`�    %��      i�     �`3+y6�}�B'�fB'�f7/8ED"g�    
q<<+��6�!?�  ?��3�1V�|1�N/�6�.�        ���,8�=�1���+<4.�,c5�E�<�>Z70�2>�:3�<,�<�+=G=E=��h=��*=���=؋1=��=��>��>#�8���1���+<4+!��+S.��6)=���@fi�@ox�@eh@e	`@d�7@d�B@d�v@d��@d�@di�@dV@dJ�{@��{@��{@��                    E���4�ZG�9�G(aFkdBE˾-E�D?&�CNz�BP9'                                                A�!C$��B��lA��dAUn�@��$?�1�>�7-=�3                                                                    EA$�@ļF�XA���1���F[�4B)M-<�D/,�#�                                                    {@�ο�����@adq7>q�%��Y$p�,�Լ|/O$�`5�?��
q<<��K�.���2�O�    >���<�{1��,@�L�}���G2�=%�    B=�eB=�eC�O�C�O�Co-u?�  >H,B3cl<C��6f�\7�3A"`7@�@�>���A��?�e]AJ�@)Di?�Q�A�	@S�A�mG����@�#f��-�����    ������6�?�u-?o�h?o�h?��6>��@,��>���>��s@�:?6JdF�*�1���/0�p,�Bz2&:6�    3o�ZFH��Fm��F�7�F�{�?Z@8��B>�@�� ?�H5?O��?+	??�>�V>���>��\8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�F&    ����{@��5s��                5��6<gb5`i�5�P                                                A�����R��R{@��C8N�@�GB(XH@_��<��{@��{@��    C���B\k�6�<    6�X27Pj�@���3�4�6���    >�4�C��7�pBLi�D���CΨ�B��B�*A�6h@��6@j0? ��                                                ?_	�A�:L@�`:?�H�?<��>�	>��=T��<<m�                                                                    D|�tF�J�F�[D�:�D��Ce?pB�W�Aǯ@ŕ"                                                @�y(C:�B���AD�@���@5>�?�'�>�[�=�F]                                                                    C�� E���Ei�TD'Cc2BB��'B	�dA+"�@(�                                                @X<8Bu��A�*5@��B@5H�?��c?��>;�=7��                                                                    4ؙ�5T�PAR�>9`xA�u�?m�?�GZ)���(4�v+�q�+�)���(qI)HN�+��#+~fg(��\7/8E                        �>3cl<1�3�            7>���*�>2�W�                        >w���                                                                 3 j8            3 j8{@��,���>L��>L��>L��>L��>L��?y�>�$�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A��cG�)p>+E�A�h                                                B�                                          B��                    A�                @�p�    Gj(-�_    3+y2?�2&:=?QA6��=���49�	;�;"9�3�O�            �X?6    4F��16�    5�A�{@��3�+5M/5�(�V�{@��{@��{@�δ��7'�            5��I    6�!{@��    6�!    5O�6wH�{@��6���{@��7��6 ڧ6�~    5g
�{@�ε�8�2I�55���5���<�Y_    B�p6s��F}t�D���>n�a=�@��/            ?�Y    : �t    >���>J�@���@�N�?�:?Us�>�"I>	{��ċ
l����ѪF���ˊ�N���
l�b�>�G�[�2p;�!G�"(?��>���>J�@���@�7�?��?Us8>�"@=��(��ċ
l����ѪF���ˊ�N���
l�b�>�G�[�2p;�!G�"(<+�a+�D�$�!p0�X{17�E;(;8=;zN6���4���=$��                                                �!���C�������)��E{�V
�Q���̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6n5
�>�Ճ7��{6�.?��?��?� �@>�=�]�7e6DJ�FH:z7;G>�B<;P>��Y>��7;GJ�P6+��7�=YCQ��D�V�D�%�DU�C��MC�B`g{A��_@�e�                                                Aѱ�CRV�CA%�BթB>�MA��@�g{@�_?#e�                                                                    E�muG
'�G��F���F@��E�%D�C���Cx�                                                C�p�EH��E<�1EC�D�[�C���C��B.~dA=�                                                                    E��oF�k�F��Fz&�F\��F.�%E�F3E��~EH�                                                D)F"D�1D��nD���D�x�D~ZeD0�=Cߢ,C�?                                                                        5GGSp\GI��G.�F��F}�SFE%E�KEP�7                                                A�1*A׬BA��A��A�u@�|                                                        >�4h@D�U@�5�AK�mA�w�B_�7B� 
B��6<�<�<�<�<�<�<�<�<�<�<�<�E�� E�jbEg��EeSD�M:DX$�C��IC��*                                                {@��{@��{@��{@��{@��A�E�B��ڧg��6�5�C�@�    ?��?��{@�ξ�R���R�C�UU{@��C��fC�tC�UU{@��@b�	{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�n>�˒D��OG@<D��&D�KQD�KQA8<;A8<;F]�B2�BC�0C��FL]�FL]�D��D��F^?B2�&                >��uC���C��~C�[#?ĭC�=FC���C���C�ڔC� _C�.�C�nC��2C���C�H)C�� C��C�qC�L�C���C�� C��C�-�C�M�C�]�C�^2C�S3C�8�C��C��C��C��C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�3�C�7�>��>� >��>��!>��_>�7�>���>�y[>�+>�_w>��Y>�e>�ܔ>�OR>��I>�U�>�>F>�/�>�8�>��8@;�M@�        @*� ƽ�ƽ��ƽ���[҄ƿ��ƿ��ƿ��sƽ=cƽ=cƽ=c�W�
        =L��A�!�    �ح��ح�{@��@���FPt_6�([7�/        7+zd                >���?@�?M&�?Z��                                                �5uA�  ?kt    <#�
%��      i�     �`                4�p65+�$4,a�3��                                                                                    3��4X�3Y��2��                                                                    G�9�G(aFkdBE˾-E�D?&�CNz�BP9'                                                                0<�U0�ז/���/R�J                                                                                    /n��/�+E.�4.���                                                                                    3�Y�46��3Zt�3�                                                                                    3��4S@m3��3��                                                                                    1���2'�1`�0��                                                                                    3ߐ�4�&3��H3,                                                                                    4�y>5Er24_�q4 W                                                                                    2�V2�y�1�j�1r��                                                                                    4�y>5Er24_�q4 W                                                                                    4�4�� 3��3Z��                                                                                    1�S2'K&1R9�0��L                                                                                    4�4�� 3��3Z��                                                                                    3Ť3�Dq2�G�2���                                                                    9(8��S;9V�;e�:�^89��c9'Y�8���                                                                                    1��p2~�Z1��1&�d                                                                                    3>�3���2�3�2Ce3                                                                                    ��h4�(�5�����                                                                                    ������A�3ӳ��X`                                                                                    2FF2�1��1G5�                                                                                    �x�
�������8�                                                                                    ��C������X�l��                                                                                    �񱊔����<��                                                                    !)qR!��6+$�R+4J�(/�� pfE$�&,�h!                                                                    &b��&�<�                                                                                                            4�Zx5>��4Tn�4��                                                                                    1F
T1ΐ�0�r0��K                                                                                    3��<4H�3h��3��                                                                                    0L�0��R0�/�Kb                                                                                    4�>5Qx�4}mf4*J�                                                                                    1S϶1��1��0��y                                                                                    4��q5j�4��3��U                                                                                    2���3�2(��1�:�                                                                                    3b�3�x�2�Q2��                                                                                    1�4O1���1O
0��                                                                                    4��n5&�144eg3���                                                                                    2�k�3>��2N*�2�                                                                                    2v}'3x73��3[��                                                                                    1G.�2H�<1�ְ21�~                                                                                    2���3���32Z\3�<_                                                                                                                                                                        %��      i�     �`6�#�    F�% <��=�    >��    6���4: @� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�X�    @�`�    16:48:09        