CDF      
      time       bnds      lndgrid       levsoi        levdcmp       cft       glc_nec    
   ltype      	   natpft        levlak     
   nvegwcs       string_length         levgrnd       hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      Sun Jan  9 16:23:30 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1961.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1961.nc
created on 12/10/21 17:07:15    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:30 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1961.nc had following "history" attribute:
created on 12/10/21 17:07:15
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�fGb @��`    @��@    +:Y      �b     �`3�q6��B'�fB'�f7+�
D"g�        +��6�L?�  ?��4���1QM�2��0��}0Mx"                9�*2�d0OH6��=�D�=1�?7�bM25�:'@<!B<�=G�=�3�=��=��=���=�9�=��e>)P>��9.02���,U;1���8k��1,:D3j-;�@^'�@m�@b�b@b�=@b��@c%@c0'@cO�@cmW@c��@c�M@c��{@��{@��{@��                    E�()5ubmG�hF�Y�F\P
E�X;E�!D9�eCJՖBO*                                                AG�C ��B��bA��AI�@�M1?��>�F!=���                                                                    EL��@�4�F)�@A��`2`dpFi�4B6x�<�D/- $                                                    {@���X�K�X�K@lL�7W�n&l�$k�Ψ��S/P%�x��>��    �j�/`�L2��8    >�y�=LӢ1郩>m����0`����<%){    BS�BS�C���C���Cs �?j�>j�3<�C�r�6?�A6� �A/��A7��>�y�Av�c@l�A���@�oV@�|A[�\@��%Aۀ������-�    ����    �������t    ?}kZ?}kZ@}8?`:�@��w?O�7?P��@�mC?��aF��2`dp/�r-�[2x`�68B�    4)9�FE�RFLr�F�tFF�M*?��@�SCB� @u�)?��?O��?+	?>��1>�N�>�Š>�|�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M���    ���|���|{@��6/��6��6�^�4ř60(�7)ߨ6}Rx5uk�4�>                                                A����3���3��{@��CB��A�B5�@i{J=4�M{@��{@��    C�sBS!�5��W    6Q�t7��@�+l3��6��    >�v�C�2�7��zA�1�D5_�C��B���B	ԄAO��@�JE@g>��C                                                ?��A1��@ذ?���?6�)>�#�>�=F��<*&?                                                                    D~FUnEaED��C�#�CFPkB�IDA�tm@�
                                                @���B�aA�S)A|�@�@F�?�\�>�Ab=���                                                                    CW�^E�B1D���C�eCF�B���B ��A!y�@hj                                                @ �HB2�A>��@��i@�?�J�?�9>2��=*�6                                                                    5��06�A!R�>O�uBB-?,�?�V�)}�h(`+�|]+��)}�h(K*�)�F�+�V�+h��(�kg7+�
                        �)��3<�21��            7)�Ĵ�MԷ)��3f�                        >uWZ��|�                                                                2��            2��{@��,�Y�>L��>L��>L��>L��?B��?7a?d�?P�.?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A���G��(>+�1A�X�                                                B�                                          B��                    A�                @�p�                3�2���2xam=K��6'�X=�V4�J�;-�;�%                ���~                    {@��4��4%ȋ6�� ��vE{@��{@��{@��1*�7���            6��    6��{@��    6��    2��%0h {@��    {@��7���6���7��    3��q{@�ε�ϸ2���        <��}    B��k6��tFx(�D��>��>�nA0V�            ?�&%    :g�    :c=>��@�`�@�~g@O}U?��@>ю>-�7h�P{�������+�ʣ�����P{�p��TI��=�f�+v��lT?Ɉ�:c;�>�
x@��@��@O|�?��@>ю>��7h�P{�������+�ʣ�����P{�p��TI��=�f�+v��lT<Nǡ    +c<2���7/�?<�F=�ރ7�j�1��2��<�6                                                �!���C�������)��q,��c���4G����̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6�! 5Ii?��78
�a6�a@�@.c"@�@�@�=�4c7�%�6&bF�iH7ճ�>���B� 5ʮ >�f�7ճ�JTO�4$7��CI.rD�ZD��D5{�C�k�C��BW�A�2@���                                                A�.rCa�ZCH�B�{�B+k�A���@��@	2?��                                                                    E�&�G��GkeF�<lF6q�E��D�k�C�a�C�K                                                C�!GER�?EC��E � D���Cђ[ChB,�6A=׳                                                                    E�nLF�F�x�F|
4F]{�F/&E�?E���EG�{                                                D)ĔDĀD�#�D�M=D�D~�	D0�
C�\�C�"�                                                                        5���G[��GOh3Gv�F̤jFyC�FdE��
EP0                                                A�i�A�`A�X�A�
�?/;D                                                            ?Ba@�TjAA�A��AB2ǯBfb�B�reB���<�<�<�<�<�<�<�<�<�<�<�<�E��iE���Eb8�E��D��DV�[C��C�`�                                                {@��{@��{@��{@��{@��A�44B��P    �f�B6�b;C�@	�    ?~�F?~�F{@�ξƾ�C��C{@��C��@C��4C��C{@��@lL�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��>�%�D�	G��D��DK|"DK|"@�~�@�~�Fj�B?[B�)�C�gxFKoPFKoPD��VD��VFk��B?                ?$2cC���C�}�C���?   C���C�U*C�U*C�_BC�n�C���C���C���C���C�4�C�huC���C���C��C�;�C�x�C��%C��hC�0�C�]VC�z�C��xC��C�S�C�.�C�)?C�,�C�cn{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C���>���>���>���>��>�-�>�G2>��C>�y$>��>��L>�W>�L�>�j&>��a>�})>�)�>�^>�,�>�� >��@F?���        @uͺ�N����^�q���,���+�����oPc���(���)���(�Z^{        =M��A�*&    ��4���4�{@��@TV�F]��6�/�7�[        0�g�=�7:=�U�<�->��?B��?H��?Q��?_�}                                                ���a5'�A�  ?k�    <#�
+:Y      �b     �`6Z1�5��@3T�5S�r63(5kDS4@D/3���                                                                    5�Ί4��2�wE4���5b#�4���3r��2Ҿ�                                                                    G�hF�Y�F\P
E�X;E�!D9�eCJՖBO*                                                1�{�1-28.��c0߇1�&60��/�$/0R�                                                                    1�x0Z�,.L0,�0���0��/ G�.^�d                                                                    4�$4P��2� 4,��5*��4xaj3r�N2��                                                                    5��4S}2�84+~5Ȕ4�W3���2�k�                                                                    2��12P�/�*1�ރ2�F�2a�1n��0��S                                                                    55�z4���2/�~4Q�5=.
4�3���2�J�                                                                    6N(�5�Z�38L5;�n6/�c5�]54q313�)�                                                                    3'�y2���0K�2rH�3��[2��v1��1<��                                                                    6N(�5�Z�38L5;�n6/�c5�]54q313�)�                                                                    5��}5��2�#�4�o5��4�<3͓.3*$�                                                                    2���2'	/��1��e3 p�2_.1cQ�0��                                                                    5��}5��2�#�4�o5��4�<3͓.3*$�                                                                    �����ግ1++@3ti4���3�h�2�}2��                                                                    ���9,�;�G^;5ry:���:)i9n]�8���                                                                    3XM�2��0<�[2X��3;��2�"&1�3W0��                                                                    4�7~4.�v1��3��S4�ZL3͇�2�1�2�                                                                    �<�γBn��U�82ϴ,�<�e�R�LR���H                                                                    ��J���
�*�������!�V����K�H�a+                                                                    3�O�3.��0��a2�m�3�x}2�'�1�[�1:�                                                                    �I���"�d���ʳ��'���^�0�f�,,p���a                                                                    ���g�����]_��u���X���
��J�@.�                                                                    ��j���IY��5.���Y��`E��?���nӰ^�                                                                                                                                                                        (��-N��                                                                                            5��k5�B&31��5>�&62�5��t4j�y3�                                                                    2%p�2��/���1�E�2�;�2��1�0��                                                                    4�K�4��@2.q4=.s58<�4�I#3��e2�ا                                                                    1"�1�
.�S0��|1��1�V0�e/��w                                                                    5��S5}��3*{a5;��6>mv5���4�#4
�                                                                    2��2 E�/�^1�{!2�v<2�1�;0�lY                                                                    5��5<42�h�5x�6&�5:J�4%�3��                                                                    3.��3��1;�3*�d4 ,M3T�2=|�1�
S                                                                    3�@�3��k1���3�_4��4��3�2z\-                                                                    2I�1��/��2

3n�2,U1�0�                                                                    5:�35�?3?�56�K6+K�5c��4J�V3�U�                                                                    3U��34��10I 3P�z4C�^3��2g�b1�a�                                                                    1��A1���0
r�2���4V�3�K33%{N36��                                                                    0���0���.��N1b>2�D2���2��2��                                                                    2f�1��0)6�2��4 ��3�?w3JA_3_R�                                                                                                                                                                        +:Y      �b     �`6D�V                    >S�x    6�1�3�FB@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��`    @��@    17:07:15        G~ @��@    @���    +:�      �~     ��3�=6���B'�fB'�f7CUJD"g�    
+p�,	��6R�7?�  ?�4���1V��2�R�0�\�0g{)        \\�@P]9��2j~�/��;6�'="�t=;�37�P�2=�~:.��<!�<�̋=G�@=�7=��=��=��=�:�=��J>)�>�'9M�2~m+런1V��7���0��i0�*�8`��@Z5T@l��@b�@b�E@b��@c@c0 @cO�@cm#@c��@c�@c��{@��{@��{@��                    E���5HU�G�!�F���F]kE�F�E{D9�CK[�BO��                                                A��C!��B�]�A��AJE@��e?�W�>�Ԍ=�UD                                                                    EL�;@�MF)��A���2?�FiBB5�!<�D/-cZ                                                    {@��@E�&@E�&@h�7h�V%���Q�ԨpAe/4_����>쏨
+p�?a�
/V��2�%�    >�y�?U��1ױ,?��K���~����<#�
    B4�B4�C���C���C�?s�V?*lX3�-C�P6-�6��A-�fA�~�>�y�B"�AF��B^��@���A'BwoA9|B4��ٖ���!��I��ٖ    ��ٖ�5&b@5Eb?}p�?}p�A�@+uA8�O?��K@z�A&_@.�
F���2?�0�-� 2�@�7&-    4�?F���F�'�G!ũG)C ?��0@�(�C9�@��H?�!�?O��?+	?>���>�Q�>�i>��*8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�*    �d��d�{@��6`36�A6~.�4f��6�6���6w|�5om�4Ȼ�                                                A�����!���!�{@��Cw�A�7BG��@b�=�8�{@��{@��    C�)EBR�35�<�    6OH7��@�D3B6�S�    >v|/C���7��9B	zDAXvC.1�B���B,��Au�8@̐�@!�>���                                                ?I0A=g�@AY�?��N?e	>��P>M=_r�<A �                                                                    D�FT�	Eh!�D���DyCK=B�@DA�v@���                                                @���B�'�A�ɭA'1�@���@�?�PL>�C=�8d                                                                    C\BE�W�D���C�,�CO��B�U�B��A&	�@!                                                @�UB! �AH�6@��V@'�<?��?��>9��=0n�                                                                    5��6�QA��>J3(BrN?(H�?�Ұ*��(��1,e+]��*��(ڌA+~+�+���)��7CUJ                        ��R3�-2Dڏ            6�R�n����R3+4�                        >�S��t�q                                                                2�8�            2�8�{@��-,Z>L��>L��>L��>L��?X�?99�?'I?>}�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Bx��G��$>+��A�C�                                                B�                                          B��                    A�                @�p�    o�P�    3�:x2���2�A>A6�ˢ>��5��;7�f;9�                5�W    3h�2)ݻ        {@��3��5�5�/E���{@��{@��{@��1���7���            5��    5��{@��    5��    4��6���{@��7f�{@��7��]6q��7I.�    5l�{@��4�ZY4�(�5��05��0<���    B���6��jFR�D�?�l�>��A���            ?� =    :w�r    :E�U>��GA�@��
@g�h?��>�j>;��7h�P{�������+�ʣ�����P{�p��TI��=�f�+v��lT?�^:E�a>���A �@��z@g��?��>�j>:k#�7h�P{�������+�ʣ�����P{�p��TI��=�f�+v��lT<���    *�\�2s�/7��E=�;�>r�
8*�F1��1ǡ>:�t)                                                �!���C�������)�Đ3�čs9�k���3�X̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7 ��5#8�?DST8k��6{�@� �@��A�KA&�M>]m84PK6H;�G'��7���?
�C��?��?��7���J��H5?L7��DCIK�D��mD��)D5�nC���C��BW""A�/Z@��@                                                A�K�Ca�mCG�)B��nB+��A���@�""@	/Z?�@                                                                    E�#�GيGOGF�QF6�0E��D�m�C�icC�                                                C�ER��EC\
E �"D��iCѡAC
B,��A=�K                                                                    E�mF��F�s�F|qF]z�F/�E�?�E��EG��                                                D)ùD�x.D�D�I�D��D~�+D0�C�^�C�#J                                                                        5���G[qGOH.G�qF̵$FyI�Fd�E���EP                                                A���AʮA�bA�h�?AQ[                                                            ?2~@� 0A�~A��B-��B^�B��B�Y�<�<�<�<�<�<�<�<�<�<�<�<�E���E�t2EbI E�D��\DV�~C�DC�a]                                                {@��{@��{@��{@��{@��A�/�B����)��fZ]6�0�B��?�t�    @jF�@jF�{@�ξ\��\�C���{@��C��<C���C���{@��@h�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G~�?1�D��G�dD�aDO�cDO�c@��@��Fi)'B=dB�nC�_�FKm�FKm�D��PD��PFj��B=^�                ?���C���C���C��k?   C��GC�yuC�yuC�|�C��kC���C���C��'C��\C��C�:�C�a�C��C��C���C�$�C�`�C��gC��vC��C�:�C�W]C�k�C�V�C�1CC�)�C�,�C�~t{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���D�T>�;=>�fN>���>��3>��>���>�o�>��[>���>���>�1>�l>�(u>�"�>�{�>�h>�wu>���>��>���@"l�?���        @�ga��"�����
�Ğ�.���$��(�� �n3����������~2$        =��A�5�    ��Em��Em{@��@r�ZF]k�6���7�n        0疕=�S�=�
<<�r>��u?;�y?@�W?GX�?Qmi                                                ��d>6���A�  ?k�    <#�
+:�      �~     ��6��5��3�57 =6N5c�s48X�3�AV                                                                    55ف4�T�2�Kc4gQ%5$�v4�Դ3h�22�*                                                                    G�!�F���F]kE�F�E{D9�CK[�BO��                                                1�1%�1/zh0�y�1��i0��/¿,/:�                                                                    0�0Q|�.F�v/�d$0��A0�).���.;<�                                                                    4'��4F��2V�'4�4��4p�3hh2�;?                                                                    5	�D4wA�2�l4?�:5y�4���3�/�2���                                                                    2���2xy01h)2�_2���2i-�1�_f0�N~                                                                    5(�4��2���4ju+5&͟4���3�:^2�                                                                    62>+5��z3�
5(;/6	�5���4l=:3�9�                                                                    3��2�;�0��2a�^3KE�2���1�n�1#�i                                                                    62>+5��z3�
5(;/6	�5���4l=:3�9�                                                                    5�^�5�	2�+4���5U�4��s3�)P33c                                                                    2�Í2'٢0H1�ו2��}2\m�1a�0�$I                                                                    5�^�5�	2�+4���5U�4��s3�)P33c                                                                    ��m�)�1K��39n�4=��3���2�i2eh                                                                    �&�9�Y;��w;5�2;�::LO�9��8�WJ                                                                    3H+>2�F0�X�2r%3%x2�g21� 0�"                                                                    4��4+Y
1��3���4SA�3ȫ�2���2T�                                                                    ��M���бX����b���>�]�(�Cl���/�                                                                    ��橯�}��m�3��vű��ݱ��	�Cz�<�                                                                    3�M3+��0�	2���3T�2�~61�@1��                                                                    ���#x��	�έ`��`��*�ʳ$�/��.                                                                    ���Ѱ�ک�����SΔ�9>���A㯺�Ư!(                                                                    �K��z	H�����Ɉ�N� ���~��f��                                                                                                                                                                        '��-
�T                                                                                            5W�5{�{3x�5$;�6*�5y3�4`�w3��                                                                    1�z�1��w/��Y1�Oy2�M2�0��Y0[æ                                                                    4R��4v��2s�4#e5��4��+3v��2Ի                                                                    0�0�s2.���0�v1�ά1��0��/p�6                                                                    5NU�5q��3m��5!��6
&O5��V4�0�3�s�                                                                    1Н�1�B�/�L1��s2��2m.1 �0�c                                                                    4Ȯr4�!#3/א5 �b5�h�54�4� 3���                                                                    2�Y�3�81H�\3"�3�w�3M�25<�1���                                                                    3�*�3��y273�4�^�4��3 %�2R�                                                                    1�U`1�NA0"d�1�˃2��2&K\1t.0pv                                                                    4�F�5i�3V�5Zx5��B5\�4Aҏ3�߆                                                                    3(�3+�a1u��33�4&3{�.2]��1���                                                                    1�\�1��0@�2p�3��+3�n33B�3S�                                                                    0iP0�s�/�1B��2���2���1�ƨ1���                                                                    1�q\1��0k��2�?>3��3��3An<3;f                                                                                                                                                                        +:�      �~     ��6Bj/                            6�fa3߷@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��@    @���    17:07:18        G� @���    @�Ӡ    +;!      ��     �p4�]6��B'�fB'�f7fD"g�    �
�,��6)�|?�  ?��4��1-J�2��0̀�0�2        ��]�$9�22��0݆�7=\�=�=3z�7�(L2Bs�:4�F<!��<��p=G��=�:=��=�	4=��C=�;�=��>*>�}9)��2��,nx2�8��81[�50���8\��@Vj@l�@b�
@b�R@b��@c�@c/�@cO}@cl�@c�_@c��@c��A�!�;���{@��                    E�s�5��IG��F��'F^�E�"�E �D:&�CK��BO��                                                A��C"�=B��/A��AJ�#@�ӭ?ĪL>�i�=���                                                                    EL�;@���F)hlA��-2��@Fh��B4�C<�D/.�N�                                                    {@��A�|	A�|	@d�i7��%V<"$x`�� _/)�e&sp�>�>��
�@���/��;2��    >�y�A��1�6|A��A}'��ּ=��    B��rB��rC�;C�;Cx�i?kb�@�[3�C�0�5��6��XA+�:B���>�y�B�	B��CH.AY%lA�HB�6iA�r�B��� r��\�#~��� r    �� r�A�.�?m,�?m,�A�T�@��AȌ@+��@��-A���@[6�F��$2��@0X -+NV2�h�7�?5 4M2�F��NF��'GsL4G���>@��A��B�/A�?�	?O��?*��?��>��>��.>�hw>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�(>    ?�|?�|{@��6\X7_�6��H5�;�6͗97P��6mlu5e�h4�8�                                                A���������{@��C1��A��B��@V�;>]��{@��{@��    C�ݛBRtQ5�)�    6L"X76
@�*�3�H76�    <��nC���7�ȬB�eDL��CMkiBłoBIZ�AxQ@�թ@"M?x�                                                ?,DpAH�7@c�O?�p?��3>�18>��=u+�<T�b                                                                    D\FT{�En��D�fDkCLlSB�8�AƇ�@���                                                @�'�B���A��A1�}@�K@v?�E�>�+�=�%�                                                                    C_�E�j*D�ѨC�#.CW��B��B��A*�s@$�%                                                @
*rB"��AR�!@�� @1�S?��?�]>@y0=6[�                                                                    5�6Oo$AT�>DQ�A���?#el?���*�^d)�o,h�0+�[*�^d)���,S�N,�s�,G�A)xd 7^x3Z#�1NF+��            7$2a3��v2Y��1N�1NF+���$2a����7$2a3q��    &�]�3MA�3M�.��    =��v7[Ha6�w�1'O�09,� �                        *at�1�1[q-�M    1'Sg2�*            3��n:g-BX�>L��>L��>L��>L��?|3'?6��?�?>7?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C�G�*�>+F�A��5                                                B�                                          B��                    A�                @�p�    �U��O3rA�4�g2�$�2�ic>��7�)&?~�6�h;M  ;H�U                7w~    5;3�uZ    8��{@��1p�$6qC�8E����{@��{@��{@��3�s�7/��            8���    8���{@��    8���    6O�`7��L{@��8��{@��6�d
6��6�J    6b�V{@��5�'�6e�6�D�6�D�<�G>    B���6���F0�DH�@��?���B���            ?�+�3Z#�:�Y�3Z#�:�C�?"$b@��F@�<�@�\'?̴�?M�>X�ċ7h�P{�������+�ʣ�����P{�p��TI��=�f�+v��lT?��L:�B�?"!�@��l@���@�\?̴�?M�>X���7h�P{�������+�ʣ�����P{�p��TI��=�f�+v��lT<��    .��3';(8*(�>d{>=I-�7@��2�w2�5�6�tg                                                �!���D;���L{�T������ĭi�ĕ8��q�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7RR�64>&N�8�6�8ݾA�� @�Z�A�KsA��_>
}Y8w��7V�IF�)�6�{�>{��B���@g�>���6�{�K�d8��B7�&�CI��D�'DǸDD6N�C�rC�BWm�A���@���                                                AɅ�Cb'CG�DB�N�B,rA��@�m�@	��?��                                                                    E�xG��G3qF�e@F6��E�+D�qpC�r�C�C                                                C��ER�EC3�E�D��CѠ>C�B,��A=�                                                                    E�k�F�6F�o<F| �F]zgF/=E�A	E��WEG��                                                D)��D�pND�(D�FYD�3D~��D0�{C�`C�$7                                                                    1'Sg5��G[OGO(�G�
F��^FyJ�Ff E��EP	�                                                A���A�	*AגzAEӘ                                                                ?]@��A*�4A�h�B(�
BX9�B��B��<�<�<�<�<�<�<�<�<�<�<�<�E�l�E�]EbYuE�TD��zDVۂC�	�C�bn                                                {@��{@��{@��{@��{@��A�/=B���U�o�fN6���BÅg?��    @��@��{@�ξ�������C�^�{@��C�&*C� MC�^�{@��@d�j{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��? �D��G�D�vDS�'DS�'@�Ņ@�ŅFh�KB;��B�FC���FKk�FKk�D��#D��#FjfB;ޫ                @&EC���C�(xC�<H?   C��C�C�C��uC���C��'C��:C�ĳC��C�
�C�-�C�M�C�q�C��/C��\C��C�(TC�a�C���C���C��C�(nC�M�C�S�C�3�C�)�C�,�C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��C�uJ>��h>��{>��>�v|>�{�>���>�y!>�V�>��>���>���>���>��>�k�>�X>��>���>��>��u>��@Vʣ@!        AJQ���,�������V���{eD�z���zk�v��ø��ø��ø�ĉ-B        >k=�B    �\�ܿ\��{@��@��cF]6�6g�?7�N        2 N>��=ރ�=�
�?+O/?5�?9G?>��?Fd@                                                ����7�,�A�  ?k�    <#�
+;!      ��     �p6�,5�g�5S�5�c6W�5X��4.�3���                                                                    5��4��X42��5��5�Ԏ4���3[�B2���                                                                    G��F��'F^�E�"�E �D:&�CK��BO��                                                2ss1U�0�M�1w;�1�:L0���/��q/>r                                                                    1E�N0� /���0�%�1�#0�'.�B�.1&t                                                                    4��4~��3��	4��5L~�4d>�3[L�2ǔ�                                                                    5Q�4��)4[�5	�5_��4�CQ3�ӿ2�sl                                                                    2�Y2gF�1���2���30w22j6�1��w0��                                                                    5�4�F�41��5+(H5��Y4�c3��\2��h                                                                    69�6 �5|�5��6X�C5wN�4c��3���                                                                    3Tr,2�p2\43dS3��2�>M1��51�_                                                                    69�6 �5|�5��6X�C5wN�4c��3���                                                                    5��5B��4J�f50O�5���4�<(3�~w3NE                                                                    2àP2bN�1���2���3#x~2U�1\�'0�D_                                                                    5��5B��4J�f50O�5���4�<(3�~w3NE                                                                    ��-�q�2�E�3���4��~3��2��2S�m                                                                    4��9�I�;�x�;)�@;�A:h�b9���8���                                                                    3���3�p2>��30��3�e'2�kT1�0��                                                                    4�R4a��3\�O46�4��~3�~�2�`%1�                                                                     �|���Bk��̇X��^�N�9�R�/�8
豗?                                                                    �6/� |�� �<�T�A�̱�Kɱ7��2f�                                                                    3�|3b��2^��37��3�yo2�t�1�	�0���                                                                    ���V�F�l�������_�B��"������bb                                                                    �	Oϰ���#�㰰.D��[��9                                                                    ��vO��E~�(���&�+��n߱������u��\�                                                                                                                                                                        *�VH0P�                                                                                            5�y�5�s�4��5�~�6Uu�5l�T4Tq�3�^�                                                                    2^��2#<�1m72V�2�F2 n0�O0O�V                                                                    4���4�)]3��4���5\��4x��3h�2��                                                                    1Z_1�0ha�1U)1� B1��0��/c�2                                                                    5�}�5��/4�35΀#6dM�5�<�4}qw3��#                                                                    2U�82��1c�V2S�72��?2'�1�0w�                                                                    5M4v5�h4��R5���6'�$5+�4X�3uw%                                                                    3j�34N.2��3�\�4?��3Ct52*�c1�D                                                                    4%�y3��13��4�zW5�`4
3#2�^2F[                                                                    2=��2��1�F�2�g>3�n2�M1	�0b�D                                                                    5z�W5@ӛ4��+5�_p6Md5Q�46��3��                                                                    3�QV3\_�2���3��[4jXs3n�]2P�@1�o�                                                                    2��1��1��E3R�4/G3��3�3ā                                                                    0��,0Ó�0���1��02�	,2��1�٘1���                                                                    24�(2�1�I�3;d�4@V3ī�36$�30�,                                                                                                                                                                        +;!      ��     �p6?vp&�]�                        6�+�3��F@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @���    @�Ӡ    17:07:20        G� @�Ӡ    @��`    +;�      ��     �59W�6��B'�fB'�f7���D"g�    '�g�+��c6Pl%?�  ?|#<4�1�3�%�2�	�1Cw�2dP�1�h�0���-��S)�M�>$��=�c�=P�<�|B=�"<���7�'2GN�:9c<!�l<��@=G��=�<�=��=�
�=��T=�<�=���>*i>��@*h�?���?A�>��->���>�[�>t�>6f@X��@io�@`:�@`�g@a	�@a}�@a��@bnw@bԫ@c$�@c[|@cw�A��<���{@��                    E���6�[G�`�F��F_yOE��8E8MD:yCLw�BPn�                                                A2�C# B�|�A�[AK۵@��?��>�n=�I                                                                    ELy�@���F)AvA���3���FhvNB3�6@���/�
�                                                    {@��B��B��@`\&8��$�~\$��ħ�b�_��&bA>�{�'�g�@|�f0�4�No    >�b�A�^�1ѝ8A��A�4�AAo�>�t    B�x�B�x�C���C���C���?+��@�
 3<MCܕ6�N6�<�A)�vC?�>�b�C /JB|4rCE�TA�TEBG' C"��A�n�B�"A���.5� ���A��    A��BX�@�Ȣ>�5.>�G�A�}sAh�A�@:�y@[�xA/�@��F��x3���1+q�.z�3�+�8_�6545O�^F�q�G1�G��G��]>"h@w@�A-t5@�J?)��?(�:?W�?�>��V>�/x>��>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���MŰ�    @�
�@�
�{@��7l�a9�8�[84:o7�+�7U�6��5��4�U                                                A���@�@�{@��BcdB-#�B��?@I��>�4{@��{@��    C�y�BQe@5�W    6H�e7qQ@�h�3
�s7$��    =�|kC�<7���B"�mDS��C`�JBܢ�B`m�Aw�@�R@-+?��                                                ?8vOAO�J@yW�@��?�v;>�	�>%��=�U(<e�                                                                    D��FT�Es�mD�GD	��CM@�B�Aʽ@�E                                                @��.B��A�xIA:VU@�h�@ W?�"b>ʭ�=��                                                                    Cb��E�dOD��C��PC^��B�|�BudA.�v@(�                                                @��B#��AZC�@��L@;[N?�ZV?��>F�
=<)�                                                                    6��J7w�%A��>?{�A���?`�?��+�
*=��,�Q�,8�+�
*|��-Y֒-G�,�_�)ڻ�7{��4�d 2�$-���            7�.�4��2�7X2�M�2�$-���.��]Ɠ7�.�3��    *���4��$4���/�|    =
�/8I�7��2�&1bq�.��b                        ,ɓ�2pU�2i��/�E�    2�2�2�)�            4�V�: �-n��?%��?�?��>�5z?n7,?Lũ?+��?Rh�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CAŷG��b>+��A�F�                                                B�                                          B��                    A�                @�p�    -�Y�)]h4�{�5��*3�U3�F�?[��8^/@��7�}�;~��;w�    /\_        7h��    5@�T6�&    7��v{@��/��5��6���ֵI{@��{@��{@��6�y�6���            7X��    7X�4{@��    7X�4    6=&�6��{@��7��{@��6��5�y�5U��    6[��{@��5��.7>�6�C66�C6<��/    B��o6�JFZC���A�"v?��,B�v�            ?Ǩ�4�d
:c�4�d
<�9>�	B@��<@p��@q2e?��o?J<>ozE�7h�P{�������,�ʣ�����P{�p��TI��=�f�+v��lU?�GY<��>��@��
@o��@o�k?Ġ�?>f[��7h�P{�������,�ʣ�����P{�p��TI��=�f�+v��lU<�`�    ,��N5�0!7��c>bfA<�p�<��<�x<r�#<��                                                ��8��,����-��Y=ļ�qĪ��ĔH�r�>̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6��a6'��=�S�9k8�.�A��A	�FA��;A��u=v�v7�=7RpE6=6M'�<}�A(�V>�m<�B6M'�Ju�;7�;q6��CJ�D�kD�TD7�C�½C-�BW�~A�68@�OH                                                A��CbkCHTB��B,½A�-�@��~@
68?OH                                                                    E��G��GuF�v�F6�E�kD�v1C�}]C��                                                C�-ER]�ECE�D��Cѝ�C$B,�DA=�`                                                                    E�kF�	yF�k�F{��F]y�F/�E�B'E���EG�0                                                D)�ED�j�D��D�ClD��D~�D0�KC�bXC�%:                                                                    2�2�6ݞwG[8�GOG��F��FyKoFg�E���EPO                                                @ݽjA.�A\�mAdG�;X��                                                            @� �AhV	A��tA���B,_B\;VB�΀B���<�<�<�<�<�<�<�<�<�<�<�<�E�]E�O�EbkE�4D��&DV�AC��C�c�                                                {@��{@��{@��{@��{@��CrO�Br��- ~��gb7��4B�l�?�֘    A�VA�V{@�ξSa�SaC��{@��C�*�C�x�C��{@��@`\&{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��>�6D�(G 1D��DVmDVm@��v@��vFh�pB;�2B7��C��FKl�FKl�D��|D��|Fj1(B;��                @Q�QC�txC�'HC��_?   C�˹C�C�C��)C�l�C�C�C�1GC�&8C�#EC�*�C�9�C�K�C�dFC���C���C��C��bC�4mC�ksC��iC��C��uC�-�C�K�C�5fC�*nC�,�C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�ޘC�cA>�`
>��>�ix>���>���>��>���>�O>���>���>��>�z�>�Q�>���>�>�Q>��>�>���>�sI@#��?��j        A��O�W�O`��O
gř���ҦP��� ��sE�CLr��&���&���&�đH         >��B2�    ?�<2?�<2{@��@p�\F\��6�G}7|�        :�Ǖ?#�?�>�?`?��?8�"?<�\?C�?M^�                                                ¤��7���A�  ?k�@;��@�M+;�      ��     �8�=�8n�7_�H6�ɐ6\@5��@4z�3��m                                                                    8ȋ7H 36�.H5��5���4���3�C02�S                                                                    G�`�F��F_yOE��8E8MD:yCLw�BPn�                                                4r(�3�Y92�'P2)�1耗1L�0\O/&ō                                                                    3��+2�c;2&N1V�1��0/�t/'1<.R��                                                                    7�6�+�6 ��5���5Qf�4��u3�2�G�                                                                    78t�76x�5�;.5c\�4��3���31�                                                                    4�1�4�+�4+I�3��34�2��b1�۞0�2                                                                    7ar)7"�]6��6N�5��4�N�4��3                                                                     8Z8-8,$�7Pߥ6���6^��5�)4� �3ȸ�                                                                    56��5$�X4z��3���3��@32'�1Cf�                                                                    8Z8-8,$�7Pߥ6���6^��5�)4� �3ȸ�                                                                    7���7�q86���5��:5���4�y<4�i3+��                                                                    4�o]4�2�3�ED3R�3(�2���1�X0���                                                                    7���7�q86���5��:5���4�y<4�i3+��                                                                    ����g/4�S4���4�C�3΃�2���2v��                                                                    5zZ�9`��;cD�;	i;	4�:cP�9���90�                                                                    5�(358:4�)4M&3�U2Ց2��1]r                                                                    6�`6��m5��r4�	4�2@3�-�2�zj2�D                                                                    ��$����!#��;��Sc*��}��'��q�                                                                    ��t��h���0xֲ0�F����q��Ѓ�S�                                                                    5�hW5���4��44�3��2��>1���1�                                                                    ���Ѷ��ֵ�+T�5r3����F��`�<���                                                                    �T:�m��g��������֪�����5'                                                                    � �ʳ�̢��Э���粮B)�� ������                                                                                                                                                                        *�B(//�                                                                                            8*a�7��w7:�6�/96Z��5�aR4�K�3�8j                                                                    4�C�4�3�*3�[2�o2{�1-�R0w0                                                                    7&��6���66:;5�']5bR�4��Q3��2�F�                                                                    3��`3z\�2�>2��1��l1$<+0>f�/��                                                                    8#p7�y\72o�6��6i�i5��P4��W4*�                                                                    4�>�4u'�3�i3��2��i2+��1O(0�q�                                                                    7�=�7u�o7�j6`XG6+� 5O�]4Vk�3��Y                                                                    5���5��5j�4�2r4D
n3m��2u[1�j�                                                                    6�Y6F�;5�K�55I�5
�[4'�m3-D�2kV8                                                                    4�#�4c3�{i3O/�3j�2?�X1F�0�zi                                                                    7�g�7�E7�6��6Q��5~�4��3��4                                                                    5�	5���56��4��o4o��3�&�2��1�e�                                                                    4c��4<��43�3ѪE4 ��3Çh3U��3+�1                                                                    384,3|�2�pO2�m&3�]2� �2,�2
��                                                                    4�M�4f��4/�4  �4Dp'3���3���3Q��                                                                                                                                                                        +;�      ��     �6<c�*���                        6�ٟ3ڇ�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�Ӡ    @��`    17:07:23        G� @��`    @��@    +;�      ��     ��5���7H?�B'�fB'�f7ĩHD"g�    (ֿ}+��7�?�  ?xa�    3�sF(��/������        /
�T*���>�e>9�.=�z=@�<$��6��,/�g2T�f:7�<!j�<��}=G��=�>c=�)=��=��&=�=Y=��_>*�>�@՝@���@�K?�A+?�6�?���?c1X?0Ҏ@hO�@h^k@^T=@^}�@^�@^�0@_*�@_y2@_�@`�@`c�@`�w@�!�=>�
@���                    E���7m�FG���F�hF_.�E�ӿED:O"CLWdBP\�                                                A��C"��B�5�A���AK��@��c?�Ֆ>��B=�6@                                                                    ELd�@�ȪF)06A��23�jwFhH�B3��A9mv/�ϛ                                                    {@��Bd��Bd��@_�8WI����    ���m/ID�#���?   (ֿ}A!�<��%Z1���    >�'Bx51��A��A-�A��        B��{B��{C�O�C�O�C� ?�@���3j7C��6�s6��eA)��C>��>�'CR�B��C`R�A���B��?C<�?BWrB��BT5�?;�"    BT5�    BT5�BQF�<U�a    6��A�(`A4�$A��@a��@)��@ͱ�?�8F��3�jw1}�.�"�4��8���6��5�k�G%�G4z�G�,G�=�=�r�            >3 >�:�>�/�>ͣ/>�Ȉ>��g>���>Ë8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    A<}A<}{@��7���9��e9hhe8�T>8IV>7��6��G6z�55��                                                A�{�A3�A3�{@��:\�GB|�gB�,'@@��>�X{@��{@��    C�u�BO�6r��    6G,7@�}�3	d�7>^    =ȸ�C�H�7���B�PDG�pC9�0B��B:�|A_6@��/@j�?	�                                                ?!7dADE�@N6?�K�?y=�>�t�>�I=fw�<Rh*                                                                    D(YFT�En��D��DIhCJ�
B��A��@ý�                                                @�a�B��nA�zA4'g@���@I�?�@�>Ƌ:=��                                                                    C`#�E��D��2C�^�C[#B�:�Bk	A+�g@'                                                @
a�B!ǋARi�@�u@6�E?�Ma?��>B��=:4                                                                    7r��7�yA�>>��A��?��?���,B�*��-�,�w�,B�*��J-��0-�F�-/O�*4)@7�Ձ5c92�`1	��            7��5Fe2�ǰ2Ҟ�2�`1	�����1���7��Z4t��    /D�<5�C5�T3=�    <��8kY8F�4���3
%�2 wa                        1CD�4Un>4%�f3>�`-���4��I2��            5�];�aA-c4�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C\�G��>+��A�Q�                                                B�                                          B��                    A�                @�p�    .�7�*�s5&��6�3��4�?�ts8F�8@X�[8Â;�M�;��5=�V            7��1    2��6���    3'�n{@��            �L�{@��{@��{@��7
j�7��p            6�ˇ    6���{@��    6���    6���    {@��3'�n{@��3'��3�1�u�    6���{@��6���7z�7�
~7�
~<��B,m�sB��6�3E���C���A�2�    C%L�            ?qP5o�9�@S5o�<ш<���@U�?�Mw?��@?9VX>��^>P�ԋ7{�P��������F�ʣ�����P��p��TI��=�}�+vӊlg?��<с<���@��?��%?���?_�>~|=�-�7{�P��������F�ʣ�����P��p��TI��=�}�+vӊlg=�+�/Aȋ*��6`Ѕ7��<#�;�Q�<��L>/��>b�=��                                                ĠH�ĝNė3Č�g�|f��T���$J��שV̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         :\�G        {@��{@��{@��{@��                /�k=                /�k=            CK�D�GD� �D9��C�R�CRNBX�2A�E�@��                                                A��Cb�GCI �B���B/R�A�RN@؉2@E�?�                                                                    E�#ZG�JG0pF�p�F6�sE�D�m#C�x�C�.                                                C�&ER��EC//EXD��kCы�C	�B,��A=�                                                                    E�lF��F�n�F{�F]zpF/TE�A�E��5EG��                                                D)��D�p�D�KD�E"D�:D~��D0��C�a�C�$�                                                                    4��I7`רG[T.GO0�G��F�ݽFyE�Ff�E��EP
�                                                                                                                                @_�)AkI8A�̸B ��B4� Bi3AB��'B�R_<�<�<�<�<�<�<�<�<�<�<�<�E�q>E�e�Eb{�E��D��DV�\C��C�c�                                                {@��{@��{@��{@��{@��D�B��.f���f�8��C/��@+��    AH��AH��{@�ξi���i��C�"?{@��C���C���C�"?{@��@_�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�<řD���G"D���DS�DS�@�u@�uFj%?B>iX    C�'}FKu�FKu�D�ւD�ւFk0B>d                @��C�	�C�]4C�T�?   C��lC���C���C�(�C���C��KC�Q�C�1C��eC�]�C�C��C��OC�}�C�[�C�E	C�=`C�G9C�aC��OC���C��CC��C�?�C�6#C�*�C�,�C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C�G7>��5>��>��>�QN>έ>ɭ�>Ċ>�^l>�.�>��q>���>�1->�B�>��}>���>�C>�}=>���>���>�d�@)�i?��        A�sǉQTǈ�ǈ��ž��������R�~������������Ğ@        >���B1-�    @�y�@�y�{@��@n�F\��75nn7l        ;�4�?<�z?=0{?>�Z?AH�?EX�?K��?U�^?f@�                                                �a��7g��A�  ?gF8A	��A z�+;�      ��     ��9G8�[�8ٻ7eV�6���5ٗ�4�l�3��7                                                                    8{��7�Ģ7*Ve6��i5��5	mb47=3>>                                                                    G���F�hF_.�E�ӿED:O"CLWdBP\�                                                4�U�4&6�3�vl2�N2)	r1e��0|�y/�f`                                                                    4��3Q�\2���2�1U�20�0�/���.�=�                                                                    7g��7Fƹ6à�6<�5��4���4�U3<��                                                                    7�c�7Zig6��b6}��5���5� 4JgQ3=>�                                                                    5$ |5��4�y4;��3m�12�7�2+'$1 l                                                                    7��7�yM7.6���5���5#� 4wa�3gL�                                                                    8��g8�jb7��7Y��6�J�5���5O-4r                                                                    5���5�<59J4���3�j3\�2�;�1�W�                                                                    8��g8�jb7��7Y��6�J�5���5O-4r                                                                    8�w7���7BG�6��6z!5O�b4�Ah3��N                                                                    5Pr5��4�*y4	�3p�(2ִy2R.1Hy                                                                    8�w7���7BG�6��6z!5O�b4�Ah3��N                                                                    �ͣϴ)��5���5^$�4�42�R3|i`2ȑ                                                                    2�_�3BT:�i:K�j:��9�٥9\�l8�lf                                                                    5ڶ�5��5ʥ4���3�y33�2i�}1Z                                                                    7(o�7��6S-�5�π5=d4@�u3l��2o                                                                    �)��´�Ec����I��}ֲ��!���                                                                    ���u��ز�O�����)"�G�����                                                                    6)y6�S5U�4���4�`3B�$2p&G1r8j                                                                    �6�7�{����9�+a�Uུ��/��W����                                                                    ���γ�3ҳ����lѱ�l�1��uԯ�i                                                                    ����x�H� �A��������.X%�PF��j�                                                                    'ހ�(s%-�=�-a��.�ۗ/�7/��P/4                                                                    +��,1Ե                                                                                            8�h8|i�7�E�7Pw6�WX5�r5��4-��                                                                    5'4�454e�$3գ�3*?h2���1�1�0�T                                                                    7�a�7wD�6ޤ86N��5��4��A4!��3>��                                                                    4�u3� �3a�2��20*�1���0�S�/��                                                                    8�[�8r Q7��7M|6�k�6`�50�4O�                                                                    5��4���4\k�3ҕ�36
2�c�1�u�0�d                                                                    8	�B7��7��<7 -:6y�'5���4��3�F�                                                                    6�K6�Y5��i574���3�#*2�Ou2(z                                                                    6�܉6��6YE6o�5I��4���3���2���                                                                    4���4�-*4���4�/3f��2�~�1�W{0Ֆ�                                                                    8(��8�I7��7Cŀ6���5��4���4U�                                                                    6@��6+�5ܱ�5_�$4�Wy3�330�2!��                                                                    4Ɖ�4���4�?4���4i�74!}s3̐�3�i                                                                    3�o3���3���3q��3=,3;2�N_2\v                                                                    4�4�4ӏj4� �4���4E`p3��3��S                                                                                                                                                                        +;�      ��     ��6:�R/D�<            G��]    = �6���3��"@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��`    @��@    17:07:25        G� @��@    @��     +<M      ��     ΀5�A�7��tB'�fB'�f8+.D"g�    )#D�+���71/�?�  ?qh�    4%��)%�I/���V�        /T�+?2�>� �>l=���<�g�:��3�a+��2[;p:+�]< �S<��c=G�?=�=c=��=��=��R=�=�=���>*�>�6@�Z�@��n@�?��?��P?�1?�B?e��@h@h�@^Z@^+�@^L�@^p�@^��@^�Y@^��@_ @_%@_)�=��=���?���                    E�/7�``G�ZF�-�F^v[E�F�E�oD9܂CK�BPC                                                A_HC!��B��A��AK&A@�h:?�\�>�be=��9                                                                    EL��@��F)T�A��;1A�Fh��B5}A=LO.�q�                                                    {@��B�kjB�kj@f��8����l�    �dV/f�&�?   )#D�A(�ʬї�1�{Q    >�y�BR�<1���AW��A�MA1�        B�B�CÂ CÂ C�=�?F�@�~�3�EC�q6��=6���A,��CF_�>�y�CZ�B���CjUA�F�B�6.CF�'B��B�H�B%\�5��    B%\    B%\B0�k            A���A<��A��@xjV@l(@���?׍�F�1A�.�m�+ӡ 1,9[8��6�65��F�KgG-��G��G��=���            >&�>��A>͏�>�#�>�Љ>�9�>�؅>��e8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M���    Aq�*Aq�*{@��8 *9�?�9���8篨8Z��7�e�6���6&��5=]�                                                A��bAq�/Aq�/{@�Π�+�B���B�<�@B��?z�{@��{@��    C���BL=6��(    6Hb#7L@�R�2ψ57h��    =���C���7�KCA��<D9l�C��B���B�[AB(k@�}�?���>�g�                                                ?C�A5[n@!y?��S?>1�>���=�n'=@�<:�                                                                    D��FT.�Egf}D�*D��CF�vB�(�A�?�@�$                                                @��VB�Aէ�A)�n@��@��?��]>�*�=��J                                                                    C[��E��ZD�3�C湇CS�\B��B ��A&7(@#q                                                @��B�OAF�u@��@,�4?��E?>l>:wW=4��                                                                    7���8X�AoE>F=(A� �?%i?�C,��g+f{�-���-+�,��g+���.���.��-�}B*��7��N4���2]n.�4�-��c-���'ڷ�7�,�4��2{�W2^O2](�.�17��,�1�{Q7�,�4���    3�2�4��4�ŝ0�3    :���8Q�8;�b6���5l�2�Py                        5]��6�R�6���4��3��87i�2�F�            4���>Q�q-S��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Cl��G�@�>+ZeA���                                                B�                                          B��                    A�                @�p�    /��+ S5M;/6(6M26�>1::G?�>L8R�@��8)�;�`�;ګd                8 �o    2=_K6�;3        {@��            �ʟ{@��{@��{@��7�Z7�#9            6�X�    6�X�{@��    6�X�    6�܃    {@��    {@��                6�܃{@��6�{7���7�9�7�9�<�Ww2	�YB���7��E��+D
�A��\    C0�=            <!ۏ4�j~6�"'4�j~<���<��=/]�<��v<U}b;��n;=�:��C� Gۋ>���SU��d����]��7���>��r�E�U�$�>��,�'�x <E�<�� <�Qc=.Y�<���<T#�;���;1�H:�L{� Gۋ>���SU��d����]��7���>��r�E�U�$�>��,�'�x 9%x�    *��\6��e7ݔ9��9Q88��a9��9?L{9�c"                                                ��,q�ǻ���'�ĳ��Ģ�ċ�g�c��)��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         ��+�        {@��{@��{@��{@��                                                    CK՛D�E_D��D:��C�˙C��BX�A�@���                                                A�՛CbE_CI�B���B0˙A���@��@?��                                                                    E�'�GݰGL9F�b�F6�$E���D�Z)C�j�Cz�                                                C�"bER��ECW�E�D���C�ldC��B,��A=�                                                                    E�mLF��F�s�F|�F]{�F/�E�@E���EG�c                                                D)��D�y�D��D�HxD� D~��D0��C�_�C�#�                                                                    7i�7~�G[y�GONHG��F�ـFy8�FcVE��#EP	                                                                                                                                @P�A_m�A���A���B.3�B`�B��(B���<�<�<�<�<�<�<�<�<�<�<�<�E��	E�{,Ebz�E�$D��DV��C�JC�b:                                                {@��{@��{@��{@��{@��D�1B �.�b��f�I8$�Ckxi@I�H    A�
+A�
+{@�ξ�������C��[{@��C��lC��lC��[{@��@f��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G.;;�~�D��/GT8D��DNm:DNm:@�6�@�6�Fk��BB!�    C���FKy�FKy�D���D���FlX�BBP                @���C��C��C���?   C�C��lC��lC�E�C��nC���C�<�C���C�t�C��C��\C�^�C��C�ŻC��C�>EC�3C�չC��lC���C���C�љC���C�2�C�5�C�+\C�,�C�98{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�<C��=?
�3?�>��+>�(�>�kS>�$�>�,>�@>Գz>υ2>�t�>ů>�Z�>�Y,>��z>�&�>��M>�o>���>�� @U�@TE        BLǪ=�ǩz�ǩo��["�#��"��"fK�|۰�ǭl�ǭ��ǭl��η        ?9�B,    A*��A*��{@��@��~F]L7o��7��        ;�L�?4�C?50�?6{0?8�?<n�?A�?INV?T~'                                                ��U7+��A�  ?k�A	��A z�+<M      ��     ΀9cM�8��8��7�[�6���5�(5�Z47{                                                                    8���7��a7?�
6�"�5���5]�4$d3'�                                                                    G�ZF�-�F^v[E�F�E�oD9܂CK�BPC                                                4��4<�3�=�3��2<	1}4�0�2�/���                                                                    4�3m��2�i32+N1m��0��y/�M�.��e                                                                    7���7a�|6ݔ�6T]15���4�`.4%G_3G��                                                                    7��37F#�6� a6Z�5�Pg4�]K46ض319�                                                                    5-��4�=4�4!?�3f�%2�i2�1�J                                                                    7�5?7r+�6� w6�<W5��(594_z�3X��                                                                    8���8��y8��7mǧ6�56��5&3;4$�                                                                    5�T�5��5�64�b�4t�3i2�2���1���                                                                    8���8��y8��7mǧ6�56��5&3;4$�                                                                    8 �t8<{7SX�6��!6ݤ5^�94�6�3���                                                                    5 �L5	?4���4�3�y2�UI2� 1i                                                                    8 �t8<{7SX�6��!6ݤ5^�94�6�3���                                                                    ����4�+�5�w�5�H�5:�4PiR3��2��{                                                                    5}3f68�#7��}7K�s6���6O�5��5k��                                                                    5�5�w�5�`4�e�3�B�3\2R��1K�@                                                                    7A�7"�6h��5�o5�4P�j3{�2y1                                                                    �B���+�`�ۢv�]{=��k���j�	RͲ��                                                                    �=H�~'����"F����[�>�	�����                                                                    6AT�6#45i��4Ƥk4��3Ri2~,T1{|3                                                                    �QX��0������{&�mĞ�����ܳ�s                                                                    �өc��xq� ����*���
�C�ɰ�p����|                                                                    ���0��gm�5����[�O��@kJ�bh$�xLY                                                                                                                                                                        +��r-�`                                                                                            8��{8�G�8 �77k�k6�#�6�+5 ��47�n                                                                    5+Ί5݈4�W,3���3=C�2���1�0�]_                                                                    7�w7�\�6���6jJ�5�Le5
.�40AJ3I�>                                                                    4(NX4��3^�2�;2C�?1��>0Ǧ�/�xX                                                                    8� �8�q7�R�7h��6�u6��5?��4[r                                                                    5$�&5
�54z�3�`�3Jl�2���1�8�0��R                                                                    85O81�7�:d73�
6�o5�i�4�ރ3�k(                                                                    64�66[5˰r5Mf�4���3�x�2��(2��                                                                    6���6���6��6;�5`�74��3��2ł�                                                                    5�5 ��4���4%�<3�v�2��,1��	0��                                                                    8A]�8* p7�՗7[�~6��5�52�4]�                                                                    6\�B6Bn65��5{�4�M�4I�3��2*�}                                                                    4�4�h4�Uz4��4�sx42KA3�k23�1�                                                                    3��3�s3��w3��T3R�y3p2��y2i
]                                                                    55j<4��4�o�4�p�4Y�34�-3�<�                                                                                                                                                                        +<M      ��     ΀6;�t3�2�            H��    > �6�bQ3���@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��@    @��     17:07:28        G @��     @���    +<�      �     �P5�I 7��4B'�fB'�f8$J�D"g�    )a��+�"[7;�6?�  ?tV�    4)�)�/����        /�8�+���>���>��=���<69��2��}+υ,2S�@: l<[V<��u=G�4=�:=��=�
n=���=�=[=���>*�>�M@�)@��@%��?�7E?�ܮ?�+r?��?��@g�$@g��@]��@]�|@]�@]�L@]�@]�x@^v@^(^@^;=@^F�=�&�=��?��x                    E���7�/�G���F�):F]��E��vEX	D9Z�CK\�BO�                                                 A�RC!H�B�-�A�0�AJ��@� �?��p>�כ=�~�                                                                    EL�7@�6�F)|�A�߼1�@FiB6��AB9�-���                                                    {@��B��fB��f@k��8��ˤ��'    ��C/5}
&6�D?   )a��A\����@1�;h    >�y�B8�2aBA<�+@��@��        BpmTBpmTC�d�C�d�C�V�?�Z@�b�3$C�
u7 �]6�(�A/�?C6�b>�y�CH�
B���C^ˤA��wB"�C;�%B�LB���B/Д��$K    B/Д    B/ДB4!�            A�BA.�2Aܔ�@lU@
�<@��?���F�$�1�@.ȕ�+��+1LB-8�6��5�#1F��EG.Q+G�فG��=��            >l�>��r>��P>�I�>�|�>���>��B>�]8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�I    A}7A}7{@��7��K9�[r9zE�8��8R��7�,�6�`6 ��56�                                                A��SA�5gA�5g{@�Ο�W�B�-IB�e\@I5I?	�U{@��{@��    C�(�BK�B6��}    6Kz�7@�h|2ڈ�7i`�    =��C��7���A�U�D,}�B�n�BRifA��NA.��@��?ڍ�>��                                                >�ŃA'	I@F�?�Fg?�+>x��=�ce='�<,�M                                                                    D��FT7sE_�fD���D ��CB�ZB�X�A��@�Ԣ                                                @�ZYB�
�Aɬ�A�@�M[@��?�&�>�=�=�                                                                     CW-PE��|D���CފCL�B��7A�m�A Ӱ@ )�                                                ?���B�oA;p@�7$@#?�?��{>��>2�9=06�                                                                    7z�7�Z	A`>L�iA�`r?*R�?���,��+k>i-���-P,��+��F.�o�.R-��q*�(7�&f4�#�2c|^-jA.��.w�(Ơ�7s��4�M2]d�2c#S2b��-i��s��1�;h7s�4�hE    3���4�04��E/��    :��89�+849�6���5Y
�2��                         5L# 6���6��>4��3�*�7G�t2��z            4��/>?�m-G��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Cc��G��m>*�9A��s                                                B�                                          B��                    A�                @�p�    /H�>+4�C5A&�6 ��2;1Y�?���8E*�@we�8�x;��;��                7���    /�F6�e�        {@��            ����{@��{@��{@��7C<8�@            6�ri    6�ri{@��    6�ri    6�d�    {@��    {@��                6�d�{@��6��7��8
708
70<��1�!�B��j7�`E܈8D.E,A���    C#�h            <,ƻ4��H6��4��H<��<��=C,�<�a[<g\�;�nS;F�2:�Iȋ �b�|����
���V��X ������|��s�V-�?M�,���J<*�z<�]<��=A��<�G<eA�;�H�;=��:��2� �b�|����
���V��X ������|��s�V-�?M�,���J9P7    +K	�6��8�#9��9_�o9��8�^�9z59#��                                                �c���c�����R��p��ш|Ķ�Ęm(̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         ��W�        {@��{@��{@��{@��                                                    CK�D�pD���D9��C���C��BVstA���@��                                                A��Ca�pCG��B���B0��A���@�st@
��?�                                                                    E�*G GhF�P�F6��E��-D�BbC�W�Cq/                                                C�%�ER�EC�
E ��D�ˇC�G�C�vB,�A=��                                                                    E�n�F�_F�y(F|�F]|�F/0E�>YE���EG��                                                D)��DăsD�$�D�LD��D~�D0�C�]�C�"�                                                                    7G�t7y3G[��GOdpG�F���Fy'oF]�E���EP�                                                                                                                                @=G�AO0�A��kA�;B$��BR�B|�B���<�<�<�<�<�<�<�<�<�<�<�<�E���E�� EbjE��D�̖DV�C�5C�`~                                                {@��{@��{@��{@��{@��D 0A�D5.�É�f%�8"s=C��@`$|    A�H�A�H�{@�ξ��c���cC��{@��C�6|C�6|C��{@��@k��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G^;q�vD�"G{�D��DJ"�DJ"�@֒�@֒�Fm�BD�~    C��eFKx_FKx_D���D���Fm��BD�6                @� C�J C��"C�;p?   C�^�C�6|C�6|C���C��yC�PC��C��C�p�C��C���C�pYC� SC���C�{�C�'5C���C��C�A�C��C���C��C��3C�'dC�4�C�+�C�,�C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�a�C��H?�o?	��?d?��?�>� *>�^`>�>�}>�:>ܑ5>��?>�\
>���>�bg>�}�>���>���>�r>��A@<�@A�        BPǜ
�Ǜcǚ�������rD��<�O��m���j���j���j��c+        ?$�Bɢ    A>{iA>{i{@��@��.F]R}7�P�7�        ;���?*/�?)��?*��?,*�?.��?2j?7w�?>�                                                ��ϓ7��A�  ?k�A	��A z�+<�      �     �P9\&w8���8p7}.�6���5�g4�� 4 pO                                                                    8�
�7�KQ7<?J6��5�G�5��4!�l3"=                                                                    G���F�):F]��E��vEX	D9Z�CK\�BO�                                                 4�47b
3�r3�z29��1zsA0��/���                                                                    4�3g�B2���2(�!1j��0�-�/��+.�w�                                                                    7��}7\6��6Q�5�b4���4"��3A�                                                                    7��7$vV6��6-[�5���4��*4ҽ3B                                                                    5�L4�'�4Z�}3���3M�2���2*�1u�                                                                    7��7Ii6�xs6S�(5��5��4?�=3B�                                                                    8ѹ�8��8�?7d��6� �6#�5��4��                                                                    5��5��95bm4��3��3]�*2��1���                                                                    8ѹ�8��8�?7d��6� �6#�5��4��                                                                    8$�8�v7H�z6���6	c5V?]4���3��                                                                    5�5�?4���4��3u�2��2�@1`[                                                                    8$�8�v7H�z6���6	c5V?]4���3��                                                                    ����4��6)�5�l�5�F4U�U3�R2�p�                                                                    5\�6(�I7�u7g�y6��6_y�5��5_m�                                                                    5���5eع4կr4[u3��3�P24z�17                                                                     7<�87!O=6_��5�"45
,�4J�3r2n��                                                                    �=���(���|C�[-�����eƳ���                                                                    �tf� �<��`
��񔲟g��Y��Li���                                                                    6<��6!25`�4��4	��3K��2te�1p��                                                                    �K$q�+⯶���}�j&���՛��uγ�d                                                                    ��c����I��Ͳ�k��-m�@�䰃���!K                                                                    �������b�2v���"��>��_uy�q��                                                                                                                                                                        ,?X�-���                                                                                            8��B8��17�L 7iDr6�q�6��5�42p�                                                                    5&�5+4��3��3:b}2���1��0� �                                                                    7��"7��i6�"�6g��5��s5�4-'�3C��                                                                    4#Ql4
0W3z�12�Y2@�p1�vr0�#�/ݡ�                                                                    8�,8�ա7��u7e�6��6z�5<dW4T�b                                                                    5��5P�4u�O3룣3GXb2�lD1�fV0�"j                                                                    8�8]D7��.71��6���5���4���3���                                                                    606�N5���5K�4���3֙$2�d)2�*                                                                    6���6�c6��(6�5^�Z4��q3��[2�"�                                                                    5I4���4�<�4$"N3~�g2�i�1�$�0ە�                                                                    8<NZ8&��7ֶ 7Y8�6�xr5�85l�4M�                                                                    6W4�6>z'5�bI5x@�4���4$�3��2&�                                                                    4�u�4��Y4���4�N�4�Zg40�3ە13�N�                                                                    3���3��K3���3�c�3Q:3�F2�p�2b�}                                                                    5V4�\44�m�4�C�4�)4X�40�3�|�                                                                                                                                                                        +<�      �     �P6>��3���                        6���3�1�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��     @���    17:07:30        G6 @���    @���    +=      �6     � 5��X7���B'�fB'�f8��D"g�    )��1+�]�7I�?�  ?zK    4��(�Щ/�AӰ�(        /�X[+�D">��`>+�D=���=-צ;��
4��g-lc�2B��:9<#6<�i=G�=�5�=��=��=��=�<�=��o>*�>�_@�2�@�`�@-a�?��?ҕ?��c?���?�F�@imM@i`@_`@^��@^�a@^��@^��@^f�@^H�@^1�@^!�@^D> �=���?s_                    E�%s7k�G�"F�I�F\��E�!iEjD8��CJ�DBOgy                                                AG�C �~B���A�x�AI�@��!?�^S>�_�=�58                                                                    EL��@�r�F)�IA��1$��Fi�:B8A�AD��-ժ�                                                    {@��Be�<Be�<@q0�8|����z�    �z3�/Lw��;?   )��1AXq���1�:    >�y�B
�2��A�<@A�m@AF�        BlQDBlQDC���C���C�U�?l�@�Gz3�GC��7
'�6�K�A3�C�p>�y�C(��B��C=��A��`BN[�CCCBGB�ߎB}i�5!2    B}i    B}iBQ�            A�s�A�A�Q@O�@?���@�'�?�m�F�1�1$��.���,��1[��8��D6�G�5��F�<�G �)GڃG���=�t�            >D�>��#>���>���>��~>�*>�z>�A�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    AS�cAS�c{@��7�Ń9�?[9a՚8�ή8<��7���6�|E6S�5'��                                                A���AV,�AV,�{@�Π3�BSǚB�#l@UQ�>��{@��{@��    C�8�BLM�6�5    6O��7`@��,2�6I7_V�    =�t�C�l�7�"UA�?�D#lfBθqB88A��A(Z�@�i�?��>�n                                                >ӷcAэ?� �?^-?)n>n��=ՌG=��<(�)                                                                    DĔFT%�EXσD��	C���C?��B���A��y@��H                                                @��B���A�WAP�@��@��?zV�>�b*=��                                                                    CSh�E��D��?C�C�CG6�B�gA��.A�B@�                                                ?���B>A1	�@��C@
�?��>��>,��=-V{                                                                    7^��7ܺ�A��>S�B��?0?�?��t,m�m+�4-^�v,ݑ.,m�m+>Z�.8@-��-�q�*�,7�x4��Y2Qyl-���.��.Z�)5p�6�ݗ4��2j|�2Q4�2P�-�Z��ݗ1�:6��84���    31�,4���4��r0�    :�Cu8�8�6}E.5!"P2��                        5a�6S�6I�g4�3Y��7�2�M)            4�hR>+U�-?�V?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CG��G�
A>+)�A���                                                B�                                          B��                    A�                @�p�    /fO�+Qx65�6b:1ޓ�1e�q?��'8.�@@O37ۉq;�}�;��                7��J        6|��        {@��            ��p{@��{@��{@��7j��7� p            6�=�    6�=�{@��    6�=�    6|��    {@��    {@��                6|��{@��6�}�7gq�7���7���<���1̞�B�)j7�FU$C�z?Aa��    C
��            <<cF4��6���4��<�Q�<�z=[1=��<�EK;�(;lD:�T.�!M��#Ԋ��\�����Kc��F芍#Ԋt"L�W�@0Ɋ-��y<8��<�Ib<�5	=YU&=}z<�8;��;Y":�|v�!M��#Ԋ��\�����Kc��F芍#Ԋt"L�W�@0Ɋ-��y9^�    +9�!6�0g8
"�9��@9�6>9@�E9FӠ9��I9^��                                                �?J��������#q�\������<̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �3�        {@��{@��{@��{@��                                                    CI��D�eD�O�D8e3C��C
ĽBTw�A���@���                                                A���CaeCGO�B�e3B/�A�Ľ@�w�@	��?��                                                                    E�*dG$8G��F�>9F6UE��dD�)�C�C!CgZ                                                C�&cESiEC��E �oD���C�#{C�]B,��A=��                                                                    E�o�F� F�}�F|!F]}�F/�E�<�E��aEG��                                                D)ŪDċ�D�+�D�O^D��D~�TD0�iC�\/C�"                                                                    7�7`�IG[�.GOz$G��F̼�Fy�FXE��BEP�                                                                                                                                @4U�AG�kA�5�A�"BBJ��Bq^�B�m�<�<�<�<�<�<�<�<�<�<�<�<�E���E��TEbRHE�
D���DV�pC�LC�^�                                                {@��{@��{@��{@��{@��D$rA���/
���f�_8�6C�j@oW    Ax0�Ax0�{@�ξ�.#��.#C��{@��C�SC�SC��{@��@q0�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G{�:��~D��G��D�rDF��DF��@΢D@΢DFm��BG<V    C�I�FKt�FKt�D��)D��)Fno�BG7                @���C��NC�C�+?   C�*:C�SC�SC��C���C��TC���C�k�C�BPC��C���C��tC�}C�BMC�1C��	C�i.C�dC�ŌC��fC�K2C�(�C�^C�!�C�2�C�,C�,�C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�+�C�l�?�1>�}�>�pD>��>��>���>���>��>�@_>��>�>��q>�n�>��>��>ʧ�>Ř*>�O�>�>�@3�|@ �        AȎ&ǂ �ǁ�Tǁ>���6"�����,l�MN�����������        ?��BtU    A΃A΃{@��@���F]��7���7�E        ;�p^?%?�?$H�?$�?%v�?&�`?(�M?,#)?0�                                                ��46�܆A�  ?k�A	��A z�+=      �6     � 9BA8�`�8�7f^�6���5٣�4��3�'�                                                                    8u_�7�D@7*g76�5��5	t�4�3i�                                                                    G�"F�I�F\��E�!iEjD8��CJ�DBOgy                                                4�5�4$!+3��2�g�2)�W1e�0yG/{��                                                                    4�V3ORR2��2��1V��0�:N/�p?.��o                                                                    7c�=7E,6Ī�6>�[5�jv4�z�4l�33�                                                                    7t�V7�6�:�6�W5d�)4��4<�3g�                                                                    5vN4�ZO4+�!3�vb34"�2�T1��0�V!                                                                    7���7$�E6���6*�j5��|4�
�4%I�30~�                                                                    8���8���7�d7Ltp6�\�5�>�5X�4�                                                                    5��F5��4��B4�3���3E�,2���1��K                                                                    8���8���7�d7Ltp6�\�5�>�5X�4�                                                                    8
Ӑ7�	#71R�6�A�5�~5@��4t�3s�r                                                                    5��504k��3��U3Z�2��(2�1��                                                                    8
Ӑ7�	#71R�6�A�5�~5@��4t�3s�r                                                                    ����1�z�5���5�~5WZ4F��3��B2ȿ�                                                                    5Ob5���7�W�7���76��6l5��1                                                                    5�~/5=�4��041W|3���2�2�1&$m                                                                    7(��7F�6Gs�5�s�4��47�<3[��2[��                                                                    �(N˶	�� �G�ߴ����ي���<��n                                                                    ���8��F}�׏��� �����H�D��(����                                                                    6(0�6Ӎ5G>4�f}3�$38W2]d�1\�                                                                    �3sӷ����9�#H�T�Y��$���A��                                                                    ��oٳ������q9��J�/�H�pm����8                                                                    ���۴y��!�����ݳ /��/vǱN�ʰ`+l                                                                                                                                                                        ,)-�O                                                                                            8���8y�H7�"�7S��6�|�5�d�5�44$�<                                                                    5E�4��`4f�^3��`3)U�2��1�	J0���                                                                    7���7t�g6�|�6R�5��h4��4�X34�P                                                                    4E�3�c3a�02�K�2/91�s�0��F/���                                                                    8���8o��7��|7P�~6��6�g5,�}4D�c                                                                    5E�4�=�4]B3վ�35B2��21ñB0��                                                                    8��7�|7�s�7"B6|]�5���4�I�3ܪ�                                                                    6�6�5�;59;�4�5�3�e�2�TC1�0�                                                                    6ܻE6�E6�ٓ6�5K��4���3�E�2�Q                                                                    4�C�4�`O4�A�4��3i�2��A1���0��c                                                                    8&�8�!7��7F�6�9`5��4�v�4�                                                                    6>�m6+2o5޹�5bep4�A�3�CR3��2�                                                                    4�-�4���4��4��t4l��4"��3���3�@F                                                                    3��3���3�Y�3u?�3?�3�2�	y2R��                                                                    4���4�a04��74�x84��<4F�w3�!3�2                                                                                                                                                                        +=      �6     � 6B��31�,                        6�{(3�W*@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @���    @���    17:07:33        GT @���    @��    +=y      �T     ��5���7��B'�fB'�f7Ώ�D"g�    )*��+�:N7DwU?�  ?~�    43V(v��/�)��N"        /^��+I��>���>0+7=�Ϛ=7K;��4�$-�j&27�*:/�<4�<�z/=G�%=�1k=��=��=��D=�< =��>*�>�p@�S�@�vZ@�?�(?�Y�?�u�?�3�?��@j�@j�W@`~"@`\�@`9�@`�@_�@_��@_�s@_c�@_C@_/�?,�=��?���                    E��7ra�G��F�F\[6E���E�rD8��CJ��BO&3                                                A�LC *'B�NZA�ωAIn�@�V�?��t>��=��                                                                    EM>�@ӱ6F)��A�F�2g3Fj'B9�IA=KO-�b�                                                    {@��B ͙B ͙@vۊ8[$l S    ��C�/!�%yl�?   )*��@��>���O1�_D    >�y�A��C2��@��~?<i�?:�W        B_cB_cC�&�C�&�C�Fu?, @b��3E�C�6�7��6�t�A6LB���>�y�B�aB�\C�/AN��A�P&BΎ*A�rBjI�@H�7��'+    @H�7    @H�7@���            A%ҙ@�X�A���@� ?��2@l��?r�F�%�2g30�m-%J2�\e82�6"�b5�{iF���G	�|G��uG�P>�            >&�>�R>��>Õ>�r>�>��t>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    A9"A9"{@��7�xD9�3�9j[/8˄L8D?7�}6��6@057�                                                A�ȜA>�RA>�R{@�Ο9@B��Bo`7@c��>�f�{@��{@��    C�M�BL�P6�k�    6S�_7�@���38�7\Ba    >�C�d�7��#A��ND��B�J�B*�|A�KkA&hd@�`\?��>���                                                >Ǯ<A��?�α?M�h>��>k��=��O=C�<'                                                                    D	��FTEQ��D�CC�C�C=}MB�^A���@�7"                                                @�P�B��rA�0�A��@��j@
Վ?tI�>��=��                                                                    CP"E�sD��C�9�CB�B���A�&A�@�                                                ?��B�'A(3�@��@|�?��:>�b>("3=+H7                                                                    7e�N7ゑA%�>Z�tB
 �?6?�֟,~�*��-5>L,�7,~�*�#�-���-���-8�*Kϳ7���4�D�2-��-a{/4t/
�*%���!��4���2{�2+�D2+��-^�07!��1�_D�!��4�>�    1�k.4��
4��/�G    :��7��X7��M5E��3��b1��                        3y�a5$xI5�2��1��5��2�*�            4���=9R<-)݌?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C�RG���>,�A��*                                                B�                                          B��                    A�                @�p�    /�)+�4�?V6q�2�J=2�c�?<��7�7?�Q�7pW�;��;ۤS                7Xe    2L$�6�n        {@��            ���8{@��{@��{@��7�z�8$�            7,|�    7,|�{@��    7,|�    6*J    {@��    {@��                6*J{@��6N��6�s-8+Ò8+Ò<�g�09J�B��M6���F-%KC��	@���    B�v            =q�4��88J�4��8<ш�<�2,>�M�>6�>=��s=9�:<ȹ~<V��� ������jϊ�` ������ߊ���s�1�V���?��-]��(�=k5�<сU<�>��->6PA=�	�=8c<��z<)u�� ������jϊ�` ������ߊ���s�1�V���?��-]��(�:�8|    *�P�6m�K7�VO:�b8:�p9̼G9��';9�;4�                                                 �����SD�ϒ����ľ��İ �ĝ!WĄ^̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �9@        {@��{@��{@��{@��                                                    CH�D���D�AD6��C�B�C	�5BR�A�z8@��                                                A��C`��CGAB���B.B�A��5@��@z8?�                                                                    E�)3G@�G��F�*QF6f1E��|D�7C�, C]                                                C�$�ESF�ECŻE ��D��iC��AC�?B,}A=��                                                                    E�p�F�%BF��yF|VF]~�F/,E�;`E��;EG��                                                D)�tDēID�2$D�RmD�UD~��D0�]C�Z�C�!.                                                                    5��7i7:G[�:GO��GtF̫Fy�FRE��uEP�                                                                                                                                @P��A]��A�R\A�zUB*�+BY��B�,�B��<�<�<�<�<�<�<�<�<�<�<�<�E��!E���Eb7RE��D���DV��C��C�]3                                                {@��{@��{@��{@��{@��D�B:/.�K˧eNs8��C�AI@rK�    Ak9Ak9{@�ξ�����C�S�{@��C�?�C�?�C�S�{@��@vۊ{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��Gz�:��D��G�hD��DC��DC��@�w�@�w�Fnc:BH�K    C�CFKp�FKp�D���D���Fn�BH�                @5*kC���C���C�Uh?   C�FbC�?�C�?�C�4�C�,�C�$�C�AC�-C��C��C���C��(C�v�C�LPC�2C��C��$C�a\C�WC��C���C�e�C�4CC�#�C�1"C�,LC�,�C�3�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�E�C�h�>��h>�{�>�uK>�,>��g>���>�8>���>�t>�G�>�z>ߞ�>�N�>؉n>�8�>ϕF>���>�6>­�>���@3��@�p        A��c�8؝�8V��8,�q��ǵ)Ǵy^ǳ�������������������        >�c�A���    A��A��{@��@�VhF]��7��7�        ;�?4j�?3��?4�q?5��?7gf?:<?>�B?E�1                                                ��R6�{0A�  ?k�A	��A z�+=y      �T     ��9D�8��8l7q��6�^5��4�\)4�F                                                                    8x�u7�9�71`6��S5���5��4 ��3%I�                                                                    G��F�F\[6E���E�rD8��CJ��BO&3                                                4��4)3�[�2�i
23��1u6�0�]�/�M�                                                                    4c�3U�+2�f�2!O�1b�S0�� /���.���                                                                    7gDV7K�p6�ދ6G�15���4�E4 _�3D��                                                                    7o��7�b6x�6�~5o��4��j43�`                                                                    4�C�4�/�4&J$3�XO3<u�2�n\1�.W1ƍ                                                                    7�_Z7 ��6���6'S�5�uF4��f4,vt3@�v                                                                    8��-8���7�Gt7Sg�6���5���5>�4 �                                                                    5���5��o4��~4�"T3�`o3NV>2��j1��<                                                                    8��-8���7�Gt7Sg�6���5���5>�4 �                                                                    8?8��74�<6� �6 ��5JW�4���3�Ӝ                                                                    5�*5�4l%Q4 #\3b��2ȠV2ˠ1Ka                                                                    8?8��74�<6� �6 ��5JW�4���3�Ӝ                                                                    ��L���6'd5��50)4T�3���2��                                                                    3�0
4��}9<,8˴W8S]�7Ч�7b}�6���                                                                    5���59V�4�^�4."�3�m�2�Nx2"�15�4                                                                    7,u7�26L��5��456�4B b3j�2p �                                                                    �+OŶ7Ե��R������]k�AԲ�R                                                                    �����̲�	|�Ļ���	�VQg�^����                                                                    6+&6�R5K�4�c�4�!3B �2kD�1q6                                                                    �6k�Ԉ��M��
�b�_�����γ�@�j                                                                    ������,��貎q����e�:鰀�I����                                                                    ��8Ĵ�ǀ�)�����č�;m��_ �v�w                                                                                                                                                                        ,-�-%9�                                                                                            8��58��Y7���7\�n6���5�J{5h44��                                                                    5s�5YT4n~�3�b32�2���1��.0̪�                                                                    7��&7|��6�l6[R-5�o�5XO4*H3F�                                                                    4h3�b�3i��2�÷28=�1�Z�0��q/�jh                                                                    8��8wV�7�E�7Y��6�;76	�`58�'4W�                                                                    5\�4��4dƊ3�%n3>nH2�1ф�0�*                                                                    8�7���7�`�7*j�6��{5��D4ۄ3�2                                                                    6��6�p5�%S5B��4��3�ۻ2��2
��                                                                    6��6��L6�rO6	��5W��4�}3�b�2�"B                                                                    5 b'4�$�4��64a�3vٺ2�c�1ʺ30�''                                                                    8)�L8�7�Y�7PIJ6�Xy5�S5&4S�                                                                    6B.61Q5�f�5n
�4��B4 ۝3P2)�                                                                    4ǎr4�04���4��k4{+�4-��3�G83�N�                                                                    3�A�3��K3���3� �3J�d3�P2�1�2g��                                                                    4���4�4�1G4�J4�~J4T��4 �3�'H                                                                                                                                                                        +=y      �T     ��6F�41�k.                        6�r3�>@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @���    @��    17:07:36        Gs @��    @��`    +=�      �s     �5��7~��B'�fB'�f7��dD"g�    (�ٯ+�[)7)/?�  ?ٴ3-�3�u-w�x/qq��ǣ        .��*��R>�S><4�=߲�=N�	<3�5}�-��j2.�:B�<�:<�^H=GԽ=�-�=�p=��=��$=�;Y=���>*�>��@α�@m�@�&?�l?�"�?I�V?,�
?
��@k\@k�@aqz@aA�@a�@`��@`�@`tJ@`J�@`*�@`X@`	�@��=.��{@��                    E�C�7u�oG�XF��UF[�aE�cEX9D7�CI�BN��                                                A[tC�rB��wA�{AHڌ@�� ?�f�>�U�=�l{                                                                    EMt@���F*�A�t�3�m�Fj��B:K.A0_-�g                                                    {@��A���A���@wW�8J-�&.�*    �mTa/E	���?1�(�ٯA
ә��Y3���    >�q@�`�2�{@ }�����Q�        B8�*B8�*C��aC��aC���?Z�?�GR3{6�C�fg6�͚6�D�A7�>B#�[>�qB8 �Ae�Bx�A>7A#tB4�xAM,�B6����I�1'    ���I    ���I���<=��    6�l@�)@ PA!�?�k>�^-?��+?��F��:3�m�1eX.BǼ3��67���5b��5�k�F��#G]G��G�A�>��+            >Omv>�(�>�'�>ܳ�>Μ�>�D>�y>�zc8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    A��A��{@��7�r9���9l��8�{98Kt27�RM7 ��6/W�5U�                                                A��A(�A(�{@��:C/�A���BNB	@mH�> �D{@��{@��    C��5BN 26��    6��36�T�@��|3��d7F�    >�M�C���7���A�XvDT\B�.�B/6A���A0�	@�p?��\>�BX                                                >��A"X?�~�?Z�?��>}�!=�J�=|�<,n                                                                    D��FTI�EKs�D�t{C� �C=WB���A���@��m                                                @}��B�8[A�11An.@���@
�$?r�>���=���                                                                    CM��E��PD�ȌC�Y�C@�2B��qA��)A @ e                                                ?��-B�lA"��@���@�J?���>��>%�=*��                                                                    7j�7羯A*9�>`	�B�?:�g?���+_<E*�,�"�,Pn�+_<E*2��,��,���,���)� [7��k3ŞZ1�YM,�!�/�B�/��+1���_=4!��2b��1�`,1�*�,��z7�_=1@з�_4���    ,�`�3�U3��/$��    <Sd96S	7 ��2�$�1��!/`j                        /t�2��2�\�05A�    2�M�3)�s            4;�<&�-%�?�  ?�  ?�  ?�  ?�  ?}H�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B���G���>+��A�m�                                                B�                                          B��                    A�                @�p�    .tiZ*^�3��6��3��}3���>Ɠ7�6�?��F6��;��1;ʤo7.��            6�[    2�f�5��    3�!{@��            �ػ{@��{@��{@��8:�k8�._            7���    7�Ѝ{@��    7�Ѝ    5��{    {@��3�!{@��3�!2�2
�`    5��{{@��6hѶ6B8�'�8�'�<�I�    B�ݫ6՝�FC�C�Y�@:91    B:�            >��U3Ţ9`�3Ţ<�_=��@�8?�$�? %>���=���=�zf� ���g��ka��`����~���A���g�s���V�/�?φ�-^�(�>«�<�M=��@k�?���?�m>��=��=2�� ���g��ka��`����~���A���g�s���V�/�?φ�-^�(�<v$/od�+��6B+�7T&!;��r;��:v�9�.�<|�!<�x                                                �\6��Zw�N�r�A��.���3��߄ZËh�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         :GE�        {@��{@��{@��{@��                                                    CG�AD�JD���D4�oC�ǫC�*BP��A�#&@�]�                                                A��AC`JCF��B��oB,ǫA��*@Ѐ�@#&?]�                                                                    E�%�G\EG��F��F6HSE��ED��KC�CM�                                                C�MESn�EC�E ��D���C���C�NB,e�A=��                                                                    E�q�F�*cF��F|�F]�F/�E�9�E���EG��                                                D)�DĚ�D�8�D�U�D�D~��D0�C�X7C��                                                                    2�M�7oG[�1GO��Ga�F̗,Fx�*FKE��EO�0                                                                                                                                @��jA��A٥kB	�oBA��B}�B���B�"�<�<�<�<�<�<�<�<�<�<�<�<�E���E���Eb�E��D���DV�[C��9C�Z�                                                {@��{@��{@��{@��{@��D��B��.קeU�8�iC���@aI    AH�AH�{@�ξ��\���\C���{@��C��!C��5C���{@��@wW�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��GK�:K�~D�pG{�D�FDA��DA��@Ə�@Ə�FnHBHb
    C�Q�FKj�FKj�D�ʫD�ʫFn�cBH\�                ?�-�C��iC�\�C���?1�C��fC��!C��!C���C��C���C���C�
C��C��C�#C�)C��C�>C��C���C��YC�{�C�C�C��C��5C�� C�Z�C�,�C�0)C�,uC�,�C��m{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C��e>��>�<q>�Y�>�Q�>�9e>��>۵h>�>� �>��<>�_B>�q�>��>�.�>ԒK>�R�>͠�>�Ğ>��>��0@a��@*��        A P�������nh��P�<ʛ�<IZ�;�L�I. �§��§��§��g۟        >x1A�Y�    @�@�{@��@�rF^*)7f�7Io        ;�O�?K�?L�?NW&?P�?UW(?\�?g�?tRX                                                �A+6�D�A�  ?-�@��RA z�+=�      �s     �9B�a8�t�8��7{�6�!5��5��4��                                                                    8v%7���75zl6��;5㈠5޿44�,3?�\                                                                    G�XF��UF[�aE�cEX9D7�CI�BN��                                                4�ݾ4*�3��e3��2>V�1��0��/�b�                                                                    4P3Wqh2���2'�51pmk0���/��R.ʗ�                                                                    7e5�7M��6ϩ�6N�5���5�43|j3d9                                                                    7i��7�6���6 f5���4��4#l�3={                                                                    4�h�4�nQ45O3��3X'�2��2
�w1 �                                                                    7��V7#�N6�ƒ62�a5���5~%4G��3g�                                                                    8�8�5�7��7ZP�6�C�6�15+-�46A#                                                                    5�'�5���4��U4�!�3�O!3]y�2�V-1�mK                                                                    8�8�5�7��7ZP�6�C�6�15+-�46A#                                                                    8ս8u�77>6��6tK5Y#�4���3�[                                                                    5'z5��4n�4��3q`�2�P�2ҁ1(��                                                                    8ս8u�77>6��6tK5Y#�4���3�[                                                                    ��Ϣ��
�6S�5���5c�4bc3��52��                                                                    1:86��:�3�:<@Q9��9b8��8�                                                                    5�k�58��4���47T�3�	�3��2;�1YM�                                                                    7*�7"`6Pu5��5
��4Pk�3��q2�%�                                                                    �*3��,۵�P��Z�⴮���h�c��-�3                                                                    ��-�����%��ꬲ�c�gxw�Ly��G�                                                                    6)ۑ65N�\4��4
"c3P.�2�]�1��X                                                                    �4_�� �춓����kĒ�������09                                                                    ��^_��X}�:Ʋ����ұG>���'����                                                                    ��N���d��-_ѳ�4��.�JN�{>ڰ�>�                                                                    )�(�*㗃/�.W�-�E,�a�/�΁0 ��                                                                    ,\,-^oy                                                                                            8�i8�U�7�k7d�6��A6��5-U4Q�                                                                    5�5��4rB�3��3;�2��1�
0���                                                                    7�mY7\6꺉6b~"5���5��4=�3e2�                                                                    4c4w3mS2�=2B2|1�K0��+0��                                                                    8�q�8zF7�ا7`ܩ6��6%�5Nr�4y^`                                                                    5�84��'4hcW3�qh3H��2�{F1��K1<@                                                                    8
S8 Z�7��Q71c�6���5�$�4�Q�4�                                                                    6�_6�5��5J�4�T3��3SZ2!&�                                                                    6�7&6�qc6��O6XQ5e?U4��3�ڹ2���                                                                    4�u4��4���4#Ҧ3���2���1�g�19                                                                    8(΂8��7ЗG7X�-6�^L5�e�5#�4,W�                                                                    6@�63J:5�c�5w�44�"�4�3,�62D�_                                                                    4�%J4�T4�m�4�H�4�eh4;�*3�,_3��J                                                                    3�"3��<3�u63�_!3W��3ȳ2Ǽl2��                                                                    4�-�4��4�{4�<�4�
*4e��4�3˘�                                                                                                                                                                        +=�      �s     �6��%,�`�F�P�    =%�T            6�*T4+�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��    @��`    17:07:38        G� @��`    @��     +>A      ��     �05�q�7C�B'�fB'�f7h�mD"g�    &�M�+�3b7�?�  ?��4`;3���1��1�r鮗��0-�    -1�(�vx>�[>e�">��=��=!r9>11�f22#;�:	�1<8�<�M�=G�2=�+�=�	�=��=��N=�:�=���>*�>��@��8@S�l?�'?r[�>�F=8a�<��>N\�@gi�@li@a�(@a��@a^�@a#�@`�e@`��@`i�@`B@`)@`@7�<��{@��                    E��7J7�G���F���F[4E��E�bD7� CI��BN�                                                A�TC<�B��A�r AHH�@��?�,�>���=�΅                                                                    EM��@�aF*-PA��3ղcFj��B:_�A �-�-�                                                    {@��@X�h@X�h@u�$8��%r��"s�J#�4/F
�%^(d?!�&�M�?�Ӡ.�(4	-v    >�N�?)�V2S�?�|�'?�$��;���    B>#�B>#�C��RC��RC��?"|�?$�3 �*C��^6�Z+6�TA8D�Ay�4>�N�A�)@�BA�#�@��F@FcHA�:�@��A������@`P7"��M����    ��������?騉=�wp=��K@1�t?a��@��?5M�>��@6�X?4áF���3ղc1H��.��;3֝�6��2[Q�5���F�"JF�=�GCSGj�?86=>�|?z^s>���>�MF>��(? PS>�˯>�\�>�l�>�(>�*8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M뾶    ?'�?'�{@��7���9P�29C�&8��83�<7���6���5�g75B��                                                A��,�Ū��Ū�{@��@�s'A.I�B<|�@q��=��{@��{@��    C�2BPR�6m͏    6Y+7�@��G3@U7-y    >�9}C�8�7��cA���D	�B��vB/�mA�T�A2�@�3�?�Ű>��                                                >�\�A��?� \?^s�?Z->�-�=�N�=#1)<(	�                                                                    D�FTvYEF�D�zC��qC<ZmB���A��(@�=#                                                @{I)B��A�>�A�a@���@
�T?vV�>�،=�ֵ                                                                    CK��E��6D���C��C>�	B�ޥA��AA�K@��                                                ?�}�Bg�A�@�f@�?��>�>&'S=)�                                                                    7?��7�vA+�<>b-B?<Z�?�R*MĤ(�
,J	r+�N*MĤ)$�j*��~,,�>,))I%97hi�0�C�/��)n.�,W.�#(��=��\i3&�K2b��.�v<.�n�(py=7�\i��6��\e4��O    (���0�)q0��*�ޘ    =��<���3�9�/n�m.���+-��                        +&�c/�/�+���    /o4?2��            2�%U<��K-	�?�  ?�  ?�  ?�V?tD?#|*?�v?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  BG�G�.>+I�A���                                                B�                                          B��                    A�                @�p�    ,��(�4t0�Q�5���3�e]3ֺv=�$�6�x�>���5x_o;{j;|P7��S            5�O�    4�|�4�r    6���{@��27�41��    �ܗ�{@��{@��{@��7�Q�8)1            7,��    8��{@��    8��    4�U5	3{@��6���{@��7�_6��6��    5�/{@��5�Z4�:�7�/7�/<��    B��6�'�FY�C�_�?��b>12�AeK�            ?v�0�D':�0�D'<��I>�>@��@V�0?��??�^>m��=�<�� ���f��ka��`����~���@���f�s���V�.�?φ�-^�(�?t��<���>�'�@�̭@Vu�?���?��>m��=�p�� ���f��ka��`����~���@���f�s���V�.�?φ�-^�(�;�
�/�\�,�$5��h92��<U�;Q:S��8 ��6���<�/�                                                ËŢ���p��a��ƑÞ�V�e���1³�<�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 3���3�a>�6�5"�m4�\�?�4�@��?�A�@N�(=1��4�KE41�Cj(i7!�H<<$?P�M=>�=�]�7!�HG��5"y�6�W�CF�ND߮�D�w"D3D�C��HC��BOw�A�:n@�̐                                                A��NC_��CFw"B�D�B+�HA���@�w�@:n?̐                                                                    E��Gr�G��F���F6)YE�tOD���C��0C<�                                                C��ES�OED�E �nD�{*CЩ.C��B,RQA=o�                                                                    E�r�F�.�F��LF|F]��F/�E�8�E���EG�M                                                D)ǾDġD�>�D�X�D��D~��D0�mC�V�C�	                                                                    /o4?7E��G\XGO�2GP�F̃lFx޸FF�E���EO��                                                @3��;�'�                                                                        @��PA��AA�B�NB]��B��B�C B�Ph<�<�<�<�<�<�<�<�<�<�<�<�E��`E��?Ea��Ew�D��VDV��C��dC�W�                                                {@��{@��{@��{@��{@��C�kB6�F,^�i�fn17�6�C[Z�@@(�    @��@��{@�ξ�Z���Z�C��{@��C�"BC�/�C��{@��@u�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��=�D�
$G<�D�	�D@kD@k@���@���Fm�BFgv@4xC���FKd�FKd�D��FD��FFn�BFb.                ?S3pC�zC�M:C��&?�wC��C�2PC�2PC���C��+C�(�C�vtC��KC�iC�h�C��tC��!C�	�C�/xC�J�C�[�C�`#C�T�C�:�C�IC���C���C�{�C�9�C�04C�,�C�,�C��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��{C�;Q>��>��r>���>���>�>���>�s >�+�>��>�3�>�a>�n=>�+Y>�DI>τ�>��,>��>ʣu>��e>���@:�@q�        @�������s�����!����Wm��1�����o\�ƿ��ƿ��ƿ���         =�SA��    �?�O�?�O{@��@��F^K�71g�7_�        ;�ш?Uw\?^{�?b_�?g�A?p�=?{;�?��?�                                                  �6`%�A�  >��d>�R`@S��+>A      ��     �09�[8��7��7^�_6���5�-�4�fU4
��                                                                    8Be7�k�7v6��5�N�4���3�fQ3/�                                                                    G���F���F[4E��E�bD7� CI��BN�                                                4�K�4��3���2�AC2&�61��0�/�v�                                                                    3�e34�2���2�1R��0Hp�/@�P.��                                                                    74��7,Mf6� `66�5�o4��3��b3PaC                                                                    78 N6�L�6bk6 ��5p�}4��_3��3(5}                                                                    4ɸ{4�C�4*3�I�3@�*2|0�1��d1e�                                                                    7`��7=6�'6�v5��L4��<3�03M��                                                                    8�O#8�A�7̻�7@��6��B5��,4��4%�i                                                                    5i��5ys�4�'H4r�~3�@t3U;2"�1�*                                                                    8�O#8�A�7̻�7@��6��B5��,4��4%�i                                                                    7��7���7X6��v5�J5��4�o3��                                                                     4�ٲ4��4KIX3�k3Sk�2���1�iK1��                                                                    7��7���7X6��v5�J5��4�o3��                                                                     ����4.5��*5�X4�-a4�3+�p2�k$                                                                    -z%9N��;BE:��B:ei�9��\9I8p�                                                                    5��/5*�4�|h4 ҡ3�h2��h1�c�1A`                                                                    7��7 ��61��5�."4�,f3���3>2|��                                                                    ��[��յ�C�BH��g���ײ�+�
�                                                                    ��U�����/e��=���`+�9籚#,����                                                                    6�
5�<50A4�:3�z�2��'2��1}�n                                                                    �N��F.�|{���ݵN �[�ٳ!k���                                                                    ���|���N��EҲ����"0��Bb����X>                                                                    �]�A�Z�$�mT��(.��ꬱ�B��-ΰ��                                                                    )��--BJ0�/K�.H�1,
k�*�0�my                                                                    +��Z.p�V                                                                                            8g�8Y�a7��7I�.6�~�5�?�4��4>��                                                                    4��4�dz4O4�3΃�3$ p2-��1E,0��                                                                    7bP�7U�R6�Ð6Hh5��4�0<3�ܱ3Q&�                                                                    3��l3���3J��2�	�2)�*16I�0X2�/���                                                                    8]��8QB7Ė�7F��6�+;5� }4Ϩn4c�N                                                                    4�4�i64F�3ˏ�3/g�2>�#1k9 1 �                                                                    7�r�7�<�7��7��6x�5ou�4z540                                                                    5���5��5��54-44�D3�՟2���2                                                                    6���6���6m��5�ˏ5I/	4A��3J02Д�                                                                    4ɽ�4��}4���4��3e�2]%�1g^0�`�                                                                    8~�8��7�ǉ7@�W6�%+5�V<4��4�;                                                                    6�6R�5�vS5\7?4��V3�=�2��}24E�                                                                    4��4��W4��D4��,4jP�3��3z(L3���                                                                    3}!�3��r3���3o �3=X[2��{2J%�2v��                                                                    4�nR4��4ōo4��4�184
�3�߽3�n                                                                                                                                                                        +>A      ��     �06K��(���                        6�4�3��F@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��`    @��     17:07:41        G� @��     @��     +a      ��     � 4��7�B'�fB'�f75hwD"g�    "�)+�P�6��}?�  ?��5�K�26�2�E�1	240*~�        (�:P$��9_��;��>-G,>��=���<��\5�/+2!):ea<^�<�Q�=G�,=�,�=�
=��=��`=�; =���>+(>��<��#=��2=��=Q�=HJ=�k7�L�=u��@b�d@k�i@a��@a��@a^�@a#�@`�c@`��@`i�@`B@`)@`�{@��{@��{@��                    E���6��G��F�1�FZ�GE�Z{E��D7ηCJ�BM��                                                A
]C��B��A�;�AH@�v5?�8�>�s�=�~$                                                                    EM� @��F*/�A��u33<PFj�IB:A�=�U-1�                                                     {@���{�e�{�e@u'�7�ԏ&J	����ɨ�v/&a�?�e>�y�"�)��O/���54/�    >���<�|B2��?���6���(�<{�    B=��B=��C��vC��vCk�:?n*�>A�G3)��C�b<6��6�@�A8�@��,>���A�@f Ab�/@5�	?�4Anf@d��A����2@ ]�    ��2    ��2��=�    ?s��?s��@��>���@Q"�?�>�e@?r�?og�F�kL33<P0�	-�#3776�    5'�FH6sFA�F��5F��?�}@�ŀBD�a@�i?~�	?G�j>��>�.�>��>�̢>��>ͯY8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��#    ���?���?{@��7*@�7�ǔ8w�8u
�8Cj7[O�6���5u�M5 eb                                                A��/�*��*�{@��B�w�@�p(B4P@sBg<��{@��{@��    C��	BQ��61:�    6]7�@�9�3cS6���    >��C�
�7�ϛA� Dc�B�aiB7T[A��A;,�@�]p?�s>ْ�                                                >�H]AX�?�
�?g��?ќ>�M�>�K=D@�<'��                                                                    D��FT�EF�cD��2C�jqC<��B��A��@�%�                                                @|��B���A�p A+f@�O�@|?zqj>���=��J                                                                    CL6,E�ێD�r�C�AC>K�B�N	A�lUA��@��                                                ?��EB�A �@��[@L?��0>�>.2�='ɐ                                                                    6�
v7CA*`�>`�B��?:�?��)|�(7e+���+��)|�(I�2)/�+� �+p{(���75hw                        ���/3)��2DG�            7��/��β���/4T5v                        >'�0�;�                                                                2�t�            2�t�{@��,�,�>��>��h?dx?ya�?t��?:�\>�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A�*�G��>+7�A��                                                B�                                          B��                    A�                @�p�    (�ê$sOs    4��43Ań378�=&b�6{=�t4#9;�[;�_6���            ��~�                    {@��3��4���7%�~����{@��{@��{@�δe�7���            7%��    7{*�{@��    7{*�    4'}4�*�{@��    {@��7���6,o�7�=�    4���{@�ζ�1���/�}/�}<���    B��6��uF|�Ds��>�N�>
�@�CD            ?�:�    :PKs    ;X?W��@���@���@��?_a?>��>	%�� ���f��ka��`����~���@���f�s���V�.�?φ�-^�(�?���;W?W��@��L@���@��?_^m>���=�@B� ���f��ka��`����~���@���f�s���V�.�?φ�-^�(�;��.�S�,F*
3z8I��:�wZ9�D�8���84�4�M<�+K                                                �"�H¨�
đ��Ą�V�g,n�=
�4�í8̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6<��5em?���7��6f�l?�Y�?���?���@M�-=��D7R��6��FOY�7�y>LH,BD�R;��>R��7�yI���4��z7�yCFIYD߶AD��bD2F,C���Co�BO�A�A@�,�                                                A�IYC_�ACE�bB�F,B*��A�o�@��@A?,�                                                                    E��GgG��F��F6_E�b�D��:C��4C/�                                                C�ES~�EC��E �,D�n�CЏ�C��B,S�A=\�                                                                    E�r�F�-8F��F|F]�{F/�E�8�E���EG��                                                D)��DĞ�D�>�D�Z;D�ZD~��D0�yC�X"C�                                                                        6�v�G[�[GO�lGF{F�w�Fx��FD�E���EO�                                                A�-�A���<�zB                                                                    ?��GA;�uA�R�B%B7�
Bn�B�LB��c<�<�<�<�<�<�<�<�<�<�<�<�E��E���Ea�zEnpD���DV��C���C�V                                                {@��{@��{@��{@��{@��A��2B��5('� �e�7s6cC.�"@#��    ?�m�?�m�{@�ξ?�8�?�8C�Y{@��C�2�C��eC�Y{@��@uE7{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�f>u�D�GYD��D@�9D@�9@���@���FlQBD��B}�C��hFK\RFK\RD��]D��]FmaBD{�                >�w1C�DEC��PC�"?[ZC��;C�qNC�qNC���C���C��*C�!0C�aC���C��C�7uC�ucC���C���C�2�C�mKC���C���C��XC��C�ڃC���C��mC�HlC�1_C�,�C�,�C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�ӭC��>���>�=>��>���>�|�>��r>�e�>�r/>�42>���>�P�>��>���>�,�>�r�>��>�xQ>ǻ�>���>�G�@m�?��0        @3�=���<������ͯ�h���ú��ø=�÷5�kX��§�§�§�h��        =L��A��    ������{@��@w��F^L�727`�        6�C�>a0�?Ks?@'�?C�L?IH?Q)Y?]G�?n�                                                ©�e4��A�  ?L�G=�Z<#�
+a      ��     � 7*73y7�Ls7$��6l��5��=4C��3��C                                                                    6I�6b��6�36PF�5�y�4�2�3v��3��                                                                    G��F�1�FZ�GE�Z{E��D7ηCJ�BM��                                                2�(�2��33�j2�<B1��1
1$/Η^/q�h                                                                    1�iq1�*2cJ�1�81��0.��/z�.�Ϡ                                                                    5;;,5�}6t��65^@+4��:3tr3+��                                                                    5��5�6"�5ůq5<�f4��3��#3	�                                                                    3�j3I�&3�a3��>3�|2h�1pxf0��                                                                    5���5ã	6F<5�5f�}4�0�3��3(��                                                                    6�Z7?zo7���7Y�6`�5�$�4pa�4�                                                                    3��4&�4��43N3���2�
m1�+1�&�                                                                    6�Z7?zo7���7Y�6`�5�$�4pa�4�                                                                    6��6���6�z96^�i5���4��3��3e|,                                                                    3��3�@�4��3��3��2h��1b^Y0���                                                                    6��6���6�z96^�i5���4��3��3e|,                                                                    ���M2}}+5��5Q�4�7�3�2�T�2��                                                                    4� `9��K;�#�;7�:�G�9�o94�68��0                                                                    3�$R3�B�4Rf�3���3a�P2��H1�1�9                                                                    5=إ5��~5�r�5o߰4��t3܇#2���2P�                                                                    �H1��e��z��h�f/v��q[�P>���>                                                                    ������^��l���Y�X�Ǳ�t��P4���{k                                                                    4<ҍ4�KT4�$24n/O3�n�2�mr1���1Pvo                                                                    �D赲$��-8k������?̳,3��W�                                                                    ���Ҳ4��"̲?�����<��r�����u�                                                                    �e�ɳi���mb����4���Ⱛ�5�X�I                                                                    &��+R�W.T�-Q`y,�Ե+�2&(:�/�n�                                                                    )BE.�Ķ                                                                                            6o7��7���7�6e�F5�Ky4j�!4'�                                                                    2��3�2�4(�3�On2�T�2��1�`0�1                                                                    5j4�6�J6��%6� 5m��4�1�3���3,Q�                                                                    1���2�8N3C2�8�1��N1sL0Ӎ/�1�                                                                    6eU�7
�7�߷7%6u��5�}4��4;{�                                                                    2��3�=�4]�3�!�31�2%��1��0�^�                                                                    5�>6�$A7J�6�6:�05Pt�4(��3�m                                                                    4+�4��%5g��5�n4Uy[3n<_2A �1��S                                                                    4��75gW	6#�5��55��4(s3�2�w#                                                                    2�3�1�4;a13���3,�2@��10�q                                                                    6
#�6��l7w�:7�H6dLs5~Ǭ4N��4mL                                                                    4��4��5���5#3�4�t�3���2li23                                                                    2�P3Z�{4_e4[G04/�(3Ŀ.3(�;3|V�                                                                    1��.20��34:�311�3 �2���2�L2K�5                                                                    2��3��h4�L>4� �4V��3�w�3N�3�5*                                                                                                                                                                        +a      ��     � 6O��    F�%     =�    ;T�    6�&3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @��     @��     17:07:43        