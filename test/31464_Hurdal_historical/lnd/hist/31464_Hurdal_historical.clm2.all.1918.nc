CDF      
      time       bnds      lndgrid       levsoi        levdcmp       cft       glc_nec    
   ltype      	   natpft        levlak     
   nvegwcs       string_length         levgrnd       hist_interval            +   CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       CF-1.0     history      Sun Jan  9 16:23:28 2022: ncks -A /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1918.nc /nird/home/ecaas/31464_Hurdal_historical/lnd/hist/31464_Hurdal_historical.clm2.all.1918.nc
created on 12/10/21 16:45:13    source        #Community Terrestrial Systems Model    title         CLM History file information   comment       :NOTE: None of the variables are weighted by land fraction!     hostname      saga   username      ecaas      version       ctsm5.1.dev043-6-g5ae72ca      revision_id       9$Id: histFileMod.F90 42903 2012-12-21 15:32:10Z muszala $      
case_title        UNSET      case_id       31464_Hurdal_hist_for_decomp   Surface_dataset       "surfdata_31464_Hurdal_simyr2000.nc     Initial_conditions_dataset        .31464_Hurdal_Spinup.clm2.r.1201-01-01-00000.nc     #PFT_physiological_constants_dataset       clm50_params.c210528.nc    ltype_vegetated_or_bare_soil            
ltype_crop              ltype_UNUSED            ltype_landice               ltype_deep_lake             ltype_wetland               ltype_urban_tbd             ltype_urban_hd              ltype_urban_md           	   ctype_vegetated_or_bare_soil            
ctype_crop              ctype_crop_noncompete         2*100+m, m=cft_lb,cft_ub   ctype_landice         4*100+m, m=1,glcnec    ctype_deep_lake             ctype_wetland               ctype_urban_roof         G   ctype_urban_sunwall          H   ctype_urban_shadewall            I   ctype_urban_impervious_road          J   ctype_urban_pervious_road            K   cft_c3_crop             cft_c3_irrigated            time_period_freq      month_1    Time_constant_3Dvars_filename         :./31464_Hurdal_hist_for_decomp.clm2.h0.1901-02-01-00000.nc     Time_constant_3Dvars      /ZSOI:DZSOI:WATSAT:SUCSAT:BSW:HKSAT:ZLAKE:DZLAKE    CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)    history_of_appended_files         �Sun Jan  9 16:23:28 2022: Appended file /nird/home/ecaas/all_sites_decomp/31464_Hurdal_hist_for_decomp/lnd/hist/31464_Hurdal_hist_for_decomp.clm2.all.1918.nc had following "history" attribute:
created on 12/10/21 16:45:13
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
>��>���?z�?L��?��?�{?ٙ�@�@   @?\)@e�@���@��@�ff@�{A z�A�RAU>�A��sA��>B'�fF�& @�=     @�D�    $��      a     3�2�u�6�\�B'�fB'�f7 �9D"g�        +��56wZ�?u�?u9�3gϭ0�F00*TV/��                8�x�1�y+<4.��5V��<v�z;E�3��):+N�< ��<��j=F��=���=��+=���=�v;=��=�Ѓ>��>�8�x�1��+<4+!��+S.��6y�+>,�	@e��@l`�@b3�@bIz@b_l@bwy@b�$@b�%@b��@b�@b��@b�"{@��{@��{@��                    E���4F]G�\�GŭFgBE��eE?HD?$�CN�|BMڿ                                                A�oC$�B�.�A�gKAP��@��?�J>�oN=ף                                                                    E@b@č�F�	A�^�1�Y�F[B*��<�D/,���                                                    {@�ο뵹�뵹@^V�7)X)%B�%������ar/<$~��6C>��8    �8\-�?2J    >�y�=f�D1۠"?z���4�q�8+z        BR�BR�C���C���Cj� ?�  >SҞ3	�QC�76)k�6��SA&k�A��>�y�AU�[@CupA��@g�F@ ��A4��@�(�A�I���@Zw�    ���    �����A:5�	?Z�?Z�@W��?8q�@z��?:�0?(i�@h~u?�3�F��1�Y�0<{5,�d�2!t6&��    3S�F#�FN�LF�� F�]�?�b    B��[@�U<?��`?O��?!�F? X@>�H>�>��>�H�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M���    ������{@��5��        2�3�}�5�25��5>`H4޾M                                                A��
�WC;�WC;{@��C_��@��;B0��@[�6=2�{@��{@��    C�\xBW�5��    6GS�7}w�@���3	�d6���    >Î�C[�>7��%B]`ED��fC�hB��B��A���@��<@�?�                                                ?ze�A�$�@�I$?͑q?V��>��>(�U=T��<W�G                                                                    D|�	F�cFŵD�_D
�C_��B��fA�~�@�s�                                                @��C
L1B���AF��@�t�@4�G?��V>���=�l�                                                                    C��RE��SEi��D
;2C[�$B���B
,bA)޽@)��                                                @\[$B}��BcF@�&Q@4�?��b?7�>?G=?��                                                                    4l�U4�]cA1�>=0A��`?��?�P<)yԄ(��+�bH+�s)yԄ(G�7)���+w��+Q�u(�G�7 �9                        ��Q�3	�Q1�b�            6�Q��5H+��Q�2^du                        >��Ҷ�!                                                                2��2            2��2{@��,�xQ>L��>L��>L��>L��>L��>��0? ��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A���G���>*��A��                                                B�                                          B��                    A�                @�p�                2�u�2��2!t=C�6��=�l4��t;K�;��                �L��        +���        {@��47H4�)�    �c��{@��{@��{@�βҸn7�X4            ,��L    ,��L{@��    ,��L    4��1�r�{@��1� {@��7��6l�w7�E    4Øl{@�ε�_*2��
-�V-�V<���    B�6vxWFt�D��>��=�uA�t            ?���    :.�.    >xڊ=���@�rY@���@r�?N�=>���>I��C��"���ɾ���Ŋ�'���Ԋv"��T���;�Z�'��w��
.�?�g�>xڊ=���@���@�?A@
f?Mr:>���=��Ջ�C��"���ɾ���Ŋ�'���Ԋv"��T���;�Z�'��w��
.�=iU&    X�15#p0J7=<�x{=,�S<�b;��s=��> ��                                                �!���C����Ï�o�ק2Ęn������̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6�l�4�6�?�~p8��6�@�k@�-?���@���=�d�7Ǥ45��F�SX7�#�>�@B��Y7O��>�z7�#�J�iX4� 57ȯCMHD�|�D�V�DM5C���C<xB]�#A�D�@��@                                                A�HCQ|�C@V�B�5B6��A�<x@ݙ#@D�? �@                                                                    E�j�G	kSGL3F�<iF?EHE�4�D���C�ܾC��                                                C��mEG��E<�EC5D�C�c�C�B.q�A<��                                                                    E諂F�WJF�ܶFzqF\��F.��E�.-E��eEH�                                                D)6�D��TD�̬D��$D�iAD~5�D0��CߙMC�x�                                                                        4�DGR��GI=EGY�F��F}8F3iE�CvEPs+                                                A���A�RUA�r�A�7A�J�A{y7=*(�                                                    >�/�@a��A�aAV0�A�rB� Bh��B��b<�<�<�<�<�<�<�<�<�<�<�<�E���E���Ef��E�ID�̰DX	yC��C���                                                {@��{@��{@��{@��{@��A�J5B��;    �f�%6��
B���@    ��P¾�P�{@�ξtċ�tċC��q{@��C���C��.C��q{@��@^W{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��GMp>��qD��mG� D��ED��D��A=A=F[�eB3��C	_3C[N�FK��FK��D��D��F]$B3�q                ?;�C�mC��2C��|?��C���C��8C��8C��NC��C�	�C�L&C���C�ǎC�,C�K C��$C���C��C�.�C�l�C���C���C�5C�A�C�YBC�aC�X�C�*�C��C�BC��C��i{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��,C��F>�4 >�4�>���>��>�>�N�>�]>�{q>���>�\�>�W�>�}�>��?>�$�>��?>�1�>�6�>��V>��>�e�@/��@e        @ed���f�������\���������J���`�����1���1���1���/P        =M0�A���    � �� �{@��@���FO�k6��27v�        8���        :C�<G� >��>�]r?#�?(oD                                                ��_55)}�A�  ?k�    <#�
$��      a     3�        1��3	Qa4��4�.<4<r3�O                                                                            0�r+2-tD3<�3��>39��2��                                                                    G�\�GŭFgBE��eE?HD?$�CN�|BMڿ                                                        -#�.��</���0ʓ/�A�/s�                                                                            ,M��-��.���/6��.���.F��                                                                            0l�31���3��3�X3:"2�;�                                                                            0_�81���3!K3�K3c�w2�/                                                                            .7�/�`�1 �1���1?B�0��/                                                                            0��"2	[L3D�3��?3�S�3��                                                                            1�F42��n46~4��`4=��32                                                                            .��%0)&�1xu�2f�1�OQ1?��                                                                            1�F42��n46~4��`4=��32                                                                            0��2Hn�3~f4ٛ3�|�3&�-                                                                            .%'1/�OS0�%1�S16��0��*                                                                            0��2Hn�3~f4ٛ3�|�3&�-                                                                            /;N�1��2W�2މ�2��-2ru                                                                    9qZ8�ٰ;A�;%�=:��9��9K%8�?                                                                            .�j@0�1A��1�o�1�\(1
L�                                                                            0	mp1R_�2xO�2��J2�P�2�E                                                                            �`@W�������S>�>;��"O                                                                            �l)g�U�� ��������L��                                                                            /�50Ra�1yyB1�q�1��E1og                                                                            �,mb������,�O�)��l���r                                                                            ��U��"��X���Vᯖ#��0~=                                                                            ����E�n<p�����}��S                                                                                                                                                                        �����G                                                                                                    1���2�< 4��4��f44�U3ԛ{                                                                            .��/�?U0��1$C�0��e0p�_                                                                            0��1�nh3�a3�t3FO�2�B                                                                            -
��.�R�/��d0,g/�/��                                                                            1�A�2���4! 4���4W�|3��                                                                            .��/~̐0��14�20�g�0��                                                                            15�u2���3�~4[�G3��3���                                                                            /O��0ް�2�2{ Q2��1��                                                                            0�1�uN2���31y�2˘�2a��                                                                            .(�/���0ؔ	1J�L0��0���                                                                            1^i�2�'�4P[4�7H4�y3��                                                                            /~/�1�2#ɍ2�c�2/�f1��                                                                            .?��0-e@1���2�2�!F3&
a                                                                            -�3/0���1���1�Q02,�                                                                            .j9�0S�1���2�V3? 3J�=                                                                                                                                                                        $��      a     3�6:�`                    >EQ�    6�S�3�
2@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�=     @�D�    16:45:13        F�^ @�D�    @�K�    $�      a/     8�2Q�6��B'�fB'�f7/VOD"g�        +�Un69�U?l]d?lZ�3 Sq0�>�0�Od07�/{O                8���1�;+<4+!��-Bv�4���<N�}4�u�:4+�<!I�<��=F�=��=���=��_=�w�=��=�ш>�>P8���1�J+<4+!��+S.�1�6�=�D�@a�N@lB@b3A@bIq@b_e@bwo@b�@b�@b��@b��@b��@b�{@��{@��{@��                    E�E�3�fG�sGVbFhLE��`E�CD?�CO-�BN/�                                                A��C%�tB�ƟA�~8AQ�e@���?�*�>�=��L                                                                    E@].@�i�F� A�A	2WFZ��B)��<�D/-��H                                                    {@�ο�v���v�@Z�u75!-"��ۤ�6���8�/N!J��>�]d    ��z-��@24�2    >�y�>��1Ɍ>���&_S�!�Q        B0��B0��C��iC��iCyo�?�  ? �3}8C�_b5���6���A$-�B��>�y�B,h~AV��Bp#�A
	A_%B,��A@q�B9����ſ��V    ���    �����p>�>Z?��?��A'{)@:kDAH��?Ԛ�@/�AD��@IN�F��Z2W/��@,��J2�m6��    2���F�E�F�?G-Gܬ?p�I    C�W@�P?��?O��?!��? �:>�ow>��~>�C�>�H38��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��    ��U���U�{@��4�[�                3�ϒ5 �n5Cƅ4��                                                A�~���J���J�{@��C�^Ag}SB?(@T*�=�S�{@��{@��    C�x#BW�5��%    6D��7|�L@�de3�36�a�    >m}CO�
7��0Bp\D�C��B�ڣBB��A�<KAKl@�d?��                                                ?���A��mA?�R.?�b�? �v>Mf�=q88<bq'                                                                    D �F���FO�D���D��Ch(�B��.Aˠ@��4                                                @�*QCZB��AR�G@�?�@@��?��2>ʞ=�1�                                                                    C��EE�\�Eln�D/CdݡB�9�Be-A.��@,t�                                                @a�B7�B�m@��)@?��?�r�?iK>F޹=Cw                                                                    4'HH4��@A�>7H\A�n?��?���)��(m��+�[+3w2)��(�^s*�)�+�+�+�/�(�4�7/VO                        ��:�3}81�d            6�:���2H��:�2;/                        >�k����                                                                2��C            2��C{@��,�>L��>L��>L��>L��>L��>L��?�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B�3G��_>,��A�r                                                B�                                          B��                    A�                @�p�                2Q�2W2�m=�4�6�y>m~�5��;+pq;/ǟ                ��w`        1�c        {@��4�4�"�    �om�{@��{@��{@�έ�g�76/            )�A�    )�A�{@��    )�A�    4%sM5]�7{@��5W��{@��7E�6slQ7��    4%y�{@�εR��41�4�4�<�O3    B��$6{YuF8AD�F?��>�8A�w'            ?���    :3z`    ?zq>��E@�=K@�'f@-l?Q��>�nQ>l Y��C��"���ɾ���Ŋ�'���Ԋv"��T���;�Z�'��w��
.�?�?zq>��E@�i@��5@

�?L	t>�\�=�R���C��"���ɾ���Ŋ�'���Ԋv"��T���;�Z�'��w��
.�=u��    ��15#h0J75<�2�=4I=H�i<��f=�p>%W                                                �!���C����"Î�î)���Si�jA"�U	$̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6���4�
P?5:j8�A�6Ҫ�A	um@�o�A�A1�>��8(��6wG�7i<�>��Cߟ<n�>�7i<�J�1�4�&7}u]CM�D�AhD��DMR�C���CXEB]�OA���@�R�                                                A��CQAhC@�B�R�B6��A�XE@��O@��?!R�                                                                    E�jkG	OVG3EF�E)F?h2E�O�D��C��GC҂                                                C���EG�7E;�|EI�D�4�C݊�C�UB.|�A<�                                                                    E��F�R1F��9FzFF\��F.�1E�1>E���EH�                                                D)6sD���D��$D��bD�h�D~7�D0��Cߛ�C�yb                                                                        4Ko1GR�zGI G]�F���F}GXF8%E�F$EPt�                                                A��XAԞ�A׋6A��A�lA��c                                                        >��@_HA	\�APY�A�P�B�B[)�B�c�<�<�<�<�<�<�<�<�<�<�<�<�E��E��Ef�YE��D���DXzC��C���                                                {@��{@��{@��{@��{@��A�3�B��D    �f6^6�BvB�,@8�    @8�@8�{@�ξ(�C�(�CC�B�{@��C�C�X	C�B�{@��@Z�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G=>ɬD���G�BD���D�>�D�>�AB��AB��FZ��B2�C-COL�FK�FK�D��D��F\w,B2l                ?�B.C�=C�2C�5�?��C�|C�X	C�X	C�]C�e-C�o�C�~4C��|C���C��9C�\C�?PC�k�C���C���C��C�D�C���C���C��C�PC�4�C�E?C�/�C�SC�3C��C�_P{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�C�9�>��o>��>�B�>���>�t�>�a�>�0\>��T>���>�J>�I�>�̛>��2>��>���>�s�>���>���>�6�>���@��?�L        @������K���!���i�"A��I���:���(�Ř� ���>���?���>��0        =]��A���    ��������{@��@b�FO�86}�d7Z�        8&��                <���>Ӫ�?;�?�c                                                ��I�6B4"A�  ?k�    <#�
$�      a/     8�                2�Q4��4l�3�                                                                                    2�y32.3;{�2�M�                                                                    G�sGVbFhLE��`E�CD?�CO-�BN/�                                                                .~��/��h/���/f�                                                                                    -��.��
.�l�.1Yj                                                                                    1�t3Fk3;P�2�JJ                                                                                    2.�3j�a3�oT2�8k                                                                                    0
%�1@n�1[�0���                                                                                    2T��3���3�k�3�%                                                                                    3�4(�.4C�*3�pw                                                                                    0Z��1�
_1��31-��                                                                                    3�4(�.4C�*3�pw                                                                                    2ZK|3���3�a[3�&                                                                                    /��h1��1?X�0��8                                                                                    2ZK|3���3�a[3�&                                                                                    1�2E�2�I�2T�h                                                                    9��I9 S;E;'m^:���9���9A�y98                                                                                    0Q�-1���1�\i1 �'                                                                                    1P� 2���2�12�~                                                                                    ��͠�	[�����B�                                                                                    ���^�zݱK(�6*H                                                                                    0S�1���1�>c1��                                                                                    ���Ʋ�0��s˲��                                                                                    �.��d乯�*ϯ{                                                                                    �@���bv����ݯ���                                                                                                                                                                        ͨ���                                                                                                            2��4�Y46h3��b                                                                                    /�3�0��0�E�0V��                                                                                    1�23!��3G��2�k                                                                                    .���/� /�,�/k�                                                                                    3U�4)44Y=�3�Pu                                                                                    /�Av0�c'0�0�-�                                                                                    2���3��3�Ls3x��                                                                                    0��G1��d2��1�"A                                                                                    1�PG2�$#2̯�2H�O                                                                                    /�7-0��(0��50e�                                                                                    2��4;T4�*3�                                                                                    1��2��20�01��                                                                                    0���2J��2�i03ғ                                                                                    /���1$�1�Z,1��                                                                                    0�D2x3d34�	                                                                                                                                                                        $�      a/     8�68��                            6��3�'�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�D�    @�K�    16:45:16        F @�K�    @�S�    $�q      aN     >�2kR�6YdB'�fB'�f7=̝D"g�    �+��@6�?r�?r�3/ޏ0�s�0�)�/�/}z         Pi�G��8���1�b�+<4+!��-���5�<�-�4�r :;RD<!�H<�4�=F�=�ڶ=��.=��c=�y6=�F=�ҵ>��>�8��$1�S�+<4+!��+S.&��5��=qX�@]�I@l�@b2�@bIj@b_a@bwf@b�@b�@b��@b��@b��@b�Ap�;�db{@��                    E��3�2�G��G��FiU�E��4E��D@�`COñBN�R                                                AB�C&i�B�^A���AR��@���?��>٫�=�e�                                                                    E@1�@�?F�YA�1}��FZh�B(}<�D/.v�j                                                    {@��@�!@�!@Vʡ7D*%(��$\g$��c�/3c&�$>�v��>S�-�b�2,I�    >�y�@�B�1�A@
Y�@X⿃L;��	    Bz1�Bz1�C�ɨC�ɨCq�?�  ?�||3@�C���5�Mr6�L3A!��B���>�y�B�%$B��C�AR�AڨBѢ A��?B���@����m����@��    @��AK�>A+$1?~�?~�A�%S@�׎AŒo@)��@��A���@�=�F���1}��/,�:x1�_�7�O4i}�2��.Fؙ4GKLGG�aGn�b>H��@Q��C1@�)P?��?O�B?!�r? ��>�U>�b�>�.�>�b8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�[�    �&(��&(�{@��4ˮ                4��5_>�5A�4���                                                A�pY��������{@��C�"A�RRBx� @IF>[�N{@��{@��    C�+EBW	5���    6A�z7|�@���3�6�\    ;q[CY�7���B�i�D�s>C���B�0�BfiA��-A!DI@/�?�                                                ?���A��A�3@8�?���?�j>m�~=�m+<m>                                                                    D��+F�K)F�OD�KD$�Cp��B��mA�_�@˸�                                                @���C��B�GrA^�d@��@L|�?�tO>��$=�c�                                                                    C�-�E�3xEo�D�Cm�Bȉ�B3?A3��@/<#                                                @gUWB�c�B�@�x�@Ky�?�@V?!��>N"=G�                                                                    4=�:4��aA
�Q>1YA�(�?�3?�'�*���)q7,��+�cF*���)M�J,��,5�+���) �87: 2���0�5�- 
�            7�3f�i1���0�u�0�5�- 
˷���+7�2ߐ    &��2�y�2�W�/�M    >cTg7
l 6�G0��/��1-�a                        *:	0���0�iK.��6    0��_2��            3&yN:7�'-��>L��>L��>L��>L��>L��>L��?��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C�G¦s>-}]A��#                                                B�                                          B��                    A�                @�p�    �>�ݛ�2���2kR�1ʁ1�_�>�J7r�e?.۸6��p;9<�;;��                61�    4;��3�    7B��{@��4
R5�m�6������{@��{@��{@��2GI�6��f            7oC    7oC{@��    7oC    5SN�7�T�{@��8qt{@��6�_�5�DU6$r    5W�7{@��3��5��4���4���<��    B�t�6~�F!%�D�@��	?��B�a�            ?��%2���:5�2���?^M�>܄`@��Y@��@�?M	�>�)�>��
��C��"���ɾ���Ŋ�'���Ԋv"��T���;�Z�'��w��
.�?���?^M�>܄`@�',@���@�v?KHC>�&=�����C��"���ɾ���Ŋ�'���Ԋv"��T���;�Z�'��w��
.�=r�    $a�0�"M0�B<�,<=4=��;��m7�5A>E�                                                �!���C�Y�ɰÎ�%î������ŁQ��w+̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 7! 5K�>F��9%Ki7��A���@�YdA���A���>-r08v�_6�DfGi�6���>��C�U@��>�G6���K},7U
26˵9CM"|D�2D���DM��C��QCv%B^K�A��x@�ρ                                                A�"|CQ2C?��B͕�B6�QA�v%@�K�@�x?!ρ                                                                    E�jG	3dGqF�M�F?�E�j�D�$�C��zC�J                                                C��[EG��E;�^EO�D�M�CݲeC�BB.��A<�                                                                    E�<F�MF�ӼFzF\��F.ǇE�46E��oEHF                                                D)5�D��D���D�ޟD�g�D~9~D0�'Cߝ�C�zb                                                                    0��_4Y��GRb#GIGa�F���F}V�F=E�IEPv�                                                A���AԇTAׄA��A��A�t�                                                        >�p@`�A	j�APi�A�żB6TABW��B�f�<�<�<�<�<�<�<�<�<�<�<�<�E�ɜE��^Ef��E��D��1DX�C���C��*                                                {@��{@��{@��{@��{@��A�/�B��U����f��6�1AB���?���    @��H@��H{@�ν�]���]�C�oZ{@��C�ޗC���C�oZ{@��@V��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G?�>��D��G��D���D��4D��4AHJAHJFZ5B0��CJqCSAgFK��FK��D�_D�_F[�UB0{�                ?���C��<C���C�F�?�*C�.C�z�C�z�C�|�C��0C���C���C���C���C��C��C�%.C�J{C�s?C���C�ϳC��C�A[C�|C���C���C�JC�)�C�.�C��C�4C��C�}�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�A�C���>�Lk>�eq>���>���>�V>�o�>���>��>���>��n>��\>��>���>���>���>��w>� m>��>��">��{?�'\?��@        A6���	�8�	o��	\i�� �~��}ڲ�}�mƯ�����������������        =�7B�&    �e��e�{@��@/��FO_�6D�76�        7��                =��?O�?��?,�                                                ���7El8A�  ?k�    <#�
$�q      aN     >�                3�4>h�43��Z                                                                                    2(��3p�(35��2���                                                                    G��G��FiU�E��4E��D@�`COñBN�R                                                                .���/Ǆ�/��/��                                                                                    -���.��.���.*x�                                                                                    1�k43G�635��2�]�                                                                                    2au�3�83�e32�-Y                                                                                    02�A1�O1h��0���                                                                                    2��3޶F3�&Z3b�                                                                                    34�4j�;4A�3���                                                                                    0}��1۶R1�K21*`�                                                                                    34�4j�;4A�3���                                                                                    2x�_3�f�3�Ou3��                                                                                    /�:�1W�1?�_0�@�                                                                                    2x�_3�f�3�Ou3��                                                                                    1F�2q�2�.q2I�)                                                                    9��9z�m;D�.;&&�:�^C9�B>9Hd;9�                                                                                    0�D1�!	1�~�1c                                                                                    1j�2���2�&@2 �                                                                                    ����8�q����M                                                                                    ��,���>D�J�.�                                                                                    0nN�1�|�1�u�1�                                                                                    ���سZt�����ț                                                                                    �?!b��IS���{���                                                                                    �Rę��d/�yگ�s�                                                                                                                                                                        &���%�E                                                                                                            3E�4O��40�[3��"                                                                                    /�d�0���0�5�0NȦ                                                                                    2	��3YҶ3A�T2�*y                                                                                    .�XI/�&/ۇ+/b��                                                                                    3�
4d�4R�M3���                                                                                    /�K�0�9�0�ز0v�~                                                                                    2�y�4�3�#�3n�                                                                                    0�,2+p�2k1�a2                                                                                    1��2�p�2��2@�2                                                                                    /��s1
��0�d;0\i^                                                                                    2���47X�4ι3�ٹ                                                                                    1�)2Q��2+5f1���                                                                                    0���2��v2�Kb3�I                                                                                    /�O�1\��1��#1�4�                                                                                    0�R2��;3�3-VY                                                                                                                                                                        $�q      aN     >�65�J&��                        6�GD3��6@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�K�    @�S�    16:45:18        F�� @�S�    @�[     $��      al     D@4���6astB'�fB'�f7ـD"g�    %R#+q��6�?�  ?~��4��3�	2��71��1�� 0��    +���'��c=�I=�f�<�K:<�{5��%<G�;e`3�r
:@�d<"cP<�R�=F�=���=��g=��=�z�=�b=�Ӯ>�>>?���>���>3�=�`|<���<c�n<�'=D��@Y�8@k!@a��@a�@b2�@bd�@b��@b�^@b�@b��@b��@b�'AIŦ<3��{@��                    E��{6���G�V+GiIFjXE���EM�DA�CPcEBN��                                                A��C'�B��A���AS��@�FJ?��Q>�R�=��8                                                                    E@ @�gF��A� 2��WFY�nB'�@&�/�                                                    {@��A�A�@RJ�7�U��E�$Jg����/i���	�>�1%R#@�!�/��4]�_    >�\�A[}<1���?�Q�B�AL�>-4�    B0��B0��C��C��C� �?"m�@W�3>C�"n5���6�� A��B��>�\�C۽BF�CM�A|��B�B���A��B|�Al�(�p:�g�rAl�(    Al�(A�A���?'?]A~��@�QfA���@2�@DRA%�@3�$F��2��W07W-�w2j�8)��5��64�NUGQ�G4��G��xG�{G>~�RAy@)B<�@qL�?K!�?CN?!��?�>�ŉ>��>�.j>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M���    @y�@y�{@��7c*8�~�8�S�7�>a6�r�5�AJ6
�Q5G�4�u                                                A�`(?�rh?�rh{@��C�:B#�B�cB@<'d>� 
{@��{@��    C��RBV&5��    6> 7{O�@��3-�7~    >&_(C���7�y�B�S.D�I�DX"CV<B�+�A�(6A2f�@=pk? Ȁ                                                ?�eA���Ap@"~�?��1?0��>���=��*<sݖ                                                                    D���F���F �D�B�D��Cx�B��fA�+�@�7                                                @���C0GB�D�Ai��@�l�@XYb?��>�=���                                                                    C�6E��
EqS�D^KCv��B���BU�A8�k@1�_                                                @lIB��B�X@�b@W �?��?)�>U[�=K0�                                                                    6� z7��AG>+��A��?Ji?�g	+6�D)�2�,�װ,I�+6�D*!�,�$M,ɞv,�b�)��B7e�45Ik2�/�            7��:4U��1���2@�2�/����:���d7��;1�ߢ    )�d�4,E`4*c�0�ʍ    =�M7���7q��2u;�1*�/H5                        ,m�P2Gv�2;�0E��    2uJ�2��E            4E��:^ߩ-Ǡ?eC>�I>�2X>��
>n]L>���?�?|��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C��G��">,��A��1                                                B�                                          B��                    A�                @�p�    +*��'�I4Iht5��U2��2j��?8��7�b�?��X7f�U;�~�;�!�    /&"S        7ܨ    5�-v5�9�    8�#l{@��0��4���8?����-�{@��{@��{@��6w�*7ly~            8��    8��{@��    8��    ���6o%a{@��8z��{@��6J��5�O�5��    58{@��6F�6�7R
�7R
�<��.    B��h6�M�F��C��AYǩ?lӓBʪ            ?��/45Ip:�@45Ip? ڬ>��@�u@p2H?��T?@[>��q>��_��D��"���ɾ���Ŋ�'���Ԋv"��T���;�Z�'��w��
.�?��J? ��>��N@��4@n�?�wD?@>���=�����D��"���ɾ���Ŋ�'���Ԋv"��T���;�Z�'��w��
.�=�>L    -��7<��6��=<�@�=	j�=�p�8Pz�8ܭ�>\��                                                ���������#D?ÆjNÝ���q�0|��*��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 6��
5��s=�HT92�w8:��AT��@��A�zA��}=��88a6�L�FG�m6~n =њ�B4A6?�)�=�S�6~n J�c�8��17G[cCMx�D�ypD��SDN`C�DC��B_>�A��w@�q�                                                A�x�CQypC?�SB�`B7DA���@�>�@�w?"q�                                                                    E�k�G	�G�F�VTF?�nE��qD�;wC��Cܔ                                                C��EGqE;�EVD�g�C���C��B.�A<��                                                                    E��F�H�F�ϷFzF\��F.��E�7	E��EH	�                                                D)5�D��TD���D���D�g6D~;|D0�5CߠC�{]                                                                    2uJ�6{��GRM�GH� Gf�F��=F}f�FB2E�LWEPx�                                                AY�Av@�A��`A��>A�. @8�                                                        @�]�A[�dA�g�A� aAѺ�BD�Be��B�f�<�<�<�<�<�<�<�<�<�<�<�<�E���E���Ef�'E�vD���DX�C���C���                                                {@��{@��{@��{@��{@��C�B���*�_�f��7X
$B��@ ¤    A
�RA
�R{@�ξ?�,�?�,C���{@��C�qkC�oC���{@��@RJ�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��Gb�>��QD���G�vD��mD��>D��>AM�AM�FZ)B/�`B���C���FK��FK��D�:D�:F[vxB/�I                @E��C�
:C�
�C���? �C�{�C�V6C�V6C�
�C��LC���C��3C��[C��]C��C���C��C�6>C�Y,C�C��C��]C�$C�H�C�C���C���C�
}C�'�C�C�DC�C� t{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���CΓ&>�f�>��>�N�>�U�>�w>��>�Q�>���>���>���>���>�U8>�K>��=>�+>�0>�!�>�,>� p>�n#@��?��        A���-��,��,�����iǞA�ǝ�kǝ��ƀ5��z�����yķ,|        >q�pBy�    @qB@qB{@��@m��FO%(6LE�7�        :�6P?	�V>���>�*�>��q>`�u? �$?�? �                                                ±�Y7���A�  ?k�>û�=0�9$��      al     D@8���7�:7
M�6J4��}4��4g�3��                                                                    7�l�7=h6.��5#�3�ͫ4O�38�&2���                                                                    G�V+GiIFjXE���EM�DA�CPcEBN��                                                4>d�3~W�2���1�O�0C��0r��/�c
/�@                                                                    3p2��G1��0��{/w>/�O�.��x.9(                                                                    6��6���5�N�4ӨU3���3���38�B2�|�                                                                    7��6���6.��5E��4;��4t�3��2�                                                                    4�4T0:3��3(22H��1�	�0��                                                                    7/nc6�e6U�-5q��4efb4���3��(3�                                                                    8"�7��7�62�4أ5�24I>�3��                                                                    5D4�
�4:�3A�h2=�[2���1�K1;ޫ                                                                    8"�7��7�62�4أ5�24I>�3��                                                                    7sh74�V6w��5U�J46(�4zq�3�6�3!3�                                                                    4x�$4K�r3��Z2�O�1�=2܏1Ix0�{�                                                                    7sh74�V6w��5U�J46(�4zq�3�6�3!3�                                                                    �(y���3��3�C2�	�3�>2�B�2Yn�                                                                    9���9,��;#ߜ;��:���9ڑ99V�!9%�b                                                                    5MR-4�Gz4c��3w�,2bP�2�t:1��1��                                                                    6�6TA�5��4U�
3)?�3_D�2��52                                                                    ��丵S`��S���iҲ������� ���g@                                                                    �T��&�y�ΌW�E_}��I�K�^�PT�=kE                                                                    5���5S�`4��"3Y?�2,q�2ck�1���1<                                                                    ��]*�X~ߵ�)	���ҳx�泭`C��t���S                                                                    ����Ĳ�˱C��:�;�߯�_��$��                                                                    ��� ��?��!׺�7Y�|��8�ݰ|v��g�                                                                                                                                                                        0��.��h                                                                                            7�z�7��6��5�O4�{�4��Y43�O3Ơ*                                                                    4u)41��3y�2pf�1F.
1���0��
0`��                                                                    6퉽6�&Z5��64��3���4�w3EM2�ʑ                                                                    3p*33.�2t�C1n��0M�0��/�}�/v��                                                                    7��7���6��5�9�4�a�5
�B4V��3���                                                                    4k+_4*n�3o�}2l��1S��1��.0�(�0�4�                                                                    7fI<7+��6�P[5�{�4�2�4��3���3�Z                                                                    5���5D�4��3�D)2��2��X2�1�ԩ                                                                    6:�6
�\5��4��F3i�3��*2���2Q�                                                                    4T��4hi3�!42��1�+a1��T0�/0n�                                                                    7���7Q��6·�5ܗ4�=�4�8f4�h3��                                                                    5�հ5o��4ވ�3�k2�k72��P2-��1���                                                                    4,��4003���3 )�2~�33%g2��(3�y                                                                    3��2�|l2�q2l�1M�>2��1�q�1�k                                                                    4S�4%:�3���3C�T2�ü3J(�3�3;��                                                                                                                                                                        $��      al     D@62R�)�d�                        6��j3�(�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�S�    @�[     16:45:20        F� @�[     @�b�    $�9      a�     J5�GQ7:�B'�fB'�f7�lD"g�    )"F+KI6�?�  ?tqe    3�D�)"~?/�d����        /X}�+Q�g>�vl>:�=��T=$vm:���8��v1(�2Yn�:B n<"��<�i�=F�>=��=��+=��w=�{�=�W=�ԅ>�r>�@��@� @#��?�=?�^?���?���?p�@e��@e��@[�b@\0�@\�@]�@]�7@^Z�@^��@_��@_�K@`;@)�7=Uq>�-=                    E�Zu7isG���G.0Fj35E��EtRDA��CP�DBO�                                                A�!C&��B��ZA�{�AS��@�nH?��z>�x@=��]                                                                    E?�.@�cFy�A���3+�dFY�oB'�A?SC/0;�                                                    {@��BntBnt@O��8f�%g�B    �w��/G�%�ׂ>���)"F@��׬���2�e    >�y�B1\�1Ş�A0QAQH�AM"�        B�ɞB�ɞC�FDC�FDC���>���@�dT3 ��CʹK6~�6��QAqCIBY>�y�C]\B�['CxP�A��B�1CR�<B��B�E�BG
����    BG
�    BG
�BLN�<T�<6���7��A�͜A>��A�ٸ@f�E@/�1@��?�~�F�V3+�d0�ٔ._�3T2�8�/V6� 
5�_�G�\G8 iG��(GݜX=q7V    8BR    > ��>�(�>Ȇ\>�S�>�(>��*>���>�m�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M���    A.��A.��{@��7��9p��9O޽8��8I#l7�<�6�&�6ғ5 ��                                                A�U�A
�A
�{@��=�gBj�aB�"@2[>�T{@��{@��    C�IBT?_6)�Z    6;�7z�a@��$2��7E��    =ly C�+	7�	5B���D��lC���B��	Bhl5A�D�A1i@0а?�                                                ?��:A���A��@�?�[�?!<>ke�=�F.<f��                                                                    D�5/F��UF[*DÅ�D��CxqB��&AԚY@�K	                                                @�&�ChB��QAc��@ې@W��?�-�>�Fu=���                                                                    C�KE���Eo��D�-Ct�B�g�B|A7�-@0��                                                @iޱB��B��@ب�@T#?�Z?'�T>T��=J�                                                                    7rD�7�sUA�c>)��Aؤ�?YD?�V	,V��+ �	-*�{,�,:,V��++�.-i�-���-VOX*T�7���4ᠵ2��U0�Ŭ            7Ѕ�4�0�2S�2�:�2��U0�Ŭ�Ѕ�2�e7ІE4(��    2�x,4�
4�I�2�5�    =�^8Y��83�5��h4���2�v(                        4���5�S�5��]3�o2�O6� �2�-            4�y�=�І-$L�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CoSG�KW>-,(A�r�                                                B�                                          B��                    A�                @�p�    /�U*��q5%͚6oh3>9�3UV�?��8D@]��8
R�;��/;�fM                7ǧ    2%�T6��    3'R�{@��    2[�    �ﵿ{@��{@��{@�ζf-�68I�            4ĂE    4ĂE{@��    4ĂE    6��[    {@��3'R�{@��3'@3�1��    6��[{@��5��Z7��(6I�6I�<���1Q�dB�� 6�8�E��C�#�A�    C,k7            >���4�x89��4�x8<��m<�l�?Z��?7�U>�9>n��>w�>
���R ������F���r̊��芑�r�u���T`n�; �'0��;�	��>f�><��j<�f�?Z�?6�P>�Q�>_�>
]=�ًR ������F���r̊��芑�r�u���T`n�; �'0��;�	��=��    'aI�6 rC6=�;q�;LP:�;<v�1<�kP=Ȫ                                                ���s���������Ķ8Ĝ������B��̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         :�ʭ        {@��{@��{@��{@��                    5�C        ;�*                CO�,D�"�D��ADQ��C�#C
�Bb��A�Ή@�+�                                                Aώ,CR"�C@�ABѪ�B:#A�
�@⠚@Ή?#+�                                                                    E�wG	12G�F�[�F?�E���D�CVC��C߈                                                C�
0EG�aE;��EY�D�lqC��%CmB.��A<�                                                                    E��F�K�F���Fz1F\�,F.�nE�7�E���EH
                                                D)6wD���D���D���D�gfD~<CD0�Cߠ�C�{�                                                                    6� �7SQ:GRh<GH��GxF�
0F}ufFF�E�N�EPy�                                                =E��;#W�        =b�N                                                            @E��AXg�A���A��B+w�B\��B�_B�6�<�<�<�<�<�<�<�<�<�<�<�<�E��xE���Ef�E�wD��VDX(�C���C��{                                                {@��{@��{@��{@��{@��D"�BUw.�QF�f�8��C��@,�    AU�AU�{@�ξ]㺾]�C�՟{@��C��C��)C�՟{@��@O��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��<I:�D��,GڿD��D�QyD�QyAI�&AI�&F[1FB1�3=�NJC�FK�+FK�+D� D� F\$fB1�                @���C��C�5�C��w?   C��C��C��C�=uC��GC�VC��7C��>C��C���C�mlC�9�C��C��C��OC��C��7C��C�,�C�Z�C��;C��C��C�tC�eC�ZC�&C�-�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�SC�4I>�[>��>�j`>��>�2	>�_�>��>��>��0>���>�7>���>��K>�X�>��H>�s>��#>�$�>�͑>�M@$/%?�^        A�h�ǘ�@ǘ3Ǘ׽��Ճ�c���#�>�Ɖ���ð"�ð8�ð!��p        ?rB2��    @�j@�j{@��@hC/FN��6�u7�        ;���?/w?0Q�?1�R?4d{?8�?>2�?E�\?P�                                                �s�7j�A�  ?k�A	��@�N$�9      a�     J9E�v8��n8�7[�n6��5ʒ�4ץ�3��                                                                    8y��7��7$��6��`5ëN4��
42�3R                                                                    G���G.0Fj35E��EtRDA��CP�DBO�                                                4��4#�43��2��2"U�1TCe0a�r/`)                                                                    4˵3N��2��b2Q=1M�0��/���.���                                                                    7Ot:72�g6� I65�
5���4֯z4�.3$o�                                                                    7�3�7Fף7��6��G6b�5=y�4P8k32�                                                                    5 ��4�gN4�ę4V�93ص{3��2/�L1&s                                                                    7���7s�7'W�6��s6&��5g��4~}�3Z��                                                                    8�d8��8��7[��6���5��5�Z4p�                                                                    5�t5��25)��4�A�4�3pm�2�'>1���                                                                    8�d8��8��7[��6���5��5�Z4p�                                                                    7�j7݈27c�W6���6��5Xs	4p.3uU�                                                                    4��m4��4�jS4��3�[Y2늳2�19                                                                    7�j7݈27c�W6���6��5Xs	4p.3uU�                                                                    �F����\5(#�5:��4��Y4�\3R�2�:+                                                                    4;{^4�x9���9��9cK 9��8�v�8���                                                                    5��5���52sz4��r4$k3_K�2p��1N��                                                                    7U�7_�6x�n5��k5�4A��3\+2Uy�                                                                    �҄�Bo���E�=L�������_��qر�a                                                                    ����ֵ/��\/����dy�2?{�޾����                                                                    6�z6Y5z��4�P�4�s3E�2_��1X��                                                                    �#��K�� }���.�Q�t���Ƴ�Ϙ�߀�                                                                    ���b��Dس�в��ձ�8��&�ɰ]ͳ�},                                                                    �A�`��u��T��=��!�M�9���H�                                                                                                                                                                        )%�)���                                                                                            8�R�8b�77�Z7Ji�6�X�5�2.5�84}                                                                    5�c4�V\4o�Y3�o�3'2r�a1�T0���                                                                    7���7^56�0�6H�t5��4�Mc4|>3'3�                                                                    4�3�3j�N2��_2,��1�0��&/�e�                                                                    8}ܜ8Y��7�\�7G�6�7�5�h�5`C45�z                                                                    5 U�4���4e�B3�x�32��2��;1��:0�f                                                                    7�D�7�$y7��#7�#6op�5�'4��^3ž�                                                                    6�5���5�75/��4�Ҙ3��v2��1��`                                                                    6��6��t6s�E5�B)5A|�4�~�3��%2��"                                                                    4�]4�:�4�)�4ܪ3] �2���1��0��                                                                    8��8$�7�,�7;��6�R�5�Y04ߘ�3�
                                                                    6/}N6r�5�| 5V��4�:H3��2��	2
s                                                                    4�T�4�:�4���4�T�4S��4��3���3j��                                                                    3�/�3�ʾ3��3\U@3+�2�đ2��2=�                                                                    4�.�4��4��?4��y4�R�40�t3��R3��O                                                                                                                                                                        $�9      a�     J60-�2�x,            G>�#    <�T`6�<3˴o@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�[     @�b�    16:45:23        F�R @�b�    @�j@    $��      a�     O�5�Rp7u�OB'�fB'�f83+0D"g�    )�3+8d78�=?�  ?q�B    3�ʷ)Fk/�!O�`�        0!I|,S>��>=�t=ȏo=I�;��^7Ҧ/���2ddq:7�/<"�G<�k�=F��=�ߓ=��2=���=�|$=��=��>��>�@��@��L@D��@��?��I?��u?�ٛ?�j�@ky@j�@`�8@`J�@`�@_�@_��@_`�@_-�@_@^�@^��=.��=~00>��                    E���7J��G���GՑFi�3EɂKE5.DAS8CP@7BN�$                                                AN�C&'�B�KKA���AS2�@�+�?ʒ2>�-�=��x                                                                    E?�c@�F�	A�Y/�z�FY�<B(�AM3.�s~                                                    {@��B��sB��s@U?"8��p�f��    ��c�/@9�$'FY?   )�3A2�����2�@    >�y�BM�1�h�A?a@�Q�@��        B� +B� +C�y)C�y)C��>�X@��3b5C��G7.6�SaA JRCZ��>�y�Cp=�B�<�C�?/A��B�rC['2B(b�B���B`s,?��    B`s,    B`s,BX�,<�    6�BA���AP��A�j�@���@-$�@�=k?씟F�#�/�z�-���+J�Y0�ٞ8�=�6�}5�%�G!t�GOw�H	��G�sn=�j-            =�8>�V>�)|>�� >�c>��>��>�(�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�+�    AHv�AHv�{@��7�`9Q�92�8���8*�i7��y6ж�6 ��50	                                                A�j&AA�AA�{@��:	g�B�CbB���@2>�?��{@��{@��    C�rOBQqI6�f�    6<L7zۖ@�RU2��7H"�    =�W^C��57��[Bk�cD�=�C��cB�eB8��A�IA3B@�P?�                                                ?�`lA�:�A��?�&�?v�s>�g]>@r1=e�<To0                                                                    D�abF�YF]�D�"DK�Cr;�B�<�A��@���                                                @�iC
[�B�9AZ�@�'@Ne�?���>�*L=�%�                                                                    C�E��fEn3�D��Cm��B�)EB�YA39y@.M=                                                @e�B~�B��@���@J��?�%�?!>Mih=Fje                                                                    7Hȑ7��MA3�>-�A�R�?3s?�@�,��+��-C�,ͮ�,��+O|�.X9*-�&U-�j*u8�7��#4\�-2��-�./=..��)��B7勢4��e1�A�2k2)�-�s]�勢2�@7��4N�    3��4TY�4S��/�HH    =F�8SP"89=6�t�5��a3>�                        5��]6�{�6���4p��3�B�7�۞2�j�            4���>�)>- @�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C~�AG��>+*�A                                                B�                                          B��                    A�                @�p�    /�5+��&5[<6[$0�/Z0ȧ�?�L�8[N�@�C�8+,|;�c-;�x^                7���    0��6�I    2��{@��            ���
{@��{@��{@��7X�7�h�            6�o0    6�o0{@��    6�o0    6�@V    {@��2��{@��2��2�o�1i�    6�@V{@��6��\7�P�7�F7�F<�Ca2-<�B�4�7 ikE�˨C��PA��    C?��            ;��C4�96
\�4�9<�qt<��4<�A?<yĦ;�;e�:��e:[�8��/�����j���Fu���;���Ҋm��M��56�!﵊]����;�SF<�m�<���<��<v$;��;[ؿ:�:H����/�����j���Fu���;���Ҋm��M��56�!﵊]����8�4    )�z5��g6K�;9"��9h& 91x�9��8��n8�]Q                                                ť�OŢ�hŠ��ŔvŅJ�j��J��2�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         :	g�        {@��{@��{@��{@��                                                    CP�iD�HHD��>DS��C��'C�sBei�A��l@���                                                A��iCRHHC@�>BӞ�B<�'A��s@�i�@�l?#��                                                                    E�~�G	G�G#?F�^F?��E��PD�B�C�yC��                                                C��EG�8E;�-E[�D�f�C���CB.�@A<�t                                                                    E��F�O�F���Fz3F\�F.�^E�7E���EH	�                                                D)7<D��BD��]D���D�hD~<,D0�ACߟ�C�{s                                                                    7�۞79^/GR��GIyG��F�1F}}hFIRE�O	EPy�                                                                                                                                @I�A*k�A���A̰(BcB;i�Ba�B��L<�<�<�<�<�<�<�<�<�<�<�<�E��oE��0Ef�E�uD��DX-�C���C���                                                {@��{@��{@��{@��{@��D+{lA�X�/�r%�g%�8 �CX&M@<<�    Aj��Aj��{@�ξ�����C��&{@��C��aC��uC��&{@��@U?"{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�g;�G�D��{G
SD��RD�c D�c AD�AD�F\�%B3ي    C�h<FK�9FK�9D�$�D�$�F]E�B3�s                @�n�C�6mC���C��f?   C��DC��aC��aC�HeC���C���C�z�C�0�C���C��8C�5.C��C��qC�k}C�*�C��C���C��0C�xGC�xC��+C��1C�ԿC��C��C�nC�@C�I�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��KC���?_>�'�>�g>���>���>ݴ�>�	�>���>̹�>ȅ�>�>�>�/>�.�>���>�R�>��>���>��>��;>�V@1|�?���        B�Hǫ� Ǫ�XǪZ��7���)���(��(MLƨf������������Ś|�        ?v�B9:B    Aw�Aw�{@��@|��FO7G�7        ;�KB?0�?��?)?�?�]?9�?��?ݳ                                                �#�7/��A�  ?k�A	��A z�$��      a�     O�9+_v8��7�q�7A�6���5��(4�VV3��                                                                    8Xx�7�^�7�6s�}5��C4�΃3�z�2                                                                    G���GՑFi�3EɂKE5.DAS8CP@7BN�$                                                4��.4�W3nf�2�M�2U�1>��0K��/F��                                                                    3��U31�a2���1���16Q10p�(/��.z��                                                                    74�H7=M6���6!��5��.4��=3�
D3�                                                                    7m��7!��6���6E^:5���5q@4 ~3�G                                                                    5�{4�Ԫ4y�4��3��2�M2��0��}                                                                    7�p�7E��6��6q:G5�E�5*m�4D(Z32e�                                                                    8�R08|�7��'79��6��[5�L5B�3���                                                                    5l!
5d�5'4�L�4V3OR2�91xޢ                                                                    8�R08|�7��'79��6��[5�L5B�3���                                                                    7ֳq7��N76�E6�ݏ5��5<�U4`�*3V�                                                                    4��4ӳ4~rq3�;u3}��2��2 �A0�F�                                                                    7ֳq7��N76�E6�ݏ5��5<�U4`�*3V�                                                                    �#��2��5�a 5Jn4�4�3M@�2���                                                                    5+�35��7	�76�yo6h�5�r5g�T4�Ks                                                                    5���5^��4���4x#�3�n|3$aZ29br1(o�                                                                    7�6�j6I{�5���4�V4*,�3B�#2;r                                                                    ��)���˵�K�&�E��a���@(��oz��'                                                                    ���A��M��z����4�hB� S���^��i�                                                                    6i5�]�5J\�4�-�3��3,ő2E�Y1=�P                                                                    ��l��O������<��=�Ĵ��ֳ���Ư�                                                                    ��~4�s������h�y���Y����IhG�aV                                                                    �^7l�A�P��ӳ�o��؛��K�'��2/�                                                                                                                                                                        +�|
,�=                                                                                            8ffh8C��7ж�74�Y6���5�d4�4�T                                                                    4��4��{4S�3��35~2].�1�i�0���                                                                    7a��7?��6�v63G�5�qN4�/o4�3��                                                                    3�3�3�̊3N��2��2wZ1h#�0��/�]�                                                                    8]�8;�@7�5a71�s6�\�5�Kz5�n4!�E                                                                    4�th4�ę4Jl*3�h`3!�62s�1��]0�.�                                                                    7ھ�7���7�|7I�6U9�5��4�73�g                                                                    5��Y5��5�H�5��4s�~3�tt2���1�u�                                                                    6��\6�t6T��5ڦ 5,M�4g%�3�fz2��P                                                                    4��4�`4r��3��n3D��2��1�u0���                                                                    8�W7��7��B7%Zo6�M�5��Q4�Ę3�a�                                                                    6�6+G5��5<��4��?3��92旉1��                                                                    4��n4�(�4�W74p/q4<x�4r3��3Ph@                                                                    3�j�3r�L3_��3B�3L�2�B�2�̑2(h�                                                                    4�G�4��4�C4�ǩ4fZ�4�3�X�3~�N                                                                                                                                                                        $��      a�     O�60��3��            G��5    =�a�6�{3�p�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�b�    @�j@    16:45:25        FÐ @�j@    @�r     $�      a�     U�5���7�|-B'�fB'�f85��D"g�    )�O�+7�k76/�?�  ?uq    4��)/*�/� A��m        0 n!+��n>��w>!�=�%�<�wy9�2�O�+��2`�d:*j�<"�<�bW=F�d=���=�ĭ=���=�{�=��=��F>��>@ݍ�@���@9��@
� ?��V?��?�fy?Í-@l+v@l$@a�P@a�-@a�
@ac@a6�@a2@`�M@`��@`|*@`e�=f;=�} >�8                    E���7s:hG� ZGw�Fi+bE�!�E��DAXCPdBO`                                                A��C%��B��A�f�AR�Z@��?�P4>��S=���                                                                    E@.@�7xF�LA�0Q�GFZapB)�0AI)�-�Wo                                                    {@��B��B��@c��8�Фԓ    ����/HI����?   )�O�A%]`��[2�o    >�y�BIT�2�A$7�A
X�A	�T        B�C�B�C�C���C���C���>���@��43��C���7�6�ITA#�xCB��>�y�CU��B���Cd��A��B��C@�B��B�[~BLN���8    BLN�    BLN�BUW�            A�Q�A<m�A���@reA@Sn@��?�<F�E�0Q�G.�W+���0以8�G6�1�5�MG�GH��H�rG��=���            >�C>��`>�>��>���>��Z>�S�>�(38��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���Mˈ    A�e�A�e�{@��7�s[9}��9W��8�I<8FO7��6��"6�f5j�                                                A���A��2A��2{@�ΠD5B� 8B�B@@.C?o�{@��{@��    C�giBRb(6�)    6>�97{m�@�>�2�Ie7|�b    =���C�T�7���B[9D�ęC�c�B�K�B�5A�<@��@@{&?��                                                ?s�A�^T@��?��M?K�>�4>&l�=OK}<M��                                                                    DI4F�>EFWD�oDs	Cm&B���Ā�@ɮ�                                                @��	C	9�B�+AQl�@ɔ�@F��?�Q�>�O=Ȅ�                                                                    C�E�d�ElYD8Ch�B���B0.A/ڒ@-.                                                @`�	B{�yB�K@�\�@B#D?�^b?ֈ>HVQ=D�X                                                                    7kd@7�]FAxC>0�cA�p#?(�?�0,�- +s��-�g-5%,�- +�wO.�-�.s�-�'�*��|7���4���2u�-���/�*�/�۽+J�7���4���1�%o2b2�-��x����2�o7���4�`    3���4�4�ɴ/�    :v*/844�8,�6�s�5��3Gq#                        5�6�6�*�6�ء4��3�L7�|2�c            4���>s�$-Fo?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  Co�G���>+��A�q�                                                B�                                          B��                    A�                @�p�    /�>�+���5UP&6�*k1=�L1
;`?Ɇd8K�@z�8�<�<�>                7��    2���6�5�        {@��            �TM�{@��{@��{@��7�и8�6            7��    7��{@��    7��    6��;    {@��    {@��                6��;{@��6���7��V8y�8y�<� %2��B��7��EӀBDԔA�h    C2f�            ;���4��n6,{~4��n<��2<�ŀ<��2<�H<Mr;��_:���:d��`����N�Җ~���r��L���dT�l�N�L�t�4���!M%�ʦ� �;�A~<��)<��	<�]}<�@J<��;���:��:Q��`����N�Җ~���r��L���dT�l�N�L�t�4���!M%�ʦ� �8��    *3�'6!�6�k�9:ڈ9`9,��9�#8��8�w                                                �k��qd�o_��mP��g�w�_���W��N7�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         �D5        {@��{@��{@��{@��                                                    CQ(�D�C�D���DS{�C�=�C�JBfb�A��@���                                                A�(�CRC�C@��B�{�B==�A��J@�b�@�?#��                                                                    E���G	`IG5|F�]7F?��E��1D�B�C��C�[                                                C� �EG��E;�E[D�_�C���C
B.�VA<�'                                                                    E�7F�S�F��>FzrF\��F.�VE�6�E���EH	�                                                D)8(D��ND��+D��SD�h�D~< D0��CߟuC�{�                                                                    7�|7]�wGR�SGI&�G��F��F}TFJE�N�EPy�                                                                                                                                @"TA6u8A� A�d�B�MB<��B^cB�5�<�<�<�<�<�<�<�<�<�<�<�<�E���E��"Ef�E�-D�pDX/=C���C���                                                {@��{@��{@��{@��{@��D%I�A��:/M�֧g&k8]�C��R@?�    A�;`A�;`{@�ξb?h�b?hC�Hn{@��C�`�C�`�C�Hn{@��@c��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G�;?.�D���G:D���D��'D��'A?�dA?�dF]��B5�B    C�2�FKފFKފD�(�D�(�F^fCB5�+                @�+C�)�C�JMC�s'?   C��uC�`�C�`�C���C�f!C��C�ƩC�k4C�C���C�:�C��OC���C�Q�C��C�� C�o~C�,[C��;C��C���C���C��,C��C�1C�vC�[C��g{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C�T�?��?	X?��?n�>�Ap>��>�\J>�1�>�s�>�"->�˫>Κh>ɲ�>���>�LH>�>�>�	�>���>���>���@$)z?�7         A���ǞbDǝ�/ǝT����6���{��ƒ����cs��c���cq�o�        ?K�B��    AJ�!AJ�!{@��@iVFO@.7e�7"m        ;�N�?�"?(�?wK?�?+�?ު?��??                                                �� m7ؑA�  ?k�A	��A z�$�      a�     U�9O� 8�JV8�7d56��5�[�4ѩn3��L                                                                    8�J�7���7-�6�5ȃ5 o�4k2�;S                                                                    G� ZGw�Fi+bE�!�E��DAXCPdBO`                                                4��|4*�3��\2�2&Y�1U50[��/PW,                                                                    4	�O3V�2�Z�2 t1R �0��!/���.��_                                                                    7[��7;-�6�6?�5�2=4ڮ�4034                                                                    7��7:�E6�CB6Cs&5�k
5!4�3{]                                                                    5��4��:4zL4�:3�.[2�_20,0��G                                                                    7��7c��6�5�6n�5�E5%E4?zu35zV                                                                    8���8���8��7X��6�v5��45
^�4�                                                                    5�W5��n5q4�KD4��3`˿2��1�]�                                                                    8���8���8��7X��6�v5��45
^�4�                                                                    8h"7杓7Y�6�GX6��5NP4m�X3_wU                                                                    5W4�5�4���4�3��2��2�50��                                                                    8h"7杓7Y�6�GX6��5NP4m�X3_wU                                                                    ���3C�Z5��5��[4ݕ�4"I	3f��2�jj                                                                    5.(5��78�37�6��6��5��t4�V6                                                                    5�)@5�64�w^4w53৸3��25#�1+�                                                                    7��7S�6r�5�(�5�|4;�E3O��2C��                                                                    �$#A����áϵE_c������Dn�۰���                                                                    ������|��ӕ��=p��_�3p�������                                                                    6�6ѵ5r�4���4
��3>&�2R��1F[                                                                    �,�y�=\����\�ٴ�X��A��i                                                                    ���_���\����aر��v�*���Y�A�lX                                                                    ��NҴk<'�!E=������*�"�T�5��:�z                                                                                                                                                                        ,	��,S�                                                                                            8�-�8mh,7��7V�6�p�5�e�5*�40�                                                                    5��4�C4~��3���3/�2x�}1��Z0��                                                                    7�R�7h��6���6U]5�(O4���4�3�6                                                                    4
�d3�$/3yc�2�[,25��1�x�0��/��m                                                                    8�wl8c�f7�C7S��6���5� 5{4)��                                                                    5�4�@4t43�Ȳ3;�2��G1�@�0�&                                                                    81�7�77��(7 )�6v5��4�F3��1                                                                    68�6Q�5�m�57 4��3��2˽�1�,�                                                                    6�C76�e6�H�6l�5F��4�[	3�I2���                                                                    4��4�٘4��	4��3cLU2���1���0�֞                                                                    8"��8�07�M7C�;6�hO5æ}4���3��;                                                                    6:a6!�[5ݿ35_�C4���3ߙ�2�{2 p�                                                                    4�pz4�E�4��4�3�4Y��4��3� �3Z~�                                                                    3�)�3�J�3���3eҬ3/�2�U)2�Sh20��                                                                    4��\4��'4�#�4�͸4��41�g3�D�3��\                                                                                                                                                                        $�      a�     U�62ʺ3���            F�YJ    =��6�0Y3��@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�j@    @�r     16:45:28        F�� @�r     @�y�    $�e      a�     [P5�M7�_)B'�fB'�f8��D"g�    )�#�+L�p76Tj?�  ?{kp    4�P(��/�}3���        /�Jq+�i>��K>'�=�<׬:�A�4��c-Ir32U��:O<!s�<�N�=F��=�ج=���=��O=�z�=�x=��5>��>/@ށy@���@90@	��?�g?ș�?��"?��p@k��@k�@ay=@aj�@a]�@aRU@aG�@a>X@a6�@a0�@a,@a)=��=���?���                    E��7m�G�UVG�Fh`^EȈ�E� D@�HCO�hBN�L                                                AO�C$�xB�^KA�AR,�@��?�ڷ>ٔ�=ز+                                                                    E@I2@�ZtF�A�7�0�VFZ�VB*��AIF-9��                                                    {@��B\QYB\QY@hҦ8~���8�    �hp�/<5��\�
?   )�#�A>Vn��.�1�"�    >�y�B�2��A�@6̑@6�        Bl��Bl��C��fC��fC��>��@��3�C�ف6��K6���A%�WC�>�y�C��Bw��C@�A�� BA�8CUA�o�B���A�ʿ׃	    A��    A��A��;y�!    6	�+AnF�A�A�d�@I ?���@�C:?��F�@�0�V.�	/+ǚ�1"pm8���6d��5�	XF��G.{�G��)G�'}=��            >YV>��b>�&�>���>��>��>���>���8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M҅�    Ag�Ag�{@��7�n9njq9S�8�#�8E�a7���66w5#��                                                A��hAr�~Ar�~{@��9�B@�BB8�@L*=>�<{@��{@��    C�#&BS��6���    6B��7|E	@�M2�z7f�t    =��C�x�7���BI�D�n�C�Bl$`A���Av�[@��?���>��                                                ?\єA��(@�8�?�>�?$aE>��B>2
=7�k<A��                                                                    D}$HF�rF�UD��}D��Cf8(B� HA��@�9�                                                @�CC�B�|dAE��@��j@<��?�S�>���=Ľ�                                                                    C��<E�6jEi�D�C`?�B��B��A+$�@*�                                                @Z��Bx�VB C�@��M@7S�?�P�?ؒ>A"=A?                                                                    7f��7��A5�>53oA�??�d?�x>,��+#}�-me,��0,��+Y��.O]�-�z�-�0*�D�7�ML4�/�2�A-��/��(/�	+t#6g
�4��1��2z12�-�pp�g
�1�"�6g 4��,    3���4�5�4���04b�    :��7�m�8�e6�b�59��37��                        5LU6��26x�4v8�3j9�7�!2�Q�            4��u>B�-
ž?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  CNȣG�e�>+{+A�+�                                                B�                                          B��                    A�                @�p�    /��+�q�5h�6zf.1���17?��H8!З@687֞9;�շ;�е                7���    1&'�6a�s    2D0�{@��            �N�{@��{@��{@��7f6d7��            6֫K    6֫K{@��    6֫K    6a��    {@��2D0�{@��2D0�2.�0��b    6a��{@��6��Y7iA7��b7��b<���1��GB���70F��C��A%A�    Cw�            <
�4��Z6~�O4��Z<�(<��`= |�<ݽ�<V�;��v;.��:����
Bʊ������������ꊏa��q��QB~�8Y�$����x���<�<�#�<���=�
<�Ѳ<R� ;�5(;"tE:�a��
Bʊ������������ꊏa��q��QB~�8Y�$����x���9A��    )�\�6�6u_A9���9��9��9j	�9A��8��8                                                �[���_<}�\���V�\�L~^�>l�.y��n�̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         9�        {@��{@��{@��{@��                                                    CQ�D�.�D���DS�C�!C��Bf)�A��\@�z�                                                A��CR.�C@��B��B=!A���@�)�@�\?#z�                                                                    E��(G	�GOF�[oF?�DE�}D�<�C��C�                                                C�-QEG��E<�EY�D�SC��/C��B.�WA<��                                                                    E讲F�YIF���FzF\�4F.�E�5�E���EH	f                                                D)9<D��<D���D��D�i�D~;�D0�3CߞMC�{3                                                                    7�!7\zZGR�GIB�G�cF�@F}z~FH�E�M�EPyE                                                                                                                                @!��A7�{A�$KA�|eB�B@��Be�B�n1<�<�<�<�<�<�<�<�<�<�<�<�E�E��Ef�iE��D�DX,�C���C��O                                                {@��{@��{@��{@��{@��D'�A�&�/A���g�8�C��W@UJ�    A�6A�6{@�ξ(T��(T�C���{@��C��	C��C���{@��@hҦ{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G2:��uD���GO�D���D���D���A:+A:+F^��B8�    C�W�FK�FK�D�,D�,F_EQB8�                @|�C���C�o.C��H?   C���C��	C��	C�[�C�;pC�C���C��tC��TC�HaC��C���C�z.C�,qC��2C���C�'�C��lC�|�C�;#C�C���C���C��C�C�mC�uC�Xe{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��5C�x�?X%?��?(�?��? >s>��y>��8>�=]>���>�^>��>�d >��k>��>�ּ>�>��>�J>�Q^>��m@��?ʊk        A��.�z���y�?�y���f��HF��b�����h�^�ҿ��ҿ��ҿ��a�        ?B$    A1��A1��{@��@N��FO�q7x�)7N         ;`{b?J?��?��?u�?��?'�?5�?"��                                                ��	 6釞A�  ?k�A	��A z�$�e      a�     [P9B� 8�/8*r7g��6��v5�b4�163��                                                                    8u�7�ϥ7+��6�C�5�*�5��4��3�                                                                    G�UVG�Fh`^EȈ�E� D@�HCO�hBN�L                                                4��4%�v3���2�i2,��1b�"0p�/i��                                                                    4 �>3Qt_2�H�2L1Z&40�2�/��!.���                                                                    7N�]77&�6�?�6CIl5���4�	q4��3+��                                                                    7~.�7-��6���6%�G5��4�5�4U3\e                                                                    5��4Ԕ$4S�w3�{[3|z2�fG2 �!1J�                                                                    7�U�7TT�6�$6J��5Ġ�5�4;3@T|                                                                    8��8�:d8�y7Z6�I�5���534;p                                                                    5�U�5�H�5t�4��4r�3e�2��e1�pM                                                                    8��8�:d8�y7Z6�I�5���534;p                                                                    7�K�7�Hn7[�$6�S6��5U�34}�3w2�                                                                    4�%F4��'4��&4Uf3�l�2��2/?1��                                                                    7�K�7�Hn7[�$6�S6��5U�34}�3w2�                                                                    �;o�i��5�;5��4�_�47�3��2��Z                                                                    4��T5'�7��$7\]G6���6I�5�=�56<G                                                                    5�_)5q(�4�24S
�3ú�3��21|�16n                                                                    7s�7QT6v��5��5-?4D?�3_>�2Y>�                                                                    �n*���� �H�ഝ�F��Wڲ𑚱�@                                                                    ��c��ň��Kr��j����3�?Zk���)��	                                                                    6�6��5vzP4��:4��3Fb�2a�C1[��                                                                    �"H��w��lն	e��e糴����v��p                                                                    ������"������8�5�y�m��P                                                                    �~�t�fLc� �6���g���-e̱FPu�Q��                                                                                                                                                                        +���,g�                                                                                            8��8hC�7���7Z]a6�H5�e�5Gy4c                                                                    5.�4��J4}w�3���37a2�qX1�Ln0��e                                                                    7�
�7c�26��H6X��5�<�5 @�4�3.�2                                                                    4x53�T3xM�2�/O2=^?1� �0��}/��)                                                                    8|��8^̢7�z�7W>'6�1!6Nb5*��4>%U                                                                    4���4�C^4s#]3ܕ�3C�2���1���0�b�                                                                    7��7�17�R�7#6�e5��54�733�QY                                                                    686��5���5:Kt4�Xs3��2��1�ʯ                                                                    6ʈ6�1 6�m6�;5N�4��3�� 2���                                                                    4�v�4�\�4�74��3l��2��\1�I%0���                                                                    8*8
��7���7G;\6��&5Оy4��3�*�                                                                    6/�6T�5�&�5c��4���3�k�3WQ2O                                                                    4��$4�V�4��s4���4bt!4GF3���3u 3                                                                    3��k3��3���3i�J36�*2��{2���2F�                                                                    4�.�4��_4˜�4��Z4�ci4=��3�j3�̑                                                                                                                                                                        $�e      a�     [P66(�3���                        6��3�~E@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�r     @�y�    16:45:30        F�
 @�y�    @؁@    $��      b     `�5��7�?B'�fB'�f7�[�D"g�    )w�+j�-7/J?�  ?g�    3ߍ(6u�/��	����        /H�$+Bc>��?>=4)=���=,�!;)P4�Y�-��52@+w:Pp<!U<�?h=F�D=��]=��=��f=�y�=��=��>�>T@�m�@�9�@�?��?���?�+O?��H?�i�@no@n�@c��@ca�@c4�@c@b�@b�3@bV�@b$9@a�"@a��?�`=:�R?���                    E��7e��G��KG��Fg��E��EE�D@>�CON�BN�R                                                A��C$p�B��wA��AQ�1@�0�?�p�>�0�=؂                                                                    E@�@ĂTFeA�Y�2zpF[N�B+��A<�-�1"                                                    {@��A��A��@m��8K� $f    ��o�/B���??   )w�A6m���1�8�    >�y�AN݈29H@��+������F        BP�VBP�VC�-?C�-?C�U>���@O�3�TC���6�\6��]A'�B�0e>�y�B�`�B	�`B���AQ�OA�;FB���A���Bh�%@�&��0.    @�&�    @�&�@��<��4st@70
|A�w@��GA�Z;@�?~,@T*?iH�F�)�2zp/�FJ,���27�z8#$6y5�1F�6G1��G��(Gȉp>�k�    6*4    >(��>��2>�%�>�~>���>�O�>�O�>�j�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M���    A��A��{@��7�f&9a%�9M�58�d�8Bz�7�~6��e6��5*�!                                                A���@�z�@�z�{@��:ۑB�Bnh5@Y��>�#�{@��{@��    C���BT�6��    6F77}9@�LR3x[7E��    >&UC��47�(&B?��D�a�C��BQDIA��AZ�+@�)
?�>��                                                ?O�PA�\@��=?z<?�>�C�=�!@=)v<8��                                                                    D{h�F��<F��D�X�D	�bC`��B�f[A���@���                                                @���C�B��A=`�@��[@4��?��1>��_=�                                                                    C�.yE��Eg��D/CZkOB�ӘB	��A'"�@(\�                                                @V��Bw�A��?@���@/b�?��4?9>:��==H�                                                                    7^#�7ۮMAǲ>:'WA�J1??�)+�t*-},�
,R��+�t*f��-}-�,֔�)�o7���4p�x2	��-3�/��c/�v�*�#�;4��1��2�2�!--�D7#�;1�8�#r4��9    2=�4h�4g��/��    :{�/7��m7�485��o4Z��2&H                        4�35�|5�3�3]5<2g6�W2��t            4�2=�y5,�cA?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  C4G�K(>)�PA�r                                                B�                                          B��                    A�                @�p�    .�n*�:4��c6o �2x!72;�;?;��7��?���7V!S;�Y�;��                 7?C�    2��5�2W    3t��{@��    -�     �ݦ{@��{@��{@��8Ҭ8O�            7Q��    7Q��{@��    7Q��    5�dN    {@��3t��{@��3t�x3C��2DzK    5�dN{@��6��h6�v<8R��8R��<�x�0��B���6׷mF$��C�z�@��    B�?�            =5��4y�U7ç�4y�U<�Zz<���>F1>~�=���=^:<��H<
�w�v2��׊�����?��(�����t׊S��9��&*Ɋ0��	�=2L�<�W*<��A>E|n>i=�-e=#<�]M;�6��v2��׊�����?��(�����t׊S��9��&*Ɋ0��	�:WP?    *�[5��/63��:�,9�x�9�s9�]N:㟪:��                                                �Ȑ�����c����������p�ю����̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         :�{r        {@��{@��{@��{@��                /|E2���        9
/|E            CP� D� �D���DRK�C�FnCFBd�!A��K@�k                                                AВ CR �C@��B�K�B<FnA�F@��!@�K?#k                                                                    E���G	��Gc�F�Y�F?��E�sBD�5lC��C�                                                C�5�EH!PE<3�EXdD�G�Cݾ�C�NB.�VA<�                                                                    E��F�]+F��fFz�F\�GF.��E�4�E��
EH�                                                D):D���D��D��{D�jbD~;xD0��Cߝ%C�z�                                                                    6�W7X��GRքGIY�G}-F���F}r�FFE�K�EPx^                                                7�b�                                                                            @R�\A]��A�c8A��B(��BT��B|vB�#�<�<�<�<�<�<�<�<�<�<�<�<�E� �E�GEf��E�=D���DX(�C��'C���                                                {@��{@��{@��{@��{@��D!�QBw(.��d�gM�8�C�mx@\�    AE,sAE,s{@�ξ��N���NC��>{@��C��HC��qC��>{@��@m��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G-�:DmD��
GOaD���D���D���A6(�A6(�F_#`B9{K7�b�C��PFK�FK�D�,�D�,�F_�1B9v4                @*	=C�H�C�*OC���?   C���C��HC��HC��>C��C�kC�"�C�.!C�5�C�7�C�2DC�(sC�uC��DC��oC���C�u.C�0C���C���C�_�C�0AC��C���C��C�RC��C��A{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C���C��>�K=>��>��>��>�ϕ>ݗ >��>�7H>�ԣ>�!�>��>�,�>�׈>��!>��n>�]�>ǟ]>� �>�O�>�u�@=e/@|{        A~���n��'d��F�Q�ǅ7�Ǆ��ǄwM��R��đ�đ,�đ��~        >w��A�ڡ    @Ղ@Ղ{@��@��VFO�7oRt7z~        ;���?5ې?4a?5�?4�?43�?4�[?6;b?:�i                                                �� 6��]A�  ?k�A	��A �$��      b     `�97I8�V�8�L7g�6�P�5�z�4���3�'                                                                    8g��7��7*86���5���5	��4��3��                                                                    G��KG��Fg��E��EE�D@>�CON�BN�R                                                4� �4 �Y3�%�2�.�2->X1d��0u	
/tҩ                                                                    3��3K
 2�J�2�81Z�|0��"/��W.���                                                                    7C#u72)?6��6B�5���4��4-�33��                                                                    7jB�7"W6���6��5���4�4��3�                                                                    5Қ4�:-48 �3�N3^b�2��Q1�1J                                                                    7�(�7F�6��f63LR5�y�5	ҋ40�V3@��                                                                    8���8�z8
gY7U`�6�-5��&5K	4s�                                                                    5{"5��5ț4�ǟ4�3_Ң2�c{1��                                                                    8���8�z8
gY7U`�6�-5��&5K	4s�                                                                    7�=�7�'#7Tf$6�j062C5RO-4}*3�                                                                    4�4�a>4�X+443�2ڂV2��1�F                                                                    7�=�7�'#7Tf$6�j062C5RO-4}*3�                                                                    �FX����5��S5�h�5^4@%�3�ɰ2�z
                                                                    3���4t�8֢�8�^;8)&97��~7"ɥ6��                                                                    5�Y^5a��4�P�4;c�3�P~3�2(Ho17��                                                                    7�7�|6o܋5��5
Fr4B�}3`dT2a�	                                                                    ��3�繵�zr�H��������0���l�\                                                                    ��KA�ն��˚X��xز�KY�A�i���u��s                                                                    6�C6��5oE�4�iN4J�3D�\2b�1d$m                                                                    �?��*ƶ��p��D�f-Q���_��v���v@                                                                    ��񑳌���<M��2}��6��q̑��t�                                                                    �p�f�`��~������&�/l��J���[�                                                                                                                                                                        +�m�,)�                                                                                            8x�*8a��7��7Yj�6�{�5���5�^4&ɓ                                                                    4��w4�ln4zc13��378�2�M21��|0��d                                                                    7s�l7]RZ6�L6W��5�rG5�4��36�                                                                    3�i`3���3uI2�8X2=��1��0��/�(;                                                                    8n��8X��7�|7VN�6�h�6-T5-�g4F�u                                                                    4�IH4�a4p/ 3۠�3C�b2��31�Ν0�c                                                                    7���7ܑA7� �7# �6��5���4Ǐ�3�o�                                                                    6_I5��5�n�5:I�4��3�Z+2��1�Z�                                                                    6�oz6�<}6}��6�5P�4���3�B�2���                                                                    4�ȋ4˲�4���4�73m��2�y�1�Lu0��                                                                    8��8ʌ7��u7G9�6�R�5��4��b4D5                                                                    6%tv65�N=5c��4��T3�5Q3`82)a                                                                    4�#�4�lO4�8 4��_4c��43�/�3���                                                                    3�$�3�#�3���3j�37�2���2�W�2O�r                                                                    4؁14��(4��4�X4�$�4@#3��3�&�                                                                                                                                                                        $��      b     `�69�J2=�                        6�=�3�*�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @�y�    @؁@    16:45:33        F�H @؁@    @؉     $�-      b$     f�5�l�7a�|B'�fB'�f7�:TD"g�    (]�^+}o�7�?�  ?�z2g��3Ր5)��/�nj����        .���*� n>��>K� =�~�=Mq�;ν�5L��-��27��:8B<!�<�;s=Fї=��u=���=���=�x�=�f=���>�>}@�~p@u�Z@��?�f"?��X?i?;��>ߣ5@m�T@nC�@c��@c�C@c�&@c�I@c~�@c]�@c<7@c.@c�@b��@���=
��{@��                    E�&�7j�=G�b�G!aFg
E�wE��D?��CN��BN6%                                                Ab�C#�"B���A�+PAQ7@��?��
>ؕd=�G                                                                    E@�@ĤyF&�A�v�3�\QF[��B,-,A2~T-��G                                                    {@��A#��A#��@oYT8<��%�2�    �L|/Dn��>�
?(]�^@�	��Q�2��Q    >���@n;2U?�v��6���KS1        BB;?BB;?C��UC��UC���>��?�^63*�<C�o6̈́[6���A)/B$g�>���B:�Aj1�B���AB/A&ÈB<�(AN��B:i~�m[������N��m[    �m[�lW�>��;:��;D7@���@"@EA(ƴ?��`?>Qv@Ti�?W�F��3�\Q1a.7Z3��[7�:�59�.5Ѽ9F�XOG+tG�9"G���>�q0    <�    >I3C>ʊ�>��>�/�>�ܹ>�i\>�/�>��8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M��]    @��@��{@��7��o9a��9Q�-8ܧV8H6n7��17 �6)��5S                                                 A��@@�K�@�K�{@��=�	AA�M}BK'@d�g>&0}{@��{@��    C�LaBV��6��    6p��7�@��63S�N7=?W    >���C��7�`TB8*}D��C���BBOpAΡ0AI>z@��Z?�;Z>嬬                                                ?G�AA�O@���?kx?Is>�G�=뵙=�<,U                                                                     Dy�BF��JFu�D�HD�C\�B��3A���@�d�                                                @〝C��B��4A6t�@�F�@.x?�c]>�
�=���                                                                    C��RE��1Ee��D~}CUP�B�sBE�A"�}@$�                                                @S5iBv	1A�E�@��O@(��?���?	��>4K=7Ԗ                                                                    7a��7�#�A��>>A�؁?T�?��C+�)�s�,��,��+�)���,0��,���,�y�)�7�9�3���1Vp�,\+ .�͏.�?)�<ҷ���3�q�1��1O9�1O|,T�7���1�������4�<�    ,d�3�E�3��.�-�    ;_B5�
16���2��81���.��                        .�O2�C2��Y0C    2І92�2            3�	8;r��,�B�?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B���G��9>+��A���                                                B�                                          B��                    A�                @�p�    .9�f**v}3�Jj6p�3�+3�_�>�37x,�?_#l6���;�2;���70�            6�J�    4Q*�5C<�    59�r{@��-��1�.    ��,�{@��{@��{@��7���8            7�\    7���{@��    7���    5��    {@��5|�;{@��5bĆ5�i4��:    5��{@��6�k5ǥ(8�
8�
<�3�    B�~6���FB~HC�a�@	y�=�g`B�            >�i�3���9C�s3���<���<��?�??��?B>�<�=���=b~��f�����؁��)A��kR���V�s�R�+�9�&�&7����>�L�<���<��C?˫l?��a?�>��G=�0=+Tl�f�����؁��)A��kR���V�s�R�+�9�&�&7����;�8�/�)�9�5�J<6�l;e�{;bd:��:��<6�<<\��                                                �c���f���Z�c�MԻ�9�^� !�⭌�l�]̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾�         <쪵        @��@�P�@�C_@� </�            50E�:1p�        <t�50E�            CO��D�)D�ƘDPտC��C�1Bc	IA�{�@�J�                                                Aϻ�CR)C@ƘB�տB:�A��1@�	I@{�?"J�                                                                    E���G	��Gz/F�T]F?m�E�d?D�'C���Cڰ                                                C�;9EHDDE<T�ET�D�8�CݩC��B.�'A<�                                                                    E谰F�aJF��WFz!3F\�|F.ȏE�3QE���EH                                                D):�D��D���D��D�kCD~:�D0߁Cߛ@C�y�                                                                    2І97`	�GR�GIr�Gu�F��F}f3FA�E�HEPu�                                                =#�                                                                            @x��A}-�A�ޤB�B=>�BtC�B��rB��<�<�<�<�<�<�<�<�<�<�<�<�E�2�E�&vEf�E�iD��JDX!C��\C���                                                {@��{@��{@��{@��{@��D��B��-�0F�f38:�Cv��@Pf�    A-�9A-�9{@�ξ�0U��0UC��{@��C���C���C��{@��@o]�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��:
qD���G)D��aD���D���A3	�A3	�F^�'B98�=#�C��,FK�FK�D�+D�+F_{aB93�                ?�O�C��|C�DC���?C�~ZC���C���C���C�� C�bC�!�C�;�C�TdC�iC�u�C�|�C�KC�|�C�t=C�d�C�LC�(�C��pC�ɌC��"C�d�C�+-C�NC��C�'C��C��y{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C�|�C��>��>�h,>��>�Dc>˽>�Q>��{>�%>��~>�g&>ђ>�d9>��>�Г>�/#>��=>��>�չ>�>��Z@=$�@�c        A�5����ݹ��ݖ�ı���R��;�D�)�ƾ/�ƾ/�ƾ/��tp�        =��oA�F�    @�C4@�C4{@��@�^�FO��7LQ�7�9        ;�?HX?HW�?JT?L��?Pxt?WM�?cm�?x��                                                �19�6���A�  ?=P@��)@�2n$�-      b$     f�97hW8��#8	͏7p��6�E5�x5�d4�                                                                    8g�77ĵ97. 6���5ި5=�4+�~38H                                                                    G�b�G!aFg
E�wE��D?��CN��BN6%                                                4�A�4#>03�p�2�X28�1|��0�K�/��(                                                                    3��~3N3�2�s�2:�1iU�0���/��".�	v                                                                    7C�/75��6�P�6J,�5�?�5"V4,��3`Z.                                                                    7eǄ7�6�?o6-�5��e4�	4�C39a�                                                                    4��4��U41��3��3[}l2��22nX1�E                                                                    7�k�7B9�6���6-��5��|5̫4@��3b��                                                                    8�y8�+[8
ً7Z��6�&�6�T5(��47�|                                                                    5z�45���55�4��4BI3n�2�i1�7�                                                                    8�y8�+[8
ً7Z��6�&�6�T5(��47�|                                                                    7��(7��7TV�6��|6.�5b�$4��3���                                                                    4�:�4�.�4��q4ˮ3��2���2�1.�                                                                    7��(7��7TV�6��|6.�5b�$4��3���                                                                    �.������5�
l5�50@4X�3���2�Q�                                                                    0��B1)&�:f��:,I*9�I+9�)8���8 ��                                                                    5��5\�`4���44��3��3
27�.1X6                                                                    7nf7v�6q%�5���5V14S��3�d�2�6�                                                                    ����
���ƒz�Q�����8��^�+�&�U                                                                    ��"=������S貖+�V3�����Ba                                                                    6j;6�|5p�4��:4�3UgB2���1��l                                                                    �����涛�?�ǋ�tC���߳�vw�6�                                                                    ��o볏b`�;��Le�}�H����믬k�                                                                    �q���dp��#������ô�By�l}���X                                                                    (�&)KS�.� �.�h�-��B-G�/��/�~�                                                                    +|q3+��                                                                                            8y�	8f9�7�u�7aS�6���6��5('
4O�f                                                                    4�V�4��|4@�3��3Bo�2�>`1�y@0�E�                                                                    7t~)7a��6�Q6_�}5�N�5��48`�3c��                                                                    3�2P3�H3z12�D�2I2$1�}0��H0 ��                                                                    8ohI8\��7�,#7^"6¢D6M5H�]4wȌ                                                                    4�
4�I4tٸ3�^3O��2���1�;Q1Vf                                                                    7��7���7�"�7*�6���5��\4�\K4c2                                                                    6��6 d5�'�5Bf�4�<*3�Hi3�2��                                                                    6�G)6���6�66	t�5^ZA4�Wf3��2���                                                                    4ۿ
4ϬD4��4�3~2�c�1֖�0��                                                                    8i8	k�7���7O��6�'t5�p�5��4%y=                                                                    6&.�6z5�X5m�4�-4e#3"H�2=�                                                                    4��W4��^4��4�3M4st�4-��3�H3���                                                                    3���3��s3��3t]A3D��3u�2��~2��                                                                    4�Vj4��4�I_4���4��f4Tr.4V�3Ģ\                                                                                                                                                                        $�-      b$     f�6d
f,d�F�P�    =%�T            6�94r$@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @؁@    @؉     16:45:35        FĄ @؉     @ؐ�    $��      bB     l`5��a7-�PB'�fB'�f7e��D"g�    'b�=+�y�6��?�  ?�2���3���)�m>1&�����-���    -�@�)��>G>�>K��>�=��=[�8�e1Tj�2/
d:�v<!G<�C�=Fи=��o=��=��Q=�x~=�\=��>�K>�@�|d@l�S@
��?��?�z?X�s? �><��@i�y@m�$@c�@@c�C@c��@c�@c�	@cu�@cV�@c8D@cb@cA+�<��f{@��                    E��P7<F�G�
�G�qFfxUE���E��D?12CN	BM�                                                AC#��B�'�A�AP��@�bm?�V>��i=�gq                                                                    E@��@ĵ�F8�A���4ޛF[��B,@�A�[-���                                                    {@��?�o?�o@nɋ8��:!�	�&��/S��#���?�'b�=�:���u�2�8�    >�`b?�2_�?����.ծ��        B8��B8��C�	AC�	AC��/>��}?�l3�C��6�Н6�%A)ovA���>�`bA�'�@���A�Q@��4@R�A���@�A�B�,��7x@__�"+ ��7x    ��7x��#i>��=E�=F7@:��?rϞ@��?G�X>��e@.�5?9?�F���4ޛ1���.��r4��6�%82.�5�.�F���Gh�GEnKGZ��>օ�    >�=�>d!�>�'M>ؗO>��K>ڷ�>�A>�%�>��k>��8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�{Q    >�~�>�~�{@��7���9#�}9-H�8�~8)�7�D6�*�6X5>�                                                A��ٿ�1��1{@��@j�TA"��B4J�@j�B=��{@��{@��    C���BX�e6U'    6T.6�n�@�o�3!�7�    >���C���7�,�B4'D�A�C��AB>��A�hiAA��@��?�?�>�ub                                                ?C�NA��3@�?i�?�	>�4=��=7<!P�                                                                    Dx��F�ѪF�)D���D!�CX��B��A���@���                                                @�?Ca\B���A2 @��i@)�?��>��o=���                                                                    C�}E��XEc��D�CQΉB��wBi&A�@!1�                                                @PԙBu�BA�ڻ@�mQ@$z�?�ԕ?�>.q=1�X                                                                    77��7�� A7>?U/A�?p%?���*:lA(�o,7��+�S#*:lA)�X*��,��, 1�)47e�^0�
.Y��(��            ��H�3�#1�B�.Yܷ.Y��(��7�H�1tg���Hw4��&    &��q0�>$0�*�+�C    <�A��U�3�ʑ.���-�	�*�Y                        )���.�/p.ɚ@+�L    .��2�J�            2��f:��,�
?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  B-��G���>,��A�'Y                                                B�                                          B��                    A�                @�p�    -= )+|�0��=6<24�4��>�36֫�>�E�5s̬;wy;y��7C��            5G�    4DȔ4�1,    5�c�{@��1͆�3�5p    �0f&{@��{@��{@��74��7���            6�sK    7�H�{@��    7�H�    4��#-Nn�{@��5��
{@��6	��5��J5�H    5�#{@�β�w#4s\�7�?�7�?�<��4    B�٤6��FU�JC��i?�?�>x��A|�            ?D/I0�9�u�0�<}�r<��@u��@5��?��?��>`=���f������~��)>��kP���T�s�R�(�9�#�&5����?@��<}��<�"@t��@5�?�AT?7~>@Q?=�r+�f������~��)>��kP���T�s�R�(�9�#�&5����<dT�/ݐ�(��F5g�5�OC<�,�<"N�;zT�;0]n<��<��                                                ���Al��S�W�@d��&e�7�î�_�/L̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 2�|b2�K)=�tK4�դ4��@k�@8@?�i�@F1S=��v4"c94z�B�Y�5��`;��>��2�oc=�W�5��`F�9�1���4�a�CN�&D��=D�£DOL	C�y�C��B`�rA�l@�2�                                                A��&CQ�=C@£B�L	B9y�A���@��r@l?!2�                                                                    E���G	��G�TF�NF?Z7E�UAD�kC��.C��                                                C�=EH^E<oEPD�*WCݓHC��B.zPA<�                                                                    E�:F�dRF��Fz$�F\��F.�ME�1�E��REH                                                D);D�ID��_D��+D�l
D~:�D0ކCߙ1C�x                                                                    .��79|GS�GI�*GmiF��IF}YHF<�E�C�EPr�                                                @��?nzB                                                                        @_�mA~|"A�}�B��B@mBz/BB���B��D<�<�<�<�<�<�<�<�<�<�<�<�E�?�E�4�Ef��E�D��|DX�C���C��J                                                {@��{@��{@��{@��{@��D
�}B!��,󳨧f��7� �CG�6@4�K    @���@���{@�ξ_�ξ_��C��{@��C��C�%
C��{@��@n�c{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��G��9�D��G�D���D��SD��SA1A1F]�TB7�u@Kt�C�[FK��FK��D�(D�(F^��B7�]                ?E�C�asC�8�C���?0�C���C�.�C�.�C�v-C���C��C�+HC�l+C���C��C�)0C�XsC���C��C�˴C��C��<C��(C��C�ȼC��-C��9C�JDC��C�C��C��C���{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C��6C�m�>���>�&%>�[�>�tm>��*>�@	>�->��>� �>���>��d>��>�$>ǂ�>�Z�>�W>�pV>�Ȳ>á_>�J�@�?�j�        @�}d�4������q�=�׳��י���~�Ğ*ƽ�1ƽ�7ƽ�1�af        =oaZA�#�    �v)��v)�{@��@u�FPq7��7��        ;i8�?0P�?H��?K�F?O2�?T}w?]K�?l�?&l                                                �jt#6p?A�  ?%yV?�bN@�v7$��      bB     l`9z�8�n�7���7L��6��e5�~�4���4�d                                                                    8'W�7�:�7�|6�xV5�$J5��4��3'~c                                                                    G�
�G�qFfxUE���E��D?12CN	BM�                                                4���4��3r�2���2 <�1_��0�A�/��X                                                                    3�js3*_2��v2��1Jg�0�I�/�E�.�s�                                                                    7�-7#6�d6+��5���4��-4%w3K��                                                                    7$�M7 \6Y��5��V5k��4Ņ�4'!3�K                                                                    4���4��4
R3�Q38�2��(1�N1�                                                                    7I�B7�X6�-�6��5��{4�jo4&h�3A�                                                                    8f�8|`�7�B78A�6�l5��K5θ4$Ml                                                                    56�5\��4�4m�3��3M�c2�f1�m�                                                                    8f�8|`�7�B78A�6�l5��K5θ4$Ml                                                                    7�
+7���7/&�6�1Y5���5E%�4��3��U                                                                    4�X�4�A�4ar"3�Y3j��2�z�2Z13J                                                                    7�
+7���7/&�6�1Y5���5E%�4��3��U                                                                    �$'賤�65�#�5���4��R4BG3�hB2�xG                                                                    2&h�-;�i;:�b�:?�X9��]8��28g�                                                                    5n�522B4��4�3�~�2�&2�19B�                                                                    6�T'6���6G��5��4�~�49��3f�2y��                                                                    ��_���@����6�2����a������e-��                                                                    ������񲯕ﲟ,���r²>8� �Ʊ���                                                                    5;k5�G�5Fg74�D�3�S.3:�v2g��1|
�                                                                    ��sj��Ѷ��z��k��R�R���#��Kĳ	϶                                                                    �`�*�m����vcT����0��{�寜�                                                                    �/@�<�j�	Y�������,A=�V��z�	                                                                    ) ~|)h/�/��Z/�.�:h.g00v��01�                                                                    *���+��                                                                                            84�~8>X7҈Z7?U6���5�m~5�4<                                                                    4���4�r�4T�\3��;3'��2��R1���0�c                                                                    70�7:w96�>&6=��5�f�4�0�4&Lu3N14                                                                    3��_3��3P�2�c�2-�}1�
0�_�/�                                                                    8-4�86�U7���7<Q6�۷6y�54��4`VN                                                                    4�4��;4L/�3���33Y2�i�1��Z0��                                                                    7�U�7�x7�Y\7@�6o1#5��84Ӡ3��y                                                                    5��5ԛ�5��E5& �4��]3��e2��p2�E                                                                    6�B�6�Tl6Z��5��i5AIS4�Z3��2�_�                                                                    4�'�4��24y��4$�3\�2��f1�p�0���                                                                    7ҡ�7�_�7�m871��6�,@5�>~5S�4Ʈ                                                                    5�6�5�5J�x4� 3�ٴ3�62,P�                                                                    4���4�X4�{�4�*�4S�4[�3�6�3��'                                                                    3P�3l(Q3fF�3P��3+g2�w�2��2l�~                                                                    4�X�4��4�%~4�޽4�a4<��3���3�5�                                                                                                                                                                        $��      bB     l`6GC�&��q                        6���3�=@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @؉     @ؐ�    16:45:38        F�� @ؐ�    @ؘ@    $�U      ba     r04�\X7��B'�fB'�f7={�D"g�    "x�+r�?6��?�  ?��5-C�2Q�2��20��X0l��        (O�$E�9f}8��=_@�=k��=�AP<���6��2-x8: �<!�}<�Z#=F�i=�ј=���=���=�x�=��=�ա>��>">h=�(>,��>
�d>W>�a>�<�>���@g�@i��@_��@`��@aKD@a�H@bd�@b�b@c8@c"9@c+�@c+�{@��{@��{@��                    E��.6D��G�OG�MFfvrE��}E]�D?CM��BMZ�                                                A;C#ٲB�F)A��APn{@�>?�#�>וK=�                                                                     E@� @ı�F6�A���3+�CF[�iB+>��-�P                                                    {@���>���>��@e7�$N�C�:$wo�Z�2/VP���f?#�"x�����/`��4��    >�y�<�s1�N@�[�h�H�ID=�[    B.�B.�C�FXC�FXCz�-?Z`�>o��3b�%C�pz6p�6��BA'U@҆�>�y�A�?ۗ:A8��@�?�� @�H@E#�A�����.V@��� Y@���.V    ��.V��ّ@3�d?9?9
�?�u�>��:@�&>�f�>�tH@p?-%�F�r�3+�C0�O�-���3Je�6<��    4�JFf�F�f	F�\F�?�LB@7�0A���@^44?u�9?La�?
I=>�S0>�>7>�R>�a�>�Y�8��8ě�8��>8�IR8���8iS8D��8*
8��8�[7���7���M�׽    ��'|��'|{@��6�7�Y�7��7��37��n7T�O6��5�)5�                                                A��������{@��B�%�@��B+)@c,�<�|�{@��{@��    C��KBX�6     6�7E�@�K�3���6��!    >���C��97��BAa�D���C��}Bl��A���Ab��@�y?楧>���                                                ?U�bA��@׎�?�g�?!�p>�`�>ω=,6�<;%^                                                                    Dy�mF��ZF+dD���D� CZ�gB�$A���@��P                                                @��C��B�#�A8@���@,��?�A�>��)=���                                                                    C��E��lEe9ED;�CU�B��AB�A!T@$H�                                                @TNBxmA�X@�y@)r�?�q�?Ȑ>1=6E;                                                                    6]��6�5A��>=��A���?�?���)���(\��+��+=G�)���(�UA)h�<+�7;+���(�j7={�                        ���3b�%1��            7����oJ���4ۉ                        =���W+                                                                3��            3��{@��,��>�'�>��>��]?2�?`��?D}q?v��?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  ?�  A�tiG���>+A�                                                B�                                          B��                    A�                @�p�    (�#�hi    5��37�T3Jl�=`�6.�=�u�4g��;8�F;>G�                ��F    4���3
Y�    7+\{{@��4X*5�w�64���e�`{@��{@��{@��4��7Ϫ            7��    7��{@��    7��    5
�6N�{@��7cY{@��7vE6b��7=՗    5dG{@�ζx�2C�J7SS!7SS!<���    Bì�6�<�Fw��DS�>�
�=��4@�6(            ?��W    :��    =��=ރ@���@� �?�'?/�j>��>�Ӌf������~��)>��kP���T�s�R�(�9�#�&5����?�d=��=�t@��=@��?�vH?-g�>bō=��\�f������~��)>��kP���T�s�R�(�9�#�&5����=�o    '��3�j3p�.<�N�;��:v�<7u=��O={̔                                                �.��T&��N�@�ēE��{�'�H�-���̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� ̾� 5	��49C?O�16�@�5���?���?؆*?��3@<'�=X6u6G�N5I�7E�x�7��h=�$�A�j%>�il=��r7��hIr6��7��=CN:�D���D��ADNc�C��dC��B_i�A���@�e�                                                A�:�CQ��C@�AB�c�B8�dA���@�i�@
��? e�                                                                    E���G	��G� F�L&F?S�E�NUD��C��C˧                                                C�5@EHQE<j]EN�D�%�C݉6CۂB.oHA<�!                                                                    E��F�b�F���Fz%?F\�F.�aE�1�E���EH                                                D):�D��D��qD��D�lpD~:�D0�PCߘC�wV                                                                        6}L2GR�bGI��Gh�F���F}RcF9�E�A EPq                                                A�ՍA��Ak��@��O                                                                ?+�@音A��0A�(�B0zSBc��B��aB��<�<�<�<�<�<�<�<�<�<�<�<�E�7�E�1eEf�E�qD��DXC�{�C���                                                {@��{@��{@��{@��{@��B�*B���'���gO�77�WC��@�X    @F;�@F;�{@�ξX�9�X�9C�Y4{@��C��CC�qQC�Y4{@��@e_{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��Gvv> �?D���GŮD���D��ID��IA4� A4� F\�UB5�Bk`]C���FKۣFKۣD�#�D�#�F^5B4��                >�P]C���C��C�J�?�hC�*C�4C�4C�Q�C��)C��-C��C�'�C�lcC��xC��iC�4�C�q>C���C��C� C�N[C�vIC���C���C��>C��C�Z�C��C��C��C��C�_�{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��{@��C� VC�Ȼ>��>�a$>���>�-�>���>���>�Y�>�I7>�ߠ>�
>�nN>��m>�>�J:>�Ik>��x>�?H>�ˀ>�fs>�G�@!U�?�b        @,dƾd_ƾdƾc��z]�ƿ�^ƿךƿ�ā%wƾƾƾ�y?4        =L��A�aU    ��R���R�{@��@v��FP
�6�V!7�#        8D�>��>w>��!?8Ip??�)?E��?N=-?[�i                                                p5^�A�  ?k�>�7<#�
$�U      ba     r07��w6�WC7A�6���6Z��5�R�4���3�zh                                                                    6��6z{64��6W5��4�h�3��r3hA                                                                    G�OG�MFfvrE��}E]�D?CM��BMZ�                                                3+�2quZ2�(�2V��1� �1�]0D�}/XJA                                                                    2,�1��1���1��21�00v�/xE�.���                                                                    5��c5��E5ېs5�-/5Tw4��{3��3#d                                                                    5��25j��5��95�I5C�4���3��3ڹ                                                                    39�3�c3C��3C�3V2e��1�z�0���                                                                    5���5�^m5�O�5���5n{4�]4t�3.��                                                                    6惢6�rh7(�6��6_�m5���4���4;                                                                    3��63ɼ�4��3�*3�`(3�'2\��1z�:                                                                    6惢6�rh7(�6��6_�m5���4���4;                                                                    6*a76,4�6X�6	M5���4�b�4D��3]�                                                                    3&ؘ3:y63�Ѩ3f��3+%2~��1�i10��j                                                                    6*a76,4�6X�6	M5���4�b�4D��3]�                                                                    �5/Ų��4�Dl4��b4���3�;3X�2�ɛ                                                                    8)?�7�S;=�;�R:���9��]9;��8�7`                                                                    3�r3��S3�{�3��3l�2�ť2G[1'7�                                                                    5O�;5M:�5w(S5��4��c3��'30� 2EPg                                                                    �P=ݴM���Դ����R����.�ƎO�� �                                                                    ���!�-����'�:�߱��ı�V����-                                                                    4N4K�E4u��4�o3���2�\~22F�1Go�                                                                    �Zz�Ra���W�oBA�ϲ�JG��������f                                                                    ��zұԵ�"���2z��Ű�<��>L��r1)                                                                    ���N�����*�D�����ŕ��T��#�E�C_�                                                                                                                                                                        )���)���                                                                                            6��6��e7&6��6_8�5�t�4�=:4��                                                                    32��3,��3���3B��2�72��1�է0�=                                                                    5�f|5�I�5��c5��5f��4�3�[�3�R                                                                    2/Q`2)#52��"2Asi1��1'��0�\�/�.�                                                                    6��*6���6���6�h6n�s5��5�94.�                                                                    3+��3%��3|su3@�2��2/�n1���0� ~                                                                    6(��6&�6�k6�5}6+~5S��4��3�R~                                                                    4A�4>�34�4���4C~G3r?�2�$1���                                                                    5�I5�5��+5j�]5
:B4+Ik3�3Z2�is                                                                    3/3.�3�1U3�~3�p2C��1���0��                                                                    6N��6L:6�cJ6�zC6Q�5��4���3�~                                                                    4l
+4i3�4꺞4���4n�3�
m2��2Q�                                                                    2�u833��^4&�4n�3���3�R�3d��                                                                    1��1��2��2к�2���2���2��28��                                                                    3G�3 A 3�<4��49K3��3º�3���                                                                                                                                                                        $�U      ba     r06���    F�%     =�    >���    6�.�446�@� @�Mk@_;�@��?�T�?�=>VT�=w��                                                                    A��@�b@>�E?���>�9=M��                                                                            A��@�b@>�E?���>�9=M��                                                                            12/10/21        @ؐ�    @ؘ@    16:45:40        