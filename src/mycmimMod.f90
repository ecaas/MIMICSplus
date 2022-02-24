!This module models the decomposition of carbon and nitrogen through the soil. The number of soil depth levels and different pools are given in the file paramMod.f90. The paramMod module contains
!all relevant parameters for the flux and balance equations used in this module. In the current setup, 7 depth levels are used, and the following 10 carbon and nitrogen pools/reservoirs are found in each level:
!LITm - metabolic litter
!LITs - Structural litter
!SAPb - bacteria saprotrophs
!SAPf - fungi saprotrophs
!EcM  - Ectomycorrhiza
!ErM  - Ericoid mycorrhiza
!AM   - Arbuscular mycorrhiza
!SOMp - Physically protected soil organic matter
!SOMa - Available soil organic matter
!SOMc - Chemically protected soil organic matter
!In addition a reservoir of inorganic nitrogen, N_IN, is found in each layer. A plant pool of carbon and nitrogen is also included (not vertically resolved).

module mycmim
  use paramMod
  use dispmodule !External module to pretty print matrices (mainly for testing purposes)
  use fluxMod
  use initMod
  use writeMod
  use testMod
  use readMod
  implicit none


  contains
    subroutine decomp(nsteps, run_name,nlevdecomp,step_frac,write_hour,pool_C_start, &
      pool_N_start,pool_C_final,pool_N_final,start_year,stop_year,clm_input_path,clm_surf_path) !Calculates the balance equations dC/dt and dN/dt for each pool at each time step based on the fluxes calculated in the same time step. Then update the pool sizes before moving on
      !INPUT
      integer,intent(in)                        :: nsteps               ! number of time steps to iterate over
      character (len=*) ,intent(in)             :: run_name             ! used for naming outputfiles
      integer,intent(in)                        :: step_frac            ! determines the size of the time step
      integer,intent(in)                        :: nlevdecomp           ! number of vertical layers
      integer,intent(in)                        :: write_hour           !How often output is written to file
      real(r8),intent(in)                      :: pool_C_start(nlevdecomp,pool_types)     ! For store and output final C pool sizes 
      real(r8),intent(in)                      :: pool_N_start(nlevdecomp,pool_types_N)   ! For storing and output final N pool sizes [gN/m3] 
      integer, intent(in)                       :: start_year !Forcing start year
      integer, intent(in)                       :: stop_year !Forcing end year, i.e. forcing loops over interval start_year-stop_year
      character (len=*) ,intent(in)             :: clm_input_path             ! file path for input
      character (len=*) ,intent(in)             :: clm_surf_path            ! file path for surface data
      
      !OUTPUT
      real(r8),intent(out)                      :: pool_C_final(nlevdecomp,pool_types)     ! For store and output final C pool sizes 
      real(r8),intent(out)                      :: pool_N_final(nlevdecomp,pool_types_N)   ! For storing and output final N pool sizes [gN/m3] 
      
      
     !Shape of pool_matrixC/change_matrixC
     !|       LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc |
     !|level1   1   2    3    4   5   6   7   8    9    10  |
     !|level2                                               |
     !| .                                                   |
     !| .                                                   |
     !|nlevdecomp __________________________________________|

     !Shape of the pool_matrixN/change_matrixN
     !|       LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc NH4 NO3|
     !|level1   1   2    3    4   5   6   7   8    9    10  11  12 |
     !|level2                                                      |
     !| .                                                          |
     !| .                                                          |
     !|nlevdecomp _________________________________________________|
     
      !LOCAL
      character (len=4)              :: year_fmt
      character (len=4)              :: year_char
      logical                        :: isVertical                           ! True if we use vertical soil layers.
      real(r8),dimension(nlevdecomp) :: HR                                   ! For storing the C  that is lost to respiration [gC/m3h]
      real(r8)                       :: HR_mass_accumulated, HR_mass
      real(r8)                       :: pool_matrixC(nlevdecomp,pool_types)     ! For storing C pool sizes [gC/m3]
      real(r8)                       :: change_matrixC(nlevdecomp,pool_types)   ! For storing dC/dt for each time step [gC/(m3*hour)]
      real(r8)                       :: pool_temporaryC(nlevdecomp,pool_types)  ! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
      real(r8)                       :: pool_matrixC_previous(nlevdecomp,pool_types) !Used for checking mass conservation
      real(r8)                       :: pool_matrixN_previous(nlevdecomp,pool_types_N) !Used for checking mass conservation
      real(r8)                       :: pool_temporaryN(nlevdecomp,pool_types_N)! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
      real(r8)                       :: pool_matrixN(nlevdecomp,pool_types_N)   ! For storing N pool sizes [gN/m3] parallell to C pools and  inorganic N
      real(r8)                       :: change_matrixN(nlevdecomp,pool_types_N) ! For storing dC/dt for each time step [gN/(m3*hour)]
      
      real(r8),dimension(nlevdecomp) :: ROI
      real(r8)                       :: sum_consN(nlevdecomp, pool_types_N) !g/m3 for calculating annual mean
      real(r8)                       :: sum_consC(nlevdecomp, pool_types) !g/m3 for calculating annual mean
      real(r8)                       :: N_DEPinput
      real(r8)                       :: C_MYCinput
      real(r8)                       :: C_leaf_litter
      real(r8)                       :: C_root_litter
      real(r8)                       :: N_leaf_litter
      real(r8)                       :: N_root_litter
      real(r8)                       :: C_CWD_litter(nlevdecomp)
      real(r8)                       :: N_CWD_litter(nlevdecomp)

      real(r8)                       :: dt                            ! size of time step
      real(r8)                       :: time                          ! t*dt
      real(r8)                       :: C_Loss, C_Gain, N_Gain, N_Loss,N_exchange
      real(r8)                       :: tot_diffC,upperC,lowerC                 ! For the call to vertical_diffusion
      real(r8)                       :: tot_diffN,upperN,lowerN                 ! For the call to vertical_diffusion
      real(r8)                       :: sum_input_step         ! Used for checking mass conservation
      real(r8)                       :: sum_input_total        ! Used for checking mass conservation
      real(r8)                       :: sum_N_input_total
      real(r8)                       :: sum_N_out_total
      real(r8)                       :: sum_N_out_step
      real(r8)                       :: sum_N_input_step       ! Used for checking mass conservation

      real(r8),allocatable           :: vertC(:,:)         !Stores the vertical change in a time step, same shape as change_matrixC
      real(r8),allocatable           :: vertN(:,:)         !Stores the vertical change in a time step, same shape as change_matrixC

      !Counters
      integer                        :: ycounter, year,write_y
      integer                        :: counter            !used for determining when to output results
      integer                        :: month_counter      !for determining when to read new input values (Temp. litfall etc.)
      integer                        :: day_counter
      integer                        :: spinup_counter
      integer                        :: current_month
      integer                        :: current_day
      logical                        :: Spinup_run
      integer                        :: j,i,t              !for iterations
      integer,parameter              ::t_init=1
      integer                        :: date
      integer                       :: input_steps

      real(r8)                       :: change_sum(nlevdecomp, pool_types)
      real(r8)                       :: vertC_change_sum(nlevdecomp, pool_types)
      real(r8),dimension(15)         :: PFT_distribution
      real(r8)                       :: nh4_frac
      real(r8)                       :: NH4_temporary
      real(r8)                       :: NO3_temporary
      
      
      !For reading soil temperature and moisture from CLM output file
      real(r8), dimension(nlevdecomp)          :: TSOIL
      real(r8), dimension(nlevdecomp)          :: SOILLIQ
      real(r8), dimension(nlevdecomp)          :: SOILICE
      real(r8), dimension(nlevdecomp)          :: WATSAT
      real(r8), dimension(nlevdecomp)          :: W_SCALAR
      real(r8), dimension(nlevdecomp)          :: T_SCALAR
      real(r8), dimension(nlevdecomp)          :: r_moist
      real(r8)                                 :: drain
      real(r8)                                 :: h2o_liq_tot
          
      integer :: ncid
      integer :: writencid
      integer :: spinupncid
      
      call system_clock(count_rate=clock_rate) !Find the time rate
      call system_clock(count=clock_start)     !Start Timer  
      
      dt= 1.0/step_frac !Setting the time step
      
      if (nlevdecomp>1) then
        soil_depth=sum(delta_z(1:nlevdecomp))
        isVertical = .True.
      else
        soil_depth=1.52_r8
        isVertical = .False.
        allocate (vertC, mold = pool_matrixC)
        allocate (vertN, mold = pool_matrixN)
                     !TODO: This can be done better
      end if
      
      allocate(CUE_bacteria_vr(nlevdecomp))
      allocate(CUE_fungi_vr(nlevdecomp))
      allocate(CUE_ecm_vr(nlevdecomp))
      allocate(CUE_erm_vr(nlevdecomp))
      allocate(CUE_am_vr(nlevdecomp))
      CUE_fungi_vr=CUE_0
      CUE_bacteria_vr=CUE_0
      CUE_ecm_vr=CUE_myc_0
      CUE_am_vr=CUE_myc_0      
      CUE_erm_vr=0.0
      c1a=0
      c1b=0
      c2=0
      c3a=0
      c3b=0
      c4a=0
      c4b=0
      
      ! Fracions of SAP that goes to different SOM pools
      fPHYS = (/ 0.3 * exp(1.3*fCLAY), 0.2 * exp(0.8*fCLAY) /)
      fCHEM = (/0.1 * exp(-3.0*fMET), 0.3 * exp(-3.0*fMET) /)
      fAVAIL = 1-(fPHYS+fCHEM)

      !Set initial concentration values:
      pool_matrixC=pool_C_start
      pool_matrixN=pool_N_start
      pool_matrixC_previous = pool_C_start
      pool_matrixN_previous = pool_N_start
      
      !Make sure things start from zero
      change_matrixC = 0.0
      change_matrixN = 0.0
      HR             = 0.0
      ROI             = 0.0
      
      vertC_change_sum=0.0
      counter  = 0
      ycounter = 0
      HR_mass_accumulated = 0
      sum_input_step   =0.0
      sum_input_total  =0.0
      sum_N_input_step =0.0
      sum_N_input_total=0.0
      sum_N_out_step   =0.0
      sum_N_out_total  =0.0

      year     = start_year
      year_fmt = '(I4)'
      write (year_char,year_fmt) year
      write_y=0
      current_month = 1
      current_day   = 1      
      month_counter = 0
      day_counter = 0
      
      !data from CLM file
      if ( start_year == 1850 ) then
        Spinup_run = .True.
        spinup_counter =1
        call check(nf90_open(trim(adjustr(clm_input_path)//'_historical.clm2.for_spinup.1850-1869.nc'), nf90_nowrite, spinupncid)) !open netcdf containing values for the next year  
        call read_time(spinupncid,input_steps) !Check if inputdata is daily or monthly:         
        call read_clm_model_input(spinupncid,nlevdecomp,Spinup_counter, &
        N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
        C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
        W_SCALAR,T_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)        
      else
        Spinup_run = .False.        
        call check(nf90_open(trim(adjustr(clm_input_path)//'_historical.clm2.all.'//year_char//'.nc'), nf90_nowrite, ncid)) !open netcdf containing values for the next year  
        print*, trim(adjustr(clm_input_path)//'_historical.clm2.all.'//year_char//'.nc')
        call read_time(ncid,input_steps) !Check if inputdata is daily or monthly: 
        call read_clm_model_input(ncid,nlevdecomp,1, &
        N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
        C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
        W_SCALAR,T_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)
      end if

      allocate(ndep_prof(nlevdecomp),leaf_prof(nlevdecomp),froot_prof(nlevdecomp))   
         
      call read_WATSAT_and_profiles(adjustr(clm_input_path)//'_historical.clm2.all.'//"1901.nc",WATSAT,ndep_prof,froot_prof,leaf_prof, nlevdecomp)        
      call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist,nlevdecomp)                   
      call read_clay(adjustr(clm_surf_path),fCLAY,nlevdecomp)
      call read_PFTs(adjustr(clm_surf_path),PFT_distribution)
      f_EcM = calc_EcMfrac(PFT_distribution)
      if ( Spinup_run ) then
        max_mining = read_maxC(spinupncid,f_EcM,input_steps)
      else
        max_mining = read_maxC(ncid,f_EcM,input_steps)        
      end if
      
      !open and prepare files to store results. Store initial values
      !call create_yearly_mean_netcdf(run_name,nlevdecomp)
      call create_netcdf(run_name, nlevdecomp)
      call check(nf90_open(output_path//trim(run_name)//".nc", nf90_write, writencid))
      
      call fill_netcdf(writencid,t_init, pool_matrixC, change_matrixC, pool_matrixN,change_matrixN, &
                       date, HR_mass_accumulated,HR, change_matrixC,change_matrixN,write_hour,current_month, &
                      TSOIL, r_moist,CUE_bacteria_vr,CUE_fungi_vr,CUE_EcM_vr,CUE_am_vr,levsoi=nlevdecomp,ROI=ROI)
            
      desorb = 1.5e-5*exp(-1.5*(fclay)) !NOTE: desorb and pscalar moved from paramMod bc fCLAY is read in decomp subroutine (13.09.2021)
      pscalar = 1.0/(2*exp(-2.0*sqrt(fCLAY)))
      Kmod = [real(r8) :: 0.125,0.5,0.25*pscalar,0.5,0.25,0.167*pscalar]
      !----------------------------------------------------------------------------------------------------------------
      do t =1,nsteps !Starting time iterations
        time = t*dt
        !Time counters:
        counter  = counter + 1
        ycounter = ycounter + 1            !Counts hours in a year 
        month_counter = month_counter + 1 !Counts hours in a month
        day_counter = day_counter + 1     !Counts hours in a day (24)

        ! !Update temp and moisture values monthly/daily 
        if (month_counter == days_in_month(current_month)*hr_pr_day*step_frac+1) then
            month_counter = 1       
                 
            if (current_month == 12) then
              current_month=1 !Update to new year
            else
              current_month = current_month + 1 
            end if   
   
            if (input_steps==12) then  
              call read_clm_model_input(ncid,nlevdecomp,current_month, &
                                N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
                                C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
                                W_SCALAR,T_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)  
              call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist,nlevdecomp)   
              max_mining = read_maxC(ncid,f_EcM,input_steps)
            end if     
            
            if (input_steps==240) then
              spinup_counter = spinup_counter+1

              call read_clm_model_input(spinupncid,nlevdecomp,spinup_counter, &
                                N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
                                C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
                                W_SCALAR,T_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)  
              call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist,nlevdecomp)   
              max_mining = read_maxC(spinupncid,f_EcM,input_steps)
            end if                             
        end if 
        
        if (day_counter == hr_pr_day*step_frac+1) then
          day_counter = 1   
          current_day = current_day +1                      
          if ( current_day == 366 ) then
            current_day = 1
          end if
          
          if (input_steps==365) then                        
            call read_clm_model_input(ncid,nlevdecomp,current_day, &
            N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
            C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
            W_SCALAR,T_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)
            call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist,nlevdecomp)        
            max_mining = read_maxC(ncid,f_EcM,input_steps)                
          end if        
        end if 

        !print initial values to terminal
        if (t == 1) then
          call disp("InitC", pool_matrixC)
          call disp("InitN", pool_matrixN)
        end if

        do j = 1, nlevdecomp !For each depth level (for the no vertical transport case, nlevdecomp = 1, so loop is only done once):

          k_mycsom  = (/1.14_r8,1.14_r8,1.14_r8/)*1e-6  ![1/h] Decay constants, mycorrhiza to SOM pools TODO: Assumed, needs revision TODO: 1/yr instead? Sulman et al. 

          !Michaelis Menten parameters:
          Km      = Km_function(TSOIL(j))
          Vmax    = Vmax_function(TSOIL(j),r_moist(j)) !  ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics.

          ![1/h] Microbial turnover rate (SAP to SOM)

          tau = (/ 5.2e-4*exp(0.3_r8*fMET)*r_moist(j), 2.4e-4*exp(0.1_r8*fMET)*r_moist(j)/)
          CUE_bacteria_vr(j) = (CUE_slope*TSOIL(j)+CUE_0)
          CUE_fungi_vr(j) = (CUE_slope*TSOIL(j)+CUE_0)
          CUE_ecm_vr(j) = CUE_myc_0
          CUE_am_vr(j) = CUE_myc_0
          
          !Determine input rates (from CLM data) in timestep
          call input_rates(j,f_EcM,C_leaf_litter,C_root_litter,N_leaf_litter,&
                                      N_root_litter,C_MYCinput, &
                                      N_CWD_litter,C_CWD_litter,&
                                      C_PlantLITm,C_PlantLITs, &
                                      N_PlantLITm,N_PlantLITs, &
                                      C_PlantEcM,C_PlantAM, &
                                      C_PlantSOMp,C_PlantSOMa,C_PlantSOMc, &
                                      N_PlantSOMp,N_PlantSOMa,N_PlantSOMc)

          !Determine deposition, Leaching and nitrification in timestep
          Deposition = set_N_dep(CLMdep = N_DEPinput*ndep_prof(j)) !NOTE: either const_dep = some_value or CLMdep = N_DEPinput*ndep_prof(j)
          Leaching = calc_Leaching(drain,h2o_liq_tot,pool_matrixN(j,12))
          nitrif_rate=nitrification((pool_matrixN(j,11)+Deposition*dt),W_SCALAR(j),T_SCALAR(j),TSOIL(j))

          !Calculate fluxes between pools in level j (file: fluxMod.f90):
          call calculate_fluxes(j,nlevdecomp, pool_matrixC, pool_matrixN,dt)
          ROI(j) = ROI_function(N_INEcM+N_SOMpEcM+N_SOMcEcM,pool_matrixC(j,5),k_mycsom(1))
!---------Calc fraction of inorg N that is NH4----------------------------------------          
          NH4_temporary = pool_matrixN(j,11) + (Deposition - nitrif_rate)*dt
          NO3_temporary = pool_matrixN(j,12) + (nitrif_rate - Leaching)*dt          
          if (NH4_temporary+NO3_temporary == 0._r8) Then
            nh4_frac = 0.5_8
          else
            nh4_frac = NH4_temporary/(NH4_temporary+NO3_temporary)
          end if
!---------------------------------------------------------------------------------------          
        
          if (counter == write_hour*step_frac .or. t==1) then !Write fluxes from calculate_fluxes to file
           call fluxes_netcdf(writencid,int(time), write_hour, j)
          end if !write fluxes
          do i = 1,pool_types_N !loop over all the pool types, i, in depth level j 
            !This if-loop calculates dC/dt and dN/dt for the different carbon pools.
            !NOTE: If pools are added/removed (i.e the actual model equations is changed), this loop needs to be updated.

            if (i==1) then !LITm
              N_Gain = N_PlantLITm
              N_Loss = N_LITmSAPb + N_LITmSAPf
              C_Gain = C_PlantLITm
              C_Loss = C_LITmSAPb + C_LITmSAPf
              
            elseif (i==2) then !LITs
              N_Gain = N_PlantLITs
              N_Loss = N_LITsSAPb + N_LITsSAPf
              C_Gain = C_PlantLITs
              C_Loss = C_LITsSAPb + C_LITsSAPf

            elseif (i==3) then !SAPb              
              C_Gain = CUE_bacteria_vr(j)*(C_LITmSAPb + C_LITsSAPb &
                + C_SOMaSAPb)
              C_Loss =  C_SAPbSOMp + C_SAPbSOMa + C_SAPbSOMc
              N_Gain = N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb
              N_Loss = N_SAPbSOMp + N_SAPbSOMa + N_SAPbSOMc
              if ( N_INSAPb>0 ) then
                N_Gain = N_Gain + N_INSAPb              
              else
                N_Loss=N_Loss-N_INSAPb     !two minus becomes +         
              end if

            elseif (i==4) then !SAPf
              C_Gain = CUE_fungi_vr(j)*(C_LITmSAPf + C_LITsSAPf &
                + C_SOMaSAPf)
              C_Loss =  C_SAPfSOMp + C_SAPfSOMa + C_SAPfSOMc
              N_Gain = N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf
              N_Loss = N_SAPfSOMp + N_SAPfSOMa + N_SAPfSOMc
              if ( N_INSAPf>0 ) then
                N_Gain = N_Gain + N_INSAPf 
              else
                N_Loss=N_Loss-N_INSAPf !two minus becomes +
              end if

            elseif (i==5) then !EcM
              C_Gain = CUE_ecm_vr(j)*C_PlantEcM + C_SOMcEcM + C_SOMpEcM !
              C_Loss = C_EcMSOMp + C_EcMSOMa + C_EcMSOMc + CUE_ecm_vr(j)*C_PlantEcM*enzyme_pct
              N_Gain = N_INEcM + N_SOMpEcM + N_SOMcEcM
              N_Loss = N_EcMPlant + N_EcMSOMa + N_EcMSOMp + N_EcMSOMc
            elseif (i==6) then !ErM
              C_Gain = CUE_erm_vr(j)*C_PlantErM
              C_Loss = C_ErMSOMp + C_ErMSOMa + C_ErMSOMc 
              N_Gain = N_INErM
              N_Loss = N_ErMPlant + N_ErMSOMa + N_ErMSOMp + N_ErMSOMc

            elseif (i==7) then !AM
              C_Gain = CUE_am_vr(j)*C_PlantAM
              C_Loss = C_AMSOMp + C_AMSOMa + C_AMSOMc
              N_Gain = N_INAM 
              N_Loss = N_AMPlant + N_AMSOMa + N_AMSOMp + N_AMSOMc

            elseif (i==8) then !SOMp
              C_Gain =  C_SAPbSOMp + C_SAPfSOMp + C_EcMSOMp + C_ErMSOMp + C_AMSOMp+ C_PlantSOMp
              C_Loss = C_SOMpSOMa+C_EcMdecompSOMp + C_SOMpEcM
              N_Gain =  N_SAPbSOMp + N_SAPfSOMp + N_EcMSOMp + N_ErMSOMp + N_AMSOMp+N_PlantSOMp
              N_Loss = N_SOMpSOMa + N_SOMpEcM

            elseif (i==9) then !SOMa
               C_Gain = C_SAPbSOMa + C_SAPfSOMa + C_EcMSOMa + C_EcMdecompSOMp + C_EcMdecompSOMc &
               + C_ErMSOMa + C_AMSOMa + C_SOMpSOMa + C_SOMcSOMa + C_PlantSOMa+ CUE_ecm_vr(j)*C_PlantEcM*enzyme_pct
               C_Loss = C_SOMaSAPb + C_SOMaSAPf 
               N_Gain = N_SAPbSOMa + N_SAPfSOMa + N_EcMSOMa + &
               N_ErMSOMa + N_AMSOMa + N_SOMpSOMa + N_SOMcSOMa +N_PlantSOMa
               N_Loss = N_SOMaSAPb + N_SOMaSAPf + N_SOMaEcM

            elseif (i==10) then !SOMc
              C_Gain =  C_SAPbSOMc + C_SAPfSOMc + C_EcMSOMc + C_ErMSOMc + C_AMSOMc+C_PlantSOMc
              C_Loss = C_SOMcSOMa+C_EcMdecompSOMc + C_SOMcEcM
              N_Gain =  N_SAPbSOMc + N_SAPfSOMc + N_EcMSOMc + N_ErMSOMc + N_AMSOMc+N_PlantSOMc
              N_Loss = N_SOMcSOMa + N_SOMcEcM
              
            elseif (i == 11) then !NH4 inorganic N
              N_Gain = Deposition
              N_exchange= nh4_frac*(N_INSAPb + N_INSAPf) !N_exchange can act both as a sink and a source, depending on the SAP demand for N
              N_Loss = nitrif_rate + nh4_frac*(N_INEcM + N_InPlant  + N_INAM )!+ N_INErM
              change_matrixN(j,i) = N_Gain - N_exchange - N_loss 
              
            elseif (i == 12) then !NO3 inorganic N
              N_Gain = nitrif_rate
              N_exchange= (1-nh4_frac)*(N_INSAPb + N_INSAPf) !N_exchange can act both as a sink and a source, depending on the SAP demand for N
              N_Loss = Leaching + (1-nh4_frac)*(N_INEcM + N_InPlant  + N_INAM )!+ N_INErM
              change_matrixN(j,i) = N_Gain - N_exchange - N_loss

            else
              print*, 'Too many pool types expected, pool_types = ',pool_types, 'i: ', i
            end if !determine total gains and losses

            if (i < 11) then !Carbon matrix does only have 10 columns (if i == 11 or 12 this is handeled inside the loop over pools)
                change_matrixC(j,i) = C_Gain - C_Loss !net change in timestep
                change_matrixN(j,i) = N_Gain - N_loss            
                !For summarize total change between each written output
                !change_sum(j,i)= change_sum(j,i) + change_matrixC(j,i)*dt                
                !Store these values as temporary so that they can be used in the vertical diffusion subroutine
                pool_temporaryC(j,i)=pool_matrixC(j,i) + change_matrixC(j,i)*dt
            end if

            pool_temporaryN(j,i) =pool_matrixN(j,i) + change_matrixN(j,i)*dt
            
            !NOTE: This is introduced to avoid very small errors that may make inorganic N pools negative (~ E-020). Overall mass balance is still within error limits
            if ( abs(pool_temporaryN(j,i)) < 1e-18 ) then
              pool_temporaryN(j,i)=0._r8
            end if
            
            if (isnan(pool_temporaryN(j,i))) then
              print*, 'NaN NITROGEN value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryN(j,i)
              stop
            end if
            if (pool_temporaryN(j,i) < 0.000) then
              print*, 'Too small pool size: NITROGEN value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryN(j,i)
              print*, pool_matrixN(j,11), pool_matrixN(j,12), pool_matrixN(j,11)+ pool_matrixN(j,12)
              print*, change_matrixN(j,11), change_matrixN(j,12), change_matrixN(j,11)+ change_matrixN(j,12)
              print*, pool_temporaryN(j,11), pool_temporaryN(j,12), pool_temporaryN(j,11)+ pool_temporaryN(j,12), nitrif_rate
              call disp(pool_temporaryN)
              call disp(pool_temporaryC)
              call disp("C:N : ",pool_matrixC/pool_matrixN(:,1:10))              
              stop
            end if
            
            if (i < 11 ) then
              if (isnan(pool_temporaryC(j,i))) then
                print*, 'NaN CARBON value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryC(j,i)
                stop
              end if
              if (pool_temporaryC(j,i) < 0.000) then
                print*, 'Too small pool size: CARBON value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryC(j,i)
                stop
              end if
              if ( abs(pool_temporaryN(j,i)) < 1e-18 ) then
                pool_temporaryN(j,i)=0._r8
              end if
            end if

          end do !i, pool_types

          !Calculate the heterotrophic respiration loss from depth level j in timestep t:
          HR(j) =(( C_LITmSAPb + C_LITsSAPb  + C_SOMaSAPb)*(1-CUE_bacteria_vr(j)) + (C_LITmSAPf &
          + C_LITsSAPf + C_SOMaSAPf)*(1-CUE_fungi_vr(j)))*dt
          if (HR(j) < 0 ) then
            print*, 'Negative HR: ', HR(j), t
          end if
          
          !Summarize in and out print timestep to check mass balance
          sum_input_step=sum_input_step+(C_PlantLITm+C_PlantLITs+CUE_ecm_vr(j)*C_PlantEcM+CUE_am_vr(j)*C_PlantAM+C_PlantSOMc+C_PlantSOMp+C_PlantSOMa)*dt*delta_z(j) !g/m2
          sum_N_input_step=sum_N_input_step+(N_PlantLITm+N_PlantLITs+N_PlantSOMc+N_PlantSOMp+N_PlantSOMa+Deposition)*dt*delta_z(j) !g/m2
          sum_N_out_step=sum_N_out_step+(N_EcMPlant+N_AMPlant+N_INPlant+Leaching)*dt*delta_z(j)

        end do !j, depth_level

        !Store accumulated HR mass
        call respired_mass(HR, HR_mass,nlevdecomp)
        HR_mass_accumulated = HR_mass_accumulated + HR_mass
        
        if (isVertical) then
          call vertical_diffusion(tot_diffC,upperC,lowerC, pool_temporaryC,vertC,D_carbon)
          call vertical_diffusion(tot_diffN,upperN,lowerN, pool_temporaryN,vertN,D_nitrogen)
          pool_matrixC =  vertC*dt + pool_temporaryC
          pool_matrixN = vertN*dt + pool_temporaryN
        else
          pool_matrixC=pool_temporaryC
          pool_matrixN=pool_temporaryN
        end if!isVertical

        if (ycounter == 365*24*step_frac) then
          ycounter = 0
          write_y =write_y+1 !For writing to annual mean file
          !call annual_mean(sum_consC,sum_consN, nlevdecomp,write_y , run_name) !calculates the annual mean and write the result to file
          if (year == stop_year) then
            year = start_year         
            spinup_counter=0            
          else 
            year = year + 1             
          end if
          
          write (year_char,year_fmt) year
          
          if ( .not. Spinup_run ) then            
            call check(nf90_close(ncid)) !Close netcdf file containing values for the past year
            call check(nf90_open(trim(adjustr(clm_input_path)//'_historical.clm2.all.'//year_char//'.nc'), nf90_nowrite, ncid)) !open netcdf containing values for the next year
            call read_time(ncid,input_steps)     
          end if 
          
          sum_consN =0
          sum_consC =0
        end if

        if (counter == write_hour*step_frac) then
          counter = 0        
          call fill_netcdf(writencid, int(time), pool_matrixC, change_matrixC, pool_matrixN,change_matrixN,&
           date, HR_mass_accumulated,HR,vertC,vertN, write_hour,current_month, &
           TSOIL, r_moist,CUE_bacteria_vr,CUE_fungi_vr,CUE_ecm_vr,CUE_am_vr,nlevdecomp,ROI=ROI)
           change_sum = 0.0
        end if!writing

        !Write end values to terminal
        if (t == nsteps) then
          call disp("pool_matrixC gC/m3 ",pool_matrixC)
          call disp("pool_matrixN gN/m3 ",pool_matrixN)          
          call disp("C:N : ",pool_matrixC/pool_matrixN(:,1:10))
          pool_C_final=pool_matrixC
          pool_N_final=pool_matrixN    
          call store_parameters(writencid)    
        end if
        
        call test_mass_conservation(sum_input_step,HR_mass,pool_matrixC_previous,pool_matrixC,nlevdecomp,pool_types)
        call test_mass_conservation(sum_N_input_step,sum_N_out_step,pool_matrixN_previous,pool_matrixN,nlevdecomp,pool_types_N)
        pool_matrixC_previous=pool_matrixC
        pool_matrixN_previous=pool_matrixN
        sum_input_total=sum_input_total+sum_input_step
        sum_input_step=0.0
        sum_N_input_total=sum_N_input_total+sum_N_input_step  
        sum_N_input_step=0.0
        sum_N_out_total=sum_N_out_total+sum_N_out_step    
        sum_N_out_step=0.0
      end do !t
      
      !Check total mass conservation 
      call total_mass_conservation(sum_input_total,HR_mass_accumulated,pool_C_start,pool_C_final,nlevdecomp,pool_types)
      call total_nitrogen_conservation(sum_N_input_total,sum_N_out_total,pool_N_start,pool_N_final,nlevdecomp,pool_types_N)
      
      call check(nf90_close(writencid))
      if ( Spinup_run ) then
        call check(nf90_close(Spinupncid))
      end if
      
      print*, c1a,c1b,c2,c3a,c3b,c4a,c4b
      !deallocation
      deallocate(ndep_prof,leaf_prof,froot_prof)
      deallocate(CUE_bacteria_vr,CUE_fungi_vr, CUE_ecm_vr,CUE_am_vr, CUE_erm_vr)
      
      !For timing
      call system_clock(count=clock_stop)      ! Stop Timer
      print*, "Total time for decomp subroutine in seconds: ", real(clock_stop-clock_start)/real(clock_rate)
      print*, "Total time for decomp subroutine in minutes: ", (real(clock_stop-clock_start)/real(clock_rate))/60
      
    end subroutine decomp

  subroutine annual_mean(yearly_sumC,yearly_sumN,nlevels, year, run_name)
    REAL(r8), DIMENSION(nlevels,pool_types)  , intent(in):: yearly_sumC
    REAL(r8), DIMENSION(nlevels,pool_types_N), intent(in):: yearly_sumN
    integer, intent(in) :: year
    integer , intent(in):: nlevels
    CHARACTER (len = *), intent(in):: run_name
    !Local
    REAL(r8), DIMENSION(nlevels,pool_types) :: yearly_meanC
    REAL(r8), DIMENSION(nlevels,pool_types_N) :: yearly_meanN
    integer, parameter                         :: hr_in_year = 24*365

    yearly_meanC=yearly_sumC/hr_in_year
    yearly_meanN=yearly_sumN/hr_in_year

    !call fill_yearly_netcdf(run_name, year, yearly_meanC,yearly_meanN,nlevels)

  end subroutine annual_mean
  
end module mycmim
