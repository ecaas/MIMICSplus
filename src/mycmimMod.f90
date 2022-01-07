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
      real(r8),intent(in)                      :: pool_N_start(nlevdecomp,pool_types+1)   ! For storing and output final N pool sizes [gN/m3] 
      integer, intent(in)                       :: start_year !Forcing start year
      integer, intent(in)                       :: stop_year !Forcing end year, i.e. forcing loops over interval start_year-stop_year
      character (len=*) ,intent(in)             :: clm_input_path             ! file path for input
      character (len=*) ,intent(in)             :: clm_surf_path            ! file path for surface data
      
      
      !OUTPUT
      real(r8),intent(out)                      :: pool_C_final(nlevdecomp,pool_types)     ! For store and output final C pool sizes 
      real(r8),intent(out)                      :: pool_N_final(nlevdecomp,pool_types+1)   ! For storing and output final N pool sizes [gN/m3] 
      
      
     !Shape of pool_matrixC/change_matrixC
     !|       LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc |
     !|level1   1   2    3    4   5   6   7   8    9    10  |
     !|level2                                               |
     !| .                                                   |
     !| .                                                   |
     !|nlevdecomp __________________________________________|

     !Shape of the pool_matrixN/change_matrixN
     !|       LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc IN|
     !|level1   1   2    3    4   5   6   7   8    9    10  11|
     !|level2                                                 |
     !| .                                                     |
     !| .                                                     |
     !|nlevdecomp ____________________________________________|
     
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
      real(r8)                       :: pool_temporaryN(nlevdecomp,pool_types+1)! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
      real(r8)                       :: pool_matrixN(nlevdecomp,pool_types+1)   ! For storing N pool sizes [gN/m3] parallell to C pools and  inorganic N
      real(r8)                       :: change_matrixN(nlevdecomp,pool_types+1) ! For storing dC/dt for each time step [gN/(m3*hour)]
      

      real(r8)                       :: sum_consN(nlevdecomp, pool_types+1) !g/m3 for calculating annual mean
      real(r8)                       :: sum_consC(nlevdecomp, pool_types) !g/m3 for calculating annual mean
      real(r8)                       :: N_DEPinput
      real(r8)                       :: C_EcMinput
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
      integer                        :: current_month
      integer                        :: current_day
            
      integer                        :: j,i,t              !for iterations
      integer,parameter              ::t_init=1
      integer                        :: date
      integer                       :: input_steps

      real(r8)                       :: change_sum(nlevdecomp, pool_types)
!      real(r8)                       :: vertN_change_sum(nlevdecomp, pool_types)
      real(r8)                       :: vertC_change_sum(nlevdecomp, pool_types)

      !For reading soil temperature and moisture from CLM output file
      real(r8), dimension(nlevdecomp)          :: TSOIL
      real(r8), dimension(nlevdecomp)          :: SOILLIQ
      real(r8), dimension(nlevdecomp)          :: SOILICE
      real(r8), dimension(nlevdecomp)          :: WATSAT
      real(r8), dimension(nlevdecomp)          :: W_SCALAR
      real(r8), dimension(nlevdecomp)          :: r_moist
      real(r8)                                 :: drain
      real(r8)                                 :: h2o_liq_tot
          
      integer :: ncid

      dt= 1.0/step_frac !Setting the time step
      
      if (nlevdecomp>1) then
        soil_depth=sum(delta_z(1:nlevdecomp))
        isVertical = .True.
      else
        soil_depth=1.52
        isVertical = .False.
        !delta_z=soil_depth !So that delta_z will not be 1st on delta_z from parametersMod
        allocate (vertC, mold = pool_matrixC)
        allocate (vertN, mold = pool_matrixN)
                     !TODO: This can be done better
      end if
      
      allocate(CUE_bacteria_vr(nlevdecomp))
      allocate(CUE_fungi_vr(nlevdecomp))
      CUE_fungi_vr=0.5
      CUE_bacteria_vr=0.5
      
      ! Fracions of SAP that goes to different SOM pools
      fPHYS = (/ 0.3 * exp(fCLAY), 0.3 * exp(fCLAY) /)
      fCHEM = (/0.1 * exp(-3.0*fMET), 0.1 * exp(-3.0*fMET) /)
      fAVAIL = 1-(fPHYS+fCHEM)

      !Set initial concentration values:
      pool_matrixC=pool_C_start
      pool_matrixN=pool_N_start
      pool_matrixC_previous = pool_C_start
      !Make sure things start from zero
      change_matrixC = 0.0
      change_matrixN = 0.0
      HR             = 0.0
      vertC_change_sum=0.0
      counter  = 0
      ycounter = 0
      HR_mass_accumulated = 0
      sum_input_step =0.0
      sum_input_total=0.0
      sum_N_input_step =0.0
      sum_N_input_total=0.0
      sum_N_out_step =0.0
      sum_N_out_total=0.0
      

      year     = start_year
      year_fmt = '(I4)'
      write (year_char,year_fmt) year
      current_month = 1
      current_day   = 1
      
      month_counter = 0
      day_counter = 0
      
      write_y=0

      call check(nf90_open(trim(adjustr(clm_input_path)//'all.'//year_char//'.nc'), nf90_nowrite, ncid)) !open netcdf containing values for the next year  

      !Check if inputdata is daily or monthly:      
      call read_time(adjustr(clm_input_path)//'all.'//year_char//'.nc',input_steps)

      !data from CLM file
      call read_clm_model_input(ncid,nlevdecomp,1, &
                                N_leaf_litter,N_root_litter,C_EcMinput,N_DEPinput, &
                                C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
                                W_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)

      allocate(ndep_prof(nlevdecomp),leaf_prof(nlevdecomp),froot_prof(nlevdecomp))   
         
      call read_WATSAT_and_profiles(adjustr(clm_input_path)//'all.'//"1901.nc",WATSAT,ndep_prof,froot_prof,leaf_prof, nlevdecomp)        
      call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist,nlevdecomp)                   
      call read_clay(adjustr(clm_surf_path),fCLAY,nlevdecomp)
  
      !open and prepare files to store results. Store initial values
      !call create_yearly_mean_netcdf(run_name,nlevdecomp)
      call create_netcdf(run_name, nlevdecomp)
      
      call fill_netcdf(run_name,t_init, pool_matrixC, change_matrixC, pool_matrixN,change_matrixN, &
                       date, HR_mass_accumulated,HR, change_matrixC,change_matrixN,write_hour,current_month, &
                      TSOIL, r_moist,CUE_bacteria_vr,CUE_fungi_vr,levsoi=nlevdecomp)
                      
      desorb = 1.5e-5*exp(-1.5*(fclay)) !NOTE: desorb and pscalar moved from paramMod bc fCLAY is read in decomp subroutine (13.09.2021)
      pscalar = 1.0/(2*exp(-2.0*sqrt(fCLAY)))
      Kmod = [real(r8) :: 0.125,0.5,0.25*pscalar,0.5,0.25,0.167*pscalar]
      !----------------------------------------------------------------------------------------------------------------
      do t =1,nsteps !Starting time iterations
        time = t*dt
        counter  = counter + 1
        ycounter = ycounter + 1
        month_counter = month_counter + 1 !Counts hours in a month
        day_counter = day_counter + 1
        ! !Update temp and moisture values monthly/daily 
        if (month_counter == days_in_month(current_month)*24*step_frac+1) then
            month_counter = 1
            
            if (current_month == 12) then
              current_month=1 !Update to new year
            else 
              current_month = current_month + 1 
            end if             
             
            if (input_steps==12) then  
              !print*, ncid,nlevdecomp,current_month,time
              call read_clm_model_input(ncid,nlevdecomp,current_month, &
                                N_leaf_litter,N_root_litter,C_EcMinput,N_DEPinput, &
                                C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
                                W_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)  
         
              call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist,nlevdecomp)   
            end if    
                                      
        end if 
        
        if (day_counter == 24*step_frac+1) then
          if ( current_day == 366 ) then
            current_day = 1
          else 
            current_day = current_day +1            
          end if          
          if (input_steps==365) then
            call read_clm_model_input(ncid,nlevdecomp,current_day, &
            N_leaf_litter,N_root_litter,C_EcMinput,N_DEPinput, &
            C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
            W_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)
            call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist,nlevdecomp)            
          end if        
          day_counter = 1   
        end if 

        !print initial values to terminal
        if (t == 1) then
          call disp("InitC", pool_matrixC)
          call disp("InitN", pool_matrixN)
        end if

        do j = 1, nlevdecomp !For each depth level (for the no vertical transport case, nlevdecomp = 1, so loop is only done once):
          !Michaelis Menten parameters:

          k_mycsom  = (/1.14,1.14,1.14/)*1e-4  ![1/h] Decay constants, mycorrhiza to SOM pools TODO: Assumed, needs revision

          Km      = Km_function(TSOIL(j))
          Vmax    = Vmax_function(TSOIL(j),r_moist(j)) !  ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?

          !NOTE FROM MIMICS_CYCLE_CN:
          ! WW also modify TAU as a function of soil moisture, so things don't
          ! colapse in frozen soils...
          !mimicsbiome%tauR(npt) = mimicsbiome%tauR(npt) * fW
          !mimicsbiome%tauK(npt) = mimicsbiome%tauK(npt) * fW
          ![1/h] Microbial turnover rate (SAP to SOM), SAPr,(/1.39E-3*exp(0.3*fMET), 2.3E-4*exp(0.1*fMET)/)

          !tau = (/ 5e-4*exp(0.1*fMET)*r_moist(j), 5e-4*exp(0.1*fMET)*r_moist(j)/)
          tau = (/ 5.2e-4*exp(0.3*fMET)*r_moist(j), 2.4e-4*exp(0.1*fMET)*r_moist(j)/)
          CUE_bacteria_vr(j) = (CUE_slope*TSOIL(j)+CUE_0)
          CUE_fungi_vr(j) = (CUE_slope*TSOIL(j)+CUE_0)

          !Calculate fluxes between pools in level j (file: fluxMod.f90):
          call input_rates(j,C_leaf_litter,C_root_litter,N_leaf_litter,&
                                      N_root_litter,C_EcMinput, &
                                      N_CWD_litter,C_CWD_litter,&
                                      C_PlantLITm,C_PlantLITs, &
                                      N_PlantLITm,N_PlantLITs, &
                                      C_PlantEcM,C_PlantAM)
          !Determine deposition NOTE: Must specify name of parameters when using the set_N_dep function 
          !,const_dep = 3d0
          Deposition = set_N_dep(CLMdep = N_DEPinput*ndep_prof(j))
          Leaching = calc_Leaching(drain,h2o_liq_tot,pool_matrixN(j,11))
          
          call calculate_fluxes(j,nlevdecomp, pool_matrixC, pool_matrixN)
          
          if (counter == write_hour*step_frac .or. t==1) then !Write fluxes from calculate_fluxes to file
           call fluxes_netcdf(int(time), write_hour, j, run_name)
          end if !write fluxes
          
          do i = 1,pool_types + 1 !loop over all the pool types, i, in depth level j (+1 bc. of the added inorganic N pool)
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
              if ( N_SAPbIN<0 ) then
                N_Gain = N_Gain - N_SAPbIN !two minus becomes +              
              else
                N_Loss=N_Loss+N_SAPbIN              
              end if

            elseif (i==4) then !SAPf
              C_Gain = CUE_fungi_vr(j)*(C_LITmSAPf + C_LITsSAPf &
                + C_SOMaSAPf)
              C_Loss =  C_SAPfSOMp + C_SAPfSOMa + C_SAPfSOMc
              N_Gain = N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf
              N_Loss = N_SAPfSOMp + N_SAPfSOMa + N_SAPfSOMc
              if ( N_SAPfIN<0 ) then
                N_Gain = N_Gain - N_SAPfIN !two minus becomes +
              else
                N_Loss=N_Loss+N_SAPfIN
              end if

            elseif (i==5) then !EcM
              C_Gain = C_PlantEcM
              C_Loss = C_EcMSOMp + C_EcMSOMa + C_EcMSOMc + (1-e_m)*C_PlantEcM + e_m*C_PlantEcM*enzyme_pct
              N_Gain = N_INEcM + N_SOMpEcM + N_SOMcEcM
              N_Loss = N_EcMPlant + N_EcMSOMa + N_EcMSOMp + N_EcMSOMc

            elseif (i==6) then !ErM
              C_Gain = C_PlantErM
              C_Loss = C_ErMSOMp + C_ErMSOMa + C_ErMSOMc + (1-e_m)*C_PlantErM
              N_Gain = N_INErM
              N_Loss = N_ErMPlant + N_ErMSOMa + N_ErMSOMp + N_ErMSOMc

            elseif (i==7) then !AM
              C_Gain = C_PlantAM
              C_Loss = C_AMSOMp + C_AMSOMa + C_AMSOMc + (1-e_m)*C_PlantAM
              N_Gain = N_INAM 
              N_Loss = N_AMPlant + N_AMSOMa + N_AMSOMp + N_AMSOMc

            elseif (i==8) then !SOMp
              C_Gain =  C_SAPbSOMp + C_SAPfSOMp + C_EcMSOMp + C_ErMSOMp + C_AMSOMp+ C_PlantSOMp
              C_Loss = C_SOMpSOMa+C_EcMdecompSOMp
              N_Gain =  N_SAPbSOMp + N_SAPfSOMp + N_EcMSOMp + N_ErMSOMp + N_AMSOMp
              N_Loss = N_SOMpSOMa + N_SOMpEcM

            elseif (i==9) then !SOMa
               C_Gain = C_SAPbSOMa + C_SAPfSOMa + C_EcMSOMa + C_EcMdecompSOMp + C_EcMdecompSOMc &
               + C_ErMSOMa + C_AMSOMa + C_SOMpSOMa + C_SOMcSOMa + C_PlantSOMa+ e_m*C_PlantEcM*enzyme_pct
               C_Loss = C_SOMaSAPb + C_SOMaSAPf 
               N_Gain = N_SAPbSOMa + N_SAPfSOMa + N_EcMSOMa + &
               N_ErMSOMa + N_AMSOMa + N_SOMpSOMa + N_SOMcSOMa
               N_Loss = N_SOMaSAPb + N_SOMaSAPf 

            elseif (i==10) then !SOMc
              C_Gain =  C_SAPbSOMc + C_SAPfSOMc + C_EcMSOMc + C_ErMSOMc + C_AMSOMc+C_PlantSOMc
              C_Loss = C_SOMcSOMa+C_EcMdecompSOMc
              N_Gain =  N_SAPbSOMc + N_SAPfSOMc + N_EcMSOMc + N_ErMSOMc + N_AMSOMc
              N_Loss = N_SOMcSOMa + N_SOMcEcM

            elseif (i == 11) then !Inorganic N
              N_Gain = Deposition
              N_exchange= N_SAPbIN + N_SAPfIN !N_exchange can act both as a sink and a source, depending on the SAP demand for N
              N_Loss = Leaching + N_INEcM + N_InPlant !+ N_INErM + N_INAM 
              change_matrixN(j,i) = N_Gain + N_exchange - N_loss 
              ! print*, "Gain: ", Deposition, j,time
              ! print*, "Loss: ", Leaching,N_INEcM,N_InPlant
              ! print*, "Exchange: ", N_SAPbIN,N_SAPfIN
              ! print*, change_matrixN(j,i), pool_matrixN(j,i)
              ! print*, "-------------------------------------------------------------------------------------------------"
            else
              print*, 'Too many pool types expected, pool_types = ',pool_types, 'i: ', i
            end if !determine total gains and losses
            

            if (i /= 11) then !Carbon matrix does only have 10 columns (if i == 11 this is handeled inside the loop over pools)
                change_matrixC(j,i) = C_Gain - C_Loss !net change in timestep
                change_matrixN(j,i) = N_Gain - N_loss            
                !For summarize total change between each written output
                !print*, change_sum(j,i) + change_matrixC(j,i)*dt
                !change_sum(j,i)= change_sum(j,i) + change_matrixC(j,i)*dt
                
                !Store these values as temporary so that they can be used in the vertical diffusion subroutine
                pool_temporaryC(j,i)=pool_matrixC(j,i) + change_matrixC(j,i)*dt

            end if

            pool_temporaryN(j,i) =pool_matrixN(j,i) + change_matrixN(j,i)*dt
            !call disp("N change matrix;",change_matrixN) 
            if (isnan(pool_temporaryN(j,i))) then
              print*, 'NaN NITROGEN value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryN(j,i)
              stop
            end if
            if (pool_temporaryN(j,i) < 0.000) then
              print*, 'Too small pool size: NITROGEN value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryN(j,i)
              call disp(pool_temporaryN)
              call disp(pool_temporaryC)
              call disp("C:N : ",pool_matrixC/pool_matrixN(:,1:10))
              
              stop
            end if

            if (i /=11 ) then
              if (isnan(pool_temporaryC(j,i))) then
                print*, 'NaN CARBON value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryC(j,i)
                stop
              end if
              if (pool_temporaryC(j,i) < 0.000) then
                print*, 'Too small pool size: CARBON value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryC(j,i)
                stop
              end if
            end if

          end do !i, pool_types

          !Calculate the heterotrophic respiration loss from depth level j in timestep t: NOTE: revise!
          HR(j) =(( C_LITmSAPb + C_LITsSAPb  + C_SOMaSAPb)*(1-CUE_bacteria_vr(j)) + (C_LITmSAPf &
          + C_LITsSAPf + C_SOMaSAPf)*(1-CUE_fungi_vr(j)))*dt
          if (HR(j) < 0 ) then
            print*, 'Negative HR: ', HR(j), t
          end if
          
          sum_input_step=sum_input_step+(C_PlantLITm+C_PlantLITs+e_m*C_PlantEcM+e_m*C_PlantAM+C_PlantSOMc+C_PlantSOMp+C_PlantSOMa)*dt*delta_z(j) !g/m2
          sum_N_input_step=sum_N_input_step+(N_PlantLITm+N_PlantLITs+Deposition)*dt*delta_z(j) !g/m2
          sum_N_out_step=sum_N_out_step+(N_EcMPlant+N_INPlant+Leaching)*dt*delta_z(j)

        end do !j, depth_level

        !Store accumulated HR mass
        call respired_mass(HR, HR_mass,nlevdecomp)
        HR_mass_accumulated = HR_mass_accumulated + HR_mass
        if (isVertical) then
          call vertical_diffusion(tot_diffC,upperC,lowerC, pool_temporaryC,vertC,D_carbon)
          call vertical_diffusion(tot_diffN,upperN,lowerN, pool_temporaryN,vertN,D_nitrogen)
          pool_matrixC =  vertC*dt + pool_temporaryC
          pool_matrixN = vertN*dt + pool_temporaryN
          !call disp("vertical change C : ", vertC)
          !call disp("vertical change N : ", vertN)
        else
          pool_matrixC=pool_temporaryC
          pool_matrixN=pool_temporaryN
        end if!isVertical


        
        !START Compute and write annual means
        sum_consN = sum_consN + pool_matrixN
        sum_consC = sum_consC + pool_matrixC
        !print*, ycounter
        if (ycounter == 365*24*step_frac) then
          call check(nf90_close(ncid)) !Close netcdf file containing values for the past year
          ! call disp("C matrix : ",pool_matrixC)
          ! call disp("N matrix : ",pool_matrixN)
          ! call disp("C:N : ",pool_matrixC/pool_matrixN(:,1:10))
          write_y =write_y+1
          !call annual_mean(sum_consC,sum_consN, nlevdecomp,write_y , run_name) !calculates the annual mean and write the result to file
          if (year == stop_year) then
            year = start_year         
          else 
            year = year + 1             
          end if          
          day_counter=0
          current_day=1
          month_counter=0
          current_month=1
          write (year_char,year_fmt) year
          call check(nf90_open(trim(adjustr(clm_input_path)//'all.'//year_char//'.nc'), nf90_nowrite, ncid)) !open netcdf containing values for the next year
          call read_time(adjustr(clm_input_path)//'all.'//year_char//'.nc',input_steps)     
          ycounter = 0
          sum_consN =0
          sum_consC =0

        end if
        !END Compute and write annual means

        if (counter == write_hour*step_frac) then
          counter = 0
          call fill_netcdf(run_name, int(time), pool_matrixC, change_matrixC, pool_matrixN,change_matrixN,&
           date, HR_mass_accumulated,HR,vertC,vertN, write_hour,current_month, &
           TSOIL, r_moist,CUE_bacteria_vr,CUE_fungi_vr,nlevdecomp)
           change_sum = 0.0
          !vertC_change_sum = 0.0
        end if!writing

        !Write end values to terminal
        if (t == nsteps) then
          call disp("pool_matrixC gC/m3 ",pool_matrixC)
          call disp("pool_matrixN gN/m3 ",pool_matrixN)          
          call disp("C:N: ",pool_matrixC/pool_matrixN(:,1:10))
          pool_C_final=pool_matrixC
          pool_N_final=pool_matrixN    
          call store_parameters(run_name)    
        end if
        
        call test_mass_conservation(sum_input_step,HR_mass,pool_matrixC_previous,pool_matrixC,nlevdecomp,pool_types)
        pool_matrixC_previous=pool_matrixC
        sum_input_total=sum_input_total+sum_input_step
        sum_input_step=0.0
        sum_N_input_total=sum_N_input_total+sum_N_input_step  
        sum_N_input_step=0.0
        sum_N_out_total=sum_N_out_total+sum_N_out_step    
        sum_N_out_step=0.0
      end do !t
      
      call total_mass_conservation(sum_input_total,HR_mass_accumulated,pool_C_start,pool_C_final,nlevdecomp,pool_types)
      call total_nitrogen_conservation(sum_N_input_total,sum_N_out_total,pool_N_start,pool_N_final,nlevdecomp,pool_types_N)
      deallocate(ndep_prof,leaf_prof,froot_prof)
      deallocate(CUE_bacteria_vr,CUE_fungi_vr)
    end subroutine decomp

  subroutine annual_mean(yearly_sumC,yearly_sumN,nlevels, year, run_name)
    REAL(r8), DIMENSION(nlevels,pool_types)  , intent(in):: yearly_sumC
    REAL(r8), DIMENSION(nlevels,pool_types+1), intent(in):: yearly_sumN
    integer, intent(in) :: year
    integer , intent(in):: nlevels
    CHARACTER (len = *), intent(in):: run_name
    !Local
    REAL(r8), DIMENSION(nlevels,pool_types) :: yearly_meanC
    REAL(r8), DIMENSION(nlevels,pool_types+1) :: yearly_meanN
    integer, parameter                         :: hr_in_year = 24*365

    yearly_meanC=yearly_sumC/hr_in_year
    yearly_meanN=yearly_sumN/hr_in_year

    !call fill_yearly_netcdf(run_name, year, yearly_meanC,yearly_meanN,nlevels)

  end subroutine annual_mean

end module mycmim
