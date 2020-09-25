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
  use fluxMod2
  use initMod
  use writeMod
  use testMod
  implicit none

  contains
    subroutine decomp(nsteps, run_name, isVertical, ecosystem,step_frac) !Calculates the balance equations dC/dt and dN/dt for each pool at each time step based on the fluxes calculated in the same time step. Then update the pool sizes before moving on
                                                                        !
      integer                        :: nsteps                               ! number of time steps to iterate over
      character (len=*)              :: run_name                             ! used for naming outputfiles
      logical                        :: isVertical                           ! True if we use vertical soil layers.
      character (len=*)              :: ecosystem                            ! 'Shrub', 'Heath' or 'Meadow' (For comparison with Sorensen et al)
      integer                        :: step_frac                            ! determines the size of the time step

      real(r8),dimension(nlevdecomp) :: HR                                   ! For storing the C  that is lost to respiration [gC/m3h]
      real(r8), dimension(nlevdecomp):: mass_HR                              ![gC/m2h]
      real(r8)                       :: HR_mass_accumulated, HR_mass, input_mass_accumulated

      real(r8)                       :: pool_matrixC(nlevdecomp,pool_types)     ! For storing C pool sizes [gC/m3]
      real(r8)                       :: mass_pool_matrixC(nlevdecomp,pool_types)! For storing C pool masses [gC/m2] (pool_matrixC*delta_z)
      !real(r8)                       :: previous_conc(nlevdecomp,pool_types)    ! stores concentration at timestep. Used for comparing in subroutine test_mass_conservation
      real(r8)                       :: change_matrixC(nlevdecomp,pool_types)   ! For storing dC/dt for each time step [gC/(m3*hour)]
      real(r8)                       :: a_matrixC(nlevdecomp, pool_types)        ! For storing the analytical solution
      real(r8)                       :: InitC(nlevdecomp, pool_types)           ! Initial C concentration, determined in initMod.f90
      real(r8)                       :: pool_temporaryC(nlevdecomp,pool_types)  ! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
                                                                                !The new value after the time step is then pool_matrixC = pool_temporaryC + vertical change
      real(r8)                       :: pool_matrixN(nlevdecomp,pool_types+1)   ! For storing N pool sizes [gN/m3] parallell to C pools and  inorganic N
      real(r8)                       :: mass_pool_matrixN(nlevdecomp,pool_types+1)! For storing N pool masses [gN/m2] (pool_matrixN*delta_z)
      real(r8)                       :: change_matrixN(nlevdecomp,pool_types+1) ! For storing dC/dt for each time step [gN/(m3*hour)]
      real(r8)                       :: a_matrixN(nlevdecomp, pool_types+1)       ! For storing the analytical solution
      real(r8)                       :: InitN(nlevdecomp, pool_types+1)         ! Initial N concentration, determined in initMod.f90
      !real(r8)                       :: pool_temporaryN(nlevdecomp,pool_types+1)

                                                                             !Shape of pool_matrixC/change_matrixC
                                                                             !|       LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc |
                                                                             !|level1   1   2    3    4   5   6   7   8    9    10  |
                                                                             !|level2                                               |
                                                                             !| .                                                   |
                                                                             !| .                                                   |
                                                                             !|nlevdecomp __________________________________________|

      real(r8)                       :: dt                                   ! size of time step
      real(r8)                       :: time                                 ! t*dt
      real(r8)                       :: Loss, Gain, Loss_term                ! Source and sink term for solving the analytical solution to the dC/dt=Gain-Loss*concentration equation.
      real(r8)                       :: tot_diff,upper,lower                 ! For the call to vertical_diffusion
      real(r8)                       :: vertC(nlevdecomp, pool_types)         !Stores the vertical change in a time step, on the same form as change_matrixC
      real(r8)                       :: vertN(nlevdecomp, pool_types)         !Stores the vertical change in a time step, on the same form as change_matrixN
!NOTE: these are not used after the plant pool was introduced:
!      real(r8)                       :: lit_inputC(nlevdecomp,no_of_litter_pools)![gC/(m3 h)] Fraction of litter input to LITm and LITs, respectively
!      real(r8)                       :: myc_inputC(nlevdecomp,no_of_myc_pools)   ![gC/(m3 h)] vector giving the input from vegetation to mycorrhiza pools
!      real(r8)                       :: som_inputC(nlevdecomp,no_of_som_pools-1) ![gC/(m3 h)] !only input to SOMp and SOMc

      integer                        :: counter                                 !used for determining when to output results
      integer                        :: month_counter
      integer                        :: j,i,t                                   !for iterations
      integer,parameter              ::t_init=1
      real(r8), dimension(nlevdecomp):: HR_sum                                  !Sums total HR between two output entries
      real(r8)                       :: change_sum(nlevdecomp, pool_types)
      real(r8)                       :: vertN_change_sum(nlevdecomp, pool_types)
      real(r8)                       :: vertC_change_sum(nlevdecomp, pool_types)
      integer,parameter              :: write_hour= 24*4
!      real(r8)                       :: ecm_frac, erm_frac, am_frac
      real(r8), dimension(nlevdecomp, 7) :: tot_input

      !Assigning values: (Had to move from paramMod to here to be able to modify them during a run)
      MGE       = (/0.3,0.3,0.3,0.4,0.4,0.4/)
!      lit_inputC = 0.0
!      myc_inputC = 0.0
!      som_inputC = 0.0
      ! Fracions of SAP that goes to different SOM pools
      fPHYS = (/ 0.3 * exp(fCLAY), 0.2 * exp(0.8*fCLAY) /)
      fCHEM =  (/0.1 * exp(-3.0*fMET), 0.3 * exp(-3*fMET) /)
      fAVAIL = 1-(fPHYS+fCHEM)

      !TODO FROM MIMICS_CYCLE_CN:
      ! WW also modify TAU as a function of soil moisture, so things don't
      ! colapse in frozen soils...
      !mimicsbiome%tauR(npt) = mimicsbiome%tauR(npt) * fW
      !mimicsbiome%tauK(npt) = mimicsbiome%tauK(npt) * fW
      tau = (/ 5e-4*exp(0.3*fMET), 5e-4*exp(0.1*fMET)/)!*0.3![1/h] Microbial turnover rate (SAP to SOM), SAPr,(/1.39E-3*exp(0.3*fMET), 2.3E-4*exp(0.1*fMET)/)

      dt= 1.0/step_frac !Setting the time step

      ! if (ecosystem == 'Heath') then
      !    GEP      = 0.281
      !    fMET     = 0.4 !Fraction of input to litter that goes to the metabolic litter pool
      !    ecm_frac = 1
      !    erm_frac = 2
      !    am_frac = 0
      ! elseif (ecosystem == 'Meadow') then
      !   GEP       = 0.385
      !   fMET      = 0.5 !Fraction of input to litter that goes to the metabolic litter pool
      !   ecm_frac = 1
      !   erm_frac = 1
      !   am_frac = 0
      ! elseif (ecosystem == 'Shrub') then
      !   GEP       = 0.491
      !   fMET      = 0.30 !Fraction of input to litter that goes to the metabolic litter pool
      !   ecm_frac = 2
      !   erm_frac = 1
      !   am_frac = 0
      ! else
      !   print*, 'Invalid ecosystem name', ecosystem
      !   stop
      ! end if
      ! print*, ecosystem

      !Set initial concentration values in pool_matrixC:
      !Set initial concentration values:
      if (isVertical) then
        call initialize_vert(InitC, InitN, P_init_C, P_init_N)
      else
        call initialize_onelayer(InitC, InitN, P_init_C, P_init_N)
      end if !isVertical

      pool_matrixC = InitC
      pool_matrixN = InitN
      CPlant = P_init_C
      NPlant = P_init_N
      a_NPlant = P_init_N
      a_CPlant = P_init_C
      change_matrixC = 0.0
      change_matrixN = 0.0
      a_matrixC      = InitC
      a_matrixN      = InitN
      HR            = 0.0
      vertC_change_sum=0.0
      vertN_change_sum=0.0

      counter = 0
      ycounter = 0
      year     = 1
      HR_sum   = 0.0 !For summing up the total respiration between two output times
      current_month = 1
      month_counter = 30


      !open and prepare files to store results. Store initial value
      call create_netcdf(run_name)
      call fill_netcdf(run_name, nlevdecomp,t_init, pool_matrixN, change_matrixN, a_matrix, HR, vertN_change_sum, write_hour,current_month)
      call fill_netcdf(run_name,t_init, pool_matrixC, change_matrixC, pool_matrixN,change_matrixN, a_matrix, HR, vertC_change_sum, write_hour,current_month)
     !stop
      counter = 0
      ycounter = 0
      year     = 1
      HR_sum   = 0.0 !For summing up the total respiration between two output times
      current_month = 1
      month_counter = 30

      call read_clmdata(clm_data_file,TSOIL,SOILLIQ,SOILICE,WATSAT,current_month)
      call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist)

      ! if (current_month == previous_month + 1) then
      !   !update temp and moisture variables
      ! end if
      !Inputs from aboveground litter to LIT and protected SOM pools: GEP/depth [gC/m3h]NOTE: No longer in use!
      ! lit_inputC(1,:) = (/fMET, (0.9-fMET)/)*GEP/(delta_z(1) + delta_z(2))!Partition between litter pools. The last 0.1 fraction goes directly to SOM
      ! som_inputC(1,:) = (/0.05, 0.05/)*GEP/(delta_z(1) + delta_z(2))
      ! myc_inputC(4,:) = (/ecm_frac, erm_frac, am_frac/)*(GEP)/(delta_z(4))
      ! myc_inputC(5,:) = (/ecm_frac, erm_frac,am_frac/)*(GEP)/(delta_z(5))

      !----------------------------------------------------------------------------------------------------------------
      do t =1,nsteps !Starting time iterations
        time = t*dt
        counter =counter + 1
        ycounter = ycounter + 1
        month_counter = month_counter + 1
        previous_conc = pool_matrixC !Store value of previous timestep. Used in subroutine test_mass_conservation

        !Update temp and moisture values monthly
        if (month_counter == days_in_month(current_month)*24) then
          previous_month = current_month
          call read_clmdata(clm_data_file,TSOIL,SOILLIQ,SOILICE,WATSAT,current_month)
          call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist)
          if (current_month == 12) then
            current_month = 1
          else
            current_month = current_month + 1
          end if
          month_counter = 0
        end if

        if (ycounter == 365*24) then
          year = year + 1
          !print*, year
          ycounter = 0
        end if
        if (t == 1) then
          call disp("InitC", pool_matrixC)
          call disp("InitN", pool_matrixN)
        end if

        !If-test used to modify something after half of the total run time
         ! if (t == nsteps/2) then
         !   ! print*, tau
         !   ! print*, 'CHANGED TAU'
         !   ! tau = (/ 5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET)/)*1.5 ![1/h] Microbial turnover rate (SAP to SOM), SAPr, SAPk
         !   ! print*, tau
         !   ! print*, MGE
         !   print*, 'Changed MGE!'
         !   MGE     = MGE*1.2
         !   print*, MGE
         !   ! print*, Km,Vmax
         !   ! Kmod    = Kmod/15!LITm, LITs, SOMa entering SAPr, LITm, LITs, SOMa entering sapk
         !   ! Vmod    = Vmod/15    !LITm, LITs, SOMa entering SAPr, LITm, LITs, SOMa entering sapk
         !   ! Km      = exp(Kslope*tsoi+3.19)*10*Kmod  ![mgC/cm3]*10e3=[gC/m3]
         !   ! Vmax    = exp(0.063*tsoi + 5.47)*8e-6*Vmod ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?
         !   ! print*, Km,Vmax
         ! end if !Change

        ! if (time >= 91*24 + (year - 1)*365*24 .and. time <= 304*24 + (year -1)*365*24  ) then
        !   tsoi = 10.0
        ! else
        !   tsoi = 0.0
        ! end if

        ! if (time >= 244*24 + (year-1)*365*24 .and. time <= 335*24 + (year-1)*365*24) then
        !   !print*, time, 335*24 + year*365*24 - 244*24 - year*365*24
        !   lit_inputC(1,:) = (/fMET, (0.9-fMET)/)*GEP/delta_z(1) !Partition between litter pools. The last 0.1 fraction goes directly to SOM
        !   som_inputC(1,:) = (/0.05, 0.05/)*GEP/delta_z(1)
        !   !print*, "TIME: ", time
        !   !print*, lit_inputC(1,:)
        ! else
        !   lit_inputC = 0
        !   som_inputC =0
        ! end if
          !Calculate fluxes between pools in level j (file: fluxMod2.f90):
          call calculate_fluxes(j, pool_matrixC, pool_matrixN, CPlant, NPlant)

          !calculate the change of N and C in the plant based on the flux equations.  TODO: Needs better way to ensure reasonable values in these pools
          !CPlant_tstep and NPlant_tstep sum up the change from each layer, and will be used to update CPlant and NPlant at the end of the timestep.
          !(This will only be one value if isVertical = False)
          possible_N_change = (N_EcMPlant+  N_ErMPlant +  N_AMPlant + N_InPlant  &
          - N_PlantLITm - N_PlantLITs)*dt*delta_z(j) !gN/m2h

          Plant_lossN = (N_PlantLITm + N_PlantLITs)*delta_z(j)
          Plant_gainN = (N_EcMPlant+  N_ErMPlant +  N_AMPlant + N_InPlant)*delta_z(j)
          !print*, "PLANTLOSS", N_PlantLITm,N_PlantLITs
          !print*, "PLANTGAIN", N_EcMPlant,N_ErMPlant, N_AMPlant, N_InPlant
          possible_C_change = (C_growth_rate - C_PlantLITm - C_PlantLITs - &
           C_PlantEcM - C_PlantErM - C_PlantAM)*dt*delta_z(j)!gC/m2h

          Plant_lossC = (C_PlantLITm + C_PlantLITs + C_PlantEcM + C_PlantErM + C_PlantAM)*delta_z(j)
          Plant_gainC = C_growth_rate*delta_z(j) !TODO:: Skal ikke ganges med delta_z her??
          growth_rate_sum = growth_rate_sum + C_growth_rate
          ! if (possible_C_change + Cplant < CPlant_min) then
          !   CPlant_tstep = CPlant_tstep + (C_growth_rate - C_PlantLITm - N_PlantLITs)*dt
          !   !Instead of the values calculated in subroutine calculate_fluxes, they are set to zero:
          !   C_PlantEcM = 0
          !   C_PlantErM = 0
          !   C_PlantAM = 0
          ! else
          ! end if
          ! if (possible_N_change + NPlant < NPlant_min) then
          !   print*,"possible_N_change + NPlant < NPlant_min", NPlant + possible_N_change
          !   !TODO: How to ensure this will not happen?
          ! end if
          NPlant_tstep = NPlant_tstep + possible_N_change
          CPlant_tstep = CPlant_tstep + possible_C_change
          !print*, NPlant_tstep, CPlant_tstep

        do j = 1, nlevdecomp !For each depth level (for the no vertical transport case, nlevdecomp = 1, so loop is only done once):
          Km      = exp(Kslope*TSOIL(j) + Kint)*a_k*Kmod                    ![mgC/cm3]*10e3=[gC/m3]
          Vmax    = exp(Vslope*TSOIL(j) + Vint)*a_v*Vmod*r_moist(j)    ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?
          k_mycsom  = (/1.4,1.4,1.4/)*10e-3*r_moist(j)   ![1/h] Decay constants, myc som

          !Calculate fluxes between pools in level j:
          !call microbial_fluxes(j, pool_matrixC,nlevdecomp)
          call som_fluxes(j, pool_matrixC)
          call litter_fluxes(j, pool_matrixC)

          do i = 1,pool_types + 1 !loop over all the pool types, i, in depth level j
            !This if-loop calculates dC/dt and dN/dt for the different carbon pools.NOTE: If pools are added/removed (i.e the actual model equations is changed), this loop needs to be updated.
            !NOTE: The "change_matrixC" values correspond to the equations A11-A17 in Wieder 2015
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
            !TODO Sap: Check if MGE/efficiency still makes sense in the new setup.
            !  C_Gain = C_LITmSAPb*MGE(1) + C_LITsSAPb*MGE(2) &
            !  + C_SOMaSAPb*MGE(3) + e_s*(Decomp_ecm + Decomp_erm + Decomp_am)
              C_Gain = e_s*(C_LITmSAPb + C_LITsSAPb &
                + C_SOMaSAPb + 0.5*(Decomp_ecm + Decomp_erm + Decomp_am))
              C_Loss =  C_SAPbSOMp + C_SAPbSOMa + C_SAPbSOMc
              !N_Gain = N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb
              !N_Loss = N_SAPbSOMp + N_SAPbSOMa + N_SAPbSOMp + N_SAPbIN
              N_Gain = e_s*U_sb/CN_ratio(3)
              N_Loss = N_SAPbSOMp + N_SAPbSOMa + N_SAPbSOMc

            elseif (i==4) then !SAPf
              !C_Gain = C_LITmSAPf*MGE(1) + C_LITsSAPf*MGE(2) &
              !+ C_SOMaSAPf*MGE(3)
              C_Gain = e_s*(C_LITmSAPf + C_LITsSAPf &
                + C_SOMaSAPf + 0.5*(Decomp_ecm + Decomp_erm + Decomp_am))
              C_Loss =  C_SAPfSOMp + C_SAPfSOMa + C_SAPfSOMc
            !  N_Gain = N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf
            !  N_Loss = N_SAPfSOMp + N_SAPfSOMa + N_SAPfSOMp + N_SAPfIN
              N_Gain = e_s*U_sf/CN_ratio(4)
              N_Loss = N_SAPfSOMp + N_SAPfSOMa + N_SAPfSOMc
            elseif (i==5) then !EcM
              C_Gain = e_m*C_PlantEcM
              C_Loss = C_EcMSOMp + C_EcMSOMa + C_EcMSOMc
              !N_Gain = N_INEcM + N_SOMaEcM
              !N_Loss = N_EcMPlant + N_EcMSOMa + N_EcMSOMp + N_EcMSOMc
              N_Gain = e_m*C_PlantEcM/CN_ratio(5)
              N_Loss = N_EcMSOMp + N_EcMSOMa + N_EcMSOMc
            elseif (i==6) then !ErM
              C_Gain = e_m*C_PlantErM
              C_Loss = C_ErMSOMp + C_ErMSOMa + C_ErMSOMc
              N_Gain = e_m*C_PlantErM/CN_ratio(6)
              N_Loss = N_ErMSOMp + N_ErMSOMa + N_ErMSOMc
              !N_Gain = N_INErM + N_SOMaErM
              !N_Loss = N_ErMPlant + N_ErMSOMa + N_ErMSOMp + N_ErMSOMc
            elseif (i==7) then !AM
              C_Gain = e_m*C_PlantAM
              C_Loss = C_AMSOMp + C_AMSOMa + C_AMSOMc
              N_Gain = e_m*C_PlantAM/CN_ratio(7)
              N_Loss = N_AMSOMp + N_AMSOMa + N_AMSOMc
            !  N_Gain = N_INAM + N_SOMaAM
            !  N_Loss = N_AMPlant + N_AMSOMa + N_AMSOMp + N_AMSOMc
              !print*, N_Loss, "am"
            elseif (i==8) then !SOMp
              C_Gain =  C_SAPbSOMp + C_SAPfSOMp + C_EcMSOMp + C_ErMSOMp + C_AMSOMp
              C_Loss = C_SOMpSOMa
              N_Gain =  N_SAPbSOMp + N_SAPfSOMp + N_EcMSOMp + N_ErMSOMp + N_AMSOMp
              N_Loss = N_SOMpSOMa
            elseif (i==9) then !SOMa
               C_Gain = C_SAPbSOMa + C_SAPfSOMa + C_EcMSOMa + C_ErMSOMa + C_AMSOMa + C_SOMpSOMa + C_SOMcSOMa
               C_Loss = C_SOMaSAPb + C_SOMaSAPf
               N_Gain = N_SAPbSOMa + N_SAPfSOMa + N_EcMSOMa + N_ErMSOMa + N_AMSOMa + N_SOMpSOMa + N_SOMcSOMa
               N_Loss = N_SOMaSAPb + N_SOMaSAPf + N_SOMaEcM + N_SOMaErM + N_SOMaAM
            elseif (i==10) then !SOMc
              C_Gain =  C_SAPbSOMc + C_SAPfSOMc + C_EcMSOMc + C_ErMSOMc + C_AMSOMc
              C_Loss = C_SOMcSOMa
              N_Gain =  N_SAPbSOMc + N_SAPfSOMc + N_EcMSOMc + N_ErMSOMc + N_AMSOMc
              N_Loss = N_SOMcSOMa
            elseif (i == 11) then !Inorganic N
                 N_Gain = N_SAPbIN + N_SAPfIN + Deposition/delta_z(j)
                 N_Loss = N_INEcM + N_INErM + N_INAM + N_InPlant + Leaching/delta_z(j)
                 change_matrixN(j,i) = N_Gain - N_loss
            else
              print*, 'Too many pool types expected, pool_types = ',pool_types, 'i: ', i
            end if !determine dC_i/dt

            if (i /= 11) then !Carbon matrix does only have 10 columns (if i == 11 this is handeled inside the loop over pools)
                change_matrixC(j,i) = C_Gain - C_Loss
                change_matrixN(j,i) = N_Gain - N_loss
                !For summarize total change between each written output
                change_sum(j,i)= change_sum(j,i) + change_matrixC(j,i)*dt
                !Store these values as temporary so that they can be used in the diffusion subroutine
                pool_temporaryC(j,i)=pool_matrixC(j,i) + change_matrixC(j,i)*dt
            end if


                !print*, "CN_ratio: ", variables(i)
                !print*, CN_ratio(i)
                !print*, a_matrixC(j,i)/a_matrixN(j,i)
                a_matrixC(j,i) = a_matrixC(j,i)*exp(-dt*Loss_termC) + (C_Gain/Loss_termC)*(1-exp(-dt*Loss_termC))
                !print*, loss_termC, C_Gain, i, pool_temporaryC(j,i)
            end if
            ! print*, "FROM IN", N_INEcM, N_INErM, N_INAM
            ! print*, "FROM SOMa" ,N_SOMaEcM, N_SOMaErM, N_SOMaAM
            ! print*, "TO plant: ",N_EcMPlant, N_ErMPlant, N_AMPlant
            ! print*, "TO SOM: ",  N_EcMSOMa, N_ErMSOMa, N_AMSOMa
            ! print*, N_EcMSOMp, N_ErMSOMp, N_AMSOMp
            ! print*, N_EcMSOMc, N_ErMSOMc, N_ErMSOMc
            Loss_termN = N_Loss/a_matrixN(j,i)
            !print*, Loss_termN, N_variables(i), a_matrixN(j,i)
            a_matrixN(j,i)= &
            a_matrixN(j,i)*exp(-dt*Loss_termN) + (N_Gain/Loss_termN)*(1-exp(-dt*Loss_termN))!
            ! if (i==11)then
            !   print*, a_matrixN(j,i)*exp(-dt*Loss_termN),(N_Gain/Loss_termN)*(1-exp(-dt*Loss_termN))
            ! end if

            !print*, "*Levetid nitrogen i ", N_variables(i),  ': ', (1/Loss_termN)/24, 'dager. '
            !print*,"N loss_term (k)",  N_Loss, a_matrixN(j,i), Loss_termN
            !print*, "------------------------------------------------------------------------------"


            pool_matrixN(j,i) =pool_matrixN(j,i) + change_matrixN(j,i)*dt      !NOTE No vertical transport of nitrogen (yet..)
            !control check
            ! if (pool_temporaryC(j,i) < 0.0) then
            ! !  print*, 'Negative concentration value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryC(j,i)
            ! !  print*, 'Value changed to: ', min_pool_value
            ! !  print*, 'Month, year: ', current_month, year
            !   !call disp('Temp: ',TSOIL)
            !   !call disp('moist: ', r_moist)
            !   !call disp('Vmax: ', Vmax)
            !   !call disp('Km: ', Km)
            !   !call disp(pool_temporaryC)
            ! !  pool_temporaryC(j,i) = min_pool_value
            !   !STOP
            ! end if

            !control check
            if (pool_temporaryC(j,i) < 0.0) then
              print*, 'Negative concentration value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryC(j,i)
              print*, 'Value changed to: ', min_pool_value
              print*, 'Month, year: ', current_month, year
              call disp('Temp: ',TSOIL)
              call disp('moist: ', r_moist)
              call disp('Vmax: ', Vmax)
              call disp('Km: ', Km)
              !call disp(pool_temporaryC)
              pool_temporaryC(j,i) = min_pool_value
              STOP
            end if

            !To calculate analytical solution:
            Loss_term = Loss/pool_matrixC(j,i)
            a_matrix(j,i) = InitC(j,i)*exp(-time*Loss_term) + Gain*(1-exp(-time*Loss_term))/Loss_term !+ vertC*dt

          end do !i, pool_types

          !Calculate the heterotrophic respiration loss from depth level j in timestep t: NOTE: revise!
          !HR(j) =( C_LITmSAPb*(1-MGE(1)) + C_LITsSAPb*(1-MGE(2))  + C_SOMaSAPb*(1-MGE(3)) + C_LITmSAPf*(1-MGE(4)) &
          !+ C_LITsSAPf*(1-MGE(5)) + C_SOMaSAPf*(1-MGE(6)) + (C_PlantEcM + C_PlantErM + C_PlantAM)*(1-e_m))*dt
          HR(j) =(( C_LITmSAPb + C_LITsSAPb  + C_SOMaSAPb + C_LITmSAPf &
          + C_LITsSAPf + C_SOMaSAPf+Decomp_ecm + Decomp_erm + Decomp_am)*(1-e_s) &
          + (C_PlantEcM + C_PlantErM + C_PlantAM)*(1-e_m))*dt
          if (HR(j) < 0 ) then
            print*, 'Negative HR: ', HR(j), t
          end if
          HR_sum(j) = HR_sum(j) + HR(j)
        end do !j, depth_level

        !Update Plant pools with the total change from all the layers
        CPlant = CPlant + CPlant_tstep!Numerial
        NPlant = NPlant + NPlant_tstep

        !"analytic":
        Loss_termNP = Plant_lossN/a_NPlant
        Loss_termCP = Plant_lossC/a_CPlant
      !  print*, "time", time, "Loss_termCP", Loss_termCP, "Loss_termNP", Loss_termNP
        a_NPlant = a_Nplant*exp(-dt*Loss_termNP) + Plant_GainN*(1-exp(-dt*Loss_termNP))/Loss_termNP!
        a_CPlant = a_Cplant*exp(-dt*Loss_termCP) + Plant_GainC*(1-exp(-dt*Loss_termCP))/Loss_termCP!

        !print*, '*Levetid nitrogen i plant pool: ', (1/Loss_termNP)/24, 'dager'
      !  print*, '*Levetid karbon i plant pool: ', (1/Loss_termCP)/24, 'dager'


        Plant_CN = CPlant/NPlant
        ! if (Plant_CN < opt_Plant_CN + 5 .and. Plant_CN > opt_Plant_CN - 5 ) then
        !   print*, "innafor"
        ! else
        !
        ! !  print*, "Diff: ", Plant_CN - opt_Plant_CN
        ! !  print*, "C: ", CPlant
        ! !  Print*, "N: ", NPlant
        !   NPlant = CPlant/opt_Plant_CN
        !
        ! end if
        !Store accumulated HR mass
        call respired_mass(HR, HR_mass)
        HR_mass_accumulated = HR_mass_accumulated + HR_mass


        if (isVertical) then
          call vertical_diffusion(tot_diff,upper,lower, pool_temporaryC,vertC,time, counter,dt)
          pool_matrixC =  vertC*dt + pool_temporaryC
          vertC_change_sum = vertC_change_sum + vertC*dt
          !call disp("pool_matrixC",pool_matrixC)
          a_matrix=a_matrix+vertC*dt
        else
          pool_matrixC=pool_temporaryC
        end if!isVertical

        if (counter == write_hour) then

          counter = 0
          !call disp("pool_matrixC ",pool_matrixC)
          call fill_netcdf(run_name, t, pool_matrixC, change_sum, pool_matrixN,change_matrixN, a_matrix,HR, vertC_change_sum, write_hour,current_month)
          HR_sum = 0.0
          change_sum = 0.0
          vertC_change_sum = 0.0
        !print*, HR_sum
        end if!writing

        if (t == nsteps) then

          do i = 1,pool_types
            mass_pool_matrixC(:,i) = pool_matrixC(:,i)*delta_z
          end do
          call disp("pool_matrixC gC/m2 ",mass_pool_matrixC)
          call disp("pool_matrixC gC/m3 ",pool_matrixC)

          !call disp("respiration", HR)
        end if

      end do !t
      !call store_parameters(run_name)

      call tot_mass_input(lit_inputC, som_inputC,myc_inputC, tot_input)!TODO: Need to multiply with dt! (current dt =1, so ok for now..09.06.20)
      input_mass_accumulated = sum(tot_input)*nsteps
      call total_mass_conservation(sum_input= input_mass_accumulated, sum_respiration = HR_mass_accumulated, old = InitC, new = pool_matrixC)
    end subroutine decomp
end module mycmim
