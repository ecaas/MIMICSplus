!This module models the decomposition of carbon through the soil. The number of soil depth levels and different pools are given in the file paramMod.f90. The paramMod module contains
!all relevant parameters for the flux and balance equations used in this module. In the current setup, 4 depth levels are used, and the following 10 carbon pools/reservoirs are found in each level:
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

!TODO: Incorporate frozen soil behavior/variable T, Variable clay content(?)
module mycmim
  use paramMod
  use dispmodule !External module to pretty print matrices (mainly for testing purposes)
  use fluxMod
  use initMod
  use writeMod
  implicit none

  contains
    subroutine decomp(nsteps, run_name, isVertical, ecosystem,step_frac) !Calculates the balance equation dC/dt for each pool at each time step based on the fluxes calculated in the same time step. Then update the pool sizes before moving on
                                                                                     !It also calculates an analytical solution to the problem.
      integer                        :: nsteps                               ! number of time steps to iterate over
      character (len=*)              :: run_name                             ! used for naming outputfiles
      logical                        :: isVertical                           ! True if we use vertical soil layers.
      character (len=*)              :: ecosystem                            ! 'Shrub', 'Heath' or 'Meadow' (For comparison with Sorensen et al)
      integer                        :: step_frac                            ! determines the size of the time step

      real(r8),dimension(nlevdecomp) :: HR                                   ! For storing the C amount that is lost to respiration
      real(r8)                       :: pool_matrixC(nlevdecomp,pool_types)   ! For storing C pool sizes [gC/m3]
      real(r8)                       :: mass_pool_matrixC(nlevdecomp,pool_types)   ! For storing C pool sizes [gC/m2] (pool_matrixC*delta_z)
      real(r8)                       :: change_matrixC(nlevdecomp,pool_types) ! For storing dC/dt for each time step [gC/(m3*hour)]
      real(r8)                       :: a_matrix(nlevdecomp, pool_types)     ! For storing the analytical solution
      real(r8)                       :: InitC(nlevdecomp, pool_types)         ! Initial C concentration, determined in initMod.f90
      real(r8)                       :: pool_temporaryC(nlevdecomp,pool_types)! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
                                                                             !The new value after the time step is then pool_matrixC = pool_temporaryC + vertical change
      real(r8)                       :: pool_matrixN(nlevdecomp,pool_types)   ! For storing C pool sizes [gC/m3]
      real(r8)                       :: mass_pool_matrixN(nlevdecomp,pool_types)   ! For storing C pool sizes [gC/m2] (pool_matrixC*delta_z)
      real(r8)                       :: change_matrixN(nlevdecomp,pool_types) ! For storing dC/dt for each time step [gC/(m3*hour)]
      !real(r8)                       :: a_matrix(nlevdecomp, pool_types)     ! For storing the analytical solution
      real(r8)                       :: InitN(nlevdecomp, pool_types)         ! Initial C concentration, determined in initMod.f90
      real(r8)                       :: pool_temporaryN(nlevdecomp,pool_types)! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
                                                                                                                                                    !The new value after the time step is then pool_matrixC = pool_temporaryC + vertical change.

                                                                             !Shape of the pool_matrixC/change_matrixC
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

      real(r8)                       :: lit_inputC(nlevdecomp,no_of_litter_pools)![gC/(m3 h)] Fraction of litter input to LITm and LITs, respectively
      real(r8)                       :: myc_inputC(nlevdecomp,no_of_myc_pools)   ![gC/(m3 h)] vector giving the input from vegetation to mycorrhiza pools
      real(r8)                       :: som_inputC(nlevdecomp,no_of_som_pools-1) ![gC/(m3 h)] !only input to SOMp and SOMc
      integer                        :: ycounter, year,ncid,varid
      integer                        :: counter                                 !used for determining when to output results
      integer                        :: month_counter
      integer                        :: j,i,t                                   !for iterations
      integer,parameter              ::t_init=1
      real(r8), dimension(nlevdecomp):: HR_sum                                  !Sums total HR between two output entries
      real(r8)                       :: change_sum(nlevdecomp, pool_types)
      real(r8)                       :: vertN_change_sum(nlevdecomp, pool_types)
      real(r8)                       :: vertC_change_sum(nlevdecomp, pool_types)
      integer,parameter              :: write_hour= 24*4
      real(r8)                       :: ecm_frac, erm_frac, am_frac
      !Assigning values: (Had to move from paramMod to here to be able to modify them during a run)
      MGE       = (/0.3,0.3,0.3,0.4,0.4,0.4/)
      lit_inputC = 0.0
      myc_inputC = 0.0
      som_inputC = 0.0
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

      if (ecosystem == 'Heath') then
         GEP      = 0.281
         fMET     = 0.4 !Fraction of input to litter that goes to the metabolic litter pool
         ecm_frac = 1
         erm_frac = 2
         am_frac = 0
      elseif (ecosystem == 'Meadow') then
        GEP       = 0.385
        fMET      = 0.5 !Fraction of input to litter that goes to the metabolic litter pool
        ecm_frac = 1
        erm_frac = 1
        am_frac = 1
      elseif (ecosystem == 'Shrub') then
        GEP       = 0.491
        fMET      = 0.30 !Fraction of input to litter that goes to the metabolic litter pool
        ecm_frac = 2
        erm_frac = 1
        am_frac = 0
      else
        print*, 'Invalid ecosystem name', ecosystem
        stop
      end if
      print*, ecosystem

      !Set initial concentration values in pool_matrixC:
      if (isVertical) then
        call initialize_vert(InitC, InitN, pool_matrixC, pool_matrixN, nlevdecomp)
      else
        call initialize_onelayer(InitC, pool_matrixC, InitN, pool_matrixN)
      end if !isVertical

      change_matrixC = 0.0
      change_matrixN = 0.0
      a_matrix      = pool_matrixC
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
      call fill_netcdf(run_name, nlevdecomp,t_init, pool_matrixC, change_matrixC, a_matrix, HR, vertC_change_sum, write_hour,current_month)
      call fill_netcdf(run_name, nlevdecomp,t_init, pool_matrixN, change_matrixN, a_matrix, HR, vertN_change_sum, write_hour,current_month)
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

      !Inputs from aboveground litter to LIT and protected SOM pools: GEP/depth [gC/m3h]
      lit_inputC(1,:) = (/fMET, (0.9-fMET)/)*GEP/(delta_z(1) + delta_z(2))!Partition between litter pools. The last 0.1 fraction goes directly to SOM
      som_inputC(1,:) = (/0.05, 0.05/)*GEP/(delta_z(1) + delta_z(2))
      myc_inputC(4,:) = (/ecm_frac, erm_frac, am_frac/)*GEP/(delta_z(6) + delta_z(6))
      myc_inputC(5,:) = (/ecm_frac, erm_frac,am_frac/)*GEP/(delta_z(5) + delta_z(5))
      !TODO Add N input to different layers
      !----------------------------------------------------------------------------------------------------------------
      do t =1,nsteps !Starting time iterations
        time = t*dt
        counter =counter + 1
        ycounter = ycounter + 1
        month_counter = month_counter + 1


        if (month_counter == days_in_month(current_month)*24) then
          previous_month = current_month
          !print*, previous_month, time
           !call disp('moist: ', r_moist)
           !call disp('Vmax: ', Vmax)
           !call disp('Km: ', Km)
          call read_clmdata(clm_data_file,TSOIL,SOILLIQ,SOILICE,WATSAT,current_month)
          call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist)
          !call disp('Temp: ',TSOIL)

          !call disp('r_moist', r_moist)

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

        do j = 1, nlevdecomp !For each depth level (for the no vertical transport case, nlevdecomp = 1, so loop is only done once):
          Km      = exp(Kslope*TSOIL(j) + Kint)*a_k*Kmod                    ![mgC/cm3]*10e3=[gC/m3]
          Vmax    = exp(Vslope*TSOIL(j) + Vint)*a_v*Vmod*r_moist(j)    ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?
          k_mycsom  = (/1.4,1.4,1.4/)*10e-3*r_moist(j)   ![1/h] Decay constants, myc som

          !Calculate fluxes between pools in level j:
          call microbial_fluxes(j, pool_matrixC,nlevdecomp)
          call som_fluxes(j, pool_matrixC,nlevdecomp)
          call litter_fluxes(j, pool_matrixC,nlevdecomp)

          do i = 1, pool_types !loop over all the pool types, i, in depth level j
            !This if-loop calculates dC/dt for the different carbon pools.NOTE: If pools are added/removed (i.e the actual model equations is changed), this loop needs to be updated.
            !The Gain and Loss variables are also used to calculate the analytical solution to dC/dt=Gain - Loss*C, a_matrix(j,i)
            !NOTE: The "change_matrixC" values correspond to the equations A11-A17 in Wieder 2015
            if (i==1) then !LITm
              Gain = lit_inputC(j,1)
              Loss = LITmSAPb + LITmSAPf

            elseif (i==2) then !LITs
              Gain = lit_inputC(j,2)
              Loss = LITsSAPb + LITsSAPf

            elseif (i==3) then !SAPb
              Gain = LITmSAPb*MGE(1) + LITsSAPb*MGE(2) &
              + SOMaSAPb*MGE(3) !+ EcMSAPb*MGE(3) + ErMSAPb*MGE(4) + AMSAPb*MGE(5)
              Loss =  SAPbSOMp + SAPbSOMa + SAPbSOMc

            elseif (i==4) then !SAPf
              Gain = LITmSAPf*MGE(4) + LITsSAPf*MGE(5) &
              + SOMaSAPf*MGE(6) !  + EcMSAPf*MGE(9) + ErMSAPf*MGE(10) + AMSAPf*MGE(11)
              Loss =  SAPfSOMp + SAPfSOMa + SAPfSOMc

            elseif (i==5) then !EcM
              Gain = myc_inputC(j,1)
              Loss = EcMSOMp + EcMSOMa + EcMSOMc !EcMSAPb + EcMSAPf

            elseif (i==6) then !ErM
              Gain = myc_inputC(j,2)
              Loss = ErMSOMp + ErMSOMa + ErMSOMc !ErMSAPb + ErMSAPf

            elseif (i==7) then !AM
              Gain = myc_inputC(j,3)
              Loss = AMSOMp + AMSOMa + AMSOMc !AMSAPb + AMSAPf

            elseif (i==8) then !SOMp
              !Use the same partitioning between the depth levels as for mycorrhiza (f_myc_levels)
              Gain = som_inputC(j,1) + SAPbSOMp + SAPfSOMp + EcMSOMp + ErMSOMp + AMSOMp
              Loss = SOMpSOMa

            elseif (i==9) then !SOMa
               Gain = SAPbSOMa + SAPfSOMa + EcMSOMa + ErMSOMa + AMSOMa +  SOMpSOMa + SOMcSOMa
               Loss = SOMaSAPb + SOMaSAPf
            elseif (i==10) then !SOMc
              Gain = som_inputC(j,2) + SAPbSOMc + SAPfSOMc + EcMSOMc + ErMSOMc + AMSOMc
              Loss = SOMcSOMa

            else
              print*, 'Too many pool types expected, pool_types = ',pool_types
            end if !determine dC_i/dt
            change_matrixC(j,i) = Gain - Loss
            change_sum(j,i)= change_sum(j,i) + change_matrixC(j,i)*dt
            !Store these values as temporary so that they can be used in the diffusion subroutine
            pool_temporaryC(j,i)=pool_matrixC(j,i) + change_matrixC(j,i)*dt

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

          !Calculate the heterotrophic respiration loss from depth level j in timestep t:
          HR(j) =( LITmSAPb*(1-MGE(1)) + LITsSAPb*(1-MGE(2))  + SOMaSAPb*(1-MGE(3)) + LITmSAPf*(1-MGE(4)) &
          + LITsSAPf*(1-MGE(5)) + SOMaSAPf*(1-MGE(6)))*dt

          HR_sum(j) = HR_sum(j) + HR(j)
        end do !j, depth_level
        !call disp('HR ', HR)
        if (isVertical) then
          call vertical_diffusion(tot_diff,upper,lower, pool_temporaryC, nlevdecomp,vertC,time, counter,dt)
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
          call fill_netcdf(run_name, nlevdecomp, t, pool_matrixC, change_sum, a_matrix,HR, vertC_change_sum, write_hour,current_month)
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

    end subroutine decomp
end module mycmim
