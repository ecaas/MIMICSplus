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
      real(r8)                       :: pool_matrix(nlevdecomp,pool_types)   ! For storing C pool sizes [gC/m3]
      real(r8)                       :: mass_pool_matrix(nlevdecomp,pool_types)   ! For storing C pool sizes [gC/m2] (pool_matrix*delta_z)
      real(r8)                       :: change_matrix(nlevdecomp,pool_types) ! For storing dC/dt for each time step [gC/(m3*hour)]
      real(r8)                       :: a_matrix(nlevdecomp, pool_types)     ! For storing the analytical solution
      real(r8)                       :: Init(nlevdecomp, pool_types)         ! Initial C concentration, determined in initMod.f90
      real(r8)                       :: pool_temporary(nlevdecomp,pool_types)! When isVertical is True, pool_temporary = pool_matrix + change_matrix*dt is used to calculate the vertical transport
                                                                             !The new value after the time step is then pool_matrix = pool_temporary + vertical change.

                                                                             !Shape of the pool_matrix/change_matrix
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
      real(r8)                       :: vert(nlevdecomp, pool_types)         !Stores the vertical change in a time step, on the same form as change_matrix

      real(r8)                       :: lit_input(nlevdecomp,no_of_litter_pools)![gC/(m3 h)] Fraction of litter input to LITm and LITs, respectively
      real(r8)                       :: myc_input(nlevdecomp,no_of_myc_pools)   ![gC/(m3 h)] vector giving the input from vegetation to mycorrhiza pools
      real(r8)                       :: som_input(nlevdecomp,no_of_som_pools-1) ![gC/(m3 h)] !only input to SOMp and SOMc
      real(r8)                       :: I_tot                                   ![gC/(m3 h)]Total input to the system
      integer                        :: ycounter, year,ncid,varid
      integer                        :: counter                                 !used for determining when to output results
      integer                        :: month_counter
      integer                        :: j,i,t                                   !for iterations
      integer,parameter              ::t_init=1
      real(r8), dimension(nlevdecomp):: HR_sum                                  !Sums total HR between two output entries
      real(r8)                       :: change_sum(nlevdecomp, pool_types)
      real(r8)                       :: vert_change_sum(nlevdecomp, pool_types)
      integer,parameter              :: write_hour= 24*4
      real(r8)                       :: ecm_frac, erm_frac, am_frac
      !Assigning values: (Had to move from paramMod to here to be able to modify them during a run)
      MGE       = (/0.3,0.3,0.3,0.3,0.3,0.3,0.4,0.4,0.3,0.3,0.3,0.4/)
      k_mycsom  = (/1.4,1.4,1.4/)*10e-4   ![1/h] Decay constants, myc som
      k_mycsap  = (/1.0,1.0,1.0/)*10e-4 ![1/h] Decay constants, myc sap
      lit_input = 0.0
      myc_input = 0.0
      som_input = 0.0
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

      !Set initial concentration values in pool_matrix:
      if (isVertical) then
        call initialize_vert(Init, pool_matrix, nlevdecomp)
      else
        call initialize_onelayer(Init, pool_matrix)
      end if !isVertical

      change_matrix = 0.0
      a_matrix      = pool_matrix
      HR            = 0.0
      vert_change_sum=0.0

      counter = 0
      ycounter = 0
      year     = 1
      HR_sum   = 0.0 !For summing up the total respiration between two output times
      current_month = 1
      month_counter = 30


      !open and prepare files to store results. Store initial value
      call create_netcdf(run_name)
      call fill_netcdf(run_name, nlevdecomp,t_init, pool_matrix, change_matrix, a_matrix, HR, vert_change_sum, write_hour,current_month)
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
      lit_input(1,:) = (/fMET, (0.9-fMET)/)*GEP/(delta_z(1) + delta_z(2))!Partition between litter pools. The last 0.1 fraction goes directly to SOM
      som_input(1,:) = (/0.05, 0.05/)*GEP/(delta_z(1) + delta_z(2))
      myc_input(4,:) = (/ecm_frac, erm_frac, am_frac/)*GEP/(delta_z(6) + delta_z(6))
      myc_input(5,:) = (/ecm_frac, erm_frac,am_frac/)*GEP/(delta_z(5) + delta_z(5))
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
          !print*,previous_month, current_month, month_counter/24
          month_counter = 0
        end if

        if (ycounter == 365*24) then
          year = year + 1
          !print*, year
          ycounter = 0
        end if
        if (t == 1) then
          call disp("Init", pool_matrix)
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
        !   lit_input(1,:) = (/fMET, (0.9-fMET)/)*GEP/delta_z(1) !Partition between litter pools. The last 0.1 fraction goes directly to SOM
        !   som_input(1,:) = (/0.05, 0.05/)*GEP/delta_z(1)
        !   !print*, "TIME: ", time
        !   !print*, lit_input(1,:)
        ! else
        !   lit_input = 0
        !   som_input =0
        ! end if
        tsoi = 1.7!8.9

        !call disp("myc",myc_input)
        do j = 1, nlevdecomp !For each depth level (for the no vertical transport case, nlevdecomp = 1, so loop is only done once):
          !call disp('moisture',r_moist)
          Km      = exp(Kslope*TSOIL(j) + Kint)*a_k*Kmod                    ![mgC/cm3]*10e3=[gC/m3]
          Vmax    = exp(Vslope*TSOIL(j) + Vint)*a_v*Vmod*r_moist(j)    ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?
          k_mycsom  = (/1.4,1.4,1.4/)*10e-3*r_moist(j)   ![1/h] Decay constants, myc som
          k_mycsap  = 0! All myc necromass goes direectly to SOMa pool !(/1.0,1.0,1.0/)*10e-3*r_moist(j) ![1/h] Decay constants, myc sap
          !print*, time, j
          !call disp(Vmax)
          !print*, r_moist(j), TSOIL(j)

          !myc_input(6,:) = (/pool_matrix(5,5)+1, pool_matrix(5,6)+1, pool_matrix(5,7)+1/)*GEP*0.001



          ! if (isVertical) then
          !   if (j==1) then !The litter input is higher in the first depth level then the rest.
          !     lit_input=(/fMET*0.25, (1-fMET)*0.25/)*I_tot*0.5
          !     !som_input = ()
          !   else
          !     lit_input=(/fMET*0.25, (1-fMET)*0.25/)*I_tot*0.5 !TODO Change this so it is not always the same
          !     !som_input()
          !   end if !j=1
          ! else
          !   som_input= (/f_som1,f_som2/)*I_tot*0.1
          !   !print*, time, 335*24 + (year-1)*365*24 -244*24 - (year-1)*365*24
          !   lit_input=(/fMET, 1-fMET/)*I_tot*0.5 ! input in gC/m³*day
          !

          ! end if !isVertical


          !Calculate fluxes between pools in level j:
          call microbial_fluxes(j, pool_matrix,nlevdecomp)
          call som_fluxes(j, pool_matrix,nlevdecomp)
          call litter_fluxes(j, pool_matrix,nlevdecomp)
        !  print*, LITsSAPf, LITmSAPb, lit_input(1,:)
          !   if (counter == write_hour) then
          !
          !     call check(nf90_open(trim(run_name)//".nc", nf90_write, ncid))
          !     call check(nf90_inq_varid(ncid, "time", varid))
          !     call check(nf90_put_var(ncid, varid, time, start = (/ t/write_hour+1/)))
          ! !    call check(nf90_inq_varid(ncid, "LITmSAPb", varid))
          ! !    call check(nf90_put_var(ncid, varid, LITmSAPb, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "LITsSAPb", varid))
          ! !    call check(nf90_put_var(ncid, varid, LITsSAPb, start = (/ t/48, j /)))
          !     !call check(nf90_inq_varid(ncid, "EcMSAPb", varid))
          !     !call check(nf90_put_var(ncid, varid, EcMSAPb, start = (/ t/48, j /)))
          !     !call check(nf90_inq_varid(ncid, "ErMSAPb", varid))
          !     !call check(nf90_put_var(ncid, varid, ErMSAPb, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "AMSAPb", varid))
          ! !    call check(nf90_put_var(ncid, varid, AMSAPb, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "LITmSAPf", varid))
          ! !    call check(nf90_put_var(ncid, varid, LITmSAPf, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "LITsSAPf", varid))
          ! !    call check(nf90_put_var(ncid, varid, LITsSAPf, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "EcMSAPf", varid))
          ! !    call check(nf90_put_var(ncid, varid, EcMSAPf, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "ErMSAPf", varid))
          ! !    call check(nf90_put_var(ncid, varid, ErMSAPf, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "AMSAPf", varid))
          ! !    call check(nf90_put_var(ncid, varid, AMSAPf, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "SOMpSOMa", varid))
          ! !    call check(nf90_put_var(ncid, varid, SOMpSOMa, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "SOMcSOMa", varid))
          ! !    call check(nf90_put_var(ncid, varid, SOMcSOMa, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "SAPbSOMa", varid))
          ! !    call check(nf90_put_var(ncid, varid, SAPbSOMa, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "SAPfSOMa", varid))
          ! !    call check(nf90_put_var(ncid, varid, SAPfSOMa, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "EcMSOMa", varid))
          ! !    call check(nf90_put_var(ncid, varid, EcMSOMa, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "ErMSOMa", varid))
          ! !    call check(nf90_put_var(ncid, varid, ErMSOMa, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "AMSOMa", varid))
          ! !    call check(nf90_put_var(ncid, varid, AMSOMa, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "SOMaSAPb", varid))
          ! !    call check(nf90_put_var(ncid, varid, SOMaSAPb, start = (/ t/48, j /)))
          ! !    call check(nf90_inq_varid(ncid, "SOMaSAPf", varid))
          ! !    call check(nf90_put_var(ncid, varid, SOMaSAPf, start = (/ t/48, j /)))
          !     call check(nf90_close(ncid))
          ! !
          !  end if !writing

          do i = 1, pool_types !loop over all the pool types, i, in depth level j
            !This if-loop calculates dC/dt for the different carbon pools.NOTE: If pools are added/removed (i.e the actual model equations is changed), this loop needs to be updated.
            !The Gain and Loss variables are also used to calculate the analytical solution to dC/dt=Gain - Loss*C, a_matrix(j,i)
            !NOTE: The "change_matrix" values correspond to the equations A11-A17 in Wieder 2015
            if (i==1) then !LITm
              Gain = lit_input(j,1)
              Loss = LITmSAPb + LITmSAPf

            elseif (i==2) then !LITs
              Gain = lit_input(j,2)
              Loss = LITsSAPb + LITsSAPf

            elseif (i==3) then !SAPb
              Gain = LITmSAPb*MGE(1) + LITsSAPb*MGE(2) &
              + EcMSAPb*MGE(3) + ErMSAPb*MGE(4) + AMSAPb*MGE(5) + SOMaSAPb*MGE(6)
              Loss =  SAPbSOMp + SAPbSOMa + SAPbSOMc


            elseif (i==4) then !SAPf
              Gain = LITmSAPf*MGE(7) + LITsSAPf*MGE(8) &
              + EcMSAPf*MGE(9) + ErMSAPf*MGE(10) + AMSAPf*MGE(11) + SOMaSAPf*MGE(12)
              Loss =  SAPfSOMp + SAPfSOMa + SAPfSOMc
            elseif (i==5) then !EcM
              Gain = myc_input(j,1)
              Loss = EcMSAPb + EcMSAPf + EcMSOMp + EcMSOMa + EcMSOMc

            elseif (i==6) then !ErM
              Gain = myc_input(j,2)
              Loss = ErMSAPb + ErMSAPf + ErMSOMp + ErMSOMa + ErMSOMc

            elseif (i==7) then !AM
              Gain = myc_input(j,3)
              Loss = AMSAPb + AMSAPf + AMSOMp + AMSOMa + AMSOMc

            elseif (i==8) then !SOMp
              !Use the same partitioning between the depth levels as for mycorrhiza (f_myc_levels)
              Gain = som_input(j,1) + SAPbSOMp + SAPfSOMp + EcMSOMp + ErMSOMp + AMSOMp
              Loss = SOMpSOMa

            elseif (i==9) then !SOMa
               Gain = SAPbSOMa + SAPfSOMa + EcMSOMa + ErMSOMa + AMSOMa +  SOMpSOMa + SOMcSOMa
               Loss = SOMaSAPb + SOMaSAPf
            elseif (i==10) then !SOMc
              Gain = som_input(j,2) + SAPbSOMc + SAPfSOMc + EcMSOMc + ErMSOMc + AMSOMc
              Loss = SOMcSOMa

            else
              print*, 'Too many pool types expected, pool_types = ',pool_types
            end if !determine dC_i/dt
            change_matrix(j,i) = Gain - Loss
            change_sum(j,i)= change_sum(j,i) + change_matrix(j,i)*dt
            !Store these values as temporary so that they can be used in the diffusion subroutine
            pool_temporary(j,i)=pool_matrix(j,i) + change_matrix(j,i)*dt

            !control check
            if (pool_temporary(j,i) < 0.0) then
              print*, 'Negative concentration value at t',t,'depth level',j,'pool number',i, ':', pool_temporary(j,i)
              print*, 'Value changed to: ', min_pool_value
              print*, 'Month, year: ', current_month, year
              call disp('Temp: ',TSOIL)
              call disp('moist: ', r_moist)
              call disp('Vmax: ', Vmax)
              call disp('Km: ', Km)
              !call disp(pool_temporary)
              pool_temporary(j,i) = min_pool_value
              STOP
            end if

            !To calculate analytical solution:
            Loss_term = Loss/pool_matrix(j,i)
            a_matrix(j,i) = Init(j,i)*exp(-time*Loss_term) + Gain*(1-exp(-time*Loss_term))/Loss_term !+ vert*dt

          end do !i, pool_types

          !Calculate the heterotrophic respiration loss from depth level j in timestep t:
          HR(j) =( LITmSAPb*(1-MGE(1)) + LITsSAPb*(1-MGE(2)) + EcMSAPb*(1-MGE(3)) + ErMSAPb*(1-MGE(4)) + AMSAPb*(1-MGE(5)) + SOMaSAPb*(1-MGE(6)) + LITmSAPf*(1-MGE(7)) &
          + LITsSAPf*(1-MGE(8)) + EcMSAPf*(1-MGE(9)) + ErMSAPf*(1-MGE(10)) + AMSAPf*(1-MGE(11)) + SOMaSAPf*(1-MGE(12)))*dt

          HR_sum(j) = HR_sum(j) + HR(j)
        end do !j, depth_level
        !call disp('HR ', HR)
        if (isVertical) then
          call vertical_diffusion(tot_diff,upper,lower, pool_temporary, nlevdecomp,vert,time, counter,dt)
          pool_matrix =  vert*dt + pool_temporary
          vert_change_sum = vert_change_sum + vert*dt
          !call disp("pool_matrix",pool_matrix)
          a_matrix=a_matrix+vert*dt
        else
          pool_matrix=pool_temporary
        end if!isVertical

        if (counter == write_hour) then

          counter = 0
          !call disp("pool_matrix ",pool_matrix)
          call fill_netcdf(run_name, nlevdecomp, t, pool_matrix, change_sum, a_matrix,HR, vert_change_sum, write_hour,current_month)
          HR_sum = 0.0
          change_sum = 0.0
          vert_change_sum = 0.0
        !print*, HR_sum
        end if!writing

        if (t == nsteps) then

          do i = 1,pool_types
            mass_pool_matrix(:,i) = pool_matrix(:,i)*delta_z
          end do
          call disp("pool_matrix gC/m2 ",mass_pool_matrix)
          call disp("pool_matrix gC/m3 ",pool_matrix)

          !call disp("respiration", HR)
        end if

      end do !t
      !call store_parameters(run_name)

    end subroutine decomp
end module mycmim
