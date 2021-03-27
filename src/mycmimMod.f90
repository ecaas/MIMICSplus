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
  implicit none


  contains
    subroutine decomp(nsteps, run_name,nlevdecomp,step_frac) !Calculates the balance equations dC/dt and dN/dt for each pool at each time step based on the fluxes calculated in the same time step. Then update the pool sizes before moving on

      integer                        :: nsteps               ! number of time steps to iterate over
      character (len=*)              :: run_name             ! used for naming outputfiles
      integer                        :: step_frac            ! determines the size of the time step
      integer                        :: nlevdecomp           ! number of vertical layers

      logical                        :: isVertical                           ! True if we use vertical soil layers.
      real(r8),dimension(nlevdecomp) :: HR                                   ! For storing the C  that is lost to respiration [gC/m3h]
      real(r8)                       :: HR_mass_accumulated, HR_mass,growth_sum
      real(r8)                       :: pool_matrixC(nlevdecomp,pool_types)     ! For storing C pool sizes [gC/m3]
      real(r8)                       :: change_matrixC(nlevdecomp,pool_types)   ! For storing dC/dt for each time step [gC/(m3*hour)]
      real(r8)                       :: a_matrixC(nlevdecomp, pool_types)       ! For storing the analytical solution
      real(r8)                       :: pool_temporaryC(nlevdecomp,pool_types)  ! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport

      real(r8)                       :: pool_temporaryN(nlevdecomp,pool_types+1)! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
      real(r8)                       :: pool_matrixN(nlevdecomp,pool_types+1)   ! For storing N pool sizes [gN/m3] parallell to C pools and  inorganic N
      real(r8)                       :: change_matrixN(nlevdecomp,pool_types+1) ! For storing dC/dt for each time step [gN/(m3*hour)]
      real(r8)                       :: a_matrixN(nlevdecomp, pool_types+1)     ! For storing the analytical solution

      real(r8)                       :: mass_N(nlevdecomp, pool_types+1)     ! gN/m2
      real(r8)                       :: mass_C(nlevdecomp, pool_types)     ! gC/m2




     !Shape of pool_matrixC/change_matrixC
     !|       LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc |
     !|level1   1   2    3    4   5   6   7   8    9    10  |
     !|level2                                               |
     !| .                                                   |
     !| .                                                   |
     !|nlevdecomp __________________________________________|

     !Shape of the pool_matrixN/change_matrixN
     !|       LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc  IN|
     !|level1   1   2    3    4   5   6   7   8    9    10   11|
     !|level2                                                  |
     !| .                                                      |
     !| .                                                      |
     !|nlevdecomp __________________________________________   |

      real(r8)                       :: dt                            ! size of time step
      real(r8)                       :: time                          ! t*dt
      real(r8)                       :: C_Loss, C_Gain, N_Gain, N_Loss
      real(r8)                       :: tot_diffC,upperC,lowerC                 ! For the call to vertical_diffusion
      real(r8)                       :: tot_diffN,upperN,lowerN                 ! For the call to vertical_diffusion

      real(r8),allocatable           :: vertC(:,:)         !Stores the vertical change in a time step, same shape as change_matrixC
      real(r8),allocatable           :: vertN(:,:)         !Stores the vertical change in a time step, same shape as change_matrixC

      !Counters
      integer                        :: ycounter, year
      integer                        :: counter            !used for determining when to output results
      integer                        :: month_counter
      integer                        :: j,i,t              !for iterations
      integer,parameter              ::t_init=1
      real(r8), allocatable:: HR_sum(:)                    !Sums total HR between two output entries

      real(r8)                       :: change_sum(nlevdecomp, pool_types)
!      real(r8)                       :: vertN_change_sum(nlevdecomp, pool_types)
      real(r8)                       :: vertC_change_sum(nlevdecomp, pool_types)
      real(r8)                      :: possible_N_change, possible_C_change


      !For reading soil temperature and moisture from CLM output file
      real(r8), dimension(nlevdecomp)         :: TSOIL!(:), SOILLIQ(:),SOILICE(:),WATSAT(:),r_moist(:)
      real(r8), dimension(nlevdecomp)          :: SOILLIQ
      real(r8), dimension(nlevdecomp)          :: SOILICE
      real(r8), dimension(nlevdecomp)          :: WATSAT
      real(r8), dimension(nlevdecomp)          :: W_SCALAR
      real(r8), dimension(nlevdecomp)          :: r_moist


      integer,parameter              :: write_hour= 1*24*365!How often output is written to file
                                      !TODO: This should be input to the subroutine!

      !ALLOCATIONS:
      allocate(HR_sum(nlevdecomp))

      dt= 1.0/step_frac !Setting the time step
      
      if (nlevdecomp>1) then
        soil_depth=sum(delta_z(1:nlevdecomp))
        isVertical = .True.
      else
        soil_depth=1.0
        isVertical = .False.
        delta_z=soil_depth !So that delta_z will not be 1st on delta_z from parametersMod
                     !TODO: This can be done better
      end if

      print*, soil_depth

      ! Fracions of SAP that goes to different SOM pools
      fPHYS = (/ 0.3 * exp(fCLAY), 0.2 * exp(0.8*fCLAY) /)
      fCHEM =  (/0.1 * exp(-3.0*fMET), 0.3 * exp(-3*fMET) /)
      fAVAIL = 1-(fPHYS+fCHEM)

      !TODO FROM MIMICS_CYCLE_CN:
      ! WW also modify TAU as a function of soil moisture, so things don't
      ! colapse in frozen soils...
      !mimicsbiome%tauR(npt) = mimicsbiome%tauR(npt) * fW
      !mimicsbiome%tauK(npt) = mimicsbiome%tauK(npt) * fW
      tau = (/ 5e-4*exp(0.3*fMET), 5e-4*exp(0.1*fMET)/)![1/h] Microbial turnover rate (SAP to SOM), SAPr,(/1.39E-3*exp(0.3*fMET), 2.3E-4*exp(0.1*fMET)/)

      !Set initial concentration values:
      call initialize(pool_matrixC,pool_matrixN,CPlant,NPlant,nlevdecomp)

      !Make sure things start from zero
      change_matrixC = 0.0
      change_matrixN = 0.0
      HR            = 0.0
      vertC_change_sum=0.0
      counter = 0
      ycounter = 0
      HR_sum   = 0.0 !For summing up the total respiration between two output times
      HR_mass_accumulated = 0
      growth_sum=0

      a_matrixC      = pool_matrixC
      a_matrixN      = pool_matrixN

      year     = 1
      current_month = 1
      month_counter = 30


      !open and prepare files to store results. Store initial values
      call create_netcdf(run_name, nlevdecomp)
      call fill_netcdf(run_name,t_init, pool_matrixC, change_matrixC, pool_matrixN,change_matrixN, &
                       HR_mass_accumulated,HR, vertC_change_sum, write_hour,current_month, &
                      NPlant, CPlant,TSOIL, r_moist, growth_sum = growth_sum,levsoi=nlevdecomp)

      !read temperature and moisture data from CLM file
      call read_clmdata(clm_data_file,TSOIL,SOILLIQ,SOILICE,WATSAT,W_SCALAR,current_month, nlevdecomp)
      call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist,nlevdecomp)

      !----------------------------------------------------------------------------------------------------------------
      do t =1,nsteps !Starting time iterations
        time = t*dt
        counter  = counter + 1
        ycounter = ycounter + 1
        month_counter = month_counter + 1
        NPlant_tstep=0
        CPlant_tstep=0

        !Update temp and moisture values monthly
        if (month_counter == days_in_month(current_month)*24) then
          previous_month = current_month
          call read_clmdata(clm_data_file,TSOIL,SOILLIQ,SOILICE,WATSAT,W_SCALAR,current_month, nlevdecomp)
          call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist,nlevdecomp)
          if (current_month == 12) then
            current_month = 1
          else
            current_month = current_month + 1
          end if
          month_counter = 0
        end if

        if (ycounter == 365*24) then
          year = year + 1
          ycounter = 0
        end if

        !print initial values to terminal
        if (t == 1) then
          call disp("InitC", pool_matrixC)
          call disp("InitN", pool_matrixN)
          print*, "C plant", CPlant
          print*, "N plant", NPlant
        end if

        do j = 1, nlevdecomp !For each depth level (for the no vertical transport case, nlevdecomp = 1, so loop is only done once):
          !Michaelis Menten parameters:
          Km      = exp(Kslope*TSOI + Kint)*a_k*Kmod               ![mgC/cm3]*10e3=[gC/m3]
          Vmax    = exp(Vslope*TSOI + Vint)*a_v*Vmod!*W_SCALAR(j)   ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?

          k_mycsom  = (/1.4,1.4,1.4/)*10e-5!*W_SCALAR(j)  ![1/h] Decay constants, mycorrhiza to SOM pools TODO: Assumed, needs revision

          !Calculate fluxes between pools in level j (file: fluxMod.f90):
          call calculate_fluxes(j,nlevdecomp, pool_matrixC, pool_matrixN, CPlant, NPlant,isVertical)

          if (counter == write_hour .or. t==1) then !Write fluxes from calculate_fluxes to file
           call fluxes_netcdf(int(time), write_hour, j, run_name)
          end if !write fluxes

          growth_sum = growth_sum + C_growth_rate*dt


          !calculate the change of N and C in the plant based on the flux equations.  TODO: Needs better way to ensure reasonable values in these pools
          !CPlant_tstep and NPlant_tstep sum up the change from each layer, and will be used to update CPlant and NPlant at the end of the timestep.
          !(This will only be one value if isVertical = False)
          possible_N_change = (N_EcMPlant+  N_ErMPlant +  N_AMPlant + N_InPlant  &
          - N_PlantLITm - N_PlantLITs)*dt*delta_z(j) !gN/m2

          possible_C_change = C_growth_rate*dt-(  C_PlantLITm + C_PlantLITs + &
           C_PlantEcM + C_PlantErM + C_PlantAM)*dt*delta_z(j)!gC/m2

          NPlant_tstep = NPlant_tstep + possible_N_change
          CPlant_tstep = CPlant_tstep + possible_C_change

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
            !TODO Sap: Check if MGE/efficiency still makes sense in the new setup.
              C_Gain = e_s*(C_LITmSAPb + C_LITsSAPb &
                + C_SOMaSAPb + 0.5*(Decomp_ecm + Decomp_erm + Decomp_am))
              C_Loss =  C_SAPbSOMp + C_SAPbSOMa + C_SAPbSOMc
              N_Gain = N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb
              N_Loss = N_SAPbSOMp + N_SAPbSOMa + N_SAPbSOMc + N_SAPbIN

            elseif (i==4) then !SAPf
              C_Gain = e_s*(C_LITmSAPf + C_LITsSAPf &
                + C_SOMaSAPf + 0.5*(Decomp_ecm + Decomp_erm + Decomp_am))
              C_Loss =  C_SAPfSOMp + C_SAPfSOMa + C_SAPfSOMc
              N_Gain = N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf
              N_Loss = N_SAPfSOMp + N_SAPfSOMa + N_SAPfSOMc + N_SAPfIN

            elseif (i==5) then !EcM
              C_Gain = e_m*C_PlantEcM
              C_Loss = C_EcMSOMp + C_EcMSOMa + C_EcMSOMc
              N_Gain = N_INEcM + N_SOMaEcM
              N_Loss = N_EcMPlant + N_EcMSOMa + N_EcMSOMp + N_EcMSOMc

            elseif (i==6) then !ErM
              C_Gain = e_m*C_PlantErM
              C_Loss = C_ErMSOMp + C_ErMSOMa + C_ErMSOMc
              N_Gain = N_INErM + N_SOMaErM
              N_Loss = N_ErMPlant + N_ErMSOMa + N_ErMSOMp + N_ErMSOMc

            elseif (i==7) then !AM
              C_Gain = e_m*C_PlantAM
              C_Loss = C_AMSOMp + C_AMSOMa + C_AMSOMc
              N_Gain = N_INAM + N_SOMaAM
              N_Loss = N_AMPlant + N_AMSOMa + N_AMSOMp + N_AMSOMc

            elseif (i==8) then !SOMp
              C_Gain =  C_SAPbSOMp + C_SAPfSOMp + C_EcMSOMp + C_ErMSOMp + C_AMSOMp
              C_Loss = C_SOMpSOMa
              N_Gain =  N_SAPbSOMp + N_SAPfSOMp + N_EcMSOMp + N_ErMSOMp + N_AMSOMp
              N_Loss = N_SOMpSOMa

            elseif (i==9) then !SOMa
               C_Gain = C_SAPbSOMa + C_SAPfSOMa + C_EcMSOMa + &
               C_ErMSOMa + C_AMSOMa + C_SOMpSOMa + C_SOMcSOMa
               C_Loss = C_SOMaSAPb + C_SOMaSAPf +Decomp_ecm + Decomp_erm + Decomp_am
               N_Gain = N_SAPbSOMa + N_SAPfSOMa + N_EcMSOMa + &
               N_ErMSOMa + N_AMSOMa + N_SOMpSOMa + N_SOMcSOMa
               N_Loss = N_SOMaSAPb + N_SOMaSAPf + N_SOMaEcM + N_SOMaErM + N_SOMaAM

            elseif (i==10) then !SOMc
              C_Gain =  C_SAPbSOMc + C_SAPfSOMc + C_EcMSOMc + C_ErMSOMc + C_AMSOMc
              C_Loss = C_SOMcSOMa
              N_Gain =  N_SAPbSOMc + N_SAPfSOMc + N_EcMSOMc + N_ErMSOMc + N_AMSOMc
              N_Loss = N_SOMcSOMa

            elseif (i == 11) then !Inorganic N
              N_Gain = Deposition+ N_SAPbIN + N_SAPfIN
              N_Loss = Leaching+N_INEcM + N_INErM + N_INAM + N_InPlant
              change_matrixN(j,i) = N_Gain - N_loss
            else
              print*, 'Too many pool types expected, pool_types = ',pool_types, 'i: ', i
            end if !determine total gains and losses


            if (i /= 11) then !Carbon matrix does only have 10 columns (if i == 11 this is handeled inside the loop over pools)
                change_matrixC(j,i) = C_Gain - C_Loss
                change_matrixN(j,i) = N_Gain - N_loss
                !For summarize total change between each written output
                change_sum(j,i)= change_sum(j,i) + change_matrixC(j,i)*dt
                !Store these values as temporary so that they can be used in the vertical diffusion subroutine
                pool_temporaryC(j,i)=pool_matrixC(j,i) + change_matrixC(j,i)*dt

            end if

            pool_temporaryN(j,i) =pool_matrixN(j,i) + change_matrixN(j,i)*dt

            if (isnan(pool_temporaryN(j,i))) then
              print*, 'NaN NITROGEN value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryN(j,i)
              stop
            end if
            if (pool_temporaryN(j,i) < 0.001) then
              print*, 'Too small pool size: NITROGEN value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryN(j,i)
              pool_temporaryN(j,i)=0.01
            end if

            if (i /=11 ) then
              if (isnan(pool_temporaryC(j,i))) then
                print*, 'NaN CARBON value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryC(j,i)
                stop
              end if
              if (pool_temporaryC(j,i) < 0.001) then
                print*, 'Too small pool size: CARBON value at t',t,'depth level',j,'pool number',i, ':', pool_temporaryC(j,i)
                stop
              end if
            end if

          end do !i, pool_types

          !Calculate the heterotrophic respiration loss from depth level j in timestep t: NOTE: revise!
          HR(j) =(( C_LITmSAPb + C_LITsSAPb  + C_SOMaSAPb + C_LITmSAPf &
          + C_LITsSAPf + C_SOMaSAPf+Decomp_ecm + Decomp_erm + Decomp_am)*(1-e_s) &
          + (C_PlantEcM + C_PlantErM + C_PlantAM)*(1-e_m))*dt
          if (HR(j) < 0 ) then
            print*, 'Negative HR: ', HR(j), t
          end if
          HR_sum(j) = HR_sum(j) + HR(j)

        end do !j, depth_level

        !Update Plant pools with the total change from all the layers
        CPlant = CPlant + CPlant_tstep
        NPlant = NPlant + NPlant_tstep

        if (CPlant < 0.001 .or. NPlant < 0.001) then
          print*, 'Too small pool size: CARBON value at t',t,'CPlant:', CPlant
          print*, 'Too small pool size: Nitrogen value at t',t,'NPlant:', NPlant

          stop
        end if

        !Store accumulated HR mass
        call respired_mass(HR, HR_mass,nlevdecomp)
        HR_mass_accumulated = HR_mass_accumulated + HR_mass

        if (isVertical) then
          call vertical_diffusion(tot_diffC,upperC,lowerC, pool_temporaryC,vertC)
          call vertical_diffusion(tot_diffN,upperN,lowerN, pool_temporaryN,vertN)
          pool_matrixC =  vertC*dt + pool_temporaryC
          pool_matrixN = vertN*dt + pool_temporaryN
        else
          pool_matrixC=pool_temporaryC
          pool_matrixN=pool_temporaryN
        end if!isVertical

        if (counter == write_hour) then
          counter = 0
          call fill_netcdf(run_name, int(time), pool_matrixC, change_matrixC, pool_matrixN,change_matrixN,&
           HR_mass_accumulated,HR,vertC_change_sum, write_hour,current_month, NPlant,&
           CPlant, TSOIL, r_moist, growth_sum,nlevdecomp)
          change_sum = 0.0
          vertC_change_sum = 0.0
        end if!writing

        !Write end values to terminal
        if (t == nsteps) then
          call store_parameters(run_name)
          call disp("pool_matrixC gC/m3 ",pool_matrixC)
          call disp("pool_matrixC gN/m3 ",pool_matrixN)
          Print*, NPlant, CPlant
        end if

      end do !t

    end subroutine decomp


end module mycmim
