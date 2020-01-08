!This module models the decomposition of carbon through the soil. The number of soil depth levels and different pools are given in the file paramMod.f90. The paramMod module contains
!all relevant parameters for the flux and balance equations used in this module. In the current setup, 4 depth levels are used, and the following 10 carbon pools/reservoirs are found in each level:
!LITm - metabolic litter
!LITs - Structural litter
!SAPr - r-strategist saprotrophs
!SAPk - K strategist saprotrophs
!EcM  - Ectomycorrhiza
!ErM  - Ericoid mycorrhiza
!AM   - Arbuscular mycorrhiza
!SOMp - Physically protected soil organic matter
!SOMa - Available soil organic matter
!SOMc - Chemically protected soil organic matter

!TODO: Incorporate frozen soil behavior/variable T, Variable clay content(?)
module mycmim
  use paramMod
  !use dispmodule !External module to pretty print matrices (mainly for testing purposes)
  use fluxMod
  use initMod
  use writeMod
  implicit none

  contains
    subroutine decomp(nsteps, run_name, isVertical, nlevdecomp, ecosystem,step_frac) !Calculates the balance equation dC/dt for each pool at each time step based on the fluxes calculated in the same time step. Then update the pool sizes before moving on
                                                                                     !It also calculates an analytical solution to the problem.
      integer                        :: nsteps                               ! number of time steps to iterate over
      character (len=*)              :: run_name                             ! used for naming outputfiles
      logical                        :: isVertical                           ! True if we use vertical soil layers.
      integer                        :: nlevdecomp                           ! number of vertical layers
      character (len=*)              :: ecosystem                            ! 'Shrub', 'Heath' or 'Meadow' (For comparison with Sorensen et al)
      integer                        :: step_frac                            ! determines the size of the time step

      real(r8)                       :: pool_matrix(nlevdecomp,pool_types)   ! For storing C pool sizes [gC/m3]
      real(r8)                       :: change_matrix(nlevdecomp,pool_types) ! For storing dC/dt for each time step [gC/(m3*day)]
      real(r8)                       :: a_matrix(nlevdecomp, pool_types)     ! For storing the analytical solution
      real(r8),dimension(nlevdecomp) :: HR                                   ! For storing the C amount that is lost to respiration
      real(r8)                       :: Init(nlevdecomp, pool_types)         ! Initial C concentration, determined in initMod.f90
      real(r8)                       :: pool_temporary(nlevdecomp,pool_types)! When isVertical is True, pool_temporary = pool_matrix + change_matrix*dt is used to calculate the vertical transport
                                                                             !The new value after the time step is then pool_matrix = pool_temporary + vertical change.

                                                                             !Shape of the pool_matrix/change_matrix
                                                                             !|       LITm LITs SAP EcM ErM AM SOMp SOMa SOMc |
                                                                             !|level1   1   2    3  4   5   6   7   8    9    |
                                                                             !|level2                                         |
                                                                             !| .                                             |
                                                                             !| .                                             |
                                                                             !|nlevdecomp ____________________________________|


      real(r8)                       :: dt                                   ! size of time step
      real(r8)                       :: time                                 ! t*dt
      real(r8)                       :: Loss, Gain, Loss_term                           ! Source and sink term for solving the analytical solution to the dC/dt=Gain-Loss*concentration equation.
      real(r8)                       :: tot_diff,upper,lower                 ! For the call to alt_vertical_diffusion
      real(r8)                       :: vert(nlevdecomp, pool_types)         !Stores the vertical change in a time step, on the same form as change_matrix
      real(kind=r8)                  :: GEP                                  ![gC/(m2 h)] Gross ecosystem productivity

      real(r8), dimension(2)         :: lit_input                            ![gC/(m3 h)] Fraction of litter input to LITm and LITs, respectively
      real(r8), dimension(3)         :: myc_input                            ![gC/(m3 h)] vector giving the input from vegetation to mycorrhiza pools
      real(r8), dimension(2)         :: som_input                            ![gC/(m3 h)]
      real(r8)                       :: I_tot                                ![gC/(m3 h)]Total input to the system
      integer                        :: ycounter, year,ncid,varid
      integer                        :: counter=0                            ! used for determining when to output results
      integer                        :: j,i,t                       ! for iterations

      !Assigning values: (Had to move from paramMod to here to be able to modify them during a run)
      MGE     = (/ 0.55,0.75,0.25,0.35, 1.0, 1.0 /)

      Kmod    = (/0.125d0, 0.50d0, 0.167d0*pscalar*10 /)!, 0.5d0, 0.25d0, 0.167d0*pscalar*10/)!LITm, LITs, SOMa entering SAP
      Vmod    = (/4.0,  5.0, 10.0 /)! 3.0, 3.0, 2.0/)                          !LITm, LITs, SOMa entering SAP

      k2 = (/7.0,7.0,1.4/)*10e-4        ![1/h] Decay constants
      k = (/5.0,5.0,0.5/)*10e-5 ![1/h] Decay constants

      dt= 1.0/step_frac
 !NOTE: Must change if isVertical is True!!
                        !NOTE: If vertical myc_input must also be changed, because different amounts go in different layers.
      if (ecosystem == 'Heath') then
        GEP       = 0.281
        I_tot = GEP/depth !NOTE: Must change if isVertical is True!!
                                      ![gC/(m2*h)] Gross ecosystem productivity

        myc_input = (/0.10,0.80,0.10/)*I_tot*0.4       ![gC/(m3*h)] For Heath, most to ErM
        fMET      = 0.5                                !Fraction of input to litter that goes to the metabolic litter pool
      elseif (ecosystem == 'Meadow') then

        GEP       = 0.385
        I_tot = GEP/depth !NOTE: Must change if isVertical is True!!

        k2 = (/7.0,7.0,1.4/)*10e-6
        myc_input = (/0.1,0.1,0.80/)*I_tot*0.4         ![gC/(m3*h)] For meadow, most to AM
        fMET      = 0.5
      elseif (ecosystem == 'Shrub') then

        GEP       = 0.491
        I_tot = GEP/depth !NOTE: Must change if isVertical is True!!

        myc_input = (/0.80,0.10,0.10/)*I_tot*0.4       ![gC/(m3*h)] For shrub, most to EcM.
        fMET      = 0.2
      else
        print*, 'Invalid ecosystem name', ecosystem
        stop
      end if

      !Assigning values. Fracions of SAP that goes to different SOM pools
      fPHYS =  0.3 * exp(1.3*fCLAY)!, 0.2 * exp(0.8*fCLAY) /)
      fCHEM =  0.1 * exp(-3.0*fMET)!, 0.3 * exp(-3*fMET) /)
      fAVAIL = 1-(fPHYS+fCHEM)
      !print*, fPHYS, fCHEM, fAVAIL
      !print*, fPHYS+fCHEM+fAVAIL
      tau = 5.2e-4*exp(0.3*fMET)*20!, 2.4e-4*exp(0.1*fMET)/)*10![1/h] Microbial turnover rate (SAP to SOM), SAPr,

      !Set initial concentration values in pool_matrix:
      if (isVertical) then
        call initialize_vert(Init, pool_matrix, nlevdecomp)
      else
        call initialize_onelayer(Init, pool_matrix)
      end if !isVertical

      !open and prepare files to store results
      call create_netcdf(run_name)

      ycounter = 0
      year = 1
      print*, step_frac
      !----------------------------------------------------------------------------------------------------------------
      do t =1,nsteps !Starting time iterations
        time = t*dt
        counter =counter + 1
        ycounter = ycounter + 1
        !print*, year, ycounter/24
        if (ycounter == 365*24) then
          year = year + 1
          ycounter = 0
        end if


        !If-test used to modify something after half of the total run time
         !if (t == nsteps/2) then
           ! print*, tau
           ! print*, 'CHANGED TAU'
           ! tau = (/ 5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET)/)*1.5 ![1/h] Microbial turnover rate (SAP to SOM), SAPr, SAPk
           ! print*, tau
           ! print*, MGE
           ! print*, 'Changed MGE!'
           ! MGE     = (/ 0.25,0.35,0.55,0.75 /)
          !  ! print*, MGE
         !   print*, Km,Vmax
         !   Kmod    = Kmod/15!LITm, LITs, SOMa entering SAPr, LITm, LITs, SOMa entering sapk
         !   Vmod    = Vmod/15    !LITm, LITs, SOMa entering SAPr, LITm, LITs, SOMa entering sapk
         !   Km      = exp(Kslope*tsoi+3.19)*10*Kmod  ![mgC/cm3]*10e3=[gC/m3]
         !   Vmax    = exp(0.063*tsoi + 5.47)*8e-6*Vmod ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?
         !   print*, Km,Vmax
         ! end if !Change

        ! if (time >= 91*24 + (year - 1)*365*24 .and. time <= 304*24 + (year -1)*365*24  ) then
        !   tsoi = 10.0
        ! else
        !   tsoi = 0.0
        ! end if
        tsoi = 1.7!8.9
        Km      = exp(Kslope*tsoi+3.19)*10*Kmod       ![mgC/cm3]*10e3=[gC/m3]
        Vmax    = exp(0.063*tsoi + 5.47)*8e-6*Vmod    ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?

        do j = 1, nlevdecomp !For each depth level (for the no vertical transport case, nlevdecomp = 1, so loop is only done once):

          if (isVertical) then
            if (j==1) then !The litter input is higher in the first depth level then the rest.
              lit_input=(/fMET*0.25, fMET*0.25/)*I_tot*0.5
              !som_input = ()
            else
              lit_input=(/fMET*0.25, fMET*0.25/)*I_tot*0.5 !TODO Change this so it is not always the same
              !som_input()
            end if !j=1
          else
            som_input= (/f_som1,f_som2/)*I_tot*0.1
            !print*, time, 335*24 + (year-1)*365*24 -244*24 - (year-1)*365*24
            lit_input=(/fMET, 1-fMET/)*I_tot*0.5 ! input in gC/m³*day

              ! if (time >= 244*24 + (year-1)*365*24 .and. time <= 335*24 + (year-1)*365*24) then
              !   !print*, time, 335*24 + year*365*24 - 244*24 - year*365*24
              !   lit_input=(/fMET, 1-fMET/)*I_tot*0.5 ! input in gC/m³*day
              ! else
              !   lit_input = 0
              ! end if
          end if !isVertical

          !Calculate fluxes between pools in level j:
          call microbial_fluxes(j, pool_matrix,nlevdecomp)
          call som_fluxes(j, pool_matrix,nlevdecomp)
          call litter_fluxes(j, pool_matrix,nlevdecomp)

          ! !TODO: This writing to file should be made much more efficient, and to binary files, not text files..
          !  if (counter == 24) then
          !    call check(nf90_open(trim(run_name)//".nc", nf90_write, ncid))
          !
          !    call check(nf90_inq_varid(ncid, "LITmSAP", varid))
          !    call check(nf90_put_var(ncid, varid, LITmSAP, start = (/ t/24, j /)))
          !    call check(nf90_inq_varid(ncid, "LITsSAP", varid))
          !    call check(nf90_put_var(ncid, varid, LITsSAP, start = (/ t/24, j /)))
          !
          !
          !    call check(nf90_inq_varid(ncid, "EcMSAP", varid))
          !    call check(nf90_put_var(ncid, varid, EcMSAP, start = (/ t/24, j /)))
          !
          !    call check(nf90_inq_varid(ncid, "ErMSAP", varid))
          !    call check(nf90_put_var(ncid, varid, ErMSAP, start = (/ t/24, j /)))
          !
          !    call check(nf90_inq_varid(ncid, "AMSAP", varid))
          !    call check(nf90_put_var(ncid, varid, AMSAP, start = (/ t/24, j /)))
          !
          !    call check(nf90_close(ncid))
          ! !   write(unit=3,fmt='(F10.0,A2,I2,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10)') &
          ! !   time,',',j,',',LITtoSAP(1),',',LITtoSAP(2),',',LITtoSAP(3),',',LITtoSAP(4)
          ! !   write(unit=4,fmt='(F10.0,A2,I2,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10)') &
          ! !   time,',',j,',',SAPtoSOM(1),',',SAPtoSOM(2),',',SAPtoSOM(3),',',SAPtoSOM(4),',',SAPtoSOM(5),',',SAPtoSOM(6)
          ! !   write(unit=7,fmt='(F10.0,A2,I2,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10)') &
          ! !   time,',',j,',',MYCtoSAP(1),',',MYCtoSAP(2),',',MYCtoSAP(3),',',MYCtoSAP(4),',',MYCtoSAP(5),',',MYCtoSAP(6)
          ! !   write(unit=8,fmt='(F10.0,A2,I2,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10)') &
          ! !   time,',',j,',',MYCtoSOM(1),',',MYCtoSOM(2),',',MYCtoSOM(3),',' &
          ! !   ,MYCtoSOM(4),',',MYCtoSOM(5),',',MYCtoSOM(6),',',MYCtoSOM(7),',',MYCtoSOM(8),',',MYCtoSOM(9)
          ! !   write(unit=9,fmt='(F10.0,A2,I2,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10)') &
          ! !   time,',',j,',',SOMtoSAP(1),',',SOMtoSAP(2),',',SOMtoSOM(1),',',SOMtoSOM(2)
          !  end if !writing

          do i = 1, pool_types !loop over all the pool types, i, in depth level j
            !This if-loop calculates dC/dt for the different carbon pools.NOTE: If pools are added/removed (i.e the actual model equations is changed), this loop needs to be updated.
            !The Gain and Loss variables are used to calculate the analytical solution to dC/dt=Gain - Loss*C, a_matrix(j,i)
            !NOTE: The "change_matrix" values correspond to the equations A11-A17 in Wieder 2015
            if (i==1) then !LITm
              Gain = lit_input(1)
              Loss = LITmSAP/pool_matrix(j, i)

            elseif (i==2) then !LITs
              Gain = lit_input(2)
              Loss = LITsSAP/pool_matrix(j, i)

            elseif (i==3) then !SAP
              Gain = LITmSAP*MGE(1) + LITsSAP*MGE(3) &
              + (EcMSAP + ErMSAP + AMSAP)*MGE(5) + SOMaSAP*MGE(1)
              Loss =  SAPSOMp + SAPSOMa + SAPSOMc

            elseif (i==4) then !EcM
              Gain = myc_input(1)
              Loss = EcMSAP + EcMSOMp + EcMSOMa + EcMSOMc

            elseif (i==5) then !ErM
              Gain = myc_input(2)
              Loss = ErMSAP + ErMSOMp + ErMSOMa + ErMSOMc

            elseif (i==6) then !AM
              Gain = myc_input(3)
              Loss = AMSAP + AMSOMp + AMSOMa + AMSOMc

            elseif (i==7) then !SOMp
              !Use the same partitioning between the depth levels as for mycorrhiza (f_myc_levels)
              Gain = som_input(1)*f_myc_levels + SAPSOMp + EcMSOMp + ErMSOMp + AMSOMp
              Loss = SOMpSOMa

            elseif (i==8) then !SOMa
               Gain = SAPSOMa + EcMSOMa + ErMSOMa + AMSOMa +  SOMpSOMa + SOMcSOMa
               Loss = SOMaSAP

            elseif (i==9) then !SOMc
              Gain = som_input(1) + SAPSOMc + EcMSOMc + ErMSOMc + AMSOMc
              Loss = SOMcSOMa

            else
              print*, 'Too many pool types expected, pool_types = ',pool_types
            end if !determine dC_i/dt
            change_matrix(j,i) = Gain - Loss
            !Store these values as temporary so that they can be used in the diffusion subroutine
            pool_temporary(j,i)=pool_matrix(j,i) + change_matrix(j,i)*dt

            !control check
            if (pool_temporary(j,i) < 0.0) then
              print*, 'Negative concentration value at t',t,'depth level',j,'pool number',i, ':', pool_temporary(j,i)
              STOP
            end if

            !To calculate analytical solution:
            Loss_term = Loss/pool_matrix(j,i)
            a_matrix(j,i) = Init(j,i)*exp(-time*Loss) + Gain*(1-exp(-time*Loss))/Loss !+ vert*dt
          end do !i, pool_types

          !Calculate the heterotrophic respiration loss from depth level j in timestep t:
          HR(j) =( LITmSAP*(1-MGE(1)) + LITsSAP*(1-MGE(2)) + SOMaSAP*(1-MGE(1)))*dt
        end do !j, depth_level

        if (isVertical) then
          call vertical_diffusion(tot_diff,upper,lower, pool_temporary, nlevdecomp,vert,time, counter,dt)
          pool_matrix =   pool_temporary+vert*dt
          a_matrix=a_matrix+vert*dt
        else
          pool_matrix=pool_temporary
        end if!isVertical

        if (counter == 24) then
          counter = 0
          call fill_netcdf(run_name, nlevdecomp, t, pool_matrix, change_matrix)
        end if!writing

      end do !t

      !call closeFiles(isVertical)
    end subroutine decomp
end module mycmim
