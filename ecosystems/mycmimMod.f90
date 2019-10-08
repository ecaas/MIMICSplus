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
  !integer, parameter                            :: nlevdecomp!Number of soil depth layers (ndecomp_cascade_transitions=7,ndecomp_pools=10)


!Shape of the pool_matrix/change_matrix
!|       LITm LITs SAPr SAPk EcM ErM AM SOMp SOMa SOMc |
!|level1   1   2    3    4    5   6   7   8    9   10  |
!|level2                                               |
!| .                                                   |
!| .                                                   |
!|nlevdecomp __________________________________________|

  contains
    subroutine decomp(nsteps, run_name, isVertical, nlevdecomp, ecosystem,step_frac) !This subroutine calculates the balance equation dC/dt for each pool at each time step based on the fluxes calculated in the same time step. Then update the pool sizes before moving on to
      !the next time step. It also calculates an analytical solution to the problem.

      logical                        :: isVertical           ! True if we use vertical soil layers.
      character (len=*)              :: run_name             ! used for naming outputfiles
      character (len=*)              :: ecosystem             ! 'Shrub', 'Heath' or 'Meadow'
      integer                        :: nlevdecomp           ! number of vertical layers

      integer                        :: nsteps               ! number of time steps to iterate over
      integer,parameter              :: step_frac!*60

      real(r8)                       :: dt= 1.0/step_frac         ! 1 sec = 1/(24*60*60) days (size of time step)
      real(r8)                       :: time                 ! t*dt

      integer                        :: counter=0              ! used for determining when to output results
      integer                        :: j,i,t                ! for iterations

      real(r8)                       :: Loss, Gain           ! Source and sink term for solving the analytical solution to the dC/dt=Gain-Loss*concentration equation.
      real(r8)                       :: tot_diff,upper,lower ! For the call to alt_vertical_diffusion

      real(r8),dimension(nlevdecomp) :: HR                   ! For storing the C amount that is lost to respiration
      real(r8)                       :: pool_matrix(nlevdecomp,pool_types)   ! For storing C pool sizes [gC/m3]
      real(r8)                       :: change_matrix(nlevdecomp,pool_types) ! For storing dC/dt for each time step [gC/(m3*day)]
      real(r8)                       :: a_matrix(nlevdecomp, pool_types)     ! For storing the analytical solution
      real(r8)                       :: Init(nlevdecomp, pool_types)         ! Initial C concentration, determined in initMod.f90
      real(r8)                       :: pool_temporary(nlevdecomp,pool_types)                  !For storing the pool(t)+flux changes. These values are used to calculate the vertical transport.
      real(r8)                       :: vert(nlevdecomp, pool_types)

      real(kind=r8)                  :: GEP       ![gC/(m2*day)] Gross ecosystem productivity
      real(r8), dimension(3)         :: myc_input !vector giving the input from vegetation to mycorrhiza pools
      real(r8)                       :: I_tot     !Total input from vegetation to litter pools

      if (ecosystem == 'Heath') then
        GEP       = 0.281*24                          ![gC/(m2*day)] Gross ecosystem productivity
        myc_input = (/0.20,0.60,0.20/)*GEP*0.5/depth  !For Heath, most to ErM ![gC/m3/day]
        fMET      = 0.6
        k = (/0.05,0.05,0.005,0.05,0.05,0.005/)
        k2 = (/0.0007,0.0007,0.00014/)
      elseif (ecosystem == 'Meadow') then
        GEP       = 0.385*24                         ![gC/(m2*day)] Gross ecosystem productivity
        myc_input = (/0.1,0.1,0.80/)*GEP*0.5/depth   !For meadow, most to AM ![gC/m3/day]
        fMET      = 0.5
        k = (/0.05,0.05,0.005,0.05,0.05,0.005/)
        k2 = (/0.0007,0.0007,0.00014/)
      elseif (ecosystem == 'Shrub') then
        GEP       = 0.491*24                 ![gC/(m2*day)] Gross ecosystem productivity
        myc_input = (/0.80,0.10,0.10/)*GEP*0.5/depth  ![gC/m3/day] For shrub, most to EcM. 0.4 bc. 0.5 goes to litter and 0.1 directly to SOM (fsom1 og fsom2)
        fMET      = 0.001
        k = (/0.05,0.05,0.005,0.05,0.05,0.005/)*10
        k2 = (/0.0007,0.0007,0.00014/)*10
      else
        print*, 'Invalid ecosystem name', ecosystem
        stop
      end if
      I_tot = GEP*0.5/depth

      fPHYS = (/ 0.3 * exp(1.3*fCLAY), 0.2 * exp(0.8*fCLAY) /)
      fCHEM = (/ 0.1 * exp(-3.0*fMET), 0.3 * exp(-3*fMET) /)
      fAVAIL = 1-(fPHYS+fCHEM)
      tau = (/ 5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET) /)*24 ![1/h]*24=[1/day]Mimics include a modification, tau_mod, not included here.

      !Set initial concentration values in pool_matrix:
      if (isVertical) then
        call initialize_vert(Init, pool_matrix, nlevdecomp)
      else
        call initialize_onelayer(Init, pool_matrix)
      end if !isVertical

      !open and prepare files to store results
      call openOutputFile(run_name, isVertical)

      !The first line in is the initial condition
      do j=1,nlevdecomp
        write(unit = 1, FMT='(F10.0,A2,I2)',advance='no') 0.0,',', j
        write(unit = 15, FMT='(F10.0,A2,I2)',advance='no') 0.0,',', j
        write(unit = 2, FMT='(F10.0,A2,I2,A2,F30.10)',advance='no') 0.0,',', j,',',0.0
        do i = 1,pool_types
          write(unit = 1,FMT='(A,F30.10)',advance='no') ',' , pool_matrix(j,i)
          write(unit = 2, FMT='(A,F30.10)', advance='no') ',' ,change_matrix(j,i)
          write(unit = 15, FMT='(A,F30.10)', advance='no') ',' ,Init(j,i)
        end do
        write(1,*) ''
        write(2,*) ''
        write(15,*)''
      end do !write
      !----------------------------------------------------------------------------------------------------------------
      do t =1,nsteps
        time = t*dt
        counter =counter +1

        do j = 1, nlevdecomp !For each depth level (for the no vertical transport case, nlevdecomp = 1, so loop is only done once):

          if (isVertical) then
            if (j==1) then !The litter input is higher in the first depth level then the rest.
              veg_input=(/fMET*0.25, fMET*0.25/)*I_tot
            else
              veg_input=(/fMET*0.25, fMET*0.25/)*I_tot
            end if !j=1
          else
            veg_input=(/fMET, 1-fMET/)*(I_tot)/0.56 !average depth = 0.56 m -> input in gC/m³*day
          end if !isVertical

          !Calculate fluxes between pools in level j:
          call microbial_fluxes(j, pool_matrix,nlevdecomp)
          call som_fluxes(j, pool_matrix,nlevdecomp)
          call litter_fluxes(j, pool_matrix,nlevdecomp)

          ! !TODO: This writing to file should be made much more efficient, and to binary files, not text files..
          if (counter == step_frac) then
            write(unit=3,fmt='(F10.0,A2,I2,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10)') &
            time,',',j,',',LITtoSAP(1),',',LITtoSAP(2),',',LITtoSAP(3),',',LITtoSAP(4)
            write(unit=4,fmt='(F10.0,A2,I2,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10)') &
            time,',',j,',',SAPtoSOM(1),',',SAPtoSOM(2),',',SAPtoSOM(3),',',SAPtoSOM(4),',',SAPtoSOM(5),',',SAPtoSOM(6)
            write(unit=7,fmt='(F10.0,A2,I2,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10)') &
            time,',',j,',',MYCtoSAP(1),',',MYCtoSAP(2),',',MYCtoSAP(3),',',MYCtoSAP(4),',',MYCtoSAP(5),',',MYCtoSAP(6)
            write(unit=8,fmt='(F10.0,A2,I2,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10)') &
            time,',',j,',',MYCtoSOM(1),',',MYCtoSOM(2),',',MYCtoSOM(3),',' &
            ,MYCtoSOM(4),',',MYCtoSOM(5),',',MYCtoSOM(6),',',MYCtoSOM(7),',',MYCtoSOM(8),',',MYCtoSOM(9)
            write(unit=9,fmt='(F10.0,A2,I2,A2,F30.10,A2,F30.10,A2,F30.10,A2,F30.10)') &
            time,',',j,',',SOMtoSAP(1),',',SOMtoSAP(2),',',SOMtoSOM(1),',',SOMtoSOM(2)
          end if !writing

          do i = 1, pool_types !loop over all the pool types, i, in depth level j

            !This if-loop calculates dC/dt for the different carbon pools.NOTE: If pools are added/removed (i.e the actual model equations is changed), this loop needs to be updated.
            !The Gain and Loss variables are used to calculate the analytical solution to dC/dt=Gain - Loss*C, a_matrix(j,i)
            !NOTE: The "change_matrix" values correspond to the equations A11-A17 in Wieder 2015
            if (i==1) then !LITm
              change_matrix(j,i) = veg_input(1)-sum(LITtoSAP(1:2))

              Gain = veg_input(1)
              Loss=sum(LITtoSAP(1:2))/pool_matrix(j, i)

            elseif (i==2) then !LITs
              change_matrix(j,i) =  veg_input(2) -sum(LITtoSAP(3:4))

              Gain = veg_input(2)
              Loss=sum(LITtoSAP(3:4))/pool_matrix(j, i)

            elseif (i==3) then !SAPr
              change_matrix(j,i) = LITtoSAP(1)*MGE(1) + LITtoSAP(3)*MGE(3) &
              + sum(MYCtoSAP(1:3)) + SOMtoSAP(1)*MGE(1) - sum(SAPtoSOM(1:3))

              Gain = LITtoSAP(1)*MGE(1) + LITtoSAP(3)*MGE(3) + sum(MYCtoSAP(1:3))+ SOMtoSAP(1)*MGE(1)
              Loss = sum(SAPtoSOM(1:3))/pool_matrix(j,i)

            elseif (i==4) then !SAPk
              change_matrix(j,i) = LITtoSAP(2)*MGE(2) + LITtoSAP(4)*MGE(4) &
               + sum(MYCtoSAP(4:6)) + SOMtoSAP(2)*MGE(2) - sum(SAPtoSOM(4:6))

               Gain= LITtoSAP(2)*MGE(2) + LITtoSAP(4)*MGE(4) + sum(MYCtoSAP(4:6))+ SOMtoSAP(2)*MGE(2)
               Loss=sum(SAPtoSOM(4:6))/pool_matrix(j,i)

            elseif (i==5) then !EcM
              change_matrix(j,i)=myc_input(1)-MYCtoSAP(1)-MYCtoSAP(4)-sum(MYCtoSOM(1:3))

              Gain = myc_input(1)
              Loss = (MYCtoSAP(1)+MYCtoSAP(4)+sum(MYCtoSOM(1:3)))/pool_matrix(j,i)

            elseif (i==6) then !ErM
              change_matrix(j,i)=myc_input(2)-MYCtoSAP(2)-MYCtoSAP(5)-sum(MYCtoSOM(4:6))

              Gain = myc_input(2)
              Loss = (MYCtoSAP(2)+MYCtoSAP(5)+sum(MYCtoSOM(4:6)))/pool_matrix(j,i)

            elseif (i==7) then !AM
              change_matrix(j,i)=myc_input(3)-MYCtoSAP(3)-MYCtoSAP(6)-sum(MYCtoSOM(7:9))

              Gain = myc_input(3)
              Loss = (MYCtoSAP(3)+MYCtoSAP(6)+sum(MYCtoSOM(7:9)))/pool_matrix(j,i)

            elseif (i==8) then !SOMp
              change_matrix(j,i)=  I_tot*f_som1*f_myc_levels+ SAPtoSOM(1) + SAPtoSOM(4) + MYCtoSOM(1) + MYCtoSOM(4) + MYCtoSOM(7)-SOMtoSOM(1)
              !Use the same partitioning between the depth levels as for mycorrhiza (f_myc_levels)
              Gain = I_tot*f_som1*f_myc_levels + SAPtoSOM(1) + SAPtoSOM(4) + MYCtoSOM(1) + MYCtoSOM(4) + MYCtoSOM(7)
              Loss = SOMtoSOM(1)/pool_matrix(j,i)

            elseif (i==9) then !SOMa
              change_matrix(j,i)=SAPtoSOM(2) + SAPtoSOM(5) +  MYCtoSOM(2) + MYCtoSOM(5) + MYCtoSOM(8) + &
               SOMtoSOM(1) + SOMtoSOM(2) - SOMtoSAP(1) - SOMtoSAP(2)

               Gain = SAPtoSOM(2) + SAPtoSOM(5) +  MYCtoSOM(2) + MYCtoSOM(5) + MYCtoSOM(8) +  SOMtoSOM(1) + SOMtoSOM(2)
               Loss = (SOMtoSAP(1) + SOMtoSAP(2))/pool_matrix(j,i)

            elseif (i==10) then !SOMc
              change_matrix(j,i)=I_tot*f_som2*f_myc_levels+SAPtoSOM(3) + SAPtoSOM(6) + MYCtoSOM(3) + MYCtoSOM(6) + MYCtoSOM(9)- SOMtoSOM(2)

              Gain = I_tot*f_som2*f_myc_levels + SAPtoSOM(3) + SAPtoSOM(6) + MYCtoSOM(3) + MYCtoSOM(6) + MYCtoSOM(9)
              Loss = SOMtoSOM(2)/pool_matrix(j,i)

            else
              print*, 'Too many pool types expected, pool_types = ',pool_types
            end if !determine dC_i/dt


            pool_temporary(j,i)=pool_matrix(j,i) + change_matrix(j,i)*dt

            if (pool_temporary(j,i) < 0.0) then
              print*, 'Negative concentration value at t',t,'depth level',j,'pool number',i
            end if

            a_matrix(j,i) = Init(j,i)*exp(-time*Loss) + Gain*(1-exp(-time*Loss))/Loss !+ vert*dt
          end do !i, pool_types

          !Calculate the heterotrophic respiration loss from depth level j in timestep t:
          HR(j) =( LITtoSAP(1)*(1-MGE(1)) + LITtoSAP(2)*(1-MGE(2)) + LITtoSAP(3)*(1-MGE(3)) &
                        + LITtoSAP(4)*(1-MGE(4)) + SOMtoSAP(1)*(1-MGE(1)) + SOMtoSAP(2)*(1-MGE(2)))*dt
        end do !j, depth_level

        if (isVertical) then
          call vertical_diffusion(tot_diff,upper,lower, pool_temporary, nlevdecomp,vert,time, counter,dt)
          pool_matrix =   pool_temporary+vert*dt
          a_matrix=a_matrix+vert*dt
        else
          pool_matrix=pool_temporary
        end if!isVertical

        !Update pool sizes

        !Check if the diffusion term should act as a source or a sink in the analytical solution.
        ! if (tot_diff < 0) then
        !   Loss = Loss + abs(tot_diff)
        ! else
        !   Gain = Gain + abs(tot_diff)
        ! end if

        !call disp('pool_matrix = ', pool_matrix)
        if (counter == step_frac) then
          counter = 0
          !Write results to file. TODO: incorporate this into the above j, i loops
          do j=1,nlevdecomp
            write(unit = 1, FMT='(F10.0,A2,I2)',advance='no') time,',', j
            write(unit = 15, FMT='(F10.0,A2,I2)',advance='no') time,',', j
            write(unit = 2, FMT='(F10.0,A2,I2,A2,F30.10)',advance='no') time,',', j,',',HR(j)
            do i = 1,pool_types
              write(unit = 1,FMT='(A,F30.10)',advance='no') ',' , pool_matrix(j,i)
              write(unit = 2, FMT='(A,F30.10)', advance='no') ',' ,change_matrix(j,i)
              write(unit = 15, FMT='(A,F30.10)', advance='no') ',' ,a_matrix(j,i)
            end do
            write(1,*) ''
            write(2,*) ''
            write(15,*)''
          end do !write
        end if!writing
      end do !t
      print*, I_tot*f_som2*f_myc_levels
      print*, veg_input
      print*, myc_input
      call closeFiles(isVertical)
    end subroutine decomp
end module mycmim
