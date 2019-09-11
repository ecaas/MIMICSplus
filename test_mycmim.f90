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
    subroutine decomp(nsteps, run_name, isVertical, nlevdecomp) !This subroutine calculates the balance equation dC/dt for each pool at each time step based on the fluxes calculated in the same time step. Then update the pool sizes before moving on to
      !the next time step. It also calculates an analytical solution to the problem.

      logical                        :: isVertical           ! True if we use vertical soil layers.
      character (len=*)              :: run_name             ! used for naming outputfiles
      integer                        :: nlevdecomp           ! number of vertical layers

      integer                        :: nsteps               ! number of time steps to iterate over
      real(r8)                       :: dt=1/(24.0)          ! 1 hour = 1/24 days (size of time step)
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
      real(r8)                       :: vert(nlevdecomp, pool_types)
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
              veg_input=(/f_litm*f_lit_1, f_lits*f_lit_1/)*I_tot
            else
              veg_input=(/f_litm*f_lit_234, f_lits*f_lit_234/)*I_tot
            end if !j=1
          else
            veg_input=(/f_litm, f_lits/)*I_tot
          end if !isVertical

          !Calculate fluxes between pools in level j:
          call microbial_fluxes(j, pool_matrix,nlevdecomp)
          call som_fluxes(j, pool_matrix,nlevdecomp)
          call litter_fluxes(j, pool_matrix,nlevdecomp)

          !TODO: This writing to file should be made much more efficient, and to binary files, not text files..
          if (counter == 24) then
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

              Gain = LITtoSAP(1)*MGE(1) + LITtoSAP(3)*MGE(3) + sum(MYCtoSAP(1:3))
              Loss = sum(SAPtoSOM(1:3))/pool_matrix(j,i)

            elseif (i==4) then !SAPk
              change_matrix(j,i) = LITtoSAP(2)*MGE(2) + LITtoSAP(4)*MGE(4) &
               + sum(MYCtoSAP(4:6)) + SOMtoSAP(2)*MGE(2) - sum(SAPtoSOM(4:6))

               Gain= LITtoSAP(3)*MGE(2) + LITtoSAP(4)*MGE(4) + sum(MYCtoSAP(4:6))
               Loss=sum(SAPtoSOM(1:3))/pool_matrix(j,i)

            elseif (i==5) then !EcM
              change_matrix(j,i)=I_tot*f_ecm*f_myc_levels-MYCtoSAP(1)-MYCtoSAP(4)-MYCtoSOM(1)-MYCtoSOM(4)-MYCtoSOM(7)

              Gain = I_tot*f_ecm*f_myc_levels
              Loss = (MYCtoSAP(1)+MYCtoSAP(4)+MYCtoSOM(1)+MYCtoSOM(4)+MYCtoSOM(7))/pool_matrix(j,i)

            elseif (i==6) then !ErM
              change_matrix(j,i)=I_tot*f_erm*f_myc_levels-MYCtoSAP(2)-MYCtoSAP(5)-MYCtoSOM(2)-MYCtoSOM(5)-MYCtoSOM(8)

              Gain = I_tot*f_erm*f_myc_levels
              Loss = (MYCtoSAP(2)+MYCtoSAP(5)+MYCtoSOM(2)+MYCtoSOM(5)+MYCtoSOM(8))/pool_matrix(j,i)

            elseif (i==7) then !AM
              change_matrix(j,i)=I_tot*f_am*f_myc_levels-MYCtoSAP(3)-MYCtoSAP(6)-MYCtoSOM(3)-MYCtoSOM(6)-MYCtoSOM(9)

              Gain = I_tot*f_am*f_myc_levels
              Loss = (MYCtoSAP(3)+MYCtoSAP(6)+MYCtoSOM(3)+MYCtoSOM(6)+MYCtoSOM(9))/pool_matrix(j,i)

            elseif (i==8) then !SOMp
              change_matrix(j,i)=I_tot*f_som1*f_myc_levels + SAPtoSOM(1) + SAPtoSOM(4) + sum(MYCtoSOM(1:3))-SOMtoSOM(1)
              !Use the same partitioning between the depth levels as for mycorrhiza (f_myc_levels)
              Gain = I_tot*f_som1*f_myc_levels + SAPtoSOM(1) + SAPtoSOM(4) + sum(MYCtoSOM(1:3))
              Loss = SOMtoSOM(1)/pool_matrix(j,i)

            elseif (i==9) then !SOMa
              change_matrix(j,i)=SAPtoSOM(2) + SAPtoSOM(5) + sum(MYCtoSOM(4:6)) + &
               SOMtoSOM(1) + SOMtoSOM(2) - SOMtoSAP(1) - SOMtoSAP(2)

               Gain = SAPtoSOM(2) + SAPtoSOM(5) + sum(MYCtoSOM(4:6)) +  SOMtoSOM(1) + SOMtoSOM(2)
               Loss = (SOMtoSAP(1)*MGE(1) + SOMtoSAP(2)*MGE(2))/pool_matrix(j,i)

            elseif (i==10) then !SOMc
              change_matrix(j,i)=I_tot*f_som2*f_myc_levels + SAPtoSOM(3) + SAPtoSOM(6) + sum(MYCtoSOM(7:9))- SOMtoSOM(2)

              Gain = I_tot*f_som2*f_myc_levels + SAPtoSOM(3) + SAPtoSOM(6) + sum(MYCtoSOM(7:9))
              Loss = SOMtoSOM(2)/pool_matrix(j,i)

            else
              print*, 'Too many pool types expected, pool_types = ',pool_types
            end if !determine dC_i/dt

            if (isVertical) then
              call alt_vertical_diffusion(j,i, tot_diff,upper,lower, pool_matrix, nlevdecomp)
              vert(j,i) = tot_diff
              !Check if the diffusion term should act as a source or a sink in the analytical solution.
                if (tot_diff < 0) then
                  Loss = Loss + abs(tot_diff)
                else
                  Gain = Gain + abs(tot_diff)
                end if

                if (counter ==24) then
                  write(unit=10,fmt='(F10.0,A2,I2,A2,I6,A2,F30.10,A2,F30.10,A2,F30.10)') &
                  time,',',j,',',i,',',vert(j,i),',', upper,',', lower
                end if !writing
            end if!isVertical

            !Analytical solution:
            a_matrix(j, i) = Init(j,i)*exp(-time*Loss) + Gain*(1-exp(-time*Loss))/Loss

          end do !i, pool_types

          !Calculate the heterotrophic respiration loss from depth level j:
          HR(j) = LITtoSAP(1)*(1-MGE(1)) + LITtoSAP(2)*(1-MGE(2)) + LITtoSAP(3)*(1-MGE(3)) &
                        + LITtoSAP(4)*(1-MGE(4)) + SOMtoSAP(1)*(1-MGE(1)) + SOMtoSAP(2)*(1-MGE(2))
        end do !j, depth_level

        !Update pool sizes
        !Assuming fluxes are given as carbon transport pr. day, and the time step t is 1 day --> multiply with dt=1 will get the units correct (I guess this is what is done in the  mimics R code..)
        pool_matrix=pool_matrix + change_matrix*dt + tot_diff

        if (counter == 24) then
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
      call closeFiles(isVertical)
    end subroutine decomp
end module mycmim
