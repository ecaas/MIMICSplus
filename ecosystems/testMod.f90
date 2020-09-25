module testMod
  use paramMod
  use dispmodule
  implicit none

  integer :: i


contains

  subroutine cons_to_mass(matrix, mass_matrix) !Convert from g/m3 to g/m2
    real(r8), intent(in), dimension(nlevdecomp, pool_types) :: matrix
    real(r8), intent(out), dimension(nlevdecomp, pool_types) :: mass_matrix
    do i = 1,pool_types
      mass_matrix(:,i) = matrix(:,i)*delta_z
    end do
  end subroutine cons_to_mass

  subroutine tot_mass_input(lit, som, myc, tot_input) !Merge input arrays and calculate the total mass of the input in the different layers

    real(r8), intent(out), dimension(nlevdecomp, 7) :: tot_input
    real(r8)                       :: lit(nlevdecomp,no_of_litter_pools)![gC/(m3 h)] Fraction of litter input to LITm and LITs, respectively
    real(r8)                       :: myc(nlevdecomp,no_of_myc_pools)   ![gC/(m3 h)] vector giving the input from vegetation to mycorrhiza pools
    real(r8)                       :: som(nlevdecomp,no_of_som_pools-1) ![gC/(m3 h)] !only input to SOMp and SOMc
    tot_input(:, 1:2) = lit
    tot_input(:, 3:4) = som
    tot_input(:, 5:7) = myc
    !call disp('tot_input', tot_input)
    do i = 1,7
      tot_input(:,i) = tot_input(:,i)*delta_z
    end do
  !  call disp('mass_tot_input', tot_input)
  !  print*, 'sum of input: ', sum(tot_input)
  end subroutine tot_mass_input
  ! subroutine mass_conservation_pool(in_layer_change, vert_change old, new) !Checking that mass is conserved during a time step in the vertical layer NOTE: Not finished!!
  !   real(r8), intent(in), dimension(nlevdecomp, pool_types) :: in_layer_change, vert_change old, new
  !   real(r8), dimension(nlevdecomp, pool_types) :: mass_old, mass_new, mass_change
  !   real(r8) :: sum_old, sum_new, sum_change, diff
  !
  !   call cons_to_mass(old, mass_old)
  !   call cons_to_mass(new, mass_new)
  !   call cons_to_mass(in_layer_change, mass_change)
  !   call cons_to_mass(vert_change, mass_vert_change)
  !
  !   sum_old = sum(mass_old)
  !   sum_new = sum(mass_new)
  !   sum_change = sum(mass_change)
  !   sum_vert_change = sum(mass_vert_change)
  !   !call disp(new)
  !   diff = sum_new - (sum_old + sum_change + mass_vert_change)
  !
  !   print*, diff, sum_new
  ! end subroutine mass_conservation_pool

  subroutine respired_mass(HR_conc_vr, HR_mass_tot) !Calculates the total mass of carbon that is respired (in a time step)
    real(r8), intent(in), dimension(nlevdecomp) :: HR_conc_vr
    real(r8), intent(out) :: HR_mass_tot
    real(r8), dimension(nlevdecomp):: HR_mass_vr

    HR_mass_vr = HR_conc_vr*delta_z


    HR_mass_tot = sum(HR_mass_vr)

  end subroutine respired_mass

  subroutine test_mass_conservation(mass_input, mass_respiration, old, new) !Checking that mass is conserved during a time step
    real(r8), intent(in), dimension(nlevdecomp, pool_types) :: old, new
    real(r8), intent(in), dimension(nlevdecomp)   :: mass_respiration
    real(r8), intent(in), dimension(nlevdecomp,7)   :: mass_input
    real(r8), dimension(nlevdecomp, pool_types) :: mass_old, mass_new
    real(r8) :: sum_old, sum_new, sum_respiration, sum_input, diff, goal

    call cons_to_mass(old, mass_old)
    call cons_to_mass(new, mass_new)

    sum_old = sum(mass_old)
    sum_new = sum(mass_new)
    sum_respiration = sum(mass_respiration)
    sum_input = sum(mass_input)

    goal = sum_old + sum_input - sum_respiration
    diff = sum_new-goal
    print*, '-----------------------------------------------------------'
    print*, '(mass at t+dt) - (mass at t + Input - respiration): ', diff
    print*, 'sum of respired mass: ', sum_respiration
    print*, 'sum of input: ', sum_input
    print*, 'diff: ', diff, 'HR-Input: ', sum_respiration - sum_input
    print*, 'diff-HR-Input: ', diff - (sum_respiration - sum_input)
  end subroutine test_mass_conservation

  subroutine total_mass_conservation(sum_input, sum_respiration, old, new) !Checking that mass is conserved over the whole time period
    real(r8), intent(in), dimension(nlevdecomp, pool_types) :: old, new
    real(r8), dimension(nlevdecomp, pool_types) :: mass_old, mass_new
    real(r8) :: sum_old, sum_new, sum_respiration, sum_input, diff, goal

    call cons_to_mass(old, mass_old)
    call cons_to_mass(new, mass_new)

    sum_old = sum(mass_old)
    sum_new = sum(mass_new)

    goal = sum_old + sum_input - sum_respiration
    diff = sum_new-goal

    if (abs(diff) > 1e-4) then
      print*, 'Mass conservation is not fulfilled: '
      print*, 'diff: ', diff
      print*, 'sum of respired mass: ', sum_respiration
      print*, 'sum of input: ', sum_input
      print*, 'sum_old: ', sum_old, 'sum_new: ', sum_new
    end if
  end subroutine total_mass_conservation


end module testMod
