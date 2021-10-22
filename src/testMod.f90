module testMod
  use paramMod
  use dispmodule
  implicit none

    integer                        :: i,j
  
contains

  subroutine cons_to_mass(matrix, mass_matrix,nlevdecomp,pool_types) !Convert from g/m3 to g/m2 !Tested and works
    !INPUT
    integer, intent(in)                                    :: nlevdecomp
    integer, intent(in)                                    :: pool_types
    real(r8),intent(in), dimension(nlevdecomp, pool_types) :: matrix
    
    !OUTPUT
    real(r8),intent(out),dimension(nlevdecomp, pool_types) :: mass_matrix


    do i = 1,pool_types
      mass_matrix(:,i) = matrix(:,i)*delta_z(1:nlevdecomp)
    end do
  end subroutine cons_to_mass

  ! subroutine tot_mass_input(change_matrix, tot_input,nlevdecomp,no_of) !Merge input arrays and calculate the total mass of the input in the different layers
  !   !INPUT
  !   integer,intent(in)  ::  nlevdecomp
  !   real(r8), intent(in) :: change_matrix()
  ! 
  !   !OUTPUT
  !   real(r8), intent(out), dimension(nlevdecomp, nlevdecomp) :: tot_input
  ! 
  !   !LOCAL
  ! 
  ! 
  !   do i = 1,nlevdecomp
  !     tot_input(:,i) = tot_input(:,i)*delta_z(i)
  !   end do
  ! !  call disp('mass_tot_input', tot_input)
  ! !  print*, 'sum of input: ', sum(tot_input)
  ! end subroutine tot_mass_input
  
  
  subroutine mass_conservation_pool(in_layer_change, vert_change ,old, new,diff,nlevdecomp) !Checking that mass is conserved during a time step in the vertical layer NOTE: Not finished!!
    
    !INPUT
    real(r8), intent(in), dimension(nlevdecomp, pool_types) :: in_layer_change
    real(r8), intent(in), dimension(nlevdecomp, pool_types) :: vert_change
    real(r8), intent(in), dimension(nlevdecomp, pool_types) :: old
    real(r8), intent(in), dimension(nlevdecomp, pool_types) :: new
    integer, intent(in)                                     :: nlevdecomp

    
    !OUTPUT
    real(r8), intent(out) :: diff
    
    !LOCAL
    real(r8), dimension(nlevdecomp, pool_types) :: mass_old
    real(r8), dimension(nlevdecomp, pool_types) :: mass_new
    real(r8), dimension(nlevdecomp, pool_types) :: mass_change
    real(r8), dimension(nlevdecomp, pool_types) :: mass_vert_change
    real(r8)                                    :: sum_old
    real(r8)                                    :: sum_new
    real(r8)                                    :: sum_change
    real(r8)                                    :: sum_vert_change
    
  
    call cons_to_mass(old, mass_old,nlevdecomp,pool_types)
    call cons_to_mass(new, mass_new,nlevdecomp,pool_types)
    call cons_to_mass(in_layer_change, mass_change,nlevdecomp,pool_types)
    call cons_to_mass(vert_change, mass_vert_change,nlevdecomp,pool_types)
  
    sum_old = sum(mass_old)
    sum_new = sum(mass_new)
    sum_change = sum(mass_change)
    sum_vert_change = sum(mass_vert_change)

    diff = sum_new - (sum_old + sum_change + sum_vert_change)
  
    print*, diff, sum_new
  end subroutine mass_conservation_pool

  subroutine respired_mass(HR_conc_vr, HR_mass_tot,nlevdecomp) !Calculates the total mass of carbon that is respired (in a time step)    
    !INPUT
    real(r8), intent(in), dimension(nlevdecomp) :: HR_conc_vr ![gC/m3]
    integer                                     :: nlevdecomp
    
    !OUTPUT
    real(r8), intent(out)             :: HR_mass_tot ![gC/m2]
    
    !LOCAL
    real(r8), dimension(nlevdecomp)   :: HR_mass_vr

    do i = 1,nlevdecomp
      HR_mass_vr(i) = HR_conc_vr(i)*delta_z(i)
    end do
    HR_mass_tot = sum(HR_mass_vr)
  end subroutine respired_mass

  subroutine test_mass_conservation(mass_input, mass_respiration, old, new,nlevdecomp,no_of_pools) !Checking that mass is conserved during a time step. Tested and works
    !INPUT
    real(r8),intent(in), dimension(nlevdecomp, no_of_pools) :: old !g/m3
    real(r8),intent(in), dimension(nlevdecomp, no_of_pools) :: new !g/m3
    real(r8),intent(in)                                    :: mass_respiration !g/m2
    real(r8),intent(in)                                    :: mass_input !g/m2
    integer, intent(in)                                    :: nlevdecomp
    integer, intent(in)                                    :: no_of_pools

    !OUTPUT
    
    !LOCAL
    real(r8) :: sum_old
    real(r8) :: sum_new
    real(r8) :: diff
    real(r8) :: goal

    real(r8), dimension(nlevdecomp, no_of_pools) :: mass_old, mass_new
    call cons_to_mass(old, mass_old,nlevdecomp,no_of_pools)
    call cons_to_mass(new, mass_new,nlevdecomp,no_of_pools)

    sum_old = sum(mass_old)
    sum_new = sum(mass_new)

    goal = sum_old + mass_input - mass_respiration
    diff = sum_new-goal
    
    if (abs(diff) > 1e-4) then
      print*, '-----------------------------------------------------------'
      print*, '(mass at t+dt) - (mass at t + Input - respiration): ', diff
      print*, 'sum of respired mass: ', mass_respiration
      print*, 'sum of input        : ', mass_input
      print*, 'mass before timestep: ', sum_old
      print*, 'mass after timestep : ', sum_new
      print*, '-----------------------------------------------------------'
    end if
  end subroutine test_mass_conservation

  subroutine total_mass_conservation(sum_input, sum_respiration, old, new,nlevdecomp,no_of_pools) !Checking that mass is conserved over the whole time period. tested and works.
    !INPUT 
    integer , intent(in)                                    :: nlevdecomp
    integer , intent(in)                                    :: no_of_pools
    real(r8),intent(in)                                     :: sum_respiration
    real(r8),intent(in)                                     :: sum_input    
    real(r8), intent(in), dimension(nlevdecomp, no_of_pools):: old
    real(r8), intent(in), dimension(nlevdecomp, no_of_pools):: new
    
    !OUTPUT
    
    !LOCAL
    real(r8), dimension(nlevdecomp, no_of_pools):: mass_old 
    real(r8), dimension(nlevdecomp, no_of_pools):: mass_new
    real(r8)                                    :: sum_old
    real(r8)                                    :: sum_new 
    real(r8)                                    :: diff
    real(r8)                                    :: goal
    

    call cons_to_mass(old, mass_old,nlevdecomp,no_of_pools)
    call cons_to_mass(new, mass_new,nlevdecomp,no_of_pools)

    sum_old = sum(mass_old)
    sum_new = sum(mass_new)

    goal = sum_old + sum_input - sum_respiration
    diff = sum_new-goal

    if (abs(diff) > 1e-4) then
      print*, 'Mass conservation is NOT fulfilled: '

    else
      print*, 'Mass conservation is fulfilled (yay!): '
    end if
    print*, 'diff                : ', diff
    print*, 'sum of respired mass: ', sum_respiration
    print*, 'sum of input        : ', sum_input
    print*, 'sum_old             : ', sum_old
    print*, 'sum_new             : ', sum_new
  end subroutine total_mass_conservation

  subroutine total_nitrogen_conservation(sum_input, sum_out, old, new,nlevdecomp,no_of_pools) !Checking that mass is conserved over the whole time period. tested and works.
    !INPUT 
    integer , intent(in)                                    :: nlevdecomp
    integer , intent(in)                                    :: no_of_pools
    real(r8),intent(in)                                     :: sum_out
    real(r8),intent(in)                                     :: sum_input    
    real(r8), intent(in), dimension(nlevdecomp, no_of_pools):: old
    real(r8), intent(in), dimension(nlevdecomp, no_of_pools):: new
    
    !OUTPUT
    
    !LOCAL
    real(r8), dimension(nlevdecomp, no_of_pools):: mass_old 
    real(r8), dimension(nlevdecomp, no_of_pools):: mass_new
    real(r8)                                    :: sum_old
    real(r8)                                    :: sum_new 
    real(r8)                                    :: diff
    real(r8)                                    :: goal
    
    print*, sum_input,sum_out
    call cons_to_mass(old, mass_old,nlevdecomp,no_of_pools)
    call cons_to_mass(new, mass_new,nlevdecomp,no_of_pools)

    sum_old = sum(mass_old)
    sum_new = sum(mass_new)

    goal = sum_old + sum_input - sum_out
    diff = sum_new-goal

    if (abs(diff) > 1e-4) then
      print*, 'Mass conservation of N is NOT fulfilled: '

    else
      print*, 'Mass conservation of N is fulfilled (yay!): '
    end if
    print*, 'diff                : ', diff
    print*, 'sum N out of system : ', sum_out
    print*, 'sum of input        : ', sum_input
    print*, 'sum_old             : ', sum_old
    print*, 'sum_new             : ', sum_new
  end subroutine total_nitrogen_conservation

end module testMod
