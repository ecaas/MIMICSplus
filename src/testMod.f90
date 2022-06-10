module testMod
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use paramMod, only : delta_z,pool_types,pool_types_N,inorg_N_pools
  use dispmodule, only: disp
  use initMod, only: nlevels
  implicit none

    integer                        :: i,j

contains

  subroutine cons_to_mass(matrix, mass_matrix,layers,pools) !Convert from g/m3 to g/m2 !Tested and works
    !INPUT
    integer, intent(in)                                    :: layers
    integer, intent(in)                                    :: pools
    real(r8),intent(in), dimension(layers, pools) :: matrix
    
    !OUTPUT
    real(r8),intent(out),dimension(layers, pools) :: mass_matrix
    do i = 1,pools
      mass_matrix(:,i) = matrix(:,i)*delta_z(1:layers)
    end do
  end subroutine cons_to_mass

  
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

  subroutine respired_mass(HR_conc_vr, HR_mass_tot) !Calculates the total mass of carbon that is respired (in a time step)    
    !INPUT
    real(r8), intent(in), dimension(nlevels) :: HR_conc_vr ![gC/m3]
    
    !OUTPUT
    real(r8), intent(out)             :: HR_mass_tot ![gC/m2]
    
    !LOCAL
    real(r8), dimension(nlevels)   :: HR_mass_vr
    do i = 1,nlevels
      HR_mass_vr(i) = HR_conc_vr(i)*delta_z(i)
    end do
    HR_mass_tot = sum(HR_mass_vr)
  end subroutine respired_mass

  subroutine test_mass_conservation_C(mass_input, mass_out, old, new,nlevdecomp) !Checking that mass is conserved during a time step. Tested and works
    !INPUT
    real(r8),intent(in), dimension(nlevdecomp, pool_types) :: old !g/m3
    real(r8),intent(in), dimension(nlevdecomp, pool_types) :: new !g/m3
    real(r8),intent(in)                                    :: mass_out !g/m2
    real(r8),intent(in)                                    :: mass_input !g/m2
    integer, intent(in)                                    :: nlevdecomp

    !OUTPUT
    
    !LOCAL
    real(r8) :: sum_old
    real(r8) :: sum_new
    real(r8) :: diff
    real(r8) :: goal
    real(r8), dimension(nlevdecomp, pool_types) :: mass_old
    real(r8), dimension(nlevdecomp, pool_types) :: mass_new

    call cons_to_mass(old, mass_old,nlevdecomp,pool_types)
    call cons_to_mass(new, mass_new,nlevdecomp,pool_types)

    sum_old = sum(mass_old)
    sum_new = sum(mass_new)

    goal = sum_old + mass_input - mass_out
    diff = sum_new-goal
    
    if (abs(diff) > 1e-10) then  
      print*, "Mass NOT conserved for CARBON"
      print*, '-----------------------------------------------------------'
      print*, '(mass at t+dt) - (mass at t + Input - out): ', diff
      print*, 'sum of mass out of system: ', mass_out
      print*, 'sum of input             : ', mass_input
      print*, 'mass before timestep     : ', sum_old
      print*, 'mass after timestep      : ', sum_new
      print*, '-----------------------------------------------------------'
    end if
  end subroutine test_mass_conservation_C
  
  subroutine test_mass_conservation_N(mass_input, mass_out, &
                                      old, old_inorg, &
                                      new,new_inorg,nlevdecomp) !Checking that mass is conserved during a time step. Tested and works
    !INPUT
    real(r8),intent(in), dimension(nlevdecomp, pool_types_N) :: old !g/m3
    real(r8),intent(in), dimension(nlevdecomp, pool_types_N) :: new !g/m3
    real(r8),intent(in), dimension(nlevdecomp, inorg_N_pools) :: old_inorg !g/m3
    real(r8),intent(in), dimension(nlevdecomp, inorg_N_pools) :: new_inorg !g/m3
    real(r8),intent(in)                                    :: mass_out !g/m2
    real(r8),intent(in)                                    :: mass_input !g/m2
    integer, intent(in)                                    :: nlevdecomp

    !OUTPUT
    
    !LOCAL
    real(r8) :: sum_old
    real(r8) :: sum_new
    real(r8) :: diff
    real(r8) :: goal
    real(r8), dimension(nlevdecomp, pool_types_N) :: mass_old, mass_new
    real(r8), dimension(nlevdecomp, inorg_N_pools) :: mass_old_inorg, mass_new_inorg

    call cons_to_mass(old, mass_old,nlevdecomp,pool_types_N)
    call cons_to_mass(new, mass_new,nlevdecomp,pool_types_N)
    
    call cons_to_mass(old_inorg,mass_old_inorg,nlevdecomp,inorg_N_pools)
    call cons_to_mass(new_inorg,mass_new_inorg,nlevdecomp,inorg_N_pools)

    sum_old = sum(mass_old) + sum(mass_old_inorg)
    sum_new = sum(mass_new) + sum(mass_new_inorg)

    goal = sum_old + mass_input - mass_out
    diff = sum_new-goal
    
    if (abs(diff) > 1e-10) then
      print*, "Mass NOT conserved for NITROGEN"
      print*, '-----------------------------------------------------------'
      print*, '(mass at t+dt) - (mass at t + Input - out): ', diff
      print*, 'sum of mass out of system: ', mass_out
      print*, 'sum of input             : ', mass_input
      print*, 'mass before timestep     : ', sum_old
      print*, 'mass after timestep      : ', sum_new
      print*, '-----------------------------------------------------------'
    end if
  end subroutine test_mass_conservation_N

  subroutine total_carbon_conservation(sum_input, sum_respiration, old, new,nlevdecomp) !Checking that mass is conserved over the whole time period. tested and works.
    !INPUT 
    integer , intent(in)                                    :: nlevdecomp
    real(r8),intent(in)                                     :: sum_respiration
    real(r8),intent(in)                                     :: sum_input    
    real(r8), intent(in), dimension(nlevdecomp, pool_types):: old
    real(r8), intent(in), dimension(nlevdecomp, pool_types):: new
    
    !OUTPUT
    
    !LOCAL
    real(r8), dimension(nlevdecomp, pool_types):: mass_old 
    real(r8), dimension(nlevdecomp, pool_types):: mass_new
    real(r8)                                    :: sum_old
    real(r8)                                    :: sum_new 
    real(r8)                                    :: diff
    real(r8)                                    :: goal
    
    call cons_to_mass(old, mass_old,nlevdecomp,pool_types)
    call cons_to_mass(new, mass_new,nlevdecomp,pool_types)

    sum_old = sum(mass_old)
    sum_new = sum(mass_new)

    goal = sum_old + sum_input - sum_respiration
    diff = sum_new-goal

    if (abs(diff) > 1e-5) then
      print*, 'Mass conservation of C is NOT fulfilled: '

    else
      print*, 'Mass conservation of C is fulfilled (yay!): '
    end if
    print*, 'diff                : ', diff
    print*, 'sum of respired mass: ', sum_respiration
    print*, 'sum of input        : ', sum_input
    print*, 'sum_old             : ', sum_old
    print*, 'sum_new             : ', sum_new
  end subroutine total_carbon_conservation

  subroutine total_nitrogen_conservation(sum_input, sum_out, &
                                        old, old_inorg, &
                                        new,new_inorg, nlevdecomp) !Checking that mass is conserved over the whole time period. tested and works.
    !INPUT 
    integer , intent(in)                                      :: nlevdecomp
    real(r8),intent(in)                                       :: sum_out
    real(r8),intent(in)                                       :: sum_input    
    real(r8), intent(in), dimension(nlevdecomp, pool_types_N) :: old
    real(r8), intent(in), dimension(nlevdecomp, inorg_N_pools):: old_inorg
    real(r8), intent(in), dimension(nlevdecomp, pool_types_N) :: new
    real(r8), intent(in), dimension(nlevdecomp, inorg_N_pools):: new_inorg
    
    !OUTPUT
    
    !LOCAL
    real(r8), dimension(nlevdecomp, pool_types_N):: mass_old 
    real(r8), dimension(nlevdecomp, pool_types_N):: mass_new
    real(r8), dimension(nlevdecomp, inorg_N_pools) :: mass_old_inorg, mass_new_inorg
    real(r8)                                    :: sum_old
    real(r8)                                    :: sum_new 
    real(r8)                                    :: diff
    real(r8)                                    :: goal
    
    call cons_to_mass(old, mass_old,nlevdecomp,pool_types_N)
    call cons_to_mass(new, mass_new,nlevdecomp,pool_types_N)

    call cons_to_mass(old_inorg,mass_old_inorg,nlevdecomp,inorg_N_pools)
    call cons_to_mass(new_inorg,mass_new_inorg,nlevdecomp,inorg_N_pools)

    sum_old = sum(mass_old) + sum(mass_old_inorg)
    sum_new = sum(mass_new) + sum(mass_new_inorg)

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
