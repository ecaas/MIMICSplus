program test
  use paramMod
  use testMod
  use dispmodule !External module to pretty print matrices (mainly for testing purposes)

  implicit none
  
  integer, parameter :: nlevdecomp = 3
  real(r8), dimension(nlevdecomp,3) ::  array=reshape( (/ 1/0.02, 1/0.02, 1/0.02, &
                                                          1/0.04, 1/0.04, 1/0.04, &
                                                          1/0.06, 1/0.06, 1/0.06/), shape(array), order=(/2,1/))
  real(r8), dimension(nlevdecomp,3) :: mass_array

  real(r8),dimension(nlevdecomp)                          :: input=10.0
  real(r8),dimension(nlevdecomp)                          :: respiration=1.0
  

  call cons_to_mass(array, mass_array,nlevdecomp, 3) !Tested & works
!  call tot_mass_input(lit, som, myc, tot_input,nlevdecomp,3)
  !call mass_conservation_pool(in_layer_change, vert_change ,old, new,diff,nlevdecomp)
  ! call respired_mass(HR_conc_vr, HR_mass_tot,nlevdecomp)
  call test_mass_conservation(15.0d0, 7.5d0, array, array*2,nlevdecomp,3)
  call total_mass_conservation(100d0, 91d0, array, array*2,nlevdecomp,3)

end program test
