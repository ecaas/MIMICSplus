program test
  use paramMod
  use testMod
  use initMod
  implicit none
  real(r8), dimension(nlevdecomp, pool_types) :: Carbon, Nitrogen
  call initialize_vert(Carbon, Nitrogen)
  call mass_conservation(Nitrogen, Carbon, Nitrogen)
end program test
