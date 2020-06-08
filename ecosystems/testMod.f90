module testing
  use paramMod

  implicit none

contains
  subroutine mass_conservation(in_layer_change, vert_change old, new) !Checking that mass is conserved during a time step
    real(r8), intent(in), dimension(nlevdecomp, pool_types) :: in_layer_change, vertical_change, old, new
  end subroutine mass_conservation

end module testing
