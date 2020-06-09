module initMod
  use paramMod
  implicit none
contains

  subroutine initialize_onelayer(initialC, initialN) !This subroutine sets the initial carbon values in the pool matrix (Carbon content in each pool in every depth level)
                                                   !/LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc/
    real(r8),intent(out) :: initialC(1, pool_types)
    real(r8),intent(out) :: initialN(1, pool_types)
    initialC(1,:) = (/500,500,5,5,5,5,5,5,5,5/)
    initialN(1,:) = (/500,500,5,5,5,5,5,5,5,5/) !TODO Multiply with N:C ratio.

  end subroutine initialize_onelayer

  subroutine initialize_vert(InitC, InitN) !This subroutine sets the initial carbon values in the pool matrix (Carbon content in each pool in every depth level)
                      !/LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc/

    integer :: levels
    real(r8),intent(out) :: InitC(nlevdecomp, pool_types), InitN(nlevdecomp, pool_types)

    InitC = 0.0
    InitC(1,:) = (/50.0,50.0,10.0,10.0, 10.0,10.0,10.0,10.0,10.0,10.0/)
    InitC(2,:) = (/30.0,30.0,20.0,20.0, 10.0,10.0,10.0,10.0,10.0,10.0/)
    InitC(3,:) = (/15.0,15.0,10.0,10.0, 10.0,10.0,10.0,20.0,20.0,20.0/)
    InitC(4,:) = (/15.0,15.0,10.0,10.0, 10.0,10.0,10.0,20.0,20.0,20.0/)
    InitC(5,:) = (/0.0, 0.0, 10.0,10.0, 20.0,20.0,20.0,30.0,30.0,30.0/)
    InitC(6,:) = (/0.0,0.0,10.0,10.0,   20.0,20.0,20.0,30.0,30.0,30.0/)
    InitC(7,:) = (/0.0,0.0,10.0,10.0,   20.0,20.0,20.0,50.0,50.0,50.0/)

    InitN = InitC !Multiply with N:C ratio

  end subroutine initialize_vert

end module initMod
