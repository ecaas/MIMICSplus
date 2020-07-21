module initMod
  use paramMod
  use dispmodule
  implicit none
contains

  subroutine initialize_onelayer(initialC, initialN) !This subroutine sets the initial carbon values in the pool matrix (Carbon content in each pool in every depth level)
                                                   !/LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc/
    real(r8),intent(out) :: initialC(1, pool_types)
    real(r8),intent(out) :: initialN(1, pool_types)
    !real(r8)             :: initialPlant
    initialC(1,:) = (/500,500,5,5,5,5,5,5,5,5/)
    initialN(1,:) = (/500,500,5,5,5,5,5,5,5,5/)/CN_ratio
  end subroutine initialize_onelayer

  subroutine initialize_vert(InitC, InitN, Init_PlantC, Init_PlantN) !This subroutine sets the initial carbon values in the pool matrix (Carbon content in each pool in every depth level)
                      !/LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc/

    real(r8),intent(out) :: InitC(nlevdecomp, pool_types), InitN(nlevdecomp, pool_types+1)
    real(r8), intent(out) :: Init_PlantC, Init_PlantN
    InitC = 0.0
    InitC(1,:) = (/50.0,50.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0/)
    InitC(2,:) = (/30.0,30.0,20.0,20.0,10.0,10.0,10.0,10.0,10.0,10.0/)
    InitC(3,:) = (/15.0,15.0,10.0,10.0,10.0,10.0,10.0,20.0,20.0,20.0/)
    InitC(4,:) = (/15.0,15.0,10.0,10.0,10.0,10.0,10.0,20.0,20.0,20.0/)
    InitC(5,:) = (/1.0, 10.0, 10.0,10.0,20.0,20.0,20.0,30.0,30.0,30.0/)
    InitC(6,:) = (/10.0,10.0,10.0,10.0, 20.0,20.0,20.0,30.0,30.0,30.0/)
    InitC(7,:) = (/10.0,10.0,10.0,10.0, 20.0,20.0,20.0,50.0,50.0,50.0/)

    InitN(:,1:10) = InitC/5
    InitN(:,11) = 5
    Init_PlantC = 100
    Init_PlantN = 30
  !  call disp(CN_ratio)
  end subroutine initialize_vert

end module initMod
