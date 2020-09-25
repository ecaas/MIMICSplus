module initMod
  use paramMod
  use dispmodule
  implicit none
contains

  subroutine initialize_onelayer(initialC, initialN, Init_PlantC, Init_PlantN) !This subroutine sets the initial carbon values in the pool matrix (Carbon content in each pool in every depth level)
                                                   !/LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc/
    real(r8),intent(out) :: initialC(pool_types)
    real(r8),intent(out) :: initialN(pool_types+1)
    real(r8), intent(out) :: Init_PlantC, Init_PlantN

    !real(r8)             :: initialPlant
    initialC = (/500,500,500,500,500,500,500,500,500,500/)*10
    !initialC = (/636.174017, 	971.97288, 	11.217038, 	15.266405, 	12.743604, 	12.743604, 	12.743604, 	820.18018, 	1407.085933, 	90.696668/)
    !initialN = (/101.394007, 	154.913711, 	0.560852, 	0.477075, 	0.63718 ,	0.63718 ,	0.63718, 	36.699793, 	30.785334, 	4.007713, 	3.682967/)
    initialN(1:10) = initialC/CN_ratio
    initialN(11) = 50
    Init_PlantC =  	200*10 !1111.376721!
    Init_PlantN = 4*10!177.142764!
  end subroutine initialize_onelayer

  subroutine initialize_vert(InitC, InitN, Init_PlantC, Init_PlantN) !This subroutine sets the initial carbon values in the pool matrix (Carbon content in each pool in every depth level)
                      !/LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc/
    !
    real(r8),intent(out) :: InitC(nlevdecomp, pool_types), InitN(nlevdecomp, pool_types+1)
    real(r8), intent(out) :: Init_PlantC, Init_PlantN
    ! InitC = 0.0
    ! InitC(1,:) = (/50.0,50.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0/)
    ! InitC(2,:) = (/30.0,30.0,20.0,20.0,10.0,10.0,10.0,10.0,10.0,10.0/)
    ! InitC(3,:) = (/15.0,15.0,10.0,10.0,10.0,10.0,10.0,20.0,20.0,20.0/)
    ! InitC(4,:) = (/15.0,15.0,10.0,10.0,10.0,10.0,10.0,20.0,20.0,20.0/)
    ! InitC(5,:) = (/1.0, 10.0, 10.0,10.0,20.0,20.0,20.0,30.0,30.0,30.0/)
    ! InitC(6,:) = (/10.0,10.0,10.0,10.0, 20.0,20.0,20.0,30.0,30.0,30.0/)
    ! InitC(7,:) = (/10.0,10.0,10.0,10.0, 20.0,20.0,20.0,50.0,50.0,50.0/)
    !
    ! InitN(:,1:10) = InitC/5
    ! InitN(:,11) = 50 !inorganic
    ! Init_PlantC = 200
    ! Init_PlantN = 20
  !  call disp(CN_ratio)
  end subroutine initialize_vert

end module initMod
