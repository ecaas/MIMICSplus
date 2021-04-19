module initMod
  use paramMod
  use dispmodule
  implicit none
contains

  subroutine initialize(pools_C, pools_N, Init_PlantC, Init_PlantN,nlevdecomp) !This subroutine sets the initial C and N values in the pool matrices (content in each soil pool in every depth level + plant)
    !Order: /LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc/
    integer :: nlevdecomp
    real(r8), intent(inout):: pools_C(:,:), pools_N(:,:)
    real(r8), intent(out)  :: Init_PlantC, Init_PlantN
    integer                :: j
    real(r8),parameter     :: C_tot=20000 !gC/m2
    pools_C=C_tot/(nlevdecomp*pool_types)
    pools_C(1,:)=pools_C(1,:)+50



    !TODO: Find better initial values!
    do j=1,nlevdecomp
      !pools_C(j,:) = C_tot/(delta_z(j)*nlevdecomp)
      pools_N(j,1:10) = pools_C(j,1:10)/CN_ratio
    end do
    Init_PlantC = 500
    Init_PlantN = 20
    pools_N(:,11) = 500

  end subroutine initialize


end module initMod
