module initMod
  use paramMod
  use dispmodule
  implicit none
contains

  subroutine initialize(pools_C, pools_N,nlevdecomp) !This subroutine sets the initial C and N values in the pool matrices (content in each soil pool in every depth level + plant)
    
    !INPUT
    
    !OUTPUT
    real(r8), intent(out):: pools_C(:,:), pools_N(:,:)
    
    !INOUT

    !LOCAL
    integer :: nlevdecomp
    integer                :: j

    do j=1,nlevdecomp
       pools_C(j,:) = (/1000.,1000.,10.,50.,100.,1.,1.,1500.,1500.,500./)
       pools_N(j,1:10) = pools_C(j,1:10)/CN_ratio
     end do
    
    pools_N(:,11) = 1.
    pools_N(:,12) = 1.
    

  end subroutine initialize


end module initMod
