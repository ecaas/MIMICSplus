module initMod
  use paramMod
  implicit none
contains

  subroutine initialize_onelayer(initial,pools) !This subroutine sets the initial carbon values in the pool matrix (Carbon content in each pool in every depth level)
                                                   !/LITm LITs SAPr SAPk EcM ErM AM SOMp SOMa SOMc/
    real(r8),intent(out) :: initial(1, pool_types)
    real(r8), intent(inout) :: pools(1,pool_types)
    initial(1,:) = (/500,500,5,5,5,5,5,5,5/)
    pools = initial
  end subroutine initialize_onelayer

  subroutine initialize_vert(Init, pools, levels) !This subroutine sets the initial carbon values in the pool matrix (Carbon content in each pool in every depth level)
                      !/LITm LITs SAPr SAPk EcM ErM AM SOMp SOMa SOMc/
            !Total:   / 220, 230, 140, 100, 100,120,80,240, 210, 220/
    integer :: levels
    real(r8),intent(out) :: Init(levels, pool_types)
    real(r8), intent(inout) :: pools(levels, pool_types)
    Init(1,:) = (/100,100,50,40,40,40,60,60,60/)
    Init(2,:) = 0.9*Init(1,:)
    Init(3,:) = (/15,15,20,30,30,30,100,100,100/)
    Init(4,:) = 1.1*Init(3,:)

    pools=Init

  end subroutine initialize_vert

end module initMod
