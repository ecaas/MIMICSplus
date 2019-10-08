program test_script
  !use TridiagonalMod , only: solve_tridiag
  !use patankar_diffusion
  !use paramOneLayerMod
  !use SoilDecompType, only : carbon_pool_type

  !use dispmodule
  use mycmim
  !use mycmimOneLayer
  implicit none
  integer :: i
  integer :: step = 24*60
  integer  :: N =20*365*step!Total number of timesteps

  character(len=12), dimension(3) :: str = (/'Heath ', 'Meadow', 'Shrub '/)
  !str = [character(len=) :: 'Heath', 'Meadow', 'Shrub' ]

  do i = 1,3
    print*, str(i)
    call decomp(nsteps=N, run_name=trim(str(i))//"5", isVertical = .False., nlevdecomp = 1, ecosystem = trim(str(i)), step_frac=)
  end do
end program test_script
