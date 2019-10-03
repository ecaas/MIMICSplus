program test_script
  !use TridiagonalMod , only: solve_tridiag
  !use patankar_diffusion
  !use paramOneLayerMod
  !use SoilDecompType, only : carbon_pool_type

  !use dispmodule
  use mycmim
  !use mycmimOneLayer
  implicit none

  integer  :: N =50*365!*60 !Total number of timesteps

  call decomp(nsteps=N, run_name='Heath2', isVertical = .False., nlevdecomp = 1, ecosystem = 'Heath')

end program test_script
