program test_script
  !use TridiagonalMod , only: solve_tridiag
  !use patankar_diffusion
  !use paramOneLayerMod
  !use SoilDecompType, only : carbon_pool_type

  !use dispmodule
  use mycmim
  !use mycmimOneLayer
  implicit none
  integer :: i,dset_id,ierr
  !call store_parameters("Heath_newtest")
  integer :: N=2*365*24!Total number of timesteps (years,days, hours)
  !
  character(len=12), dimension(3) :: str = (/ 'Heath ','Meadow', 'Shrub '/)
  !call decomp(nsteps=N, run_name='Heath'//"_test", isVertical = .False., nlevdecomp = 1, ecosystem = 'Heath', step_frac=1)


  call decomp(nsteps=N, run_name=trim(str(1))//"_tautest", isVertical = .True., nlevdecomp = 4, ecosystem = trim(str(1)), step_frac=1)
  !call decomp(nsteps=N, run_name=trim(str(2))//"_tautest", isVertical = .False., nlevdecomp = 1, ecosystem = trim(str(2)), step_frac=1)
  !call decomp(nsteps=N, run_name=trim(str(3))//"_tautest", isVertical = .False., nlevdecomp = 1, ecosystem = trim(str(3)), step_frac=1)

end program test_script
