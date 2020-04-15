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
  integer :: N=1*365*24!Total number of timesteps (years,days, hours)
  !
  character(len=12), dimension(3) :: str = (/ 'Heath ','Meadow', 'Shrub '/)
  !do i = 1,3
  !  call decomp(nsteps=N, run_name=trim(str(i))//"_input", isVertical = .True., ecosystem = trim(str(i)), step_frac=1)
  !end do

   call decomp(nsteps=N, run_name=trim(str(3))//"_input", isVertical = .True., ecosystem = trim(str(3)), step_frac=1)
  ! call decomp(nsteps=N, run_name=trim(str(2))//"_tautest", isVertical = .True., ecosystem = trim(str(2)), step_frac=1)
  ! call decomp(nsteps=N, run_name=trim(str(3))//"_tautest", isVertical = .True., ecosystem = trim(str(3)), step_frac=1)

end program test_script
