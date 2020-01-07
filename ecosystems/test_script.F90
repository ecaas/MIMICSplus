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

  integer :: N=1*365*24!Total number of timesteps (years,days, hours)

  character(len=12), dimension(3) :: str = (/ 'Heath ','Meadow', 'Shrub '/)
  call decomp(nsteps=N, run_name='Heath'//"_nctest", isVertical = .False., nlevdecomp = 1, ecosystem = 'Heath', step_frac=1)

   ! do i = 1,3
   !   print*, str(i)
   !   call decomp(nsteps=N, run_name=trim(str(i))//"_newtest", isVertical = .False., nlevdecomp = 1, ecosystem = trim(str(i)), step_frac=1)
   ! end do
end program test_script
