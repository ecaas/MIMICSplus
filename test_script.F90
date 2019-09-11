program test_script
  !use TridiagonalMod , only: solve_tridiag
  !use patankar_diffusion
  !use paramOneLayerMod
  !use SoilDecompType, only : carbon_pool_type
  !use paramMod
  !use dispmodule
  use mycmim
  !use mycmimOneLayer
  implicit none
  !real(r8) :: dt =1/24.0!timestep
  integer  :: N =365*24 !Total number of timesteps
  !real :: test(3)=(/1.0,3.9,3.4/)
  !real,dimension(4)                :: node_z = (/0.05,0.15,0.5,1.10/)
  !print*, test, size(test)
  !call decomp_onelayer(nsteps=N)
  !call decomp(nsteps=N)
  !call disp('pool_matrix = ', pool_matrix)
  !call disp('change_matrix =', change_matrix)
  !print*, size(pool_matrix), size(change_matrix)
  call decomp(nsteps=N, run_name='vert', isVertical = .True., nlevdecomp = 4)
!   logical   :: isVertical !True if we use vertical soil layers.
!   isVertical = .True.
!   if (isVertical) then
!     call test
!   end if
!
! contains
!   subroutine test
!     print*, 'yo'
!   end subroutine test
end program test_script
