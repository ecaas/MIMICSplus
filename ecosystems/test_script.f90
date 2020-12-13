program test_script
  !use TridiagonalMod , only: solve_tridiag
  !use patankar_diffusion
  !use paramOneLayerMod
  !use SoilDecompType, only : carbon_pool_type

  !use dispmodule
  use mycmim
  !use mycmimOneLayer
  implicit none

  !call store_parameters("Heath_newtest")
  integer :: N=100*365*24!Total number of timesteps (years,days, hours) *10 when step_frac 10
  !real(r8), dimension(nlevdecomp) :: TSOIL
  !real(r8), dimension(nlevdecomp) ::  SOILLIQ
  !real(r8), dimension(nlevdecomp) ::  SOILICE
  !real(r8), dimension(nlevdecomp) :: WATSAT
  !call read_clmdata('/home/ecaas/clm/cruncep_iso_hist/Dovre/clm50_clm50d001_1deg_CRUNCEPV7_iso_hist.clm2.h0.SOILLIQ_SOILICE_TSOI.185001-201412.nc_Dovre.nc', TSOIL, SOILLIQ,SOILICE,WATSAT, 100)
  character(len=12), dimension(3) :: str = (/ 'Heath ','Meadow', 'Shrub '/)
  !do i = 1,3
  !  call decomp(nsteps=N, run_name=trim(str(i))//"_input", isVertical = .True., ecosystem = trim(str(i)), step_frac=1)
  !end do


  call decomp(nsteps=N, run_name="high_dep", isVertical = .False., ecosystem = trim(str(1)), step_frac=1)


end program test_script