program main
  use mycmim

  implicit none


  integer :: N=1*24*365!Total number of timesteps (years,days, hours) *10 when step_frac 10
  integer,parameter:: levels=10
  call decomp(nsteps=N, run_name="test_mmk",nlevdecomp=levels, step_frac=1, write_hour=1*24*30)


end program main
