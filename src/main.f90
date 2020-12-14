program main
  use mycmim

  implicit none

  integer :: N=100*365*24!Total number of timesteps (years,days, hours) *10 when step_frac 10
  integer,parameter:: levels=1

  ! call decomp(nsteps=N, run_name="make_test2",nlevdecomp=levels, step_frac=1)

end program main
