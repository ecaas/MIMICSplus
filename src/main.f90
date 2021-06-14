program main
  use mycmim

  implicit none


  integer :: N=600*365*24!Total number of timesteps (years,days, hours) *10 when step_frac 10
  integer,parameter:: levels=10

  call decomp(nsteps=N, run_name="default_desorption_tv80_includeEcM_mod4",nlevdecomp=levels, step_frac=1)


end program main
