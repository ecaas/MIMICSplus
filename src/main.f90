program main
  use mycmim
  use paramMod
  implicit none


  integer :: N=1*24*365!Total number of timesteps (years,days, hours) *10 when step_frac 10
  integer,parameter:: levels=10
  integer,parameter:: Cpools=10
  integer,parameter:: Npools=11
  
  real(r8)                       :: C_matrix_init(levels,Cpools)     ! For storing C pool sizes [gC/m3]
  real(r8)                       :: N_matrix_init(levels,Npools)     ! For storing C pool sizes [gC/m3]
  
  real(r8)                       :: C_matrix_Spunup(levels,Cpools)     ! For storing C pool sizes [gC/m3]
  real(r8)                       :: N_matrix_Spunup(levels,Npools)     ! For storing C pool sizes [gC/m3]
    
  real(r8)                       :: C_matrix_1970(levels,Cpools)     ! For storing C pool sizes [gC/m3]
  real(r8)                       :: N_matrix_1970(levels,Npools)     ! For storing C pool sizes [gC/m3]
  
  real(r8)                       :: C_matrix_final(levels,Cpools)     ! For storing C pool sizes [gC/m3]
  real(r8)                       :: N_matrix_final(levels,Npools)     ! For storing C pool sizes [gC/m3]

  
  !1: INITIALIZE
  call initialize(C_matrix_init,N_matrix_init,levels)
  !2: SPINUP
  call decomp(nsteps=500*24*365, run_name="Kongsvinger_Spinup",nlevdecomp=levels, step_frac=1, write_hour=1*24*5*365,&
  pool_C_start=C_matrix_init,pool_N_start=N_matrix_init, pool_C_final=C_matrix_Spunup,pool_N_final=N_matrix_Spunup,start_year=1850,stop_year=1870)
  !|
  !| Use output of last timestep to initialize step 2
  !V
  !3: RUN WITH MONTHLY UPDATES of inputdata:
  call decomp(nsteps=70*24*365, run_name="Kongsvinger_1900-1970",nlevdecomp=levels, step_frac=1, write_hour=1*24*20*365, &
  pool_C_start=C_matrix_Spunup,pool_N_start=N_matrix_Spunup, pool_C_final=C_matrix_1970,pool_N_final=N_matrix_1970,start_year=1900,stop_year=1970)
  !|
  !| Use output of last timestep to initialize step 4
  !V
  !4: RUN WITH DAILY UPDATES of inputdata, output every year:
  call decomp(nsteps=17*24*365, run_name="Kongsvinger_1971-87",nlevdecomp=levels, step_frac=1, write_hour=1*24*10*365, &
  pool_C_start=C_matrix_1970,pool_N_start=N_matrix_1970, pool_C_final=C_matrix_final,pool_N_final=N_matrix_final,start_year=1971,stop_year=1987)
  
  !5: RUN WITH DAILY UPDATES of inputdata,output every day:
  call decomp(nsteps=5*24*365, run_name="Kongsvinger_1988-1992",nlevdecomp=levels, step_frac=1, write_hour=1*24, &
  pool_C_start=C_matrix_final,pool_N_start=N_matrix_final, pool_C_final=C_matrix_final,pool_N_final=N_matrix_final,start_year=1988,stop_year=1992)
  ! 
  ! 
  ! 

end program main
