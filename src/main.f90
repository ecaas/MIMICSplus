program main
  use mycmim
  use paramMod
  implicit none

  integer :: levels
  integer,parameter:: Cpools=10
  integer,parameter:: Npools=11
  character (len=*),parameter :: site = "32374_Saltdal" !"31767_Kongsvinger" "32374_Saltdal" "32087_Dovre"
  character (len=*),parameter :: description = "highN_all_4"
  character (len=*),parameter :: clm_data_file = '/home/ecaas/nird/'//site//'_hist_for_decomp/lnd/hist/'//site//'_hist_for_decomp.clm2.'
  character (len=*),parameter :: clm_surface_file = '/home/ecaas/nird/surface_files/'//site//'/surfdata_'//site//'_simyr2000.nc'
  
  
  real(r8), dimension (:,:), allocatable   :: C_matrix_init    ! For storing C pool sizes [gC/m3]
  real(r8), dimension (:,:), allocatable   :: N_matrix_init   ! For storing N pool sizes [gN/m3]
  
  real(r8), dimension (:,:), allocatable   :: C_matrix_Spunup     ! For storing C pool sizes [gC/m3]
  real(r8), dimension (:,:), allocatable   :: N_matrix_Spunup     ! For storing N pool sizes [gN/m3]
    
  real(r8), dimension (:,:), allocatable   :: C_matrix_1970    
  real(r8), dimension (:,:), allocatable   :: N_matrix_1970   
  
  real(r8), dimension (:,:), allocatable   :: C_matrix_1987    
  real(r8), dimension (:,:), allocatable   :: N_matrix_1987    
  
  real(r8), dimension (:,:), allocatable   :: C_matrix_final   
  real(r8), dimension (:,:), allocatable   :: N_matrix_final   

  call read_nlayers(trim(clm_data_file//'all.1901.nc'), levels)
  
  !ALLOCATE:
  allocate(C_matrix_init(levels,Cpools),N_matrix_init(levels,Npools))
  allocate(C_matrix_Spunup(levels,Cpools),N_matrix_Spunup(levels,Npools))
  allocate(C_matrix_1970(levels,Cpools),N_matrix_1970(levels,Npools))
  allocate(C_matrix_1987(levels,Cpools),N_matrix_1987(levels,Npools))
  allocate(C_matrix_final(levels,Cpools),N_matrix_final(levels,Npools))

  !1: INITIALIZE
  call initialize(C_matrix_init,N_matrix_init,levels)
  !2: SPINUP
  call decomp(nsteps=5000*24*365, run_name=trim(site//"_"//description//"_"//"Spunup"),nlevdecomp=levels, step_frac=1, write_hour=24*365*100,&
  pool_C_start=C_matrix_init,pool_N_start=N_matrix_init, pool_C_final=C_matrix_Spunup,pool_N_final=N_matrix_Spunup,&
  start_year=1850,stop_year=1870,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file)
  !|
  !| Use output of last timestep to initialize step 2
  ! 
  !3: RUN WITH MONTHLY UPDATES of inputdata:
  call decomp(nsteps=71*24*365, run_name=trim(site//"_"//description//"_"//"to1970"),nlevdecomp=levels, step_frac=1, write_hour=1*24*365, &
  pool_C_start=C_matrix_Spunup,pool_N_start=N_matrix_Spunup, pool_C_final=C_matrix_1970,pool_N_final=N_matrix_1970,&
  start_year=1900,stop_year=1970,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file)
  !|
  !| Use output of last timestep to initialize step 4
  !V
  !4: RUN WITH DAILY UPDATES of inputdata, output every 20 year:
  call decomp(nsteps=17*24*365, run_name=trim(site//"_"//description//"_"//"to1987"),nlevdecomp=levels, step_frac=1, write_hour=1*24*365, &
  pool_C_start=C_matrix_1970,pool_N_start=N_matrix_1970, pool_C_final=C_matrix_1987,pool_N_final=N_matrix_1987,&
  start_year=1971,stop_year=1987,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file)
  
  !5: RUN WITH DAILY UPDATES of inputdata,output every day:
  call decomp(nsteps=5*24*365, run_name=trim(site//"_"//description//"_"//"to1992"),nlevdecomp=levels, step_frac=1, write_hour=1*24, &
  pool_C_start=C_matrix_1987,pool_N_start=N_matrix_1987, pool_C_final=C_matrix_final,pool_N_final=N_matrix_final,&
  start_year=1988,stop_year=1992,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file)
  ! 
  ! 
  deallocate (C_matrix_init,N_matrix_init,C_matrix_Spunup,N_matrix_Spunup)
  deallocate (C_matrix_1970,N_matrix_1970,C_matrix_1987,N_matrix_1987)
  deallocate (C_matrix_final,N_matrix_final)
  

end program main
