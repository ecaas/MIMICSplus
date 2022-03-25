program main
  use mycmim
  use paramMod
  implicit none

  integer :: levels,name

  !character (len=*),parameter :: site = "32087_Dovre" !"31767_Kongsvinger" "32374_Saltdal" "32087_Dovre"
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

  character (len=*),parameter :: description = "test"
  character (len=200)  :: clm_data_file 
  character (len=200)  :: clm_surface_file
  do name = 1,1, 1    
    call system_clock(count_rate=full_clock_rate) !Find the time rate
    call system_clock(count=full_clock_start)     !Start Timer 
    clm_data_file='/home/ecaas/nird/'//trim(site_names(name))//'_historical/lnd/hist/'//trim(site_names(name))
    print*, trim(adjustr(clm_data_file)//'.clm2.all.1901.nc')
    print*, trim(adjustr(clm_data_file))
    clm_surface_file = '/home/ecaas/nird/surface_files/'//trim(site_names(name))//'/surfdata_'//trim(site_names(name))//'_simyr2000.nc'
    call read_nlayers(trim(adjustr(clm_data_file)//'_historical.clm2.all.1901.nc'), levels)

    !ALLOCATE:
    allocate(C_matrix_init(levels,pool_types),N_matrix_init(levels,pool_types_N))
    allocate(C_matrix_Spunup(levels,pool_types),N_matrix_Spunup(levels,pool_types_N))
    allocate(C_matrix_1970(levels,pool_types),N_matrix_1970(levels,pool_types_N))
    allocate(C_matrix_1987(levels,pool_types),N_matrix_1987(levels,pool_types_N))
    allocate(C_matrix_final(levels,pool_types),N_matrix_final(levels,pool_types_N))

    !1: INITIALIZE
    call initialize(C_matrix_init,N_matrix_init,levels)
    !2: SPINUP
    call decomp(nsteps=1500*24*365, run_name=trim(trim(site_names(name))//"_"//description//"_"//"Spunup"),nlevdecomp=levels, step_frac=1, write_hour=200+24*365*100,&
    pool_C_start=C_matrix_init,pool_N_start=N_matrix_init, pool_C_final=C_matrix_Spunup,pool_N_final=N_matrix_Spunup,&
    start_year=1850,stop_year=1869,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file,use_ROI=.False.,use_Sulman=.False., use_ENZ=.False.)
   !|
   !!| Use output of last timestep to initialize step 2
  !
    !!3: RUN WITH MONTHLY UPDATES of inputdata:
    call decomp(nsteps=71*24*365, run_name=trim(trim(site_names(name))//"_"//description//"_"//"to1970"),nlevdecomp=levels, step_frac=1, write_hour=1*24*365*10, &
    pool_C_start=C_matrix_Spunup,pool_N_start=N_matrix_Spunup, pool_C_final=C_matrix_1970,pool_N_final=N_matrix_1970,&
    start_year=1900,stop_year=1970,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file,use_ROI=.False.,use_Sulman=.False., use_ENZ=.False.)
    !|
    !| Use output of last timestep to initialize step 4
    !V
    !4: RUN WITH DAILY UPDATES of inputdata, output every year:
    call decomp(nsteps=17*24*365, run_name=trim(trim(site_names(name))//"_"//description//"_"//"to1987"),nlevdecomp=levels, step_frac=1, write_hour=1*24*365, &
    pool_C_start=C_matrix_1970,pool_N_start=N_matrix_1970, pool_C_final=C_matrix_1987,pool_N_final=N_matrix_1987,&
    start_year=1971,stop_year=1987,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file,use_ROI=.False.,use_Sulman=.False., use_ENZ=.False.)
    ! 
    !5: RUN WITH DAILY UPDATES of inputdata,output every day:
    call decomp(nsteps=5*24*365, run_name=trim(trim(site_names(name))//"_"//description//"_"//"to1992"),nlevdecomp=levels, step_frac=1, write_hour=1*24, &
    pool_C_start=C_matrix_1987,pool_N_start=N_matrix_1987, pool_C_final=C_matrix_final,pool_N_final=N_matrix_final,&
    start_year=1988,stop_year=1992,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file,use_ROI=.False.,use_Sulman=.False., use_ENZ=.False.)
  ! 

    deallocate (C_matrix_init,N_matrix_init,C_matrix_Spunup,N_matrix_Spunup)
    deallocate (C_matrix_1970,N_matrix_1970,C_matrix_1987,N_matrix_1987)
    deallocate (C_matrix_final,N_matrix_final)
  call system_clock(count=full_clock_stop)      ! Stop Timer
  print*, "Total time for simulation: ", real(full_clock_stop-full_clock_start)/real(full_clock_rate)
  print*, "Total time for simulation in minutes: ", (real(full_clock_stop-full_clock_start)/real(full_clock_rate))/60
end do
end program main









