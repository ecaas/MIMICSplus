program main
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use mycmimMod,  only: decomp
  use initMod, only: read_nlayers,nlevels, initialize
  use paramMod, only: full_clock_rate,full_clock_stop,full_clock_start, &
                      pool_types, pool_types_N,inorg_N_pools, read_some_parameters, &
                      use_ROI, use_ENZ, use_Sulman, dt
  implicit none

  !character (len=*),parameter :: site = "32087_Dovre" !"31767_Kongsvinger" "32374_Saltdal" "32087_Dovre"
  real(r8), dimension (:,:), allocatable   :: C_matrix_init    ! For storing C pool sizes [gC/m3]
  real(r8), dimension (:,:), allocatable   :: N_matrix_init   ! For storing N pool sizes [gN/m3]
  real(r8), dimension (:,:), allocatable   :: N_inorg_matrix_init   ! For storing N pool sizes [gN/m3]
  

  real(r8), dimension (:,:), allocatable   :: C_matrix_Spunup     ! For storing C pool sizes [gC/m3]
  real(r8), dimension (:,:), allocatable   :: N_matrix_Spunup     ! For storing N pool sizes [gN/m3]
  real(r8), dimension (:,:), allocatable   :: N_inorg_matrix_Spunup     ! For storing N pool sizes [gN/m3]
  

  real(r8), dimension (:,:), allocatable   :: C_matrix_1970    
  real(r8), dimension (:,:), allocatable   :: N_matrix_1970 
  real(r8), dimension (:,:), allocatable   :: N_inorg_matrix_1970   
    

  real(r8), dimension (:,:), allocatable   :: C_matrix_1987    
  real(r8), dimension (:,:), allocatable   :: N_matrix_1987    
  real(r8), dimension (:,:), allocatable   :: N_inorg_matrix_1987    

  real(r8), dimension (:,:), allocatable   :: C_matrix_final   
  real(r8), dimension (:,:), allocatable   :: N_matrix_final   
  real(r8), dimension (:,:), allocatable   :: N_inorg_matrix_final   

  character (len=200)  :: clm_data_file 
  character (len=200)  :: clm_surface_file
  character (len=200)  :: namelist_file
  character (len=200)  :: description
  character (len=200)  :: site

  call system_clock(count_rate=full_clock_rate) !Find the time rate
  call system_clock(count=full_clock_start)     !Start Timer 
  call GET_COMMAND_ARGUMENT(number=1,value=site)
  call GET_COMMAND_ARGUMENT(number=2,value=description)
  call GET_COMMAND_ARGUMENT(number=3,value=clm_data_file)
  call GET_COMMAND_ARGUMENT(number=4,value=clm_surface_file)
  call GET_COMMAND_ARGUMENT(number=5,value=namelist_file)
  
  call read_nlayers(trim(adjustr(clm_data_file)//'_historical.clm2.all.1901.nc'))
  call read_some_parameters(trim(namelist_file),use_ROI, use_Sulman, use_ENZ, dt)
  print*, use_ROI,use_Sulman,use_ENZ, dt
    
  !ALLOCATE:
  allocate(C_matrix_init(nlevels,pool_types),N_matrix_init(nlevels,pool_types_N),N_inorg_matrix_init(nlevels,inorg_N_pools))
  allocate(C_matrix_Spunup(nlevels,pool_types),N_matrix_Spunup(nlevels,pool_types_N),N_inorg_matrix_Spunup(nlevels,inorg_N_pools))
  allocate(C_matrix_1970(nlevels,pool_types),N_matrix_1970(nlevels,pool_types_N),N_inorg_matrix_1970(nlevels,inorg_N_pools))
  allocate(C_matrix_1987(nlevels,pool_types),N_matrix_1987(nlevels,pool_types_N),N_inorg_matrix_1987(nlevels,inorg_N_pools))
  allocate(C_matrix_final(nlevels,pool_types),N_matrix_final(nlevels,pool_types_N),N_inorg_matrix_final(nlevels,inorg_N_pools))

  
!   !1: INITIALIZE

  call initialize(C_matrix_init,N_matrix_init,N_inorg_matrix_init)

  !2: SPINUP
  call decomp(nsteps=20*24*365, &
              run_name=trim(trim(site)//"_"//trim(description)//"_"//"Spunup"), &
              write_hour=24*365*1+100*24,&
              pool_C_start=C_matrix_init,pool_N_start=N_matrix_init,inorg_N_start=N_inorg_matrix_init,&
              pool_C_final=C_matrix_Spunup,pool_N_final=N_matrix_Spunup,inorg_N_final=N_inorg_matrix_Spunup,&
              start_year=1850,stop_year=1869,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file)
 ! ! |
 ! ! | Use output of last timestep to initialize step 2
 ! ! ! V
 !  !!3: RUN WITH MONTHLY UPDATES of inputdata:
 !  call decomp(nsteps=71*24*365, &
 !              run_name=trim(trim(site)//"_"//trim(description)//"_"//"to1970"), &
 !              write_hour=1*24*365*10+1+100*24,&
 !              pool_C_start=C_matrix_Spunup,pool_N_start=N_matrix_Spunup,inorg_N_start=N_inorg_matrix_Spunup,&
 !              pool_C_final=C_matrix_1970,pool_N_final=N_matrix_1970,inorg_N_final=N_inorg_matrix_1970,&
 !              start_year=1900,stop_year=1970,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file)
 !  !|
 !  !| Use output of last timestep to initialize step 4
 !  !V
 !  !4: RUN WITH DAILY UPDATES of inputdata, output every year:
 !  call decomp(nsteps=17*24*365, &
 !              run_name=trim(trim(site)//"_"//trim(description)//"_"//"to1987"), &
 !              write_hour=1*24*365,&
 !              pool_C_start=C_matrix_1970,pool_N_start=N_matrix_1970,inorg_N_start=N_inorg_matrix_1970,&
 !              pool_C_final=C_matrix_1987,pool_N_final=N_matrix_1987,inorg_N_final=N_inorg_matrix_1987,&
 !              start_year=1971,stop_year=1987,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file)
 ! 
 !  !5: RUN WITH DAILY UPDATES of inputdata,output every day:
 !  call decomp(nsteps=5*24*365, &
 !              run_name=trim(trim(site)//"_"//trim(description)//"_"//"to1992"), &
 !              write_hour=1*24,&
 !              pool_C_start=C_matrix_1987,pool_N_start=N_matrix_1987,inorg_N_start=N_inorg_matrix_1987,&
 !              pool_C_final=C_matrix_final,pool_N_final=N_matrix_final,inorg_N_final=N_inorg_matrix_final,&
 !              start_year=1988,stop_year=1992,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file)

  deallocate (C_matrix_init,N_matrix_init,N_inorg_matrix_init,C_matrix_Spunup,N_matrix_Spunup,N_inorg_matrix_Spunup)
  deallocate (C_matrix_1970,N_matrix_1970,N_inorg_matrix_1970,C_matrix_1987,N_matrix_1987,N_inorg_matrix_1987)
  deallocate (C_matrix_final,N_matrix_final,N_inorg_matrix_final)
  call system_clock(count=full_clock_stop)      ! Stop Timer
  print*, "Total time for simulation in minutes: ", (real(full_clock_stop-full_clock_start)/real(full_clock_rate))/60

end program main









