program main
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mycmimMod,    only: decomp
  use initMod,      only: read_nlayers,nlevels, initialize
  use paramMod,     only: full_clock_rate,full_clock_stop,full_clock_start, &
                          pool_types, pool_types_N,inorg_N_pools, &
                          use_ROI, use_ENZ, use_Sulman, dt
  use readMod,      only: read_namelist  
  use dispmodule,   only: disp !External module to pretty print matrices (mainly for testing purposes)
                      
  implicit none

  real(r8), dimension (:,:), allocatable   :: C_matrix_init   ! For storing C pool sizes [gC/m3]
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

  real(r8), dimension (:,:), allocatable   :: C_matrix_2009
  real(r8), dimension (:,:), allocatable   :: N_matrix_2009    
  real(r8), dimension (:,:), allocatable   :: N_inorg_matrix_2009   
  
  real(r8), dimension (:,:), allocatable   :: C_matrix_final   
  real(r8), dimension (:,:), allocatable   :: N_matrix_final   
  real(r8), dimension (:,:), allocatable   :: N_inorg_matrix_final   

  character (len=200)  :: clm_data_file 
  character (len=200)  :: clm_surface_file
  character (len=200)  :: namelist_file
  character (len=200)  :: output_path
  character (len=200)  :: description
  character (len=200)  :: site
  character (len=20)   :: spinup_char
  integer              :: spinup_years
  character (len=20)   :: spinup_only_char
  logical              :: spinup_only
  
  call system_clock(count_rate=full_clock_rate) !Find the time rate
  call system_clock(count=full_clock_start)     !Start Timer 
  
  !Information from run.bash script: 
  call GET_COMMAND_ARGUMENT(number=1,value=site)
  call GET_COMMAND_ARGUMENT(number=2,value=description)
  call GET_COMMAND_ARGUMENT(number=3,value=clm_data_file)
  call GET_COMMAND_ARGUMENT(number=4,value=clm_surface_file)
  call GET_COMMAND_ARGUMENT(number=5,value=namelist_file)
  call get_command_argument(number=6,value=output_path)
  call get_command_argument(number=7,value=spinup_char)
  call get_command_argument(number=8,value=spinup_only_char)
  spinup_char=trim(spinup_char)
  read(spinup_char,*) spinup_years 
  read(spinup_only_char,*) spinup_only

  !Get number of active soil layers from CLM file:
  call read_nlayers(trim(adjustr(clm_data_file)//'all.1901.nc'))

  !Read namelist options:
  call read_namelist(trim(namelist_file),use_ROI, use_Sulman, use_ENZ, dt)
  print*, use_ROI,use_Sulman,use_ENZ, dt,nlevels, spinup_years
  
  !ALLOCATE:
  allocate(C_matrix_init(nlevels,pool_types),N_matrix_init(nlevels,pool_types_N),N_inorg_matrix_init(nlevels,inorg_N_pools))
  allocate(C_matrix_Spunup(nlevels,pool_types),N_matrix_Spunup(nlevels,pool_types_N),N_inorg_matrix_Spunup(nlevels,inorg_N_pools))
  
  !1: INITIALIZE
  call initialize(C_matrix_init,N_matrix_init,N_inorg_matrix_init)
  print*,trim(trim(site)//"_"//trim(description)//"_"//"Spunup")

  !2: SPINUP
  call decomp(nsteps=spinup_years*24*365, &
  run_name=trim(trim(site)//"_"//trim(description)//"_"//"Spunup"), &
  write_hour=24*365*100+100*24,&
  pool_C_start=C_matrix_init,pool_N_start=N_matrix_init,inorg_N_start=N_inorg_matrix_init,&
  pool_C_final=C_matrix_Spunup,pool_N_final=N_matrix_Spunup,inorg_N_final=N_inorg_matrix_Spunup,&
  start_year=1850,stop_year=1869,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file, out_path = output_path)
  
  
  if ( .not. spinup_only ) then
    allocate(C_matrix_1970(nlevels,pool_types),N_matrix_1970(nlevels,pool_types_N),N_inorg_matrix_1970(nlevels,inorg_N_pools))
    allocate(C_matrix_2009(nlevels,pool_types),N_matrix_2009(nlevels,pool_types_N),N_inorg_matrix_2009(nlevels,inorg_N_pools))
    allocate(C_matrix_final(nlevels,pool_types),N_matrix_final(nlevels,pool_types_N),N_inorg_matrix_final(nlevels,inorg_N_pools))  

    ! | Use output of last timestep of spinup to initialize step 2
    ! ! V
    !!3: MONTHLY inputdata:
    call decomp(nsteps=71*24*365, &
                run_name=trim(trim(site)//"_"//trim(description)//"_"//"to1970"), &
                write_hour=1*24*365*10+1+100*24,&
                pool_C_start=C_matrix_Spunup,pool_N_start=N_matrix_Spunup,inorg_N_start=N_inorg_matrix_Spunup,&
                pool_C_final=C_matrix_1970,pool_N_final=N_matrix_1970,inorg_N_final=N_inorg_matrix_1970,&
                start_year=1900,stop_year=1970,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file, out_path = output_path)
    !|
    !| Use output of last timestep to initialize step 4
    !V
    !4: DAILY inputdata, output every year:
    call decomp(nsteps=38*24*365, &
                run_name=trim(trim(site)//"_"//trim(description)//"_"//"to2010"), &
                write_hour=1*24*365,&
                pool_C_start=C_matrix_1970,pool_N_start=N_matrix_1970,inorg_N_start=N_inorg_matrix_1970,&
                pool_C_final=C_matrix_2009,pool_N_final=N_matrix_2009,inorg_N_final=N_inorg_matrix_2009,&
                start_year=1971,stop_year=2009,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file, out_path = output_path)
  
    !5: RUN WITH DAILY UPDATES of inputdata,output every day:
    call decomp(nsteps=4*24*365, &
                run_name=trim(trim(site)//"_"//trim(description)//"_"//"to2014"), &
                write_hour=1*24,&
                pool_C_start=C_matrix_2009,pool_N_start=N_matrix_2009,inorg_N_start=N_inorg_matrix_2009,&
                pool_C_final=C_matrix_final,pool_N_final=N_matrix_final,inorg_N_final=N_inorg_matrix_final,&
                start_year=2010,stop_year=2014,clm_input_path=clm_data_file,clm_surf_path=clm_surface_file, out_path = output_path)
    deallocate (C_matrix_1970,N_matrix_1970,N_inorg_matrix_1970,C_matrix_2009,N_matrix_2009,N_inorg_matrix_2009)
    deallocate (C_matrix_final,N_matrix_final,N_inorg_matrix_final)
  endif

  deallocate (C_matrix_init,N_matrix_init,N_inorg_matrix_init,C_matrix_Spunup,N_matrix_Spunup,N_inorg_matrix_Spunup)

  call system_clock(count=full_clock_stop)      ! Stop Timer
  print*, "Total time for simulation in minutes: ", (real(full_clock_stop-full_clock_start)/real(full_clock_rate))/60

end program main









