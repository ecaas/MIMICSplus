!This module models the decomposition of carbon and nitrogen through the soil. The number of soil depth levels and different pools are given in the file paramMod.f90. The paramMod module contains
!all relevant parameters for the flux and balance equations used in this module. In the current setup, 7 depth levels are used, and the following 10 carbon and nitrogen pools/reservoirs are found in each level:
!LITm - metabolic litter
!LITs - Structural litter
!SAPb - bacteria saprotrophs
!SAPf - fungi saprotrophs
!EcM  - Ectomycorrhiza
!ErM  - Ericoid mycorrhiza
!AM   - Arbuscular mycorrhiza
!SOMp - Physically protected soil organic matter
!SOMa - Available soil organic matter
!SOMc - Chemically protected soil organic matter
!In addition a reservoir of inorganic nitrogen, N_IN, is found in each layer. A plant pool of carbon and nitrogen is also included (not vertically resolved).

module mycmimMod
  use shr_kind_mod, only : r8 => shr_kind_r8
  use paramMod
  use initMod,    only: nlevels,calc_init_NH4
  use readMod,    only: read_maxC, read_time,read_clm_model_input,read_clay,&
                        read_WATSAT_and_profiles,read_PFTs
  use fluxMod,    only: calc_nitrification,calc_Leaching,set_N_dep,calc_desorp,MMK_flux,input_rates,&
                        calculate_fluxes,vertical_diffusion,myc_to_plant, NH4_sol_final, NH4_sorp_final,NO3_final
  use writeMod,   only: create_netcdf,create_yearly_mean_netcdf,fill_netcdf, &
                        fill_yearly_netcdf,fluxes_netcdf,store_parameters,check
  use testMod,    only: respired_mass, test_mass_conservation_C,test_mass_conservation_N, &
                        total_carbon_conservation,total_nitrogen_conservation
  use dispmodule, only: disp !External module to pretty print matrices (mainly for testing purposes)
  use netcdf,     only: nf90_nowrite,nf90_write,nf90_close,nf90_open
  
  implicit none
    private 
    public :: decomp
    
contains
  
  subroutine annual_mean(yearly_sumC,yearly_sumN,yearly_sumNinorg,nlevels, year, run_name)
    !Input 
    REAL(r8), DIMENSION(nlevels,pool_types)  , intent(in):: yearly_sumC
    REAL(r8), DIMENSION(nlevels,pool_types_N), intent(in):: yearly_sumN
    REAL(r8), DIMENSION(nlevels,inorg_N_pools), intent(in):: yearly_sumNinorg
    
    integer,  intent(in) :: year
    integer , intent(in) :: nlevels
    CHARACTER (len = *), intent(in):: run_name
    !Local
    REAL(r8), DIMENSION(nlevels,pool_types) :: yearly_meanC
    REAL(r8), DIMENSION(nlevels,pool_types_N) :: yearly_meanN
    REAL(r8), DIMENSION(nlevels,inorg_N_pools) :: yearly_meanNinorg
    integer, parameter                         :: hr_in_year = 24*365
    
    yearly_meanC=yearly_sumC/hr_in_year
    yearly_meanN=yearly_sumN/hr_in_year
    yearly_meanNinorg=yearly_sumNinorg/hr_in_year
    
    call fill_yearly_netcdf(run_name, year, yearly_meanC,yearly_meanN,yearly_meanNinorg)

  end subroutine annual_mean

  
  subroutine decomp(nsteps,   &
                    run_name, &
                    step_frac, write_hour,&
                    pool_C_start,pool_N_start,inorg_N_start, &
                    pool_C_final,pool_N_final,inorg_N_final,  &
                    start_year,stop_year,clm_input_path,clm_surf_path) !Calculates the balance equations dC/dt and dN/dt for each pool at each time step based on the fluxes calculated in the same time step. Then update the pool sizes before moving on
      !INPUT
      integer,intent(in)                        :: nsteps               ! number of time steps to iterate over
      character (len=*) ,intent(in)             :: run_name             ! used for naming outputfiles
      integer,intent(in)                        :: step_frac            ! determines the size of the time step
    !  integer,intent(in)                        :: nlevels           ! number of vertical layers
      integer,intent(in)                        :: write_hour           !How often output is written to file
      real(r8),intent(in)                       :: pool_C_start(nlevels,pool_types)     ! For store and output final C pool sizes 
      real(r8),intent(in)                       :: pool_N_start(nlevels,pool_types_N)   ! For storing and output final N pool sizes [gN/m3] 
      real(r8),intent(in)                       :: inorg_N_start(nlevels,inorg_N_pools)   ! For storing and output final N pool sizes [gN/m3] 
      
      integer, intent(in)                       :: start_year !Forcing start year
      integer, intent(in)                       :: stop_year !Forcing end year, i.e. forcing loops over interval start_year-stop_year
      character (len=*) ,intent(in)             :: clm_input_path             ! file path for input
      character (len=*) ,intent(in)             :: clm_surf_path            ! file path for surface data
  
      !OUTPUT
      real(r8),intent(out)                      :: pool_C_final(nlevels,pool_types)     ! For store and output final C pool sizes 
      real(r8),intent(out)                      :: pool_N_final(nlevels,pool_types_N)   ! For storing and output final N pool sizes [gN/m3] 
      real(r8),intent(out)                      :: inorg_N_final(nlevels,inorg_N_pools)   ! For storing and output final N pool sizes [gN/m3] 
      
      real(r8)                    :: pool_C_start_for_mass_cons(nlevels,pool_types)     ! To use in total mass conservation subroutine
      real(r8)                    :: pool_N_start_for_mass_cons(nlevels,pool_types_N)   ! 
      real(r8)                    :: pool_Ninorg_start_for_mass_cons(nlevels,inorg_N_pools)   ! 
     !Shape of pool_matrixC/change_matrixC
     !|       LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc |
     !|level1   1   2    3    4   5   6   7   8    9    10  |
     !|level2                                               |
     !| .                                                   |
     !| .                                                   |
     !|nlevels _____________________________________________|

     !Shape of the pool_matrixN/change_matrixN
     !|       LITm LITs SAPb SAPf EcM ErM AM SOMp SOMa SOMc|
     !|level1   1   2    3    4   5   6   7   8    9    10 |
     !|level2                                              |
     !| .                                                  |
     !| .                                                  |
     !|nlevels ____________________________________________|
     
     !Shape of inorg_N_matrix
     !|       NH4_sol NH4_sorp   NO3|
     !|level1    1       2       3  |
     !|level2                       |
     !| .                           |
     !| .                           |
     !|nlevels _____________________|
     
      !LOCAL
      character (len=4)              :: year_fmt
      character (len=4)              :: year_char
      logical                        :: isVertical                           ! True if we use vertical soil layers.
      real(r8),dimension(nlevels) :: HR                                   ! For storing the C  that is lost to respiration [gC/m3h]
      real(r8)                       :: HR_mass_accumulated, HR_mass
      real(r8),dimension(nlevels) :: HRb,HRf !For storing respiration separately for bacteria and fungi
      real(r8)                       :: pool_matrixC(nlevels,pool_types)     ! For storing C pool sizes [gC/m3]
      real(r8)                       :: change_matrixC(nlevels,pool_types)   ! For storing dC/dt for each time step [gC/(m3*hour)]
      real(r8)                       :: pool_temporaryC(nlevels,pool_types)  ! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
      real(r8)                       :: pool_matrixC_previous(nlevels,pool_types) !Used for checking mass conservation
      real(r8)                       :: pool_matrixN_previous(nlevels,pool_types_N) !Used for checking mass conservation
      real(r8)                       :: inorg_N_matrix_previous(nlevels,inorg_N_pools) !Used for checking mass conservation
      real(r8)                       :: pool_temporaryN(nlevels,pool_types_N)! When isVertical is True, pool_temporaryC = pool_matrixC + change_matrixC*dt is used to calculate the vertical transport
      real(r8)                       :: pool_matrixN(nlevels,pool_types_N)   ! For storing N pool sizes [gN/m3] parallell to C pools and  inorganic N
      real(r8)                       :: change_matrixN(nlevels,pool_types_N) ! For storing dC/dt for each time step [gN/(m3*hour)]
      real(r8)                       :: inorg_N_matrix(nlevels, inorg_N_pools) ! Inorganic N pool sizes (gN/m3)
      real(r8)                       :: inorg_Ntemporary_matrix(nlevels, inorg_N_pools) ! Inorganic N pool sizes (gN/m3)
      
      real(r8),dimension(nlevels) :: ROI_EcM
      real(r8),dimension(nlevels) :: ROI_AM
      real(r8),dimension(nlevels,2)          :: f_alloc
      real(r8)                       :: sum_consN(nlevels, pool_types_N) !g/m3 for calculating annual mean
      real(r8)                       :: sum_consNinorg(nlevels, inorg_N_pools) !g/m3 for calculating annual mean
      real(r8)                       :: sum_consC(nlevels, pool_types) !g/m3 for calculating annual mean
      real(r8)                       :: N_DEPinput
      real(r8)                       :: C_MYCinput
      real(r8)                       :: C_leaf_litter
      real(r8)                       :: C_root_litter
      real(r8)                       :: N_leaf_litter
      real(r8)                       :: N_root_litter
      real(r8)                       :: C_CWD_litter(nlevels)
      real(r8)                       :: N_CWD_litter(nlevels)

!      real(r8)                       :: dt                            ! size of time step
      real(r8)                       :: time                          ! t*dt
      real(r8)                       :: C_Loss, C_Gain, N_Gain, N_Loss
      real(r8)                       :: tot_diffC,upperC,lowerC                 ! For the call to vertical_diffusion
      real(r8)                       :: tot_diffN,upperN,lowerN                 ! For the call to vertical_diffusion
      real(r8)                       :: tot_diffNinorg,upperNinorg,lowerNinorg                 ! For the call to vertical_diffusion
      
      real(r8)                       :: sum_input_step         ! Used for checking mass conservation
      real(r8)                       :: sum_input_total        ! Used for checking mass conservation
      real(r8)                       :: sum_N_input_total
      real(r8)                       :: sum_N_out_total
      real(r8)                       :: sum_N_out_step
      real(r8)                       :: sum_N_input_step       ! Used for checking mass conservation

      real(r8),allocatable           :: vertC(:,:)         !Stores the vertical change in a time step, same shape as change_matrixC
      real(r8),allocatable           :: vertN(:,:)         !Stores the vertical change in a time step, same shape as change_matrixC
      real(r8),allocatable           :: vert_inorgN(:,:)         !Stores the vertical change in a time step, same shape as change_matrixC
      

      !Counters
      integer                        :: ycounter, year,write_y
      integer                        :: counter            !used for determining when to output results
      integer                        :: month_counter      !for determining when to read new input values (Temp. litfall etc.)
      integer                        :: day_counter
      integer                        :: spinup_counter
      integer                        :: current_month
      integer                        :: current_day
      logical                        :: Spinup_run
      integer                        :: j,i,t              !for iterations
      integer,parameter              ::t_init=1
      integer                        :: date
      integer                       :: input_steps
      real(r8)                       :: change_sum(nlevels, pool_types)
      real(r8)                       :: vertC_change_sum(nlevels, pool_types)
      real(r8),dimension(15)         :: PFT_distribution
      real(r8)                       :: C_EcMenz_prod
      real(r8),dimension(no_of_som_pools,no_of_sap_pools)        :: f_partitions
      
      !For reading soil temperature and moisture from CLM output file
      real(r8), dimension(nlevels)          :: TSOIL
      real(r8), dimension(nlevels)          :: SOILLIQ
      real(r8), dimension(nlevels)          :: SOILICE
      real(r8), dimension(nlevels)          :: WATSAT
      real(r8), dimension(nlevels)          :: W_SCALAR
      real(r8), dimension(nlevels)          :: T_SCALAR
      real(r8)                                 :: drain
      real(r8)                                 :: h2o_liq_tot
      real(r8)                              :: H2OSOI
      
      real(r8) :: save_N,save_C
          
      integer :: ncid
      integer :: writencid
      integer :: spinupncid
      
      call system_clock(count_rate=clock_rate) !Find the time rate
      call system_clock(count=clock_start)     !Start Timer  
      
!      dt = calc_timestep(step_frac)
      
      if (nlevels>1) then 
        soil_depth=sum(delta_z(1:nlevels))
        isVertical = .True.
      else
        soil_depth=1.52_r8
        isVertical = .False.
        allocate (vertC, mold = pool_matrixC)
        allocate (vertN, mold = pool_matrixN)
        allocate (vert_inorgN, mold = inorg_N_matrix)
                     !TODO: This can be done better
      end if
          
      !Allocate and initialize
      allocate(CUE_bacteria_vr(nlevels))
      CUE_bacteria_vr=CUE_0
      allocate(CUE_fungi_vr(nlevels))
      CUE_fungi_vr=CUE_0
      allocate(CUE_ecm_vr(nlevels))
      CUE_ecm_vr=CUE_myc_0
      allocate(CUE_erm_vr(nlevels))
      CUE_erm_vr=0.0
      allocate(CUE_am_vr(nlevels))
      CUE_am_vr=CUE_myc_0      
      allocate(r_moist(nlevels))
      allocate(f_enzprod(nlevels))
      f_enzprod=f_enzprod_0
      
      input_mod=1.0 !Initialize input modifier (r_input)
      !For counting mineralization/immobilization occurences:
      c1a=0
      c1b=0
      c2=0
      c3a=0
      c3b=0
      c4a=0
      c4b=0
      save_N=0._r8
      save_C=0._r8

      !Set initial concentration values:
      pool_matrixC=pool_C_start
      pool_matrixN=pool_N_start
      inorg_N_matrix=inorg_N_start
      pool_matrixC_previous = pool_C_start
      pool_matrixN_previous = pool_N_start
      inorg_N_matrix_previous=inorg_N_start
      pool_C_start_for_mass_cons=pool_C_start
      pool_N_start_for_mass_cons=pool_N_start
      pool_Ninorg_start_for_mass_cons=inorg_N_start
      
      !Make sure things start from zero
      sum_consN      = 0
      sum_consNinorg = 0
      sum_consC      = 0
      change_matrixC = 0.0
      change_matrixN = 0.0
      HR             = 0.0
      HRb            = 0.0
      HRf            = 0.0
      ROI_EcM        = 0.0
      ROI_AM         = 0.0
      f_alloc        = 0.0
      
      vertC_change_sum=0.0
      HR_mass_accumulated = 0
      sum_input_step   =0.0
      sum_input_total  =0.0
      sum_N_input_step =0.0
      sum_N_input_total=0.0
      sum_N_out_step   =0.0
      sum_N_out_total  =0.0

      !For counting:
      year     = start_year
      year_fmt = '(I4)'
      write (year_char,year_fmt) year
      write_y=0
      current_month = 1
      current_day   = 1      
      month_counter = 0
      day_counter = 0
      counter  = 0
      ycounter = 0
      
      !read data from CLM file
      if ( start_year == 1850 ) then
        Spinup_run = .True.
        spinup_counter =1
        call check(nf90_open(trim(adjustr(clm_input_path)//'_historical.clm2.for_spinup.1850-1869.nc'), nf90_nowrite, spinupncid)) !open netcdf containing values for the next year  
        call read_time(spinupncid,input_steps) !Check if inputdata is daily or monthly:         
        call read_clm_model_input(spinupncid,Spinup_counter, &
        N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
        C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
        W_SCALAR,T_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)        
      else
        Spinup_run = .False.        
        call check(nf90_open(trim(adjustr(clm_input_path)//'_historical.clm2.all.'//year_char//'.nc'), nf90_nowrite, ncid)) !open netcdf containing values for the next year  
        call read_time(ncid,input_steps) !Check if inputdata is daily or monthly: 
        call read_clm_model_input(ncid,1, &
        N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
        C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
        W_SCALAR,T_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)
      end if
      allocate(ndep_prof(nlevels),leaf_prof(nlevels),froot_prof(nlevels))   
         
      call read_WATSAT_and_profiles(adjustr(clm_input_path)//'_historical.clm2.all.'//"1901.nc",WATSAT,ndep_prof,froot_prof,leaf_prof)         
      call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist)                   
      call read_clay(adjustr(clm_surf_path),fCLAY)
      
      if ( .not. use_ROI ) then !use static PFT determined fractionation between EcM and AM C input
        call read_PFTs(adjustr(clm_surf_path),PFT_distribution)
        f_EcM = calc_EcMfrac(PFT_distribution)
        f_alloc(:,1) = f_EcM
        f_alloc(:,2) = 1-f_EcM   
      else
        f_EcM = 999_r8
      end if
      if (f_alloc(1,1)==1.0 ) then !To avoid writing errors when AM = 0
        pool_matrixC(:,7)=0.0
        pool_matrixN(:,7)=0.0
        pool_matrixC_previous(:,7)=0.0
        pool_matrixN_previous(:,7)=0.0
        pool_C_start_for_mass_cons=pool_matrixC
        pool_N_start_for_mass_cons=pool_matrixN
      end if
      if ( Spinup_run ) then
        max_mining = read_maxC(spinupncid,input_steps)
      else
        max_mining = read_maxC(ncid,input_steps)        
      end if
            
      !open and prepare files to store results. Store initial values
      if ( Spinup_run ) then
        call create_yearly_mean_netcdf(run_name)  
      end if
      call create_netcdf(run_name)
      call check(nf90_open(output_path//trim(run_name)//".nc", nf90_write, writencid))
      
      call fill_netcdf(writencid,t_init, pool_matrixC, pool_matrixN,inorg_N_matrix, &
                       date, HR_mass_accumulated,HR,HRb,HRf, change_matrixC,change_matrixN,write_hour,current_month, &
                      TSOIL, r_moist,CUE_bacteria_vr,CUE_fungi_vr,CUE_EcM_vr,CUE_am_vr,ROI_EcM=ROI_EcM,ROI_AM=ROI_AM,enz_frac=f_enzprod,f_alloc=f_alloc)
            
      desorp= calc_desorp(fCLAY)
      Kmod  = calc_Kmod(fCLAY)
      k_mycsom  = calc_myc_mortality()  
      ! Fracions of SAP that goes to different SOM pools
      f_partitions = calc_sap_to_som_fractions(fCLAY,fMET)
      fPHYS = f_partitions(1,:)
      fCHEM = f_partitions(2,:)
      fAVAIL  = f_partitions(3,:)

      !print initial values to terminal      
      call disp("InitC", pool_matrixC)
      call disp("InitN", pool_matrixN)
      call disp("InitN inorganic", inorg_N_matrix)
      !----------------------------------------------------------------------------------------------------------------
      do t =1,nsteps !Starting time iterations
        time = t*dt
        !Time counters:
        counter  = counter + 1
        ycounter = ycounter + 1            !Counts hours in a year 
        month_counter = month_counter + 1 !Counts hours in a month
        day_counter = day_counter + 1     !Counts hours in a day (24)

        ! !Update temp and moisture values monthly/daily 
        if (month_counter == days_in_month(current_month)*hr_pr_day*step_frac+1) then
            month_counter = 1       
                 
            if (current_month == 12) then
              current_month=1 !Update to new year
            else
              current_month = current_month + 1 
            end if   
   
            if (input_steps==12) then  
              call read_clm_model_input(ncid,current_month, &
                                N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
                                C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
                                W_SCALAR,T_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)  
              call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist)   
              max_mining = read_maxC(ncid,input_steps)
            end if     
            
            if (input_steps==240) then
              spinup_counter = spinup_counter+1

              call read_clm_model_input(spinupncid,spinup_counter, &
                                N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
                                C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
                                W_SCALAR,T_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)  
              call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist)   
              max_mining = read_maxC(spinupncid,input_steps)
            end if                             
        end if 
        
        if (day_counter == hr_pr_day*step_frac+1) then
          day_counter = 1   
          current_day = current_day +1                      
          if ( current_day == 366 ) then
            current_day = 1
          end if
          
          if (input_steps==365) then                        
            call read_clm_model_input(ncid,current_day, &
            N_leaf_litter,N_root_litter,C_MYCinput,N_DEPinput, &
            C_leaf_litter,C_root_litter,date,TSOIL,SOILLIQ,SOILICE, &
            W_SCALAR,T_SCALAR,drain,h2o_liq_tot,C_CWD_litter,N_CWD_litter)
            call moisture_func(SOILLIQ,WATSAT, SOILICE,r_moist)        
            max_mining = read_maxC(ncid,input_steps)                
          end if        
        end if         
        

        
        input_mod = r_input(C_MYCinput,max_mining) !calculate factor that scales mychorrizal activity based on C payment from plant
        if ( abs(input_mod) > 1.0 ) then
          print*, input_mod, time
        end if
        
        
        do j = 1, nlevels !For each depth level (for the no vertical transport case, nlevels = 1, so loop is only done once):
          
          H2OSOI=SOILLIQ(j)+SOILICE(j)

          
          !Michaelis Menten parameters:
          Km      = Km_function(TSOIL(j))
          Vmax    = Vmax_function(TSOIL(j),r_moist(j)) !  ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics.

          ![1/h] Microbial turnover rate (SAP to SOM)
          k_sapsom = calc_sap_turnover_rate(fMET,r_moist(j)) 
          
          CUE_bacteria_vr(j) = (CUE_slope*TSOIL(j)+CUE_0)
          CUE_fungi_vr(j) = (CUE_slope*TSOIL(j)+CUE_0)
          CUE_ecm_vr(j) = CUE_myc_0
          CUE_am_vr(j) = CUE_myc_0
          f_enzprod=0.1_r8
          
          !Determine input rates (from CLM data) in timestep
          call input_rates(j,C_leaf_litter,C_root_litter,N_leaf_litter,&
                                      N_root_litter, &
                                      N_CWD_litter,C_CWD_litter,&
                                      C_PlantLITm,C_PlantLITs, &
                                      N_PlantLITm,N_PlantLITs, &
                                      C_PlantSOMp,C_PlantSOMa,C_PlantSOMc, &
                                      N_PlantSOMp,N_PlantSOMa,N_PlantSOMc)
                                      
          Deposition = set_N_dep(CLMdep = N_DEPinput*ndep_prof(j)) !NOTE: either const_dep = some_value or CLMdep = N_DEPinput*ndep_prof(j)
          Leaching = calc_Leaching(drain,h2o_liq_tot,inorg_N_matrix(j,3))
          nitrif_rate=calc_nitrification((inorg_N_matrix(j,1)+Deposition*dt),W_SCALAR(j),T_SCALAR(j),TSOIL(j)) !NOTE: Uses NH4 + Deposiiton in timestep

          !Calculate fluxes between pools in level j (file: fluxMod.f90):
          call calculate_fluxes(j,TSOIL(j),H2OSOI, pool_matrixC, pool_matrixN,inorg_N_matrix,Deposition, Leaching, nitrif_rate,dt)
          inorg_Ntemporary_matrix(j,1)=NH4_sol_final
          inorg_Ntemporary_matrix(j,2)=NH4_sorp_final
          inorg_Ntemporary_matrix(j,3)=NO3_final
          
          if ( inorg_N_matrix(j,1) < 2.2204460492503131E-016 ) then
            save_N=save_N + inorg_N_matrix(j,1)*dt*delta_z(j)          
            inorg_N_matrix(j,1)=0._r8
          end if
          if ( inorg_N_matrix(j,2) < 2.2204460492503131E-016 ) then
            save_N=save_N + inorg_N_matrix(j,2)*dt*delta_z(j)          
            inorg_N_matrix(j,2)=0._r8
          end if
          if ( inorg_N_matrix(j,3) < 2.2204460492503131E-016 ) then
            save_N=save_N + inorg_N_matrix(j,3)*dt*delta_z(j)          
            inorg_N_matrix(j,3)=0._r8
          end if
          if ( use_ROI ) then
            ROI_EcM(j) = ROI_function(N_INEcM+N_SOMpEcM+N_SOMcEcM,pool_matrixC(j,5),k_mycsom(1))
            ROI_AM(j) = ROI_function(N_INAM, pool_matrixC(j,7),k_mycsom(3))
            if ( C_MYCinput .NE. 0.0  ) then !To avoid NaN when both ROI is zero
              if ( ROI_EcM(j) +ROI_AM(j) ==0.0 ) then !Too little myc in layer to contribute
                f_alloc(j,:)=0.0 !Value does not really matter bc. C_MYCinput is zero
              else
                f_alloc(j,1) = ROI_EcM(j)/(ROI_EcM(j)+ROI_AM(j))
                f_alloc(j,2) = ROI_AM(j)/(ROI_EcM(j)+ROI_AM(j))
              end if
            else
              f_alloc(j,:)=0.5 !Value does not really matter bc. C_MYCinput is zero
            end if
          end if
          
          C_PlantEcM = f_alloc(j,1)*C_MYCinput*froot_prof(j)
          C_PlantAM = f_alloc(j, 2)*C_MYCinput*froot_prof(j)          
          call myc_to_plant(j,C_PlantEcM,C_PlantAM,N_AMPlant,N_EcMPlant,N_ErMPlant,CUE_EcM_vr(j),CUE_AM_vr(j),f_enzprod(j))
      !---------------------------------------------------------------------------------------          
          if (counter == write_hour*step_frac .or. t==1) then !Write fluxes from calculate_fluxes to file            
            call fluxes_netcdf(writencid,int(time), write_hour, j)
          end if !write fluxes

          do i = 1,pool_types !loop over all the pool types, i, in depth level j 
            !This if-loop calculates dC/dt and dN/dt for the different carbon pools.
            !NOTE: If pools are added/removed (i.e the actual model equations is changed), this loop needs to be updated.

            if (i==1) then !LITm
              N_Gain = N_PlantLITm
              N_Loss = N_LITmSAPb + N_LITmSAPf
              C_Gain = C_PlantLITm
              C_Loss = C_LITmSAPb + C_LITmSAPf
              
            elseif (i==2) then !LITs
              N_Gain = N_PlantLITs
              N_Loss = N_LITsSAPb + N_LITsSAPf
              C_Gain = C_PlantLITs
              C_Loss = C_LITsSAPb + C_LITsSAPf

            elseif (i==3) then !SAPb              
              C_Gain = CUE_bacteria_vr(j)*(C_LITmSAPb + C_LITsSAPb &
                + C_SOMaSAPb)
              C_Loss =  C_SAPbSOMp + C_SAPbSOMa + C_SAPbSOMc
              N_Gain = (N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb)*NUE
              N_Loss = N_SAPbSOMp + N_SAPbSOMa + N_SAPbSOMc
              if ( N_INSAPb>0 ) then
                N_Gain = N_Gain + N_INSAPb              
              else
                N_Loss=N_Loss-N_INSAPb     !two minus becomes +         
              end if

            elseif (i==4) then !SAPf
              C_Gain = CUE_fungi_vr(j)*(C_LITmSAPf + C_LITsSAPf &
                + C_SOMaSAPf)
              C_Loss =  C_SAPfSOMp + C_SAPfSOMa + C_SAPfSOMc
              N_Gain = (N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf)*NUE
              N_Loss = N_SAPfSOMp + N_SAPfSOMa + N_SAPfSOMc
              if ( N_INSAPf>0 ) then
                N_Gain = N_Gain + N_INSAPf 
              else
                N_Loss=N_Loss-N_INSAPf !two minus becomes +
              end if

            elseif (i==5) then !EcM
              C_EcMenz_prod=CUE_ecm_vr(j)*C_PlantEcM*f_enzprod(j)
              
              C_Gain = CUE_ecm_vr(j)*C_PlantEcM + C_SOMcEcM + C_SOMpEcM !
              C_Loss = C_EcMSOMp + C_EcMSOMa + C_EcMSOMc + CUE_ecm_vr(j)*C_PlantEcM*f_enzprod(j)
              N_Gain = N_INEcM + N_SOMpEcM + N_SOMcEcM
              N_Loss = N_EcMPlant + N_EcMSOMa + N_EcMSOMp + N_EcMSOMc
            elseif (i==6) then !ErM
              C_Gain = CUE_erm_vr(j)*C_PlantErM
              C_Loss = C_ErMSOMp + C_ErMSOMa + C_ErMSOMc 
              N_Gain = N_INErM
              N_Loss = N_ErMPlant + N_ErMSOMa + N_ErMSOMp + N_ErMSOMc

            elseif (i==7) then !AM
              C_Gain = CUE_am_vr(j)*C_PlantAM
              C_Loss = C_AMSOMp + C_AMSOMa + C_AMSOMc
              N_Gain = N_INAM 
              N_Loss = N_AMPlant + N_AMSOMa + N_AMSOMp + N_AMSOMc

            elseif (i==8) then !SOMp
              C_Gain =  C_SAPbSOMp + C_SAPfSOMp + C_EcMSOMp + C_ErMSOMp + C_AMSOMp+ C_PlantSOMp
              C_Loss = C_SOMpSOMa+C_EcMdecompSOMp + C_SOMpEcM
              N_Gain =  N_SAPbSOMp + N_SAPfSOMp + N_EcMSOMp + N_ErMSOMp + N_AMSOMp+N_PlantSOMp
              N_Loss = N_SOMpSOMa + N_SOMpEcM

            elseif (i==9) then !SOMa
               C_Gain = C_SAPbSOMa + C_SAPfSOMa + C_EcMSOMa + C_EcMdecompSOMp + C_EcMdecompSOMc &
               + C_ErMSOMa + C_AMSOMa + C_SOMpSOMa + C_SOMcSOMa + C_PlantSOMa+ CUE_ecm_vr(j)*C_PlantEcM*f_enzprod(j)
               C_Loss = C_SOMaSAPb + C_SOMaSAPf 
               N_Gain = N_SAPbSOMa + N_SAPfSOMa + N_EcMSOMa + &
               N_ErMSOMa + N_AMSOMa + N_SOMpSOMa + N_SOMcSOMa +N_PlantSOMa
               N_Loss = N_SOMaSAPb + N_SOMaSAPf
            elseif (i==10) then !SOMc
              C_Gain =  C_SAPbSOMc + C_SAPfSOMc + C_EcMSOMc + C_ErMSOMc + C_AMSOMc+C_PlantSOMc
              C_Loss = C_SOMcSOMa+C_EcMdecompSOMc + C_SOMcEcM
              N_Gain =  N_SAPbSOMc + N_SAPfSOMc + N_EcMSOMc + N_ErMSOMc + N_AMSOMc+N_PlantSOMc
              N_Loss = N_SOMcSOMa + N_SOMcEcM
              
            else
              print*, 'Too many pool types expected, pool_types = ',pool_types, 'i: ', i
            end if !determine total gains and losses

            change_matrixC(j,i) = C_Gain - C_Loss !net change in timestep
            change_matrixN(j,i) = N_Gain - N_loss            
            !Store these values as temporary so that they can be used in the vertical diffusion subroutine
            pool_temporaryC(j,i)=pool_matrixC(j,i) + change_matrixC(j,i)*dt
            pool_temporaryN(j,i)=pool_matrixN(j,i) + change_matrixN(j,i)*dt
            
            if ( pool_temporaryC(j,i) < 2.2204460492503131E-016 ) then
              save_C=save_C+pool_temporaryC(j,i)
              pool_temporaryC(j,i)=0._r8
            end if   
            if ( pool_temporaryN(j,i) < 2.2204460492503131E-016 ) then  
              save_N=save_N + pool_temporaryN(j,i)*dt*delta_z(j)
              pool_temporaryN(j,i)=0._r8                
            end if                                    
          end do !i, pool_types
          
          !Calculate the heterotrophic respiration loss from depth level j in timestep t:
          HR(j) =(( C_LITmSAPb + C_LITsSAPb  + C_SOMaSAPb)*(1-CUE_bacteria_vr(j)) + (C_LITmSAPf &
          + C_LITsSAPf + C_SOMaSAPf)*(1-CUE_fungi_vr(j)))*dt
          HRb(j) = ( C_LITmSAPb + C_LITsSAPb  + C_SOMaSAPb)*(1-CUE_bacteria_vr(j))*dt
          HRf(j)=(C_LITmSAPf + C_LITsSAPf + C_SOMaSAPf)*(1-CUE_fungi_vr(j))*dt

          if (HR(j) < 0 ) then
            print*, 'Negative HR: ', HR(j), t,j
            print*, "Pools C", pool_matrixC(j,1),pool_matrixC(j,2),pool_matrixC(j,3),pool_matrixC(j,4), pool_matrixC(j,9)
            print*, "pools N", inorg_N_matrix(j,1), inorg_N_matrix(j,3), N_INSAPb, N_INSAPf
            print*, "Fluxes", C_LITmSAPb,C_LITsSAPb,C_SOMaSAPb,C_LITmSAPf,C_LITsSAPf,C_SOMaSAPf
            stop
          end if
          
          !Summarize in and out print timestep to check mass balance
          sum_input_step=sum_input_step+(C_PlantLITm+C_PlantLITs+CUE_ecm_vr(j)*C_PlantEcM+CUE_am_vr(j)*C_PlantAM+C_PlantSOMc+C_PlantSOMp+C_PlantSOMa)*dt*delta_z(j) !g/m2
          sum_N_input_step=sum_N_input_step+(N_PlantLITm+N_PlantLITs+N_PlantSOMc+N_PlantSOMp+N_PlantSOMa+Deposition)*dt*delta_z(j) !g/m2
          sum_N_out_step=sum_N_out_step+(N_EcMPlant+N_AMPlant+N_INPlant+Leaching)*dt*delta_z(j)

        end do !j, depth_level
        !Store accumulated HR mass
        call respired_mass(HR, HR_mass)
        HR_mass_accumulated = HR_mass_accumulated + HR_mass
        !TODO: tot_diffC, upperC, lowerC is not used and can be removed!
        if (isVertical) then
          call vertical_diffusion(tot_diffC,upperC,lowerC, pool_temporaryC,vertC,D_carbon)
          call vertical_diffusion(tot_diffN,upperN,lowerN, pool_temporaryN,vertN,D_nitrogen)
          call vertical_diffusion(tot_diffNinorg,upperNinorg,lowerNinorg,inorg_Ntemporary_matrix,vert_inorgN,D_nitrogen)

          pool_matrixC    = vertC*dt        + pool_temporaryC
          pool_matrixN    = vertN*dt        + pool_temporaryN
          inorg_N_matrix  = vert_inorgN*dt  + inorg_Ntemporary_matrix
          
        else
          pool_matrixC    = pool_temporaryC
          pool_matrixN    = pool_temporaryN
          inorg_N_matrix  = inorg_Ntemporary_matrix
        end if!isVertical

        sum_consN       = sum_consN       + pool_matrixN
        sum_consNinorg  = sum_consNinorg  + inorg_N_matrix
        sum_consC       = sum_consC       + pool_matrixC

        if (ycounter == 365*24*step_frac) then
          ycounter = 0
          write_y =write_y+1 !For writing to annual mean file
          
          if ( Spinup_run ) then
             call annual_mean(sum_consC,sum_consN,sum_consNinorg, nlevels,write_y , run_name) !calculates the annual mean and write the result to file
          end if
          if (year == stop_year) then
            year = start_year         
            spinup_counter=0            
          else 
            year = year + 1             
          end if          
          write (year_char,year_fmt) year          
          if ( .not. Spinup_run ) then            
            call check(nf90_close(ncid)) !Close netcdf file containing values for the past year
            call check(nf90_open(trim(adjustr(clm_input_path)//'_historical.clm2.all.'//year_char//'.nc'), nf90_nowrite, ncid)) !open netcdf containing values for the next year
            call read_time(ncid,input_steps)     
          end if           
          sum_consN =0
          sum_consC =0
        end if
        if (counter == write_hour*step_frac) then
          counter = 0        
          call fill_netcdf(writencid, int(time), pool_matrixC, pool_matrixN,inorg_N_matrix,&
           date, HR_mass_accumulated,HR,HRb,HRf,vertC,vertN, write_hour,current_month, &
           TSOIL, r_moist,CUE_bacteria_vr,CUE_fungi_vr,CUE_ecm_vr,CUE_am_vr,ROI_EcM=ROI_EcM,ROI_AM=ROI_AM,enz_frac=f_enzprod,f_alloc=f_alloc)
          change_sum = 0.0
        end if!writing

        !Write end values to terminal
        if (t == nsteps) then
          call disp("pool_matrixC gC/m3 ",pool_matrixC)
          call disp("pool_matrixN gN/m3 ",pool_matrixN)          
          call disp("Inorganic N gN/m3 ",inorg_N_matrix)          
          call disp("C:N : ",pool_matrixC/pool_matrixN)
          pool_C_final  = pool_matrixC
          pool_N_final  = pool_matrixN    
          inorg_N_final = inorg_N_matrix              
          call store_parameters(writencid)    
        end if
        
        call test_mass_conservation_C(sum_input_step,HR_mass, &
                                      pool_matrixC_previous,pool_matrixC, &
                                      nlevels)
        call test_mass_conservation_N(sum_N_input_step,sum_N_out_step, &
                                      pool_matrixN_previous,inorg_N_matrix_previous,&
                                      pool_matrixN,inorg_N_matrix,nlevels)
        pool_matrixC_previous=pool_matrixC
        pool_matrixN_previous=pool_matrixN
        inorg_N_matrix_previous=inorg_N_matrix
        sum_input_total=sum_input_total+sum_input_step
        sum_input_step=0.0
        sum_N_input_total=sum_N_input_total+sum_N_input_step  
        sum_N_input_step=0.0
        sum_N_out_total=sum_N_out_total+sum_N_out_step    
        sum_N_out_step=0.0
      end do !t
      
      !Check total mass conservation 
      call total_carbon_conservation(sum_input_total,HR_mass_accumulated, &
                                      pool_C_start_for_mass_cons,pool_C_final,&
                                      nlevels)
      call total_nitrogen_conservation(sum_N_input_total,sum_N_out_total, &
                                      pool_N_start_for_mass_cons,pool_Ninorg_start_for_mass_cons, &
                                      pool_N_final,inorg_N_final,nlevels)
      
      call check(nf90_close(writencid))
      if ( Spinup_run ) then
        call check(nf90_close(Spinupncid))
      end if
      print*, "AMOUNT OF N DISCARDED: ", save_N
      print*, "AMOUNT OF C DISCARDED: ", save_C
      
      print*, "Immobilization, not enough N:",c1a
      print*, "Immobilization enough N: ",c1b
      print*, "Mineralization: ",c2
      print*, "Bacteria needs, fungi mineralize, not enough N: ",c3a
      print*, "Bacteria needs, fungi mineralize, enough N: ",c3b
      print*, "Fungi needs, bacteria mineralize, not enough N: ",c4a
      print*, "Fungi needs, bacteria mineralize, enough N: ",c4b
      !deallocation
      deallocate(ndep_prof,leaf_prof,froot_prof,r_moist)
      deallocate(CUE_bacteria_vr,CUE_fungi_vr, CUE_ecm_vr,CUE_am_vr, CUE_erm_vr,f_enzprod)
      
      !For timing
      call system_clock(count=clock_stop)      ! Stop Timer
      print*, "Total time for decomp subroutine in minutes: ", (real(clock_stop-clock_start)/real(clock_rate))/60
  end subroutine decomp



end module mycmimMod
