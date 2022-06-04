module initMod
  use shr_kind_mod, only : r8 => shr_kind_r8
  use netcdf, only : nf90_noerr,nf90_strerror,nf90_nowrite,nf90_inq_varid,nf90_get_var,nf90_close,nf90_open

  implicit none
  private 
  public  :: nlevels, initialize, read_nlayers,calc_init_NH4
  integer :: nlevels
  
contains
  subroutine initialize(pools_C, pools_N) !This subroutine sets the initial C and N values in the pool matrices (content in each soil pool in every depth level + plant)
    !OUTPUT
    real(r8), intent(out):: pools_C(:,:), pools_N(:,:)

    !LOCAL
    integer                :: j
    real(r8), dimension(10), parameter   :: CN_ratio = (/15,15,5,8,20,20,20,11,8,11/) !Fungi/bacteria: Tang, Riley, Maggi 2019 as in Mouginot et al. 2014
    
    do j=1,nlevels
       pools_C(j,:) = (/4000.,4000.,10.,50.,100.,1.,1.,4000.,4000.,4000./)
       pools_N(j,1:10) = pools_C(j,1:10)/CN_ratio
     end do
         
    pools_N(:,11) = 1.
    pools_N(:,12) = 1.
  end subroutine initialize

  subroutine read_nlayers(clm_history_file) 
    !INPUT
    character (len = *),intent(in):: clm_history_file
    
    !OUTPUT
    !integer, intent(out)         :: nlevels
    
    !LOCAL
    integer            :: ncid, nid
    
    call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))
    call check(nf90_inq_varid(ncid, 'nbedrock', nid))
    call check(nf90_get_var(ncid, nid, nlevels))
    call check(nf90_close(ncid))
  end subroutine read_nlayers   
  
  subroutine check(status)
    integer, intent ( in) :: status
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check

  subroutine calc_init_NH4(NH4_tot,water_content,NH4_sorp_eq,NH4_sol_eq)
    !IN:
    real(r8), intent(in)  :: NH4_tot   !g/m3, total NH4, both in soil solution and adsorbed
    real(r8), intent(in)  :: water_content   !m3water/m3soil (input from CLM data)
    
    !Out:
    real(r8),intent(out)            :: NH4_sol_eq  !g/m3, NH4 in soil solution at eq
    real(r8),intent(out)            :: NH4_sorp_eq !g/m3, NH4 sorbed to particles at eq
    
    
    
    real(r8), parameter :: mg_pr_g = 1e3 !mg/g
    real(r8), parameter :: m3_pr_L = 1e-3 !m3/L
    real(r8), parameter :: BD_soil=1.6e6  !g/m3 (loam) soil from DOI: 10.3390/APP6100269 Table 1
    real(r8), parameter :: NH4_sorp_max = 0.09*BD_soil/mg_pr_g    !mg NH4 /g soil
    real(r8), parameter :: KL = 0.4      !L/mg
    real(r8)            :: KL_prime       !m3/g
    !1) Calculate NH4_sorp_eq 
    KL_prime = KL*mg_pr_g*m3_pr_L/water_content 
    NH4_sorp_eq=(1+KL_prime*NH4_tot+NH4_sorp_max*KL_prime)/(2*KL_prime) - sqrt((1+KL_prime*NH4_tot+NH4_sorp_max*KL_prime)**2-4*KL_prime**2*NH4_sorp_max*NH4_tot)/(2*KL_prime)
    NH4_sol_eq=NH4_tot - NH4_sorp_eq
  end subroutine calc_init_NH4
end module initMod
