module initMod
  use shr_kind_mod, only : r8 => shr_kind_r8
  use netcdf, only : nf90_noerr,nf90_strerror,nf90_nowrite,nf90_inq_varid,nf90_get_var,nf90_close,nf90_open

  implicit none
  private 
  public  :: nlevels, initialize, read_nlayers,calc_init_NH4
  integer :: nlevels
  
contains
  subroutine initialize(pools_C, pools_N,inorg_N) !This subroutine sets the initial C and N values in the pool matrices (content in each soil pool in every depth level + plant)
    !OUTPUT
    real(r8), intent(out):: pools_C(:,:)
    real(r8), intent(out):: pools_N(:,:)
    real(r8), intent(out):: inorg_N(:,:)

    !LOCAL
    integer                :: j
    real(r8), dimension(9), parameter   :: CN_ratio = (/15,15,5,8,20,20,11,8,11/) !Fungi/bacteria: Tang, Riley, Maggi 2019 as in Mouginot et al. 2014
    real(r8)               :: NH4_sorp
    real(r8)               :: NH4_sol 
    
    call calc_init_NH4(tot=10._r8,water_content=0.5_r8,sorp_eq=NH4_sorp,sol_eq=NH4_sol)
    
    do j=1,nlevels
      pools_C(j,:) = (/100.,100.,50.,50.,50.,50.,200.,200.,500./)
      pools_N(j,:) = pools_C(j,:)/CN_ratio
      inorg_N(j,:)=(/NH4_sol,NH4_sorp,10._r8/)
    end do

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

  subroutine calc_init_NH4(tot,water_content,sorp_eq,sol_eq)
    !IN:
    real(r8), intent(in)  :: tot   !g/m3, total NH4, both in soil solution and adsorbed
    real(r8), intent(in)  :: water_content   !m3water/m3soil (input from CLM data)
    
    !Out:
    real(r8),intent(out)            :: sol_eq  !g/m3, NH4 in soil solution at eq
    real(r8),intent(out)            :: sorp_eq !g/m3, NH4 sorbed to particles at eq
    
    
    
    real(r8), parameter :: mg_pr_g = 1e3 !mg/g
    real(r8), parameter :: m3_pr_L = 1e-3 !m3/L
    real(r8), parameter :: BD_soil=1.6e6  !g/m3 (loam) soil from DOI: 10.3390/APP6100269 Table 1
    real(r8), parameter :: NH4_sorp_max = 0.09*BD_soil/mg_pr_g    !mg NH4 /g soil
    real(r8), parameter :: KL = 0.4      !L/mg
    real(r8)            :: KL_prime       !m3/g
    !1) Calculate sorp_eq 
    KL_prime = KL*mg_pr_g*m3_pr_L/water_content 
    sorp_eq=(1+KL_prime*tot+NH4_sorp_max*KL_prime)/(2*KL_prime) - sqrt((1+KL_prime*tot+NH4_sorp_max*KL_prime)**2-4*KL_prime**2*NH4_sorp_max*tot)/(2*KL_prime)
    sol_eq=tot - sorp_eq
  end subroutine calc_init_NH4
end module initMod
