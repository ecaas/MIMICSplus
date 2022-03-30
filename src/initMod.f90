module initMod
  use shr_kind_mod, only : r8 => shr_kind_r8
  use netcdf, only : nf90_noerr,nf90_strerror,nf90_nowrite,nf90_inq_varid,nf90_get_var,nf90_close,nf90_open
  
  implicit none
  private 
  public  :: nlevels, initialize, read_nlayers
  integer :: nlevels
  
contains
  subroutine initialize(pools_C, pools_N) !This subroutine sets the initial C and N values in the pool matrices (content in each soil pool in every depth level + plant)
    !OUTPUT
    real(r8), intent(out):: pools_C(:,:), pools_N(:,:)

    !LOCAL
    integer                :: j
    real(r8), dimension(10), parameter   :: CN_ratio = (/15,15,5,8,20,20,20,11,8,11/) !Fungi/bacteria: Tang, Riley, Maggi 2019 as in Mouginot et al. 2014
    
    do j=1,nlevels
       pools_C(j,:) =  (/1000.,1000.,10.,50.,100.,1.,1.,1500.,1500.,500./)
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

end module initMod
