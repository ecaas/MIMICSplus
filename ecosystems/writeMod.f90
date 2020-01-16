module writeMod
  use paramMod
  use netcdf
  implicit none
  integer :: ncid, varid, varidan
  integer :: grid_dimid, col_dimid, t_dimid, lev_dimid,mmk_dimid, MGE_dimid


  contains

    subroutine check(status)
      integer, intent ( in) :: status
      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 2
      end if
    end subroutine check

    subroutine create_netcdf(run_name)
      character (len = *):: run_name

      integer, parameter :: NDIMS = 4 !gridcell, column, levsoi
      integer, parameter :: gridcell = 1, column = 1, levsoi = 1 !TODO change levsoi to match input to decomp subroutine!
      integer :: x,v

      call check(nf90_create(trim(run_name)//".nc", NF90_HDF5, ncid))

      call check(nf90_def_dim(ncid, "time", nf90_unlimited, t_dimid))
      call check(nf90_def_dim(ncid, "gridcell", gridcell, grid_dimid))
      call check(nf90_def_dim(ncid, "column", column, col_dimid))
      call check(nf90_def_dim(ncid, "levsoi", levsoi, lev_dimid))
      call check(nf90_def_dim(ncid, "NoMMKeqs", MM_eqs, mmk_dimid))
      call check(nf90_def_dim(ncid, "NoMGEvalues", size(MGE), MGE_dimid))

      do v = 1, size(variables)
        call check(nf90_def_var(ncid, trim(variables(v)), NF90_DOUBLE, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "change"//trim(variables(v)), NF90_DOUBLE, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "an"//trim(variables(v)), NF90_DOUBLE, (/t_dimid, lev_dimid/), varid))
      end do
      call check(nf90_def_var(ncid,"HR", NF90_DOUBLE, (/t_dimid, lev_dimid /), varid ))

      do v = 1, size(name_fluxes)
        call check(nf90_def_var(ncid, trim(name_fluxes(v)), NF90_double, (/t_dimid, lev_dimid /), varid))
      end do



      call check(nf90_enddef(ncid))
      call check( nf90_close(ncid) )
    end subroutine create_netcdf

    subroutine fill_netcdf(run_name, soil_levels, time, pool_matrix, change_matrix, a_matrix, HR)
      character (len = *):: run_name

      integer :: soil_levels, time, i , j, varidchange
      real(r8), intent(in)                      :: pool_matrix(soil_levels,pool_types)   ! For storing C pool sizes [gC/m3]
      real(r8),intent(in)                       :: change_matrix(soil_levels,pool_types) ! For storing dC/dt for each time step [gC/(m3*day)]
      real(r8),intent(in)                       :: a_matrix(soil_levels,pool_types) ! For storing analytical solution
      real(r8), dimension(1)         :: HR
      call check(nf90_open(trim(run_name)//".nc", nf90_write, ncid))

      do j=1,soil_levels
        call check(nf90_inq_varid(ncid, "HR", varid))
        call check(nf90_put_var(ncid, varid, HR(j), start = (/time/48, j/)))
        do i = 1,pool_types

          call check(nf90_inq_varid(ncid, trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, pool_matrix(j,i), start = (/ time/48, j /)))

          call check(nf90_inq_varid(ncid, trim(change_variables(i)), varidchange))
          call check(nf90_put_var(ncid, varidchange, change_matrix(j,i), start = (/ time/48, j /)))

          call check(nf90_inq_varid(ncid, trim(an_variables(i)), varidan))
          call check(nf90_put_var(ncid, varidan, a_matrix(j,i), start = (/ time/48, j /)))
        end do
      end do

      call check(nf90_close(ncid))
    end subroutine fill_netcdf

   subroutine store_parameters(run_name)
     character (len = *):: run_name
     integer :: tsoiID, clayID, desorbID, MgeID, kmID, vmID, fmetID, tauID, gepID, depthID, fphysID, fchemID, favailID
     call check(nf90_open(trim(run_name)//".nc", nf90_write, ncid))

     call check(nf90_def_var(ncid, "tsoi", NF90_double, tsoiID))
     call check(nf90_def_var(ncid, "f_clay", NF90_double, clayID))
     call check(nf90_def_var(ncid, "desorb", NF90_double, desorbID))
     call check(nf90_def_var(ncid, "MGE",NF90_double,MGE_dimid, mgeID))
     call check(nf90_def_var(ncid, "Km",NF90_double,mmk_dimid, kmID))
     call check(nf90_def_var(ncid, "Vmax", NF90_double,mmk_dimid,vmID))
     call check(nf90_def_var(ncid, "GEP", NF90_double,gepID))
     call check(nf90_def_var(ncid, "tau", NF90_double,tauID))
     call check(nf90_def_var(ncid, "f_phys", NF90_double,fphysID))
     call check(nf90_def_var(ncid, "f_avail", NF90_double,favailID))
     call check(nf90_def_var(ncid, "f_chem", NF90_double,fchemID))
     call check(nf90_def_var(ncid, "depth", NF90_double,depthID))
     call check(nf90_def_var(ncid, "f_met", NF90_double,fmetID))

     call check(nf90_enddef(ncid))

     call check(nf90_put_var(ncid, tsoiID, tsoi))
     call check(nf90_put_var(ncid, clayID, fCLAY))
     call check(nf90_put_var(ncid, desorbID, desorb))
     call check(nf90_put_var(ncid, mgeID, MGE))
     call check(nf90_put_var(ncid, gepID, GEP))
     call check(nf90_put_var(ncid, fphysID, fPHYS))
     call check(nf90_put_var(ncid, fchemID, fCHEM))
     call check(nf90_put_var(ncid, favailID, fAVAIL))
     call check(nf90_put_var(ncid, depthID, depth))
     call check(nf90_put_var(ncid, fmetID, fMET))
     call check(nf90_put_var(ncid, vmID, Vmax))
     call check(nf90_put_var(ncid, kmID, Km))
     call check(nf90_put_var(ncid, tauID, tau))







     call check(nf90_close(ncid))
   end subroutine store_parameters

    ! subroutine fluxes_netcdf()
    !   do i = 1, size(name_fluxes)
    !     call check(nf90_inq_varid(ncid, trim(name_fluxes(i)), varid))
    !     call check(nf90_put_var(ncid, varid, ))
    !   end do
    ! end subroutine fluxes_netcdf

    subroutine openOutputFile(name_ad, isVertical)
      character (len=*) :: name_ad
      character (len=37),parameter :: path='/home/ecaas/decomposition/ecosystems/'
      logical, intent(in) :: isVertical
      open(unit=1,file = trim(path)//trim(name_ad)//"_pool.txt",   form="formatted", action="write", status='replace', iostat=ios)
      open(unit=15,file = trim(path)//trim(name_ad)//"_an.txt",   form="formatted", action="write", status='replace', iostat=ios)
      open(unit=2,file = trim(path)//trim(name_ad)//"_change.txt", form="formatted", action="write", status='replace', iostat=ios)
      open(unit=3,file = trim(path)//trim(name_ad)//"_LITflux.txt", form="formatted", action="write", status='replace', iostat=ios)
      open(unit=4,file = trim(path)//trim(name_ad)//"_SAPSOMflux.txt", form="formatted", action="write", status='replace', iostat=ios)
      open(unit=7,file = trim(path)//trim(name_ad)//"_MYCSAPflux.txt", form="formatted", action="write", status='replace', iostat=ios)
      open(unit=8,file = trim(path)//trim(name_ad)//"_MYCSOMflux.txt", form="formatted", action="write", status='replace', iostat=ios)
      open(unit=9,file = trim(path)//trim(name_ad)//"_SOMflux.txt", form="formatted", action="write", status='replace', iostat=ios)
      write(unit=3, fmt=*) 'time, depth_level,LITmSAPr,LITmSAPk,LITsSAPr,LITsSAPk'
      write(unit=4, fmt=*) 'time, depth_level,SAPrSOMp,SAPrSOMa,SAPrSOMc,SAPkSOMp,SAPkSOMa,SAPkSOMc'
      write(unit=7, fmt=*) 'time, depth_level,EcM_SAPr,EcM_SAPk,ErM_SAPr,ErM_SAPk,AM_SAPr,AM_SAPk'
      write(unit=8, fmt=*) 'time,depth_level,EcM_SOMp, EcMSOMa,EcMSOMc, ErMSOMp,ErMSOMa,ErMSOMc,AMSOMp,AMSOMa,AMSOMc'
      Write(unit=9, fmt=*) 'time,depth_level,SOMaSAPr,SOMaSAPk,SOMpSOMa,SOMcSOMa'!header


      if (isVertical) then
        open(unit=10,file = trim(path)//trim(name_ad)//"_vertical.txt", form="formatted", action="write", status='replace', iostat=ios)
        Write(unit=10, fmt=*) 'time,depth_level,pool_nr,net_vertical_transport, transport_upper, transport_lower'!header
      end if

      if( ios/=0) then
        write(6,*) 'Error opening file for writing'
        stop
      endif
      write(unit=1, fmt=*) 'time, depth_level, LITm,  LITs,  SAPr,  SAPk, EcM, ErM, AM,  SOMp,  SOMa,  SOMc'!header
      write(unit=15, fmt=*) 'time, depth_level, LITm,  LITs,  SAPr,  SAPk, EcM, ErM, AM,  SOMp,  SOMa,  SOMc'!header
      write(unit=2, fmt=*) 'time, depth_level, HR, LITm,  LITs,  SAPr,  SAPk, EcM, ErM, AM,  SOMp,  SOMa,  SOMc'!header
    end subroutine openOutputFile

    subroutine closeFiles(isVertical)
      logical, intent(in) :: isVertical
      close(unit=1)
      close(unit=2)
      close(unit=3)
      close(unit=4)
      close(unit=7)
      close(unit=8)
      close(unit=9)
      close(unit=15)
      if (isVertical) then
        close(unit=10)
      end if
    end subroutine closeFiles

end module writeMod
