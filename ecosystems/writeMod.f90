module writeMod
  use paramMod
  use netcdf
  use dispmodule
  implicit none
  ! varidan
  integer :: grid_dimid, col_dimid, t_dimid, lev_dimid,mmk_dimid, MGE_dimid,fracid,i


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
      integer :: ncid, varid
      integer, parameter :: NDIMS = 4 !gridcell, column, levsoi, time, nommeqs, nomgevalues,sappools
      integer, parameter :: gridcell = 1, column = 1, levsoi = nlevdecomp
      integer :: x,v

      ! if (inquire(file = trim(run_name)//".nc", Exist = file_exists)) then
      !   call check(nf90_create(trim(run_name)//".nc", NF90_CLOBBER,NF90_NETCDF4 ncid))
      ! else
      !   call check(nf90_create(trim(run_name)//".nc", NF90_CLOBBER, ncid))
      ! end if
      call check(nf90_create(output_path//trim(run_name)//".nc",NF90_NETCDF4,ncid))


      call check(nf90_def_dim(ncid, "time", nf90_unlimited, t_dimid))
      call check(nf90_def_dim(ncid, "gridcell", gridcell, grid_dimid))
      call check(nf90_def_dim(ncid, "column", column, col_dimid))
      call check(nf90_def_dim(ncid, "levsoi", levsoi, lev_dimid))
      call check(nf90_def_dim(ncid, "NoMMKeqs", MM_eqs, mmk_dimid))
      call check(nf90_def_dim(ncid, "NoMGEvalues", size(MGE), MGE_dimid))
      call check(nf90_def_dim(ncid,"SAPpools",size(fPHYS),fracid))

      do v = 1, size(variables)
        call check(nf90_def_var(ncid, trim(variables(v)), NF90_DOUBLE, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "change"//trim(variables(v)), NF90_DOUBLE, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "an"//trim(variables(v)), NF90_DOUBLE, (/t_dimid, lev_dimid/), varid))
      end do
      call check(nf90_def_var(ncid,"HR", NF90_DOUBLE, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_def_var(ncid, "vert_change", NF90_DOUBLE, (/t_dimid, lev_dimid/), varid))
      do v = 1, size(name_fluxes)
        call check(nf90_def_var(ncid, trim(name_fluxes(v)), NF90_double, (/t_dimid, lev_dimid /), varid))
      end do
      call check(nf90_def_var(ncid, "time", NF90_DOUBLE, (/t_dimid /), varid))
      call check(nf90_def_var(ncid, "month", NF90_DOUBLE, (/t_dimid /), varid))


      call check(nf90_enddef(ncid))
      call check( nf90_close(ncid) )
    end subroutine create_netcdf

    subroutine fill_netcdf(run_name, soil_levels, time, pool_matrix, change_matrix, a_matrix, HR, vert_sum, write_hour,month)
      character (len = *):: run_name

      integer :: soil_levels, time, i , j, varidchange,varid,ncid,varidan, timestep,vertid,write_hour
      real(r8), intent(in)                       :: pool_matrix(soil_levels,pool_types)   ! For storing C pool sizes [gC/m3]
      real(r8), intent(in)                       :: change_matrix(soil_levels,pool_types) ! For storing dC/dt for each time step [gC/(m3*day)]
      real(r8), intent(in)                       :: a_matrix(soil_levels,pool_types) ! For storing analytical solution
      real(r8), intent(in)                       :: vert_sum(soil_levels, pool_types)
      integer, intent(in)                       :: month
      real(r8), dimension(soil_levels)           :: HR,HR_sum
      if (time == 1) then
        !print*, "INSIDE"
        timestep = 1
      else
        timestep = time/write_hour+1
      end if
      !  print*, HR_sum
      call check(nf90_open(output_path//trim(run_name)//".nc", nf90_write, ncid))

      call check(nf90_inq_varid(ncid, "time", varid))
      call check(nf90_put_var(ncid, varid, time, start = (/ timestep /)))

      call check(nf90_inq_varid(ncid, "month", varid))
      call check(nf90_put_var(ncid, varid, month, start = (/ timestep /)))
      do j=1,soil_levels
        call check(nf90_inq_varid(ncid, "HR", varid))
      !  print*,  HR_sum(j)
        call check(nf90_put_var(ncid, varid, HR(j), start = (/timestep, j/)))

        do i = 1,pool_types
          !print*, "INSIDE again"
          call check(nf90_inq_varid(ncid, trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, pool_matrix(j,i), start = (/ timestep, j /)))
          call check(nf90_inq_varid(ncid, trim(change_variables(i)), varidchange))
          call check(nf90_put_var(ncid, varidchange, change_matrix(j,i), start = (/ timestep, j /)))
          call check(nf90_inq_varid(ncid, trim(an_variables(i)), varidan))
          call check(nf90_put_var(ncid, varidan, a_matrix(j,i), start = (/ timestep, j /)))
          call check(nf90_inq_varid(ncid, "vert_change", vertid))
          call check(nf90_put_var(ncid, vertid, vert_sum(j,i), start = (/timestep,j/)))
        end do
      end do
      call check(nf90_close(ncid))
    end subroutine fill_netcdf

   subroutine store_parameters(run_name)
     character (len = *):: run_name
     integer :: tsoiID, clayID, desorbID, MgeID, kmID, vmID, fmetID, tauID, gepID, depthID, fphysID, fchemID, favailID
     integer:: ncid
     call check(nf90_open(output_path//trim(run_name)//".nc", nf90_write, ncid))

     call check(nf90_def_var(ncid, "tsoi", NF90_double, tsoiID))
     call check(nf90_def_var(ncid, "f_clay", NF90_double, clayID))
     call check(nf90_def_var(ncid, "desorb", NF90_double, desorbID))
     call check(nf90_def_var(ncid, "MGE",NF90_double,MGE_dimid, mgeID))
     call check(nf90_def_var(ncid, "Km",NF90_double,mmk_dimid, kmID))
     call check(nf90_def_var(ncid, "Vmax", NF90_double,mmk_dimid,vmID))
     call check(nf90_def_var(ncid, "GEP", NF90_double,gepID))
     call check(nf90_def_var(ncid, "tau", NF90_double,fracid,tauID))

     call check(nf90_def_var(ncid, "f_phys", NF90_double,fracid,fphysID))
     call check(nf90_def_var(ncid, "f_avail", NF90_double,fracid,favailID))
     call check(nf90_def_var(ncid, "f_chem", NF90_double,fracid,fchemID))
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

    subroutine openOutputFile(name_ad, isVertical)!FOR TEXTFILES
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

    subroutine closeFiles(isVertical)!FOR TEXTFILES
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

    subroutine read_clmdata(clm_history_file, TSOI, SOILLIQ,SOILICE,WATSAT,month)
      character (len = *):: clm_history_file
      real(r8),intent(out), dimension(nlevdecomp)           :: TSOI
      real(r8),intent(out), dimension(nlevdecomp)           :: SOILLIQ
      real(r8), intent(out), dimension(nlevdecomp)          :: SOILICE
      real(r8), intent(out), dimension(nlevdecomp)          :: WATSAT
      integer            :: ncid, WATSATid, TSOIid, SOILICEid, SOILLIQid, month
      WATSAT=0.0
      TSOI= 0.0
      call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))
      call check(nf90_inq_varid(ncid, 'WATSAT', WATSATid))
      call check(nf90_get_var(ncid, WATSATid, WATSAT, count=(/1,1,nlevdecomp/)))

      call check(nf90_inq_varid(ncid, 'TSOI', TSOIid))
      call check(nf90_get_var(ncid, TSOIid, TSOI, start=(/1,1,1, month/), count=(/1,1,nlevdecomp,1/)))

      call check(nf90_inq_varid(ncid, 'SOILLIQ', SOILLIQid))
      call check(nf90_get_var(ncid, SOILLIQid, SOILLIQ, start=(/1,1,1, month/), count=(/1,1,nlevdecomp,1/)))

      call check(nf90_inq_varid(ncid, 'SOILICE', SOILICEid))
      call check(nf90_get_var(ncid, SOILICEid, SOILICE, start=(/1,1,1, month/), count=(/1,1,nlevdecomp,1/)))

      !Unit conversions:
      TSOI = TSOI - 273.15 !Kelvin to Celcius
      !PRINT*, month
      !print*, TSOI
      !call disp(TSOI)
      do i = 1, nlevdecomp
        SOILICE(i) = SOILICE(i)/(delta_z(i)*917) !kg/m2 to m3/m3 rho_ice=917kg/m3
        SOILLIQ(i) = SOILLIQ(i)/(delta_z(i)*1000) !kg/m2 to m3/m3 rho_liq=1000kg/m3
      end do

      !call disp(SOILLIQ)

      call check(nf90_close(ncid))
    end subroutine read_clmdata

end module writeMod
