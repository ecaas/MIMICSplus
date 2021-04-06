module writeMod
  use paramMod
  use netcdf
  use dispmodule
  implicit none

  integer,private :: grid_dimid, col_dimid, t_dimid, lev_dimid,mmk_dimid,fracid,i

  contains

    subroutine check(status)
      integer, intent ( in) :: status
      if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 2
      end if
    end subroutine check

    subroutine create_netcdf(run_name, levsoi)
      character (len = *):: run_name
      integer :: ncid, varid
      integer, parameter :: gridcell = 1, column = 1
      integer :: levsoi
      integer :: v
      call check(nf90_create(output_path//trim(run_name)//".nc",NF90_NETCDF4,ncid))

      call check(nf90_def_dim(ncid, "time", nf90_unlimited, t_dimid))
      call check(nf90_def_dim(ncid, "gridcell", gridcell, grid_dimid))
      call check(nf90_def_dim(ncid, "column", column, col_dimid))
      call check(nf90_def_dim(ncid, "levsoi", levsoi, lev_dimid))
      call check(nf90_def_dim(ncid, "NoMMKeqs", MM_eqs, mmk_dimid))
      call check(nf90_def_dim(ncid,"SAPpools",size(fPHYS),fracid))
      do v = 1, size(variables)
        call check(nf90_def_var(ncid, trim(variables(v)), NF90_DOUBLE, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "change"//trim(variables(v)), NF90_DOUBLE, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "vert_change"//trim(variables(v)), NF90_DOUBLE, (/t_dimid, lev_dimid/), varid))
        call check(nf90_def_var(ncid, "N_"//trim(variables(v)), NF90_DOUBLE, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "N_change"//trim(variables(v)), NF90_DOUBLE, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "N_vert_change"//trim(variables(v)), NF90_DOUBLE, (/t_dimid, lev_dimid/), varid))

      end do
      call check(nf90_def_var(ncid, "N_inorganic", NF90_DOUBLE, (/t_dimid, lev_dimid /), varid))
      call check(nf90_def_var(ncid, "N_plant", NF90_DOUBle, (/t_dimid/),varid))
      call check(nf90_def_var(ncid, "C_plant", NF90_DOUBle, (/t_dimid/),varid))
      call check(nf90_def_var(ncid, "C_Growth_sum", NF90_DOUBle, (/t_dimid/),varid))
      call check(nf90_def_var(ncid, "C_Growth_flux", NF90_DOUBle, (/t_dimid/),varid))
      call check(nf90_def_var(ncid,"HR_sum", NF90_DOUBLE, (/t_dimid /), varid ))
      call check(nf90_def_var(ncid,"HR_flux", NF90_DOUBLE, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_def_var(ncid, "Temp", NF90_DOUBLE, (/t_dimid, lev_dimid/),varid))
      call check(nf90_def_var(ncid, "Moisture", NF90_DOUBLE, (/t_dimid, lev_dimid/),varid))
      call check(nf90_def_var(ncid, "N_changeinorganic", NF90_DOUBLE,(/t_dimid, lev_dimid/), varid))

      do v = 1, size(C_name_fluxes)
         call check(nf90_def_var(ncid, "C_"//trim(C_name_fluxes(v)), NF90_double, (/t_dimid, lev_dimid /), varid))
      end do

      do v = 1, size(N_name_fluxes)
         call check(nf90_def_var(ncid, "N_"//trim(N_name_fluxes(v)), NF90_double, (/t_dimid, lev_dimid /), varid))
      end do

      call check(nf90_def_var(ncid, "time", NF90_DOUBLE, (/t_dimid /), varid))
      call check(nf90_def_var(ncid, "month", NF90_DOUBLE, (/t_dimid /), varid))
      call check(nf90_enddef(ncid))

      call check( nf90_close(ncid) )
    end subroutine create_netcdf

    subroutine fill_netcdf(run_name, time, pool_matrix, change_matrix, Npool_matrix, Nchange_matrix, &
      HR_sum, HR_flux, vert_sum,Nvert_sum, write_hour,month, N_plant, C_plant, TSOIL, MOIST,growth_sum,levsoi)
      character (len = *):: run_name
      integer:: levsoi
      integer :: time, i , j, varidchange,varid,ncid, timestep,vertid,write_hour
      real(r8), intent(in)          :: pool_matrix(levsoi,pool_types), Npool_matrix(levsoi,pool_types_N)   ! For storing C pool sizes [gC/m3]
      real(r8), intent(in)                       :: N_plant, C_plant,growth_sum
      real(r8), intent(in)                       :: change_matrix(levsoi,pool_types), Nchange_matrix(levsoi,pool_types_N) ! For storing dC/dt for each time step [gC/(m3*day)]
      real(r8), intent(in)                       :: vert_sum(levsoi,pool_types)
      real(r8), intent(in)                       :: Nvert_sum(levsoi,pool_types)

      integer, intent(in)                        :: month
      real(r8)                                   :: HR_sum
      real(r8) ,dimension(levsoi)                :: HR_flux!(levsoi) !HR_mass_accumulated
      real(r8),dimension(levsoi)         :: TSOIL, MOIST


      call get_timestep(time, write_hour, timestep)
      call check(nf90_open(output_path//trim(run_name)//".nc", nf90_write, ncid))
      !print*, time, timestep, "write"
      call check(nf90_inq_varid(ncid, "time", varid))
      call check(nf90_put_var(ncid, varid, time, start = (/ timestep /)))
      call check(nf90_inq_varid(ncid, "month", varid))
      call check(nf90_put_var(ncid, varid, month, start = (/ timestep /)))
      call check(nf90_inq_varid(ncid, "N_plant", varid))
      call check(nf90_put_var(ncid, varid, N_Plant, start = (/ timestep /)))
      call check(nf90_inq_varid(ncid, "C_plant", varid))
      call check(nf90_put_var(ncid, varid, C_Plant, start = (/ timestep /)))

      call check(nf90_inq_varid(ncid, "C_Growth_sum", varid))
      call check(nf90_put_var(ncid, varid, growth_sum, start = (/ timestep /)))

      call check(nf90_inq_varid(ncid, "HR_sum", varid))
      call check(nf90_put_var(ncid, varid, HR_sum, start = (/ timestep /)))

      do j=1,levsoi
        call check(nf90_inq_varid(ncid, "Temp",varid))
        call check(nf90_put_var(ncid, varid, TSOIL(j), start = (/timestep,j/)))

        call check(nf90_inq_varid(ncid, "Moisture",varid))
        call check(nf90_put_var(ncid, varid, MOIST(j), start = (/timestep,j/)))

        call check(nf90_inq_varid(ncid, "HR_flux", varid))
        call check(nf90_put_var(ncid, varid, HR_flux(j), start = (/timestep, j/)))

        call check(nf90_inq_varid(ncid, "N_inorganic", varid))
        call check(nf90_put_var(ncid, varid, Npool_matrix(j,11), start = (/timestep, j/)))

        call check(nf90_inq_varid(ncid, "N_changeinorganic", varid))
        call check(nf90_put_var(ncid, varid, Nchange_matrix(j,11), start = (/timestep, j/)))
        do i = 1,pool_types
          !C:
          call check(nf90_inq_varid(ncid, trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, pool_matrix(j,i), start = (/ timestep, j /)))
          !N:
          call check(nf90_inq_varid(ncid, "N_"//trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, Npool_matrix(j,i), start = (/ timestep, j /)))

          !C change:
          call check(nf90_inq_varid(ncid, trim(change_variables(i)), varidchange))
          call check(nf90_put_var(ncid, varidchange, change_matrix(j,i), start = (/ timestep, j /)))
          call check(nf90_inq_varid(ncid, "vert_change"//trim(variables(i)), vertid))
          call check(nf90_put_var(ncid, vertid, vert_sum(j,i), start = (/timestep,j/)))
          !N change:
          call check(nf90_inq_varid(ncid, "N_"//trim(change_variables(i)), varidchange))
          call check(nf90_put_var(ncid, varidchange, Nchange_matrix(j,i), start = (/ timestep, j /)))
          call check(nf90_inq_varid(ncid, "N_vert_change"//trim(variables(i)), vertid))
          call check(nf90_put_var(ncid, vertid, Nvert_sum(j,i), start = (/timestep,j/)))
        end do !pool_types
      end do ! levels
      call check(nf90_close(ncid))
    end subroutine fill_netcdf

   subroutine store_parameters(run_name)
     character (len = *):: run_name
     integer :: tsoiID, clayID, desorbID, kmID, vmID, fmetID, tauID, gepID, depthID, fphysID, fchemID, favailID
     integer:: ncid
     call check(nf90_open(output_path//trim(run_name)//".nc", nf90_write, ncid))

     call check(nf90_def_var(ncid, "tsoi", NF90_double, tsoiID))
     call check(nf90_def_var(ncid, "f_clay", NF90_double, clayID))
     call check(nf90_def_var(ncid, "desorb", NF90_double, desorbID))
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
     call check(nf90_put_var(ncid, gepID, GEP))
     call check(nf90_put_var(ncid, fphysID, fPHYS))
     call check(nf90_put_var(ncid, fchemID, fCHEM))
     call check(nf90_put_var(ncid, favailID, fAVAIL))
     call check(nf90_put_var(ncid, depthID, soil_depth))
     call check(nf90_put_var(ncid, fmetID, fMET))
     call check(nf90_put_var(ncid, vmID, Vmax))
     call check(nf90_put_var(ncid, kmID, Km))
     call check(nf90_put_var(ncid, tauID, tau))

     call check(nf90_close(ncid))
   end subroutine store_parameters


   !NOTE: This should maybe be somwhere else?
    subroutine read_clmdata(clm_history_file, TSOI, SOILLIQ,SOILICE,WATSAT,W_SCALAR,month, nlevdecomp)
      integer,intent(in)            :: nlevdecomp
      character (len = *),intent(in):: clm_history_file
      real(r8),intent(out), dimension(nlevdecomp)          :: TSOI
      real(r8),intent(out), dimension(nlevdecomp)          :: SOILLIQ
      real(r8), intent(out),dimension(nlevdecomp)          :: SOILICE
      real(r8), intent(out),dimension(nlevdecomp)          :: WATSAT
      real(r8), intent(out),dimension(nlevdecomp)          :: W_SCALAR

      integer            :: ncid, WATSATid, TSOIid, SOILICEid, SOILLIQid,W_SCALARid, month
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

      call check(nf90_inq_varid(ncid, 'W_SCALAR', W_SCALARid))
      call check(nf90_get_var(ncid, W_SCALARid, W_SCALAR, start=(/1,1,1, month/), count=(/1,1,nlevdecomp,1/)))
      !Unit conversions:
      TSOI = TSOI - 273.15 !Kelvin to Celcius
      do i = 1, nlevdecomp
        SOILICE(i) = SOILICE(i)/(delta_z(i)*917) !kg/m2 to m3/m3 rho_ice=917kg/m3
        SOILLIQ(i) = SOILLIQ(i)/(delta_z(i)*1000) !kg/m2 to m3/m3 rho_liq=1000kg/m3
      end do
      call check(nf90_close(ncid))
    end subroutine read_clmdata

    subroutine get_timestep(time, write_hour, timestep)
      integer, intent(in) :: time
      integer, intent(in) :: write_hour !hours between every output-writing.
      integer, intent(out):: timestep
      if (time == 1) then
        timestep = 1
      else
        timestep = time/write_hour+1 !NOTE: Del write_hour på step_frac hvis step_frac ikke er lik 1!
      end if
    end subroutine get_timestep

    subroutine fluxes_netcdf(time, write_hour, depth_level, run_name)
      integer, intent(in) :: time
      integer, intent(in) :: write_hour
      integer, intent(in) :: depth_level
      character (len = *), intent(in):: run_name

      !Local:
      integer :: timestep
      integer :: ncid
     !---------------------------------------------
      call get_timestep(time, write_hour, timestep)
      call check(nf90_open(output_path//trim(run_name)//".nc", nf90_write, ncid))
      call write_Nfluxes(ncid, timestep, depth_level)
      call write_Cfluxes(ncid, timestep, depth_level)
      call check(nf90_close(ncid))
    end subroutine fluxes_netcdf

    subroutine write_Cfluxes(ncid,timestep,depth_level)
      integer, intent(in) :: ncid
      integer, intent(in) :: timestep
      integer, intent(in) :: depth_level

      !Local:
      integer :: varid

      call check(nf90_inq_varid(ncid, "C_Growth_flux", varid))
      call check(nf90_put_var(ncid, varid, C_growth_rate, start = (/ timestep /)))
      call check(nf90_inq_varid(ncid, "C_PlantLITm", varid))
      call check(nf90_put_var(ncid, varid, C_PlantLITm, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_PlantLITs", varid))
      call check(nf90_put_var(ncid, varid, C_PlantLITs, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_PlantEcM", varid))
      call check(nf90_put_var(ncid, varid, C_PlantEcM, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_PlantErM", varid))
      call check(nf90_put_var(ncid, varid, C_PlantErM, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_PlantAM", varid))
      call check(nf90_put_var(ncid, varid, C_PlantAM, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_LITmSAPb", varid))
      call check(nf90_put_var(ncid, varid, C_LITmSAPb, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_LITsSAPb", varid))
      call check(nf90_put_var(ncid, varid, C_LITsSAPb, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_LITmSAPf", varid))
      call check(nf90_put_var(ncid, varid, C_LITmSAPf, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_LITsSAPf", varid))
      call check(nf90_put_var(ncid, varid, C_LITsSAPf, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_SOMpSOMa", varid))
      call check(nf90_put_var(ncid, varid, C_SOMpSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_SOMcSOMa", varid))
      call check(nf90_put_var(ncid, varid, C_SOMcSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_SAPbSOMa", varid))
      call check(nf90_put_var(ncid, varid, C_SAPbSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_SAPfSOMa", varid))
      call check(nf90_put_var(ncid, varid, C_SAPfSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_EcMSOMa", varid))
      call check(nf90_put_var(ncid, varid, C_EcMSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_EcMSOMc", varid))
      call check(nf90_put_var(ncid, varid, C_EcMSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_EcMSOMp", varid))
      call check(nf90_put_var(ncid, varid, C_EcMSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_ErMSOMa", varid))
      call check(nf90_put_var(ncid, varid, C_ErMSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_ErMSOMc", varid))
      call check(nf90_put_var(ncid, varid, C_ErMSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_ErMSOMp", varid))
      call check(nf90_put_var(ncid, varid, C_ErMSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_AMSOMa", varid))
      call check(nf90_put_var(ncid, varid, C_AMSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_AMSOMc", varid))
      call check(nf90_put_var(ncid, varid, C_AMSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_SOMaSAPb", varid))
      call check(nf90_put_var(ncid, varid, C_SOMaSAPb, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_SOMaSAPf", varid))
      call check(nf90_put_var(ncid, varid, C_SOMaSAPf, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_SAPbSOMp", varid))
      call check(nf90_put_var(ncid, varid, C_SAPbSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_SAPfSOMp", varid))
      call check(nf90_put_var(ncid, varid, C_SAPfSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_SAPbSOMc", varid))
      call check(nf90_put_var(ncid, varid, C_SAPbSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_SAPfSOMc", varid))
      call check(nf90_put_var(ncid, varid, C_SAPfSOMc, start = (/ timestep, depth_level /)))
    end subroutine write_Cfluxes
    subroutine write_Nfluxes(ncid,timestep,depth_level)
      integer, intent(in) :: ncid
      integer, intent(in) :: timestep
      integer, intent(in) :: depth_level
      !Local:
      integer :: varid
      call check(nf90_inq_varid(ncid, "N_LITmSAPb", varid))
      call check(nf90_put_var(ncid, varid, N_LITmSAPb, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_LITsSAPb", varid))
      call check(nf90_put_var(ncid, varid, N_LITsSAPb, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_LITmSAPf", varid))
      call check(nf90_put_var(ncid, varid, N_LITmSAPf, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_LITsSAPf", varid))
      call check(nf90_put_var(ncid, varid, N_LITsSAPf, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SOMpSOMa", varid))
      call check(nf90_put_var(ncid, varid, N_SOMpSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SOMcSOMa", varid))
      call check(nf90_put_var(ncid, varid, N_SOMcSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPbSOMa", varid))
      call check(nf90_put_var(ncid, varid, N_SAPbSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPfSOMa", varid))
      call check(nf90_put_var(ncid, varid, N_SAPfSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_EcMSOMa", varid))
      call check(nf90_put_var(ncid, varid, N_EcMSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_EcMSOMc", varid))
      call check(nf90_put_var(ncid, varid, N_EcMSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_EcMSOMp", varid))
      call check(nf90_put_var(ncid, varid, N_EcMSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_ErMSOMa", varid))
      call check(nf90_put_var(ncid, varid, N_ErMSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_ErMSOMc", varid))
      call check(nf90_put_var(ncid, varid, N_ErMSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_ErMSOMp", varid))
      call check(nf90_put_var(ncid, varid, N_ErMSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_AMSOMa", varid))
      call check(nf90_put_var(ncid, varid, N_AMSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_AMSOMc", varid))
      call check(nf90_put_var(ncid, varid, N_AMSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SOMaSAPb", varid))
      call check(nf90_put_var(ncid, varid, N_SOMaSAPb, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SOMaSAPf", varid))
      call check(nf90_put_var(ncid, varid, N_SOMaSAPf, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPbSOMp", varid))
      call check(nf90_put_var(ncid, varid, N_SAPbSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPfSOMp", varid))
      call check(nf90_put_var(ncid, varid, N_SAPfSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPbSOMc", varid))
      call check(nf90_put_var(ncid, varid, N_SAPbSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPfSOMc", varid))
      call check(nf90_put_var(ncid, varid, N_SAPfSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPfIN", varid))
      call check(nf90_put_var(ncid, varid, N_SAPfIN, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPbIN", varid))
      call check(nf90_put_var(ncid, varid, N_SAPbIN, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_INAM", varid))
      call check(nf90_put_var(ncid, varid, N_INAM, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_INEcM", varid))
      call check(nf90_put_var(ncid, varid, N_INEcM, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_INErM", varid))
      call check(nf90_put_var(ncid, varid, N_INErM, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_Deposition", varid))
      call check(nf90_put_var(ncid, varid, Deposition, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_Leaching", varid))
      call check(nf90_put_var(ncid, varid, Leaching, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_PlantLITm", varid))
      call check(nf90_put_var(ncid, varid, N_PlantLITm, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_PlantLITs", varid))
      call check(nf90_put_var(ncid, varid, N_PlantLITs, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_EcMPlant", varid))
      call check(nf90_put_var(ncid, varid, N_EcMPlant, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_ErMPlant", varid))
      call check(nf90_put_var(ncid, varid, N_ErMPlant, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_AMPlant", varid))
      call check(nf90_put_var(ncid, varid, N_AMPlant, start = (/ timestep, depth_level /)))



    end subroutine write_Nfluxes
end module writeMod
