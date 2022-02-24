module writeMod
  use paramMod
  use netcdf
  use dispmodule
  implicit none

  integer,private :: grid_dimid, col_dimid, t_dimid, lev_dimid,mmk_dimid,fracid

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
        call check(nf90_def_var(ncid, trim(variables(v)), NF90_FLOAT, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "vert_change"//trim(variables(v)), NF90_FLOAT, (/t_dimid, lev_dimid/), varid))
        call check(nf90_def_var(ncid, "N_"//trim(variables(v)), NF90_FLOAT, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "N_vert_change"//trim(variables(v)), NF90_FLOAT, (/t_dimid, lev_dimid/), varid))

      end do
      call check(nf90_def_var(ncid, "NH4", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      call check(nf90_def_var(ncid, "NO3", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      call check(nf90_def_var(ncid, "N_SMIN", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      call check(nf90_def_var(ncid,"HR_sum", NF90_FLOAT, (/t_dimid /), varid ))
      call check(nf90_def_var(ncid,"HR_flux", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_def_var(ncid,"CUEb", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_def_var(ncid,"CUEf", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_def_var(ncid,"CUE_ecm", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_def_var(ncid,"CUE_am", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_def_var(ncid,"ROI", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_def_var(ncid, "Temp", NF90_FLOAT, (/t_dimid, lev_dimid/),varid))
      call check(nf90_def_var(ncid, "Moisture", NF90_FLOAT, (/t_dimid, lev_dimid/),varid))
      !call check(nf90_def_var(ncid, "N_changeinorganic", NF90_FLOAT,(/t_dimid, lev_dimid/), varid))
      call check(nf90_def_var(ncid, "N_InPlant", NF90_FLOAT,(/t_dimid, lev_dimid/), varid))
      


      do v = 1, size(C_name_fluxes)
         call check(nf90_def_var(ncid, "C_"//trim(C_name_fluxes(v)), NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      end do

      do v = 1, size(N_name_fluxes)
         call check(nf90_def_var(ncid, "N_"//trim(N_name_fluxes(v)), NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      end do
      
      call check(nf90_def_var(ncid, "mcdate", NF90_INT,(/t_dimid/), varid))
      call check(nf90_def_var(ncid, "time", NF90_INT, (/t_dimid /), varid))
      call check(nf90_def_var(ncid, "month", NF90_INT, (/t_dimid /), varid))
      call check(nf90_enddef(ncid))

      call check( nf90_close(ncid) )
    end subroutine create_netcdf

    subroutine fill_netcdf(ncid, time, pool_matrix, change_matrix, Npool_matrix, Nchange_matrix, &
      mcdate,HR_sum, HR_flux, vert_sum,Nvert_sum, write_hour,month, TSOIL, MOIST,CUE_bacteria,CUE_fungi,CUE_ecm,CUE_am,levsoi,ROI)
      !INPUT:
      integer,intent(in)               :: ncid 
      real(r8), intent(in)             :: pool_matrix(levsoi,pool_types), Npool_matrix(levsoi,pool_types_N)   ! For storing C pool sizes [gC/m3]
      real(r8), intent(in)             :: change_matrix(levsoi,pool_types), Nchange_matrix(levsoi,pool_types_N) ! For storing dC/dt for each time step [gC/(m3*day)]
      real(r8), intent(in)             :: vert_sum(levsoi,pool_types)
      real(r8), intent(in)             :: Nvert_sum(levsoi,pool_types)
      integer, intent(in)              :: mcdate
      integer, intent(in)              :: write_hour
      integer, intent(in)              :: levsoi
      integer, intent(in)              :: month
      integer, intent(in)              :: time
      real(r8),intent(in)              :: HR_sum
      real(r8),dimension(levsoi), intent(in)       :: HR_flux
      real(r8),dimension(levsoi), intent(in)       :: ROI      
      real(r8),dimension(levsoi), intent(in)       :: TSOIL, MOIST  
      real(r8),dimension(levsoi), intent(in)       :: CUE_bacteria,CUE_fungi, CUE_ecm,CUE_am
          
      !OUTPUT:
      
      !LOCAL:
      integer                          :: i , j !for looping
      integer                          :: varidchange,varid, timestep,vertid 
      real(r8)                         :: N_SMIN


      call get_timestep(time, write_hour, timestep)

      call check(nf90_inq_varid(ncid, "time", varid))
      call check(nf90_put_var(ncid, varid, time, start = (/ timestep /)))
      call check(nf90_inq_varid(ncid, "mcdate", varid))
      call check(nf90_put_var(ncid, varid, mcdate, start = (/ timestep /)))      
      call check(nf90_inq_varid(ncid, "month", varid))
      call check(nf90_put_var(ncid, varid, month, start = (/ timestep /)))
      call check(nf90_inq_varid(ncid, "HR_sum", varid))
      call check(nf90_put_var(ncid, varid, HR_sum, start = (/ timestep /)))

      do j=1,levsoi
        call check(nf90_inq_varid(ncid, "Temp",varid))
        call check(nf90_put_var(ncid, varid, TSOIL(j), start = (/timestep,j/)))

        call check(nf90_inq_varid(ncid, "Moisture",varid))
        call check(nf90_put_var(ncid, varid, MOIST(j), start = (/timestep,j/)))

        call check(nf90_inq_varid(ncid, "HR_flux", varid))
        call check(nf90_put_var(ncid, varid, HR_flux(j), start = (/timestep, j/)))
        
        call check(nf90_inq_varid(ncid, "ROI", varid))
        call check(nf90_put_var(ncid, varid, ROI, start = (/timestep, j/)))
        
        call check(nf90_inq_varid(ncid, "CUEb", varid))
        call check(nf90_put_var(ncid, varid, CUE_bacteria(j), start = (/timestep, j/)))
        
        call check(nf90_inq_varid(ncid, "CUEf", varid))
        call check(nf90_put_var(ncid, varid, CUE_fungi(j), start = (/timestep, j/)))
        
        call check(nf90_inq_varid(ncid, "CUE_ecm", varid))
        call check(nf90_put_var(ncid, varid, CUE_ecm(j), start = (/timestep, j/)))
        
        call check(nf90_inq_varid(ncid, "CUE_am", varid))
        call check(nf90_put_var(ncid, varid, CUE_am(j), start = (/timestep, j/)))
        
        call check(nf90_inq_varid(ncid, "NH4", varid))
        call check(nf90_put_var(ncid, varid, Npool_matrix(j,11), start = (/timestep, j/)))
        call check(nf90_inq_varid(ncid, "NO3", varid))
        call check(nf90_put_var(ncid, varid, Npool_matrix(j,12), start = (/timestep, j/)))

        N_SMIN = Npool_matrix(j,11)+Npool_matrix(j,12)
        call check(nf90_inq_varid(ncid, "N_SMIN", varid))
        call check(nf90_put_var(ncid, varid, N_SMIN, start = (/timestep, j/)))        

        do i = 1,pool_types
          !C:
          call check(nf90_inq_varid(ncid, trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, pool_matrix(j,i), start = (/ timestep, j /)))
          !N:
          call check(nf90_inq_varid(ncid, "N_"//trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, Npool_matrix(j,i), start = (/ timestep, j /)))

          !C change:
          call check(nf90_inq_varid(ncid, "vert_change"//trim(variables(i)), vertid))
          call check(nf90_put_var(ncid, vertid, vert_sum(j,i), start = (/timestep,j/)))
          !N change:
          call check(nf90_inq_varid(ncid, "N_vert_change"//trim(variables(i)), vertid))
          call check(nf90_put_var(ncid, vertid, Nvert_sum(j,i), start = (/timestep,j/)))
        end do !pool_types
      end do ! levels
    end subroutine fill_netcdf

   subroutine store_parameters(ncid)
     integer :: clayID, desorbID, fmetID, tauID, depthID, fphysID, fchemID, favailID,fecmID
     integer,intent(in):: ncid

     call check(nf90_def_var(ncid, "f_clay", NF90_FLOAT, clayID))
     call check(nf90_def_var(ncid, "desorb", NF90_FLOAT, desorbID))
     call check(nf90_def_var(ncid, "tau", NF90_FLOAT,fracid,tauID))
     call check(nf90_def_var(ncid, "f_phys", NF90_FLOAT,fracid,fphysID))
     call check(nf90_def_var(ncid, "f_avail", NF90_FLOAT,fracid,favailID))
     call check(nf90_def_var(ncid, "f_chem", NF90_FLOAT,fracid,fchemID))
     call check(nf90_def_var(ncid, "depth", NF90_FLOAT,depthID))
     call check(nf90_def_var(ncid, "f_met", NF90_FLOAT,fmetID))
     call check(nf90_def_var(ncid, "f_EcM", NF90_FLOAT,fecmID))
     

     call check(nf90_enddef(ncid))

     call check(nf90_put_var(ncid, clayID, fCLAY))
     call check(nf90_put_var(ncid, desorbID, desorb))
     call check(nf90_put_var(ncid, fphysID, fPHYS))
     call check(nf90_put_var(ncid, fchemID, fCHEM))
     call check(nf90_put_var(ncid, favailID, fAVAIL))
     call check(nf90_put_var(ncid, depthID, soil_depth))
     call check(nf90_put_var(ncid, fmetID, fMET))
     call check(nf90_put_var(ncid, tauID, tau))
     call check(nf90_put_var(ncid, fecmID, f_ecm))
     

   end subroutine store_parameters

    subroutine get_timestep(time, write_hour, timestep)
      integer, intent(in) :: time
      integer, intent(in) :: write_hour !hours between every output-writing.
      
      integer, intent(out):: timestep
      
      integer,parameter:: step_frac=10

      if (write_hour == 1) then
        timestep = time/write_hour
      else
        if (time == 1) then
          timestep = 1
        else
          timestep = time/(write_hour)+1 
        end if
      end if
    end subroutine get_timestep

    subroutine fluxes_netcdf(ncid,time, write_hour, depth_level)
      integer,intent(in)  :: ncid
      integer, intent(in) :: time
      integer, intent(in) :: write_hour
      integer, intent(in) :: depth_level

      !Local:
      integer :: timestep
     !---------------------------------------------
      call get_timestep(time, write_hour, timestep)
      call write_Nfluxes(ncid, timestep, depth_level)
      call write_Cfluxes(ncid, timestep, depth_level)
    end subroutine fluxes_netcdf

    subroutine write_Cfluxes(ncid,timestep,depth_level)
      
      !INPUT:
      integer, intent(in) :: ncid
      integer, intent(in) :: timestep
      integer, intent(in) :: depth_level

      !Local:
      integer :: varid

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
      call check(nf90_inq_varid(ncid, "C_PlantSOMc", varid))
      call check(nf90_put_var(ncid, varid, C_PlantSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_PlantSOMp", varid))
      call check(nf90_put_var(ncid, varid, C_PlantSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_PlantSOMa", varid))
      call check(nf90_put_var(ncid, varid, C_PlantSOMa, start = (/ timestep, depth_level /)))
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
      call check(nf90_inq_varid(ncid, "C_EcMdecoSOMp", varid))
      call check(nf90_put_var(ncid, varid, C_EcMdecompSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "C_EcMdecoSOMc", varid))
      call check(nf90_put_var(ncid, varid, C_EcMdecompSOMc, start = (/ timestep, depth_level /)))

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
      call check(nf90_inq_varid(ncid, "N_AMSOMp", varid))
      call check(nf90_put_var(ncid, varid, N_AMSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SOMaSAPb", varid))
      call check(nf90_put_var(ncid, varid, N_SOMaSAPb, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SOMaSAPf", varid))
      call check(nf90_put_var(ncid, varid, N_SOMaSAPf, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SOMaEcM", varid))
      call check(nf90_put_var(ncid, varid, N_SOMaEcM, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPbSOMp", varid))
      call check(nf90_put_var(ncid, varid, N_SAPbSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPfSOMp", varid))
      call check(nf90_put_var(ncid, varid, N_SAPfSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPbSOMc", varid))
      call check(nf90_put_var(ncid, varid, N_SAPbSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SAPfSOMc", varid))
      call check(nf90_put_var(ncid, varid, N_SAPfSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_INSAPf", varid))
      call check(nf90_put_var(ncid, varid, N_INSAPf, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_INSAPb", varid))
      call check(nf90_put_var(ncid, varid, N_INSAPb, start = (/ timestep, depth_level /)))
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
      call check(nf90_inq_varid(ncid, "N_PlantSOMp", varid))
      call check(nf90_put_var(ncid, varid, N_PlantSOMp, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_PlantSOMa", varid))
      call check(nf90_put_var(ncid, varid, N_PlantSOMa, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_PlantSOMc", varid))
      call check(nf90_put_var(ncid, varid, N_PlantSOMc, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_EcMPlant", varid))
      call check(nf90_put_var(ncid, varid, N_EcMPlant, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_ErMPlant", varid))
      call check(nf90_put_var(ncid, varid, N_ErMPlant, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_AMPlant", varid))
      call check(nf90_put_var(ncid, varid, N_AMPlant, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_InPlant", varid))
      call check(nf90_put_var(ncid, varid, N_InPlant, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SOMpEcM", varid))
      call check(nf90_put_var(ncid, varid, N_SOMpEcM, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_SOMcEcM", varid))
      call check(nf90_put_var(ncid, varid, N_SOMcEcM, start = (/ timestep, depth_level /)))
      call check(nf90_inq_varid(ncid, "N_nitrif_rate", varid))
      call check(nf90_put_var(ncid, varid, nitrif_rate, start = (/ timestep, depth_level /)))

    end subroutine write_Nfluxes

    subroutine create_yearly_mean_netcdf(run_name, levsoi)
      character (len = *), intent(in):: run_name
      integer :: ncid, varid
      integer, parameter :: gridcell = 1, column = 1
      integer :: levsoi
      integer :: v
      call check(nf90_create(output_path//trim(run_name)//"_yearly_mean.nc",NF90_NETCDF4,ncid))

      call check(nf90_def_dim(ncid, "time", nf90_unlimited, t_dimid))
      call check(nf90_def_dim(ncid, "gridcell", gridcell, grid_dimid))
      call check(nf90_def_dim(ncid, "column", column, col_dimid))
      call check(nf90_def_dim(ncid, "levsoi", levsoi, lev_dimid))

      do v = 1, size(variables)
        call check(nf90_def_var(ncid, trim(variables(v)), NF90_FLOAT, (/ t_dimid, lev_dimid /), varid))
        call check(nf90_def_var(ncid, "N_"//trim(variables(v)), NF90_FLOAT, (/ t_dimid, lev_dimid /), varid))
      end do

      call check(nf90_def_var(ncid, "N_inorganic", NF90_FLOAT, (/t_dimid, lev_dimid /), varid))
      call check(nf90_def_var(ncid,"HR_sum", NF90_FLOAT, (/t_dimid /), varid ))
      call check(nf90_def_var(ncid,"HR_flux", NF90_FLOAT, (/t_dimid, lev_dimid /), varid ))
      call check(nf90_def_var(ncid, "Temp", NF90_FLOAT, (/t_dimid, lev_dimid/),varid))
      call check(nf90_def_var(ncid, "Moisture", NF90_FLOAT, (/t_dimid, lev_dimid/),varid))

      call check(nf90_def_var(ncid, "year_since_start", NF90_FLOAT, (/t_dimid /), varid))
      call check(nf90_enddef(ncid))

      call check( nf90_close(ncid) )
    end subroutine create_yearly_mean_netcdf

    subroutine fill_yearly_netcdf(run_name, year, Cpool_yearly, Npool_yearly, levsoi) !TODO: yearly HR and climate variables (if needed?)
      !INPUT
      character (len = *),intent(in):: run_name
      integer,intent(in)            :: year
      real(r8), intent(in)          :: Cpool_yearly(levsoi,pool_types)  ! For storing C pool sizes [gC/m3]
      real(r8),intent(in)           :: Npool_yearly(levsoi,pool_types_N)  
      
      !OUTPUT
      !LOCAL
      integer :: levsoi
      integer :: i,j,varid,ncid
      
    !  real(r8)                                   :: HR_sum
    !  real(r8) ,dimension(levsoi)                :: HR_flux!(levsoi) !HR_mass_accumulated
    !  real(r8),dimension(levsoi)         :: TSOIL, MOIST

      call check(nf90_open(output_path//trim(run_name)//"_yearly_mean.nc", nf90_write, ncid))
      call check(nf90_inq_varid(ncid, "year_since_start", varid))
      call check(nf90_put_var(ncid, varid, year , start = (/ year /)))

      ! call check(nf90_inq_varid(ncid, "HR_sum", varid))
      ! call check(nf90_put_var(ncid, varid, HR_sum, start = (/ year /)))

      do j=1,levsoi
        ! call check(nf90_inq_varid(ncid, "Temp",varid))
        ! call check(nf90_put_var(ncid, varid, TSOIL(j), start = (/year,j/)))
        ! 
        ! call check(nf90_inq_varid(ncid, "Moisture",varid))
        ! call check(nf90_put_var(ncid, varid, MOIST(j), start = (/year,j/)))

        ! call check(nf90_inq_varid(ncid, "HR_flux", varid))
        ! call check(nf90_put_var(ncid, varid, HR_flux(j), start = (/year, j/)))

        call check(nf90_inq_varid(ncid, "N_inorganic", varid))
        call check(nf90_put_var(ncid, varid, Npool_yearly(j,11), start = (/year, j/)))

        do i = 1,pool_types
          !C:
          call check(nf90_inq_varid(ncid, trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, Cpool_yearly(j,i), start = (/ year, j /)))
          !N:
          call check(nf90_inq_varid(ncid, "N_"//trim(variables(i)), varid))
          call check(nf90_put_var(ncid, varid, Npool_yearly(j,i), start = (/ year, j /)))

        end do !pool_types
      end do ! levels
      call check(nf90_close(ncid))
    end subroutine fill_yearly_netcdf
end module writeMod
