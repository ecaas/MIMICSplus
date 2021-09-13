module readMod
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use netcdf
  use dispmodule
  use paramMod
  implicit none

    contains

      subroutine check(status)
        integer, intent ( in) :: status
        if(status /= nf90_noerr) then
          print *, trim(nf90_strerror(status))
          stop 2
        end if
      end subroutine check

      subroutine read_WATSAT(clm_history_file,WATSAT, nlevdecomp)
        integer,intent(in)            :: nlevdecomp
        character (len = *),intent(in):: clm_history_file
        real(r8), intent(out),dimension(nlevdecomp)          :: WATSAT
        integer            :: ncid, WATSATid

       call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))

       call check(nf90_inq_varid(ncid, 'WATSAT', WATSATid))
       call check(nf90_get_var(ncid, WATSATid, WATSAT,count=(/1,nlevdecomp/)))
       call check(nf90_close(ncid))
      end subroutine read_WATSAT

       subroutine read_clmdata(clm_history_file, TSOI, SOILLIQ,SOILICE,W_SCALAR,month, nlevdecomp)
      
      subroutine read_clay(clm_surface_file,mean_clay_content) !Needed bc. WATSAT is only given in the first outputfile of the simulation.
        !INPUT
        character (len = *),intent(in):: clm_surface_file
        
        !OUTPUT
        real(r8), intent(out)         :: mean_clay_content
        
        !LOCAL
        integer            :: ncid, clayid
        integer,dimension(1)            :: pct_clay
        
        
        call check(nf90_open(trim(clm_surface_file), nf90_nowrite, ncid))
        call check(nf90_inq_varid(ncid, 'PCT_CLAY', clayid))
        call check(nf90_get_var(ncid, clayid, pct_clay))
        call check(nf90_close(ncid))
        
        mean_clay_content = sum(pct_clay)/size(pct_clay)
      end subroutine read_clay
      
      subroutine read_LEACHING(clm_history_file,nlevdecomp,time_entry,N_leach)!Needed bc. WATSAT is only given in the first outputfile of the simulation.
        !INPUT
        integer,intent(in)            :: nlevdecomp
        integer,intent(in)            :: time_entry
        character (len = *),intent(in):: clm_history_file
        
        !OUTPUT
        real(r8),intent(out)          :: N_leach(:)
        
        !LOCAL
        integer            :: ncid, leachid
        
        call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))
        call check(nf90_inq_varid(ncid, 'SMIN_NO3_LEACHED_vr', leachid))

        call check(nf90_get_var(ncid, leachid, N_leach, start=(/1,1,time_entry/),count=(/1,nlevdecomp,1/)))
        N_leach=N_leach*sec_pr_hr !gN/(m3 s) to gN/(m3 h)
        call check(nf90_close(ncid))
      end subroutine read_LEACHING      

         integer,intent(in)            :: nlevdecomp
         character (len = *),intent(in):: clm_history_file
         real(r8),intent(out), dimension(nlevdecomp)          :: TSOI    !Celcius
         real(r8),intent(out), dimension(nlevdecomp)          :: SOILLIQ !m3/m3
         real(r8), intent(out),dimension(nlevdecomp)          :: SOILICE !m3/m3
         real(r8), intent(out),dimension(nlevdecomp)          :: W_SCALAR

         integer            :: ncid, TSOIid, SOILICEid, SOILLIQid,W_SCALARid, month,i
         
         TSOI= 0.0
         call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))
         
         call check(nf90_inq_varid(ncid, 'TSOI', TSOIid))
         call check(nf90_get_var(ncid, TSOIid, TSOI, start=(/1,1,month/), count=(/1,nlevdecomp,1/)))
         
         call check(nf90_inq_varid(ncid, 'SOILLIQ', SOILLIQid))
         call check(nf90_get_var(ncid, SOILLIQid, SOILLIQ, start=(/1,1,month/), count=(/1,nlevdecomp,1/)))

         call check(nf90_inq_varid(ncid, 'SOILICE', SOILICEid))
         call check(nf90_get_var(ncid, SOILICEid, SOILICE, start=(/1,1,month/), count=(/1,nlevdecomp,1/)))

         call check(nf90_inq_varid(ncid, 'W_SCALAR', W_SCALARid))
         call check(nf90_get_var(ncid, W_SCALARid, W_SCALAR, start=(/1,1,month/), count=(/1,nlevdecomp,1/)))
         !Unit conversions:
         TSOI = TSOI - 273.15 !Kelvin to Celcius 
         do i = 1, nlevdecomp ! Unit conversion
           SOILICE(i) = SOILICE(i)/(delta_z(i)*917) !kg/m2 to m3/m3 rho_ice=917kg/m3
           SOILLIQ(i) = SOILLIQ(i)/(delta_z(i)*1000) !kg/m2 to m3/m3 rho_liq=1000kg/m3
         end do
         call check(nf90_close(ncid))
       end subroutine read_clmdata

       subroutine read_clm_model_input(clm_history_file, LITFALL, NPP_NACTIVE,NDEP_TO_SMINN,time_entry)
          !INPUT
          character (len = *),intent(in):: clm_history_file
          integer,intent(in)            :: time_entry
          
          !OUTPUT 
          real(r8), intent(out)        :: LITFALL      !litterfall (leaves and fine roots)  [gC/m^2/hour] (converted from [gC/m^2/s])  
          real(r8), intent(out)        :: NPP_NACTIVE  !Mycorrhizal N uptake used C        [gC/m^2/hour] (converted from [gC/m^2/s])  
          real(r8), intent(out)        :: NDEP_TO_SMINN   !atmospheric N deposition to soil mineral N [gN/m^2/hour] (converted from [gN/m^2/s]) 

          !Local
          integer            :: ncid, varid

          call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))

          call check(nf90_inq_varid(ncid, 'LITFALL', varid))
          call check(nf90_get_var(ncid, varid, LITFALL,start=(/1, time_entry/)))
          LITFALL = LITFALL*sec_pr_hr  ![gC/m^2/h]

          call check(nf90_inq_varid(ncid, 'NPP_NACTIVE', varid))
          call check(nf90_get_var(ncid, varid, NPP_NACTIVE,start=(/1, time_entry/)))
          NPP_NACTIVE = NPP_NACTIVE*sec_pr_hr ![gC/m^2/h]

          call check(nf90_inq_varid(ncid, 'NDEP_TO_SMINN', varid))
          call check(nf90_get_var(ncid, varid, NDEP_TO_SMINN,start=(/1, time_entry/)))
          NDEP_TO_SMINN = NDEP_TO_SMINN*sec_pr_hr ![gN/m^2/h]

          call check(nf90_close(ncid))
       end subroutine read_clm_model_input

end module readMod
