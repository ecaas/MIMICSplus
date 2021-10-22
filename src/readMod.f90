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
      
      subroutine read_time(clm_history_file,steps)
        !INPUT
        character (len = *),intent(in):: clm_history_file
        
        !OUTPUT
        integer, intent(out):: steps
        
        !LOCAL
        integer            :: ncid, timeid
        
        call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))
        call check(nf90_inquire(ncid, unlimiteddimid = timeid))
        call check(nf90_inquire_dimension(ncid, timeid, len = steps))
        call check(nf90_close(ncid))
      end subroutine read_time

      subroutine read_WATSAT(clm_history_file,WATSAT, nlevdecomp) !Needed bc. WATSAT is only given in the first outputfile of the simulation.
        !INPUT
        integer,intent(in)            :: nlevdecomp
        character (len = *),intent(in):: clm_history_file
        
        !OUTPUT
        real(r8), intent(out),dimension(nlevdecomp)          :: WATSAT
        
        !LOCAL
        integer            :: ncid, WATSATid
        
        call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))
        call check(nf90_inq_varid(ncid, 'WATSAT', WATSATid))
        call check(nf90_get_var(ncid, WATSATid, WATSAT,count=(/1,nlevdecomp/)))
        call check(nf90_close(ncid))
      end subroutine read_WATSAT
      
      subroutine read_clay(clm_surface_file,mean_clay_content, nlevdecomp) !Needed bc. WATSAT is only given in the first outputfile of the simulation.
        !INPUT
        character (len = *),intent(in):: clm_surface_file
        integer,intent(in)                       :: nlevdecomp
        
        !OUTPUT
        real(r8), intent(out)         :: mean_clay_content
        
        !LOCAL
        integer            :: ncid, clayid
        real(r8),dimension(nlevdecomp)           :: pct_clay

        call check(nf90_open(trim(clm_surface_file), nf90_nowrite, ncid))
        call check(nf90_inq_varid(ncid, 'PCT_CLAY', clayid))
        call check(nf90_get_var(ncid, clayid, pct_clay, count=(/1,1,nlevdecomp/)))
        call check(nf90_close(ncid))

        mean_clay_content = (sum(pct_clay)/size(pct_clay))/100.0 !output as fraction
        

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
      subroutine read_clm_model_input(ncid, nlevdecomp,time_entry, &
                                      LITFALL, LEAFN_TO_LITTER, NPP_NACTIVE,NDEP_TO_SMINN, &
                                      LEAFC_TO_LITTER,mcdate,TSOI,SOILLIQ,SOILICE, &
                                      W_SCALAR,C_CWD,N_CWD,N_leach)
        !INPUT
        integer,intent(in)            :: ncid
        integer,intent(in)            :: nlevdecomp          
        integer,intent(in)            :: time_entry
          
        !OUTPUT 
        real(r8),intent(out)                                :: LITFALL      !litterfall (leaves and fine roots)  [gC/m^2/hour] (converted from [gC/m^2/s])  
        real(r8),intent(out)                                :: LEAFN_TO_LITTER !litterfall N from leaves  [gN/m^2/hour] (converted from [gN/m^2/s])             
        real(r8),intent(out)                                :: NPP_NACTIVE  !Mycorrhizal N uptake used C        [gC/m^2/hour] (converted from [gC/m^2/s])  
        real(r8),intent(out)                                :: NDEP_TO_SMINN   !atmospheric N deposition to soil mineral N [gN/m^2/hour] (converted from [gN/m^2/s]) 
        real(r8),intent(out)                                :: LEAFC_TO_LITTER
        integer ,intent(out)                                :: mcdate  !"Current date" , end of averaging interval
        real(r8),intent(out), dimension(nlevdecomp)         :: TSOI    !Celcius
        real(r8),intent(out), dimension(nlevdecomp)         :: SOILLIQ !m3/m3
        real(r8),intent(out),dimension(nlevdecomp)          :: SOILICE !m3/m3
        real(r8),intent(out),dimension(nlevdecomp)          :: W_SCALAR          
        real(r8),intent(out)                                :: C_CWD(:) !gC/(m3 h)
        real(r8),intent(out)                                :: N_CWD(:) !gN/(m3 h)
        real(r8),intent(out)                                :: N_leach(:)


          !Local
        integer                   :: varid
        integer                   :: i
        real(r8)                    :: C_CWD2(nlevdecomp)
        real(r8)                    :: C_CWD3(nlevdecomp)        
        real(r8)                    :: N_CWD2(nlevdecomp)
        real(r8)                    :: N_CWD3(nlevdecomp)                  
          
        call check(nf90_inq_varid(ncid, 'CWDC_TO_LITR2C_vr', varid))
        call check(nf90_get_var(ncid, varid, C_CWD2, start=(/1,1,time_entry/),count=(/1,nlevdecomp,1/)))
        C_CWD2=C_CWD2*sec_pr_hr !gC/(m3 s) to gC/(m3 h)
        call check(nf90_inq_varid(ncid, 'CWDC_TO_LITR3C_vr', varid))
        call check(nf90_get_var(ncid, varid, C_CWD3, start=(/1,1,time_entry/),count=(/1,nlevdecomp,1/)))
        C_CWD3=C_CWD3*sec_pr_hr !gC/(m3 s) to gC/(m3 h)
        C_CWD=C_CWD2+C_CWD3
        
        call check(nf90_inq_varid(ncid, 'CWDN_TO_LITR2N_vr', varid))
        call check(nf90_get_var(ncid, varid, N_CWD2, start=(/1,1,time_entry/),count=(/1,nlevdecomp,1/)))
        
        N_CWD2=N_CWD2*sec_pr_hr !gN/(m3 s) to gN/(m3 h)
        call check(nf90_inq_varid(ncid, 'CWDN_TO_LITR3N_vr', varid))
        call check(nf90_get_var(ncid, varid, N_CWD3, start=(/1,1,time_entry/),count=(/1,nlevdecomp,1/)))
        N_CWD3=N_CWD3*sec_pr_hr !gN/(m3 s) to gN/(m3 h)        
        N_CWD=N_CWD2+N_CWD3
          call check(nf90_inq_varid(ncid, 'LITFALL', varid))
          call check(nf90_get_var(ncid, varid, LITFALL,start=(/1, time_entry/)))
          LITFALL = LITFALL*sec_pr_hr  ![gC/m^2/h]
          
          call check(nf90_inq_varid(ncid, 'LEAFN_TO_LITTER', varid))
          call check(nf90_get_var(ncid, varid, LEAFN_TO_LITTER,start=(/1, time_entry/)))
          LEAFN_TO_LITTER = LEAFN_TO_LITTER*sec_pr_hr  ![gC/m^2/h]

          call check(nf90_inq_varid(ncid, 'NPP_NACTIVE', varid))
          call check(nf90_get_var(ncid, varid, NPP_NACTIVE,start=(/1, time_entry/)))
          NPP_NACTIVE = NPP_NACTIVE*sec_pr_hr ![gC/m^2/h]

          call check(nf90_inq_varid(ncid, 'NDEP_TO_SMINN', varid))
          call check(nf90_get_var(ncid, varid, NDEP_TO_SMINN,start=(/1, time_entry/)))
          NDEP_TO_SMINN = NDEP_TO_SMINN*sec_pr_hr ![gN/m^2/h]

          call check(nf90_inq_varid(ncid, 'LEAFC_TO_LITTER', varid))
          call check(nf90_get_var(ncid, varid, LEAFC_TO_LITTER,start=(/1, time_entry/)))
          
          LEAFC_TO_LITTER = LEAFC_TO_LITTER*sec_pr_hr  ![gC/m^2/h]
          

         call check(nf90_inq_varid(ncid, 'mcdate', varid))
         call check(nf90_get_var(ncid, varid, mcdate, start=(/time_entry/)))
         
         call check(nf90_inq_varid(ncid, 'TSOI', varid))
         call check(nf90_get_var(ncid, varid, TSOI, start=(/1,1,time_entry/), count=(/1,nlevdecomp,1/)))

         call check(nf90_inq_varid(ncid, 'SOILLIQ', varid))
         call check(nf90_get_var(ncid, varid, SOILLIQ, start=(/1,1,time_entry/), count=(/1,nlevdecomp,1/)))

         call check(nf90_inq_varid(ncid, 'SOILICE', varid))
         call check(nf90_get_var(ncid, varid, SOILICE, start=(/1,1,time_entry/), count=(/1,nlevdecomp,1/)))

         call check(nf90_inq_varid(ncid, 'W_SCALAR', varid))
         call check(nf90_get_var(ncid, varid, W_SCALAR, start=(/1,1,time_entry/), count=(/1,nlevdecomp,1/)))
         !Unit conversions:
         TSOI = TSOI - 273.15 !Kelvin to Celcius 
         do i = 1, nlevdecomp ! Unit conversion
           SOILICE(i) = SOILICE(i)/(delta_z(i)*917) !kg/m2 to m3/m3 rho_ice=917kg/m3
           SOILLIQ(i) = SOILLIQ(i)/(delta_z(i)*1000) !kg/m2 to m3/m3 rho_liq=1000kg/m3
         end do          
        call check(nf90_inq_varid(ncid, 'SMIN_NO3_LEACHED_vr', varid))

        call check(nf90_get_var(ncid, varid, N_leach, start=(/1,1,time_entry/),count=(/1,nlevdecomp,1/)))
        N_leach=N_leach*sec_pr_hr !gN/(m3 s) to gN/(m3 h)         
         
       end subroutine read_clm_model_input
       

end module readMod
