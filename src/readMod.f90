module readMod
  use shr_kind_mod, only : r8 => shr_kind_r8
  use netcdf
  use dispmodule , only: disp
  use paramMod,    only: sec_pr_hr,delta_z
  use writeMod,    only: check
  use initMod,     only: nlevels
  implicit none
    private
    public :: read_time,read_maxC,read_WATSAT_and_profiles,read_clay,read_PFTs,read_clm_model_input

  contains
      
    subroutine read_time(ncid,steps)
      !INPUT
      integer,intent(in):: ncid
      
      !OUTPUT
      integer, intent(out):: steps
      
      !LOCAL
      integer            :: timeid

      call check(nf90_inquire(ncid, unlimiteddimid = timeid))
      call check(nf90_inquire_dimension(ncid, timeid, len = steps))
    end subroutine read_time
      
    function read_maxC(ncid,time_steps) result(max_Cpay)
      !INput
      real(r8)              :: max_Cpay
      integer, intent(in)   :: ncid
      integer, intent(in)   :: time_steps
      
      !Local
      integer    :: varid
      real(r8), dimension(time_steps) :: NPP_tot
      real(r8), dimension(time_steps) :: NPP_nonmyc
      real(r8), dimension(time_steps) :: NPP_myc
      
      call check(nf90_inq_varid(ncid, 'NPP_NACTIVE', varid))
      call check(nf90_get_var(ncid, varid, NPP_tot(:),start = (/1,1/),count = (/1,time_steps/)))
      call check(nf90_inq_varid(ncid, 'NPP_NNONMYC', varid))
      call check(nf90_get_var(ncid, varid, NPP_nonmyc(:),start = (/1,1/),count = (/1,time_steps/)))
      
      NPP_myc = (NPP_tot-NPP_nonmyc)*sec_pr_hr

      max_Cpay = maxval(NPP_myc)        
    end function read_maxC

    subroutine read_WATSAT_and_profiles(clm_history_file,WATSAT,NDEP_PROF,FROOT_PROF,LEAF_PROF) !Needed bc. WATSAT is only given in the first outputfile of the simulation.
      !INPUT
      character (len = *),intent(in):: clm_history_file       
      !OUTPUT
      real(r8), intent(out),dimension(nlevels)          :: WATSAT
      real(r8), intent(out),dimension(nlevels)          :: NDEP_PROF
      real(r8), intent(out),dimension(nlevels)          :: FROOT_PROF
      real(r8), intent(out),dimension(nlevels)          :: LEAF_PROF
   
      !LOCAL
      integer            :: ncid, varid
      
      call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))
      call check(nf90_inq_varid(ncid, 'WATSAT', varid))
      call check(nf90_get_var(ncid, varid, WATSAT,count=(/1,nlevels/)))
      
      call check(nf90_inq_varid(ncid, 'LEAF_PROF', varid))
      call check(nf90_get_var(ncid, varid, LEAF_PROF,count=(/1,nlevels/)))       
      call check(nf90_inq_varid(ncid, 'FROOT_PROF', varid))
      call check(nf90_get_var(ncid, varid, FROOT_PROF,count=(/1,nlevels/))) 
      call check(nf90_inq_varid(ncid, 'NDEP_PROF', varid))
      call check(nf90_get_var(ncid, varid, NDEP_PROF,count=(/1,nlevels/))) 
      
      call check(nf90_close(ncid))
    end subroutine read_WATSAT_and_profiles
    
    subroutine read_clay(clm_surface_file,mean_clay_content)
      !INPUT
      character (len = *),intent(in):: clm_surface_file
      
      !OUTPUT
      real(r8), intent(out)         :: mean_clay_content
      
      !LOCAL
      integer            :: ncid, clayid
      real(r8),dimension(nlevels)           :: pct_clay

      call check(nf90_open(trim(clm_surface_file), nf90_nowrite, ncid))
      call check(nf90_inq_varid(ncid, 'PCT_CLAY', clayid))
      call check(nf90_get_var(ncid, clayid, pct_clay, count=(/1,1,nlevels/)))
      call check(nf90_close(ncid))

      mean_clay_content = (sum(pct_clay)/size(pct_clay))/100.0 !output as fraction
    end subroutine read_clay
    
    subroutine read_PFTs(clm_surface_file,PFT_dist)
      !INPUT
      character (len = *),intent(in):: clm_surface_file
      
      !OUTPUT
      real(r8), dimension(15),intent(out)         :: PFT_dist
      
      !LOCAL
      integer            :: ncid, pftid

      call check(nf90_open(trim(clm_surface_file), nf90_nowrite, ncid))
      call check(nf90_inq_varid(ncid, 'PCT_NAT_PFT', pftid))
      call check(nf90_get_var(ncid, pftid, PFT_dist, count=(/1,1,15/)))
      call check(nf90_close(ncid))
    end subroutine read_PFTs
    
    ! function read_nlayers(clm_history_file) result(nlevels)
    !   !INPUT
    !   character (len = *),intent(in):: clm_history_file
    ! 
    !   !OUTPUT
    !   integer       :: nlevels
    ! 
    !   !LOCAL
    !   integer            :: ncid, nid
    ! 
    !   call check(nf90_open(trim(clm_history_file), nf90_nowrite, ncid))
    !   call check(nf90_inq_varid(ncid, 'nbedrock', nid))
    !   call check(nf90_get_var(ncid, nid, nlevels))
    !   call check(nf90_close(ncid))
    ! end function read_nlayers      
    

    
                
    subroutine read_clm_model_input(ncid,time_entry, &
                                    LEAFN_TO_LITTER,FROOTN_TO_LITTER, NPP_MYC,NDEP_TO_SMINN, &
                                    LEAFC_TO_LITTER,FROOTC_TO_LITTER,mcdate,TSOI,SOILLIQ,SOILICE, &
                                    W_SCALAR,T_SCALAR,QDRAI,h2o_liq_tot,C_CWD,N_CWD)
      !INPUT
      integer,intent(in)            :: ncid
      integer,intent(in)            :: time_entry
        
      !OUTPUT 
      real(r8),intent(out)                                :: LEAFN_TO_LITTER !litterfall N from leaves  [gN/m^2/hour] (converted from [gN/m^2/s])      
      real(r8),intent(out)                                :: FROOTN_TO_LITTER !litterfall N from leaves  [gN/m^2/hour] (converted from [gN/m^2/s])
      real(r8),intent(out)                                :: NPP_MYC                    
      real(r8),intent(out)                                :: NDEP_TO_SMINN   !atmospheric N deposition to soil mineral N [gN/m^2/hour] (converted from [gN/m^2/s]) 
      real(r8),intent(out)                                :: LEAFC_TO_LITTER
      real(r8),intent(out)                                :: FROOTC_TO_LITTER
      integer ,intent(out)                                :: mcdate  !"Current date" , end of averaging interval
      real(r8),intent(out),dimension(nlevels)          :: TSOI    !Celcius
      real(r8),intent(out),dimension(nlevels)          :: SOILLIQ !m3/m3 (converted from kg/m2)
      real(r8),intent(out),dimension(nlevels)          :: SOILICE !m3/m3 (converted from kg/m2)
      real(r8),intent(out),dimension(nlevels)          :: W_SCALAR    
      real(r8),intent(out),dimension(nlevels)          :: T_SCALAR      
      real(r8),intent(out)                                :: QDRAI   !kgH2O/(m2 h) converted from kgH2O/(m2 s)
      real(r8),intent(out)                                :: h2o_liq_tot
      real(r8),intent(out)                                :: C_CWD(:) !gC/(m3 h)
      real(r8),intent(out)                                :: N_CWD(:) !gN/(m3 h)

        !Local
      integer                   :: varid
      integer                   :: i
      real(r8)                    :: C_CWD2(nlevels)
      real(r8)                    :: C_CWD3(nlevels)        
      real(r8)                    :: N_CWD2(nlevels)
      real(r8)                    :: N_CWD3(nlevels)       
      real(r8)                    :: NPP_NACTIVE  !Mycorrhizal N uptake used C        [gC/m^2/hour] (converted from [gC/m^2/s]) 
      real(r8)                    :: NPP_NNONMYC  !NONMycorrhizal N uptake used C        [gC/m^2/hour] (converted from [gC/m^2/s]) 
                 
        
      ! C in Coarse Woody Debris
      call check(nf90_inq_varid(ncid, 'CWDC_TO_LITR2C_vr', varid))
      call check(nf90_get_var(ncid, varid, C_CWD2, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))
      C_CWD2=C_CWD2*sec_pr_hr !gC/(m3 s) to gC/(m3 h)
      call check(nf90_inq_varid(ncid, 'CWDC_TO_LITR3C_vr', varid))
      call check(nf90_get_var(ncid, varid, C_CWD3, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))
      C_CWD3=C_CWD3*sec_pr_hr !gC/(m3 s) to gC/(m3 h)
      
      C_CWD=C_CWD2+C_CWD3
      
      ! N in Coarse Woody Debris
      call check(nf90_inq_varid(ncid, 'CWDN_TO_LITR2N_vr', varid))
      call check(nf90_get_var(ncid, varid, N_CWD2, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))
      
      N_CWD2=N_CWD2*sec_pr_hr !gN/(m3 s) to gN/(m3 h)
      call check(nf90_inq_varid(ncid, 'CWDN_TO_LITR3N_vr', varid))
      call check(nf90_get_var(ncid, varid, N_CWD3, start=(/1,1,time_entry/),count=(/1,nlevels,1/)))
      N_CWD3=N_CWD3*sec_pr_hr !gN/(m3 s) to gN/(m3 h)        
      N_CWD=N_CWD2+N_CWD3
      
      !C and N litter from leafs and fine roots:
      call check(nf90_inq_varid(ncid, 'LEAFN_TO_LITTER', varid))
      call check(nf90_get_var(ncid, varid, LEAFN_TO_LITTER,start=(/1, time_entry/)))
      LEAFN_TO_LITTER = LEAFN_TO_LITTER*sec_pr_hr  ![gC/m^2/h]
      
      call check(nf90_inq_varid(ncid, 'FROOTN_TO_LITTER', varid))
      call check(nf90_get_var(ncid, varid, FROOTN_TO_LITTER,start=(/1, time_entry/)))
      FROOTN_TO_LITTER = FROOTN_TO_LITTER*sec_pr_hr  ![gC/m^2/h]    

      call check(nf90_inq_varid(ncid, 'LEAFC_TO_LITTER', varid))
      call check(nf90_get_var(ncid, varid, LEAFC_TO_LITTER,start=(/1, time_entry/)))          
      LEAFC_TO_LITTER = LEAFC_TO_LITTER*sec_pr_hr  ![gC/m^2/h]

      call check(nf90_inq_varid(ncid, 'FROOTC_TO_LITTER', varid))
      call check(nf90_get_var(ncid, varid, FROOTC_TO_LITTER,start=(/1, time_entry/)))          
      FROOTC_TO_LITTER = FROOTC_TO_LITTER*sec_pr_hr  ![gC/m^2/h]

      !C payment for total payed for N (NPP_NACTIVE) and non-mycorrhizal (NPP_NNONMYC)
      call check(nf90_inq_varid(ncid, 'NPP_NACTIVE', varid))
      call check(nf90_get_var(ncid, varid, NPP_NACTIVE,start=(/1, time_entry/)))
      NPP_NACTIVE = NPP_NACTIVE*sec_pr_hr ![gC/m^2/h]

      call check(nf90_inq_varid(ncid, 'NPP_NNONMYC', varid))
      call check(nf90_get_var(ncid, varid, NPP_NNONMYC,start=(/1, time_entry/)))
      NPP_NNONMYC = NPP_NNONMYC*sec_pr_hr ![gC/m^2/h]
      
      !Part of total NPP_NACTIVE going to mycorrhizal associations
      NPP_MYC = NPP_NACTIVE - NPP_NNONMYC
      
      !N deposition:
      call check(nf90_inq_varid(ncid, 'NDEP_TO_SMINN', varid))
      call check(nf90_get_var(ncid, varid, NDEP_TO_SMINN,start=(/1, time_entry/)))
      NDEP_TO_SMINN = NDEP_TO_SMINN*sec_pr_hr ![gN/m^2/h]

      !Environmental variables:      
      call check(nf90_inq_varid(ncid, 'TSOI', varid))
      call check(nf90_get_var(ncid, varid, TSOI, start=(/1,1,time_entry/), count=(/1,nlevels,1/)))
      TSOI = TSOI - 273.15 !UNIT CONVERSION: Kelvin to Celcius 
      
      call check(nf90_inq_varid(ncid, 'SOILLIQ', varid))
      call check(nf90_get_var(ncid, varid, SOILLIQ, start=(/1,1,time_entry/), count=(/1,nlevels,1/)))

      call check(nf90_inq_varid(ncid, 'SOILICE', varid))
      call check(nf90_get_var(ncid, varid, SOILICE, start=(/1,1,time_entry/), count=(/1,nlevels,1/)))

      call check(nf90_inq_varid(ncid, 'W_SCALAR', varid))
      call check(nf90_get_var(ncid, varid, W_SCALAR, start=(/1,1,time_entry/), count=(/1,nlevels,1/)))
      
      call check(nf90_inq_varid(ncid, 'T_SCALAR', varid))
      call check(nf90_get_var(ncid, varid, T_SCALAR, start=(/1,1,time_entry/), count=(/1,nlevels,1/)))

      call check(nf90_inq_varid(ncid, 'QDRAI', varid)) !mmH2O/s = kg H2O/(m2 s)
      call check(nf90_get_var(ncid, varid, QDRAI,start=(/1, time_entry/)))    
      QDRAI = QDRAI*sec_pr_hr !kgH2O/(m2 h)

       h2o_liq_tot=0.0
       do i = 1, nlevels ! Total liquid water in column, used for calculating leaching
         h2o_liq_tot=h2o_liq_tot+SOILLIQ(i) !kg/m2
       end do          
       do i = 1, nlevels ! Unit conversion
         SOILICE(i) = SOILICE(i)/(delta_z(i)*917) !kg/m2 to m3/m3 rho_ice=917kg/m3
         SOILLIQ(i) = SOILLIQ(i)/(delta_z(i)*1000) !kg/m2 to m3/m3 rho_liq=1000kg/m3
       end do          
      
      !CLM dates
      call check(nf90_inq_varid(ncid, 'mcdate', varid))
      call check(nf90_get_var(ncid, varid, mcdate, start=(/time_entry/)))
    end subroutine read_clm_model_input

end module readMod
