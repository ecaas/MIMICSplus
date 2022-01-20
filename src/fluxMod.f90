module fluxMod
  use paramMod
  use dispmodule !External module to pretty print matrices (mainly for testing purposes)
  implicit none

contains 
  function nitrification(nh4,t_scalar,w_scalar,soil_temp) result(f_nit)
    real(r8) :: f_nit
    !IN:
    real(r8),intent(in):: nh4
    real(r8),intent(in) :: t_scalar
    real(r8),intent(in) :: w_scalar
    real(r8),intent(in) :: soil_temp
    
    !local
    real(r8) :: anaerobic_frac
    real(r8),parameter :: pH = 6.5_r8
    real(r8),parameter :: rpi = 3.14159265358979323846_R8
    real(R8),parameter :: SHR_CONST_TKFRZ   = 0.0_R8! freezing T of fresh water          ~ degC
    real(r8),parameter :: k_nitr_max = 0.1_r8/24._r8 !from paramfile ctsm51_params.c210528.nc = 0.1/day, converted to /hour
    real(r8) :: k_nitr_t_vr,k_nitr_ph_vr,k_nitr_h2o_vr,k_nitr
    ! follows CENTURY nitrification scheme (Parton et al., (2001, 1996))

    ! assume nitrification temp function equal to the HR scalar
    k_nitr_t_vr = min(t_scalar, 1._r8)

    ! ph function from Parton et al., (2001, 1996)
    k_nitr_ph_vr = 0.56_r8 + atan(rpi * 0.45_r8 * (-5._r8+pH)/rpi)

    ! moisture function-- assume the same moisture function as limits heterotrophic respiration
    ! Parton et al. base their nitrification- soil moisture rate constants based on heterotrophic rates-- can we do the same?
    k_nitr_h2o_vr = w_scalar

    ! nitrification constant is a set scalar * temp, moisture, and ph scalars
    ! note that k_nitr_max_perday is converted from 1/day to 1/s
    k_nitr = k_nitr_max * k_nitr_t_vr * k_nitr_h2o_vr * k_nitr_ph_vr

    ! first-order decay of ammonium pool with scalar defined above
    f_nit = max(nh4 * k_nitr, 0._r8) !g/m3 h
    anaerobic_frac=0._r8
    ! limit to oxic fraction of soils
    f_nit  = f_nit* (1._r8 - anaerobic_frac)

    !limit to non-frozen soil layers
    if ( soil_temp <= SHR_CONST_TKFRZ ) then
       f_nit = 0._r8
     end if 
  end function nitrification
  
  function set_N_dep(CLMdep,const_dep) result(Dep)
    real(r8)           :: Dep
    real(r8), optional :: CLMdep
    real(r8), optional :: const_dep
    Deposition = -999
    if (present(CLMdep) .and. .not. present(const_dep)) then
      Dep = CLMdep
    elseif (present(const_dep) .and. .not. present(CLMdep)) then
      Dep = const_dep/(hr_pr_yr*soil_depth)
    elseif (.not. present(CLMdep) .and. .not. present(const_dep)) then
      print*, "N dep not set correctly, stopping"
      Dep = -999
      stop
    end if
  end function set_N_dep
  
  function calc_Leaching(drain,h2o_tot, N_inorganic) result(Leach)
    real(r8)           :: Leach       !gN/m3 h
    real(r8)           :: drain       !mmH20/h = kgH20/m2 h
    real(r8)           :: h2o_tot     !kgH20/m2
    real(r8)           :: N_inorganic !gN/m3
    Leach = N_inorganic*drain/h2o_tot
  end function calc_Leaching


  function MMK_flux(C_SAP,C_SUBSTRATE,MMK_nr) result(flux)
    !Compute C flux from substrate pool to saprotroph pool by using Michaelis Menten Kinetics.
    !NOTE: On the way, a fraction 1-MGE is lost as respiration. This is handeled in the "decomp" subroutine.
    real(r8):: flux ![gC/(m3 hr)]
    real(r8), intent(in) :: C_SAP
    real(r8), intent(in) :: C_SUBSTRATE
    integer, intent (in) :: MMK_nr
    !TODO: this works, but should not depend on Vmax & Km from mycmim mod
    flux = C_SAP*Vmax(MMK_nr)*C_SUBSTRATE/(Km(MMK_nr)+C_SUBSTRATE)
  end function MMK_flux
  
    function reverse_MMK_flux(C_SAP,C_SUBSTRATE,MMK_nr) result(flux)
      !Compute C flux from substrate pool to saprotroph pool by using Michaelis Menten Kinetics.
      !NOTE: On the way, a fraction 1-MGE is lost as respiration. This is handeled in the "decomp" subroutine.
      real(r8):: flux ![gC/(m3 hr)]
      real(r8), intent(in) :: C_SAP
      real(r8), intent(in) :: C_SUBSTRATE
      integer, intent (in) :: MMK_nr
      !TODO: this works, but should not depend on Vmax & Km from mycmim mod
      flux = C_SUBSTRATE*Vmax(MMK_nr)*C_SAP/(Km(MMK_nr)+C_SAP)
    end function reverse_MMK_flux

  function Km_function(temperature) result(K_m)
    real(r8),dimension(MM_eqs)             :: K_m
    real(r8), intent(in)                   :: temperature
    K_m      = exp(Kslope*temperature + Kint)*a_k*Kmod               ![mgC/cm3]*1e4=[gC/m3]
    
  end function Km_function

  function Vmax_function(temperature, moisture) result(V_max)
    real(r8),dimension(MM_eqs)             :: V_max
    real(r8), intent(in)                   :: temperature
    real(r8), intent(in)                   :: moisture
    V_max    = exp(Vslope*temperature + Vint)*a_v*Vmod*moisture   ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?
  end function Vmax_function

  subroutine calculate_fluxes(depth,nlevdecomp,C_pool_matrix,N_pool_matrix,dt) !This subroutine calculates the fluxes in and out of the SOM pools.
    integer         :: depth !depth level
    integer         :: nlevdecomp
    real(r8)        :: dt ! timestep
    real(r8),target :: C_pool_matrix(nlevdecomp, pool_types)
    real(r8),target :: N_pool_matrix(nlevdecomp, pool_types_N)
    !LOCAL:
    real(r8)  :: N_for_sap
    real(r8)  :: N_IN
    

    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: C_LITm, C_LITs, C_SOMp,C_SOMa,C_SOMc,C_EcM,C_ErM,C_AM, &
    C_SAPb, C_SAPf, N_LITm, N_LITs, N_SOMp,N_SOMa,N_SOMc,N_EcM,N_ErM,N_AM, N_SAPb, N_SAPf, N_NH4,N_NO3
    C_LITm => C_pool_matrix(depth, 1)
    C_LITs => C_pool_matrix(depth, 2)
    C_SAPb => C_pool_matrix(depth, 3)
    C_SAPf => C_pool_matrix(depth, 4)
    C_ErM =>  C_pool_matrix(depth, 6)
    C_EcM =>  C_pool_matrix(depth, 5)
    C_AM =>   C_pool_matrix(depth, 7)
    C_SOMp => C_pool_matrix(depth, 8)
    C_SOMa => C_pool_matrix(depth, 9)
    C_SOMc => C_pool_matrix(depth, 10)

    N_LITm => N_pool_matrix(depth, 1)
    N_LITs => N_pool_matrix(depth, 2)
    N_SAPb => N_pool_matrix(depth, 3)
    N_SAPf => N_pool_matrix(depth, 4)
    N_EcM =>  N_pool_matrix(depth, 5)
    N_ErM =>  N_pool_matrix(depth, 6)
    N_AM =>   N_pool_matrix(depth, 7)
    N_SOMp => N_pool_matrix(depth, 8)
    N_SOMa => N_pool_matrix(depth, 9)
    N_SOMc => N_pool_matrix(depth, 10)
    N_NH4 => N_pool_matrix(depth, 11)
    N_NO3 => N_pool_matrix(depth, 12)
    

    !------------------CARBON FLUXES----------------------------:
    !Decomposition of LIT by SAP:
    !On the way, a fraction 1-MGE is lost as respiration. This is handeled in the "decomp" subroutine.
    C_LITmSAPb=MMK_flux(C_SAPb,C_LITm,1)
    C_LITsSAPb=MMK_flux(C_SAPb,C_LITs,2)
    C_LITmSAPf=MMK_flux(C_SAPf,C_LITm,4)
    C_LITsSAPf=MMK_flux(C_SAPf,C_LITs,5)
    !Decomposition of SOMa by SAP. Based on the equations from SOMa to microbial pools in mimics.
    !On the way, a fraction 1-MGE is lost as respiration. This is handeled in the "decomp" subroutine.
    C_SOMaSAPb=MMK_flux(C_SAPb,C_SOMa,3)
    C_SOMaSAPf=MMK_flux(C_SAPf,C_SOMa,6)


    !Dead mycorrhizal biomass enters the SOM pools:  gC/m3h
    C_EcMSOMp=C_EcM*k_mycsom(1)*fEcMSOM(1)!somp
    C_EcMSOMa=C_EcM*k_mycsom(1)*fEcMSOM(2)!soma
    C_EcMSOMc=C_EcM*k_mycsom(1)*fEcMSOM(3)!somc

    C_ErMSOMp=0.0!C_ErM*k_mycsom(2)*fErMSOM(1)
    C_ErMSOMa=0.0!C_ErM*k_mycsom(2)*fErMSOM(2)
    C_ErMSOMc=0.0!C_ErM*k_mycsom(2)*fErMSOM(3)

    C_AMSOMp=C_AM*k_mycsom(3)*fAMSOM(1)
    C_AMSOMa=C_AM*k_mycsom(3)*fAMSOM(2)
    C_AMSOMc=C_AM*k_mycsom(3)*fAMSOM(3)

    !Turnover from SAP to SOM. Based on the turnover equations used in mimics for flux from microbial pools to SOM pools (correspond to eq A4,A8 in Wieder 2015)
    C_SAPbSOMp=C_SAPb*tau(1)*fPHYS(1)   !gC/m3h
    C_SAPbSOMa=C_SAPb*tau(1)*fAVAIL(1)
    C_SAPbSOMc=C_SAPb*tau(1)*fCHEM(1)

    C_SAPfSOMp=C_SAPf*tau(2)*fPHYS(2)
    C_SAPfSOMa=C_SAPf*tau(2)*fAVAIL(2)
    C_SAPfSOMc=C_SAPf*tau(2)*fCHEM(2)

    !Desorbtion controls transport from physically protected to available SOM
    C_SOMpSOMa=C_SOMp*desorb

    !Oxidation from SOMc to SOMa
    !From equations for decomposing structural litter in mimics,eq. A10
    !KO modifies Km which is used in the litter->SAP equations.
    C_SOMcSOMa    = ( C_SAPb * Vmax(2) * C_SOMc / (KO(1)*Km(2) + C_SOMc)) + &
                   (C_SAPf * Vmax(5) * C_SOMc / (KO(2)*Km(5) + C_SOMc))

    !Baskaran et al: Rates of decomposition of available SOM mediated by mycorrhizal enzymes:
    C_EcMdecompSOMp = K_MO*soil_depth*C_EcM*C_SOMp*(C_PlantEcM/(max_mining*froot_prof(depth)))   ![gC/m3h]
    C_EcMdecompSOMc = K_MO*soil_depth*C_EcM*C_SOMc*(C_PlantEcM/(max_mining*froot_prof(depth)))   ![gC/m3h]

    !-----------------------------------NITROGEN FLUXES----------------------------:
    N_IN = N_NH4+ N_NO3
    !Nitrogen aquired bymycorrhiza via oxidation of protected SOM pools.  gN/m3h
    !TODO: Should we also include SOMa?
    N_SOMpEcM = C_EcMdecompSOMp*N_SOMp/C_SOMp
    N_SOMcEcM = C_EcMdecompSOMc*N_SOMc/C_SOMc
    
    !Inorganic N taken up directly by plant roots
    N_InPlant =  5E-7*N_IN
    
    N_INEcM = V_max_myc*N_IN*(C_EcM/(C_EcM + Km_myc/soil_depth))*(C_PlantEcM/(max_mining*froot_prof(depth)))  !NOTE: MMK parameters should maybe be specific to mycorrhizal type?
    N_INErM = 0.0!V_max_myc*N_IN*(C_ErM/(C_ErM + Km_myc/delta_z(depth)))   !Unsure about units
    N_INAM = V_max_myc*N_IN*(C_AM/(C_AM + Km_myc/soil_depth))*(C_PlantEcM/(max_mining*froot_prof(depth)))

    !Decomposition of LIT and SOMa by SAP
    N_LITmSAPb = C_LITmSAPb*N_LITm/C_LITm
    N_LITsSAPb = C_LITsSAPb*N_LITs/C_LITs

    N_LITmSAPf = C_LITmSAPf*N_LITm/C_LITm
    N_LITsSAPf = C_LITsSAPf*N_LITs/C_LITs

    N_SOMaSAPb = C_SOMaSAPb*N_SOMa/C_SOMa
    N_SOMaSAPf = C_SOMaSAPf*N_SOMa/C_SOMa

    !Dead mycorrhizal biomass enters SOM pools
    N_EcMSOMp = C_EcMSOMp*(N_EcM/C_EcM)
    N_EcMSOMa = C_EcMSOMa*(N_EcM/C_EcM)
    N_EcMSOMc = C_EcMSOMc*(N_EcM/C_EcM)
    N_ErMSOMp = 0.0!C_ErMSOMp*(N_ErM/C_ErM)
    N_ErMSOMa = 0.0!C_ErMSOMa*(N_ErM/C_ErM)
    N_ErMSOMc = 0.0!C_ErMSOMc*(N_ErM/C_ErM)
    if ( C_AM .LT. 1.175494351E-38 ) then
      N_AMSOMp=0.0
      N_AMSOMa=0.0
      N_AMSOMc=0.0
    else
      N_AMSOMp = max(C_AMSOMp*(N_AM/C_AM),1.175494351E-38) !Hack to ensure that the numbers can be written as floats to netcdf files
      N_AMSOMa = max(C_AMSOMa*(N_AM/C_AM),1.175494351E-38)
      N_AMSOMc = max(C_AMSOMc*(N_AM/C_AM),1.175494351E-38)
  end if
    
    !Dead saphrotroph biomass enters SOM pools
    N_SAPbSOMp = C_SAPbSOMp*N_SAPb/C_SAPb
    N_SAPbSOMa = C_SAPbSOMa*N_SAPb/C_SAPb
    N_SAPbSOMc = C_SAPbSOMc*N_SAPb/C_SAPb
    N_SAPfSOMp = C_SAPfSOMp*N_SAPf/C_SAPf
    N_SAPfSOMa = C_SAPfSOMa*N_SAPf/C_SAPf
    N_SAPfSOMc = C_SAPfSOMc*N_SAPf/C_SAPf

    !Desorption of SOMp to SOMa
    N_SOMpSOMa = C_SOMpSOMa*N_SOMp/C_SOMp

    !Transport from SOMc to SOMa:
    N_SOMcSOMa = C_SOMcSOMa*N_SOMc/C_SOMc

    !All N the Mycorrhiza dont need for its own, it gives to the plant:
    N_EcMPlant = N_INEcM + N_SOMpEcM + N_SOMcEcM - e_m*(1-enzyme_pct)*C_PlantEcM/CN_ratio(5)  !gN/m3h
    if ( N_EcMPlant .LT. 0.) then
      N_EcMPlant = 0.0
    end if
    

    if ( C_PlantAM == 0.0 ) then
      N_AMPlant = 0.0
    else 
      N_AMPlant = N_INAM  - e_m*C_PlantAM/CN_ratio(7)  !gN/m3h
      if ( N_AMPlant .LT. 0.) then
        N_AMPlant = 0.0
      end if
    end if    

    N_ErMPlant = 0.0


    !Calculate amount of inorganic N saprotrophs have access to: 
    N_for_sap  = N_IN + (Deposition - Leaching - N_INPlant - N_INEcM - N_INAM)*dt

    !total C uptake (growth + respiration) of saprotrophs
    U_sb = C_LITmSAPb + C_LITsSAPb + C_SOMaSAPb  
    U_sf = C_LITmSAPf + C_LITsSAPf + C_SOMaSAPf

    !Calculate SAP demand/excess of N (to ensure constant C:N ratio) 
    N_SAPbIN = (N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb) - CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
    N_SAPfIN = (N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf) - CUE_fungi_vr(depth)*U_sf/CN_ratio(4)

    !If there is not enough inorganic N to fill SAP demand, decrease CUE:
    do while (N_for_sap + (N_SAPbIN+N_SAPfIN)*dt <0)
      CUEmod_bacteria=0.9
      CUEmod_fungi=0.9
      CUE_bacteria_vr(depth)=CUE_bacteria_vr(depth)*CUEmod_bacteria
      CUE_fungi_vr(depth)=CUE_fungi_vr(depth)*CUEmod_fungi

      N_SAPbIN = (N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb) - CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
      N_SAPfIN = N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf - CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
    end do 
    
    nullify( C_LITm,C_LITs,C_SOMp,C_SOMa,C_SOMc,C_EcM,C_ErM,C_AM, C_SAPb,C_SAPf)
    nullify( N_LITm,N_LITs,N_SOMp,N_SOMa,N_SOMc,N_EcM,N_ErM,N_AM, N_SAPb,N_SAPf,N_NH4,N_NO3)

  end subroutine calculate_fluxes


  subroutine vertical_diffusion(tot_diffusion_dummy,upper_diffusion_flux,lower_diffusion_flux,pool_matrix,vert,D) !This subroutine calculates the vertical transport of carbon through the soil layers.

      real(r8), intent(in)   :: pool_matrix(:,:)
      real(r8), intent(out)  :: upper_diffusion_flux, lower_diffusion_flux
      real(r8), intent(out)  :: tot_diffusion_dummy ![gC/h]
      real(r8), allocatable, intent(out)  :: vert(:,:)
      real(r8), intent(in)   :: D

      !Local
      integer                :: depth, pool !For iteration
      integer,dimension(1)   :: max_pool, max_depth !For iteration

      allocate (vert, mold = pool_matrix)

      !Get how many depth levels and pools we will loop over.
      max_depth=shape(pool_matrix(:,1)) !TODO: Easier way to do this?
      max_pool=shape(pool_matrix(1,:))
      !In a timestep, the fluxes between pools in the same layer is calculated before the vertical diffusion. Therefore, a loop over all the entries in
      !pool_matrix is used here to calculate the diffusion.

      do depth = 1,max_depth(1)
        do pool =1, max_pool(1)
          !eq. 6.18 and 6.20 from Soetaert & Herman, A practical guide to ecological modelling.
          if (depth == 1) then
            upper_diffusion_flux= 0.0
            lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
          elseif (depth==max_depth(1)) then
            upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
            lower_diffusion_flux= 0.0
          else
            upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
            lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
          end if

          tot_diffusion_dummy=(upper_diffusion_flux-lower_diffusion_flux)/delta_z(depth)
          vert(depth,pool) = tot_diffusion_dummy
        end do !pool
      end do !depth


  end subroutine vertical_diffusion

  subroutine moisture_func(theta_l,theta_sat, theta_f,r_moist,nlevdecomp) !NOTE: Should maybe be placed somewhere else?
    integer :: nlevdecomp
    real(r8), intent(out), dimension(nlevdecomp) :: r_moist
    real(r8), intent(in), dimension(nlevdecomp)  :: theta_l, theta_sat, theta_f
    real(r8), dimension(nlevdecomp)  :: theta_frzn, theta_liq, air_filled_porosity
      !FROM mimics_cycle.f90 in testbed:
      ! ! Read in soil moisture data as in CORPSE
      !  theta_liq  = min(1.0, casamet%moistavg(npt)/soil%ssat(npt))     ! fraction of liquid water-filled pore space (0.0 - 1.0)
      !  theta_frzn = min(1.0, casamet%frznmoistavg(npt)/soil%ssat(npt)) ! fraction of frozen water-filled pore space (0.0 - 1.0)
      !  air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
      !
      !  if (mimicsbiome%fWFunction .eq. CORPSE) then
      !    ! CORPSE water scalar, adjusted to give maximum values of 1
      !    fW = (theta_liq**3 * air_filled_porosity**2.5)/0.022600567942709
      !    fW = max(0.05, fW)


    theta_liq  = min(1.0, theta_l/theta_sat)     ! fraction of liquid water-filled pore space (0.0 - 1.0)
    theta_frzn = min(1.0, theta_f/theta_sat)     ! fraction of frozen water-filled pore space (0.0 - 1.0)
    air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)

    r_moist = ((theta_liq**3)*air_filled_porosity**2.5)/0.022600567942709
    r_moist = max(0.05, r_moist)

  end subroutine moisture_func

  subroutine input_rates(layer_nr,f_EcM,&
                        LEAFC_TO_LIT,FROOTC_TO_LIT,LEAFN_TO_LIT,FROOTN_TO_LIT, &
                        C_MYCinput,N_CWD,C_CWD, &
                        C_inLITm,C_inLITs,N_inLITm,N_inLITs, C_inEcM,C_inAM, &
                        C_inSOMp,C_inSOMa,C_inSOMc)

    !NOTE: Which and how many layers that receives input from the "outside" (CLM history file) is hardcoded here. This may change in the future.
    !in:
    integer,  intent(in) :: layer_nr
    real(r8), intent(in) :: f_EcM
    real(r8), intent(in) :: LEAFC_TO_LIT
    real(r8), intent(in) :: FROOTC_TO_LIT
    real(r8), intent(in) :: LEAFN_TO_LIT
    real(r8), intent(in) :: FROOTN_TO_LIT    
    real(r8), intent(in) :: C_MYCinput
    real(r8), intent(in) :: N_CWD(:)
    real(r8), intent(in) :: C_CWD(:)
    
    !out:
    real(r8), intent(out) :: C_inLITm
    real(r8), intent(out) :: C_inLITs
    real(r8), intent(out) :: N_inLITm
    real(r8), intent(out) :: N_inLITs
    real(r8), intent(out) :: C_inEcM
    real(r8), intent(out) :: C_inAM
    real(r8), intent(out) :: C_inSOMp
    real(r8), intent(out) :: C_inSOMa
    real(r8), intent(out) :: C_inSOMc
    
    !local:
    real(r8)           :: totC_LIT_input
    real(r8)           :: totN_LIT_input

    totC_LIT_input= FROOTC_TO_LIT*froot_prof(layer_nr) + LEAFC_TO_LIT*leaf_prof(layer_nr) !gC/m3h 
    C_inLITm   = fMET*totC_LIT_input*(1-f_met_to_som)
    C_inLITs   = (1-fMET)*totC_LIT_input + C_CWD(layer_nr)
    C_inSOMp = fMET*totC_LIT_input*f_met_to_som*fPHYS(1)
    C_inSOMc = fMET*totC_LIT_input*f_met_to_som*fCHEM(1)
    C_inSOMa = fMET*totC_LIT_input*f_met_to_som*fAVAIL(1)
        
    totN_LIT_input = FROOTN_TO_LIT*froot_prof(layer_nr) + LEAFN_TO_LIT*leaf_prof(layer_nr)!gN/m3h 
    N_inLITm    = fMET*totN_LIT_input
    N_inLITs    = (1-fMET)*totN_LIT_input + N_CWD(layer_nr)
            
    C_inEcM = f_EcM*C_MYCinput*froot_prof(layer_nr)
    C_inAM = (1-f_EcM)*C_MYCinput*froot_prof(layer_nr)

  end subroutine input_rates
  
  function calc_EcMfrac(PFT_dist) result(EcM_frac)
    real(r8)                             :: EcM_frac
    real(r8),dimension(15)               :: PFT_dist
    real(r8),dimension(15),parameter     :: EcM_fraction=(/1.,1.,1.,1.,0.,0.,0.,0.5,1.,1.,1.,1.,1.,0.,0./)
    integer                              :: i
    EcM_frac = 0.0
    do i = 1, 15, 1
        EcM_frac = EcM_frac + PFT_dist(i)*EcM_fraction(i)
    end do
    EcM_frac = EcM_frac/100.
  end function calc_EcMfrac

end module fluxMod
