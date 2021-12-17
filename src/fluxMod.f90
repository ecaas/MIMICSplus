module fluxMod
  use paramMod
  use dispmodule !External module to pretty print matrices (mainly for testing purposes)
  implicit none

  contains

  function set_N_dep(CLMdep,const_dep) result(Deposition)
    real(r8)           :: Deposition
    real(r8), optional :: CLMdep
    real(r8), optional :: const_dep
    Deposition = -999
    if (present(CLMdep) .and. .not. present(const_dep)) then
      Deposition = CLMdep
    elseif (present(const_dep) .and. .not. present(CLMdep)) then
      Deposition = const_dep/(hr_pr_yr*soil_depth)
    elseif (.not. present(CLMdep) .and. .not. present(const_dep)) then
      print*, "N dep not set correctly, stopping"
      Deposition = -999
      stop
    end if
  end function set_N_dep

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
    K_m      = exp(Kslope*temperature + Kint)*a_k*Kmod               ![mgC/cm3]*10e3=[gC/m3]
    
  end function Km_function

  function Vmax_function(temperature, moisture) result(V_max)
    real(r8),dimension(MM_eqs)             :: V_max
    real(r8), intent(in)                   :: temperature
    real(r8), intent(in)                   :: moisture
    V_max    = exp(Vslope*temperature + Vint)*a_v*Vmod*moisture   ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?
  end function Vmax_function

  subroutine calculate_fluxes(depth,nlevdecomp,C_pool_matrix,N_pool_matrix) !This subroutine calculates the fluxes in and out of the SOM pools.
    integer         :: depth !depth level
    integer         :: nlevdecomp
    real(r8),target :: C_pool_matrix(nlevdecomp, pool_types)
    real(r8),target :: N_pool_matrix(nlevdecomp, pool_types_N)
    !LOCAL:
    real(r8)  :: N_for_sap

    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: C_LITm, C_LITs, C_SOMp,C_SOMa,C_SOMc,C_EcM,C_ErM,C_AM, &
    C_SAPb, C_SAPf, N_LITm, N_LITs, N_SOMp,N_SOMa,N_SOMc,N_EcM,N_ErM,N_AM, N_SAPb, N_SAPf, N_IN
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
    N_IN => N_pool_matrix(depth, 11)
    

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
    C_EcMdecompSOMp = K_MO*soil_depth*C_EcM*C_SOMp   ![gC/m3h]
    C_EcMdecompSOMc = K_MO*soil_depth*C_EcM*C_SOMc   ![gC/m3h]
    

    !-----------------------------------NITROGEN FLUXES----------------------------:
    !Nitrogen aquired bymycorrhiza via oxidation of protected SOM pools.  gN/m3h
    !TODO: Should we also include SOMa?
    N_SOMpEcM = C_EcMdecompSOMp*N_SOMp/C_SOMp
    N_SOMcEcM = C_EcMdecompSOMc*N_SOMc/C_SOMc
    

    !Inorganic N taken up directly by plant roots
    N_InPlant = 4E-7*N_IN
    
    N_INEcM = V_max_myc*N_IN*(C_EcM/(C_EcM + Km_myc/soil_depth))   !NOTE: MMK parameters should maybe be specific to mycorrhizal type?
    N_INErM = 0.0!V_max_myc*N_IN*(C_ErM/(C_ErM + Km_myc/delta_z(depth)))   !Unsure about units
    N_INAM = V_max_myc*N_IN*(C_AM/(C_AM + Km_myc/delta_z(depth)))

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
    N_ErMSOMp = C_ErMSOMp*(N_ErM/C_ErM)
    N_ErMSOMa = C_ErMSOMa*(N_ErM/C_ErM)
    N_ErMSOMc = C_ErMSOMc*(N_ErM/C_ErM)
    N_AMSOMp = C_AMSOMp*(N_AM/C_AM)
    N_AMSOMa = C_AMSOMa*(N_AM/C_AM)
    N_AMSOMc = C_AMSOMc*(N_AM/C_AM)
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

    !Leaching based on Baskaran et al leaching rate:
    Leaching=L_rate*N_IN

    !All N the Mycorrhiza dont need for its own, it gives to the plant:
    N_EcMPlant = N_INEcM + N_SOMpEcM + N_SOMcEcM - e_m*C_PlantEcM/CN_ratio(5)  !gN/m3h
    if ( N_EcMPlant .LT. 0.) then
      N_EcMPlant = 0.0
    end if
    
    N_AMPlant = N_INAM  - e_m*C_PlantAM/CN_ratio(7)  !gN/m3h
    if ( N_AMPlant .LT. 0.) then
      N_AMPlant = 0.0
    end if    


    N_ErMPlant = 0.0


    !Calculate amount of inorganic N saprotrophs have access to: 
    N_for_sap  = N_IN + Deposition - Leaching - N_INPlant - N_INEcM - N_INAM
    ! print*, N_IN,Deposition,Leaching,N_INPlant,N_INEcM, N_for_sap
    ! print*, N_for_sap + N_SAPbIN+N_SAPfIN, N_SAPbIN,N_SAPfIN
    ! print*,"*************"
    !total C uptake (growth + respiration) of saprotrophs
    U_sb = C_LITmSAPb + C_LITsSAPb + C_SOMaSAPb  
    U_sf = C_LITmSAPf + C_LITsSAPf + C_SOMaSAPf

    !Calculate SAP demand/excess of N (to ensure constant C:N ratio) 
    N_SAPbIN = (N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb) - CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
    N_SAPfIN = N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf - CUE_fungi_vr(depth)*U_sf/CN_ratio(4)

    !If there is not enough inorganic N to fill SAP demand, decrease CUE:
    do while (N_for_sap + N_SAPbIN+N_SAPfIN <0)
      CUEmod_bacteria=0.9
      CUEmod_fungi=0.9
      CUE_bacteria_vr(depth)=CUE_bacteria_vr(depth)*CUEmod_bacteria
      CUE_fungi_vr(depth)=CUE_fungi_vr(depth)*CUEmod_fungi

      N_SAPbIN = (N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb) - CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
      N_SAPfIN = N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf - CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
    end do 
    
    nullify( C_LITm,C_LITs,C_SOMp,C_SOMa,C_SOMc,C_EcM,C_ErM,C_AM, C_SAPb,C_SAPf)
    nullify( N_LITm,N_LITs,N_SOMp,N_SOMa,N_SOMc,N_EcM,N_ErM,N_AM, N_SAPb,N_SAPf,N_IN)
    
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

  subroutine input_rates(layer_nr,&
                                  TOTC_LITFALL,LEAFC_TO_LIT,FROOTC_TO_LIT,LEAFN_TO_LIT,FROOTN_TO_LIT, &
                                  C_EcMinput,N_LEACHinput,&
                                  N_CWD,C_CWD, &
                                  C_PlantLITm,C_PlantLITs,&
                                  N_PlantLITm,N_PlantLITs,&
                                  N_LEACH,C_PlantEcM,C_PlantAM)

    !NOTE: Which and how many layers that receives input from the "outside" (CLM history file) is hardcoded here. This may change in the future.
    !in:
    integer,  intent(in) :: layer_nr
    real(r8), intent(in) :: TOTC_LITFALL
    real(r8), intent(in) :: LEAFC_TO_LIT
    real(r8), intent(in) :: FROOTC_TO_LIT
    real(r8), intent(in) :: LEAFN_TO_LIT
    real(r8), intent(in) :: FROOTN_TO_LIT    
    real(r8), intent(in) :: C_EcMinput
    real(r8), intent(in) :: N_LEACHinput(:)
    real(r8), intent(in) :: N_CWD(:)
    real(r8), intent(in) :: C_CWD(:)
    
    !out:
    real(r8), intent(out) :: C_PlantLITm
    real(r8), intent(out) :: C_PlantLITs
    real(r8), intent(out) :: N_PlantLITm
    real(r8), intent(out) :: N_PlantLITs
    real(r8), intent(out) :: N_LEACH
    real(r8), intent(out) :: C_PlantEcM
    real(r8), intent(out) :: C_PlantAM
    

    !local:
    real(r8)           :: totC_LIT_input
    real(r8)           :: totN_LIT_input

    totC_LIT_input= FROOTC_TO_LIT*froot_prof(layer_nr) + LEAFC_TO_LIT*leaf_prof(layer_nr) !+ mortality*froot_prof(layer_nr) !gC/m3h 
    C_PlantLITm   = fMET*totC_LIT_input*(1-f_met_to_som)
    C_PlantLITs   = (1-fMET)*totC_LIT_input + C_CWD(layer_nr)
    C_PlantSOMp = fMET*totC_LIT_input*f_met_to_som*fPHYS(1)
    C_PlantSOMc = fMET*totC_LIT_input*f_met_to_som*fCHEM(1)
    C_PlantSOMa = fMET*totC_LIT_input*f_met_to_som*fAVAIL(1)
    
    
    totN_LIT_input = FROOTN_TO_LIT*froot_prof(layer_nr) + LEAFN_TO_LIT*leaf_prof(layer_nr)!gN/m3h 
    N_PlantLITm    = fMET*totN_LIT_input
    N_PlantLITs    = (1-fMET)*totN_LIT_input + N_CWD(layer_nr)


    C_PlantEcM = (C_EcMinput*froot_prof(layer_nr))
    C_PlantAM = 0.0!(C_EcMinput*froot_prof(layer_nr))
    
    !TODO: Figure out how to do mycorrhizal input vertically
  end subroutine input_rates

end module fluxMod
