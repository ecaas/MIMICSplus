module fluxMod
  use paramMod
  use dispmodule, only: disp !External module to pretty print matrices (mainly for testing purposes)
  use initMod, only: nlevels
  implicit none
  PRIVATE
  public :: nitrification,calc_Leaching,set_N_dep,calc_desorp,MMK_flux,input_rates,calculate_fluxes,vertical_diffusion,myc_to_plant
contains 
  function nitrification(nh4,t_scalar,w_scalar,soil_temp) result(f_nit)
    real(r8) :: f_nit
    !IN:
    real(r8),intent(in):: nh4 !gN/m3
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
    f_nit  = f_nit*(1._r8 - anaerobic_frac)

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
      Dep = const_dep
    elseif (.not. present(CLMdep) .and. .not. present(const_dep)) then
      print*, "N dep not set correctly, stopping"
      Dep = -999
      stop
    end if
  end function set_N_dep
  
  function calc_Leaching(drain,h2o_tot, N_NO3) result(Leach)
    real(r8)           :: Leach       !gN/m3 h
    real(r8)           :: drain       !mmH20/h = kgH20/m2 h
    real(r8)           :: h2o_tot     !kgH20/m2
    real(r8)           :: N_NO3 !gN/m3
    Leach = N_NO3*drain/h2o_tot
  end function calc_Leaching

  function calc_desorp(clay_fraction) result(d)
    real(r8)             :: d
    real(r8), intent(in) :: clay_fraction
    d = 1.5e-5*exp(-1.5*(clay_fraction))
  end function calc_desorp
  
  function MMK_flux(C_SAP,C_SUBSTRATE,MMK_nr) result(flux)
    !Compute C flux from substrate pool to saprotroph pool by using Michaelis Menten Kinetics.
    !NOTE: On the way, a fraction 1-CUE is lost as respiration. This is handeled in the "decomp" subroutine.
    real(r8):: flux ![gC/(m3 hr)]
    real(r8), intent(in) :: C_SAP
    real(r8), intent(in) :: C_SUBSTRATE
    integer, intent (in) :: MMK_nr
    !TODO: this works, but should not depend on Vmax & Km from mycmim mod
    flux = C_SAP*Vmax(MMK_nr)*C_SUBSTRATE/(Km(MMK_nr)+C_SUBSTRATE)
  end function MMK_flux
  
  function reverse_MMK_flux(C_SAP,C_SUBSTRATE,MMK_nr) result(flux)
     !Compute C flux from substrate pool to saprotroph pool by using Michaelis Menten Kinetics.
     !NOTE: On the way, a fraction 1-CUE is lost as respiration. This is handeled in the "decomp" subroutine.
     real(r8):: flux ![gC/(m3 hr)]
     real(r8), intent(in) :: C_SAP
     real(r8), intent(in) :: C_SUBSTRATE
     integer, intent (in) :: MMK_nr
     !TODO: this works, but should not depend on Vmax & Km from mycmim mod
     flux = C_SUBSTRATE*Vmax(MMK_nr)*C_SAP/(Km(MMK_nr)+C_SAP)
  end function reverse_MMK_flux

  subroutine mining_rates_Sulman(C_EcM,C_substrate,N_substrate,moisture_function,T,mining_mod, D_Cmine,D_Nmine) !Sulman et al 2019 eq 34-35 + max_mining modifier
    !NOTE: T dependence (Arrhenius) seems a bit weird, makes flux very low...
    !INPUT
    real(r8),intent(in) :: C_EcM
    real(r8),intent(in) :: C_substrate
    real(r8),intent(in) :: N_substrate
    real(r8),intent(in) :: moisture_function
    real(r8),intent(in) :: mining_mod
    real(r8), intent(in) :: T !Kelvin
    
    !OUTPUT
    real(r8),intent(out) :: D_Cmine
    real(r8),intent(out):: D_Nmine
    
    !LOCAL
    !NOTE: V_max(T) in article, but not sure how this temperature dependence is?
    real(r8),parameter :: V_max = 0.3/hr_pr_yr !Sulman 2019 supplement page 7, Assumed SOMp,SOMc ~ slow SOM
    real(r8),parameter :: K_m = 0.015 
    real(r8),parameter :: E_a = 54000 !J/mol
    real(r8),parameter :: R   = 8.31 !J/(K mol)
    D_Cmine = V_max*exp(-E_a/(R*T))*moisture_function*C_substrate*((C_EcM/C_substrate)/(C_EcM/C_substrate+K_m))*mining_mod
    D_Nmine = V_max*exp(-E_a/(R*T))*moisture_function*N_substrate*((C_EcM/C_substrate)/(C_EcM/C_substrate+K_m))*mining_mod    
  end subroutine mining_rates_Sulman
  
  subroutine mining_rates_Baskaran(C_EcM,C_substrate,N_substrate,mining_mod,D_Cmine,D_Nmine) !Baskaran + max_mining modifier
    !INPUT
    real(r8),intent(in) :: C_EcM
    real(r8),intent(in) :: C_substrate
    real(r8),intent(in) :: N_substrate
    real(r8),intent(in) :: mining_mod
    
    !OUTPUT
    real(r8),intent(out) :: D_Cmine
    real(r8),intent(out):: D_Nmine

    D_Cmine = K_MO*soil_depth*C_EcM*C_substrate*mining_mod
    D_Nmine = D_Cmine*N_substrate/C_substrate
  end subroutine mining_rates_Baskaran
  
  subroutine myc_to_plant(layer_nr,use_enz,C_PlantEcM,C_PlantAM,N_AMPlant,N_EcMPlant,N_ErMPlant,CUE_EcM,CUE_AM,enzyme_prod)
    !INPUT
    integer,intent(in)  :: layer_nr 
    logical,intent(in)  :: use_enz
    real(r8),intent(in) :: C_PlantEcM
    real(r8),intent(in) :: C_PlantAM
    
    
    !INOUT: 
    real(r8), intent(inout) :: CUE_EcM
    real(r8), intent(inout) :: CUE_AM
    real(r8), intent(inout) :: enzyme_prod
    
    !OUTPUT
    real(r8),intent(out)  ::  N_AMPlant
    real(r8), intent(out) ::  N_EcMPlant
    real(r8), intent(out) ::  N_ErMPlant
    
    !LOCAL
    real(r8) ::     AM_N_demand
    real(r8) ::     AM_N_uptake
    real(r8) ::     EcM_N_demand
    real(r8) ::     EcM_N_uptake
    !----------------------------------------------------------------------------------------------------------------------------------
    !All N the Mycorrhiza dont need for its own, it gives to the plant:
    AM_N_demand = CUE_AM*C_PlantAM/CN_ratio(7)
    AM_N_uptake = N_INAM     
    if ( AM_N_uptake >= AM_N_demand ) then   
      N_AMPlant = AM_N_uptake - AM_N_demand
    else
      N_AMPlant = (1-f_growth)*AM_N_uptake
      CUE_AM = f_growth*AM_N_uptake*CN_ratio(7)/(C_PlantAM)
    end if
    if ( N_AMPlant .ne. 0 ) then      
      if ( abs(N_AMPlant) < 1e-18 ) then
        !print*, N_AMPlant, "N_AMPlant", layer_nr
        N_AMPlant=0.0
      end if
    end if
    !All N the Mycorrhiza dont need for its own, it gives to the plant:
    EcM_N_demand = (CUE_EcM*(1-enzyme_prod)*C_PlantEcM+C_SOMcEcM+C_SOMpEcM)/CN_ratio(5)
    EcM_N_uptake = N_INEcM + N_SOMpEcM + N_SOMcEcM 
    if ( EcM_N_uptake >= EcM_N_demand ) then   
        N_EcMPlant=EcM_N_uptake-EcM_N_demand      
    else
        N_EcMPlant = (1-f_growth)*EcM_N_uptake
        if ( use_ENZ ) then
          enzyme_prod = 1 - (f_growth*EcM_N_uptake*CN_ratio(5)-C_SOMpEcM-C_SOMcEcM)/(CUE_EcM*C_PlantEcM)
        else
          CUE_EcM = (f_growth*EcM_N_uptake*CN_ratio(5)-(C_SOMcEcM+C_SOMpEcM))/((1-f_enzprod(layer_nr))*C_PlantEcM)
        end if
      end if
    if ( N_EcMPlant .ne. 0 ) then          
      if ( abs(N_EcMPlant) < 1e-18 ) then
        !print*, N_EcMPlant, "N_EcMPlant",layer_nr
        N_EcMPlant=0.0
      end if
    end if 
    N_ErMPlant = 0.0
  end subroutine myc_to_plant 
  
  subroutine calculate_fluxes(depth,sulman_mining,Temp_Celsius,C_pool_matrix,N_pool_matrix,dt) !This subroutine calculates the fluxes in and out of the SOM pools.
    integer,intent(in)         :: depth !depth level
    logical,intent(in)         :: sulman_mining
    real(r8), intent(in)       :: Temp_Celsius
    real(r8),intent(in)       :: dt ! timestep
    real(r8),target :: C_pool_matrix(nlevels, pool_types)
    real(r8),target :: N_pool_matrix(nlevels, pool_types_N)
    
    !LOCAL:
    real(r8)  :: N_for_sap
    real(r8)  :: N_IN
    real(r8)  :: f_b
    real(r8)  :: minedSOMp
    real(r8)  :: minedSOMc
    real(r8)  :: Temp_Kelvin
    
    
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
    
    Temp_Kelvin = Temp_Celsius+abs_zero
    
    !------------------CARBON FLUXES----------------------------:
    !Decomposition of LIT and SOMa by SAP:
    !On the way, a fraction 1-CUE is lost as respiration. This is handeled in the "decomp" subroutine.
    C_LITmSAPb=MMK_flux(C_SAPb,C_LITm,1)
    C_LITsSAPb=MMK_flux(C_SAPb,C_LITs,2)
    C_SOMaSAPb=MMK_flux(C_SAPb,C_SOMa,3)
    C_LITmSAPf=MMK_flux(C_SAPf,C_LITm,4)
    C_LITsSAPf=MMK_flux(C_SAPf,C_LITs,5)
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
    C_SAPbSOMp=C_SAPb*k_sapsom(1)*fPHYS(1)   !gC/m3h
    C_SAPbSOMa=C_SAPb*k_sapsom(1)*fAVAIL(1)
    C_SAPbSOMc=C_SAPb*k_sapsom(1)*fCHEM(1)

    C_SAPfSOMp=C_SAPf*k_sapsom(2)*fPHYS(2)
    C_SAPfSOMa=C_SAPf*k_sapsom(2)*fAVAIL(2)
    C_SAPfSOMc=C_SAPf*k_sapsom(2)*fCHEM(2)

    !Desorbtion controls transport from physically protected to available SOM
    C_SOMpSOMa=C_SOMp*desorp

    !Oxidation from SOMc to SOMa
    !From equations for decomposing structural litter in mimics,eq. A10
    !KO modifies Km which is used in the litter->SAP equations.
    C_SOMcSOMa    = ( C_SAPb * Vmax(2) * C_SOMc / (KO(1)*Km(2) + C_SOMc)) + &
                   (C_SAPf * Vmax(5) * C_SOMc / (KO(2)*Km(5) + C_SOMc))

    !Ectomycorrhizal mining options:
    if ( sulman_mining ) then
      call mining_rates_Sulman(C_EcM,C_SOMc,N_SOMc,r_moist(depth),Temp_Kelvin,input_mod, minedSOMc,N_SOMcEcM)
      call mining_rates_Sulman(C_EcM,C_SOMp,N_SOMp,r_moist(depth),Temp_Kelvin, input_mod, minedSOMp,N_SOMpEcM)
    else
      call mining_rates_Baskaran(C_EcM,C_SOMp,N_SOMp,input_mod,minedSOMp,N_SOMpEcM) 
      call mining_rates_Baskaran(C_EcM,C_SOMc,N_SOMc,input_mod,minedSOMc,N_SOMcEcM)               
    end if
                  
    C_SOMpEcM = minedSOMp*f_use
    C_EcMdecompSOMp = (1_r8-f_use)*minedSOMp   ![gC/m3h]

    C_SOMcEcM = minedSOMc*f_use
    C_EcMdecompSOMc = (1_r8-f_use)*minedSOMc   ![gC/m3h]

    
    !-----------------------------------NITROGEN FLUXES----------------------------:
    N_IN = N_NH4+ N_NO3+(Deposition - Leaching)*dt !
    if ( N_IN < 0._r8 ) then
      print*, "Negative inorganic N pool at layer", depth, N_IN
      print*, N_NH4,N_NO3,Deposition,Leaching
    end if
    !Nitrogen aquired bymycorrhiza via oxidation of protected SOM pools.  gN/m3h
        
    !Inorganic N taken up directly by plant roots
    N_InPlant = 5E-7*N_IN
    
    N_INEcM = V_max_myc*N_IN*(C_EcM/(C_EcM + Km_myc/soil_depth))*input_mod !NOTE: MMK parameters should maybe be specific to mycorrhizal type?
    if ( N_INEcM .NE. 0.0 ) then
        N_INEcM=max(N_INEcM,1.175494351E-38)
    end if
    N_INErM = 0.0
    N_INAM = V_max_myc*N_IN*(C_AM/(C_AM + Km_myc/soil_depth))*input_mod
    if ( N_INAM .NE. 0.0 ) then
        N_INAM=max(N_INEcM,1.175494351E-38)
    end if

    !Decomposition of LIT and SOMa by SAP
    N_LITmSAPb = calc_parallel_Nrates(C_LITmSAPb,N_LITm,C_LITm)
    N_LITsSAPb = calc_parallel_Nrates(C_LITsSAPb,N_LITs,C_LITs)

    N_LITmSAPf = calc_parallel_Nrates(C_LITmSAPf,N_LITm,C_LITm)
    N_LITsSAPf = calc_parallel_Nrates(C_LITsSAPf,N_LITs,C_LITs)

    N_SOMaSAPb = calc_parallel_Nrates(C_SOMaSAPb,N_SOMa,C_SOMa)
    N_SOMaSAPf = calc_parallel_Nrates(C_SOMaSAPf,N_SOMa,C_SOMa)

    !Dead mycorrhizal biomass enters SOM pools
    N_EcMSOMp = calc_parallel_Nrates(C_EcMSOMp,N_EcM,C_EcM)
    N_EcMSOMa = calc_parallel_Nrates(C_EcMSOMa,N_EcM,C_EcM)
    N_EcMSOMc = calc_parallel_Nrates(C_EcMSOMc,N_EcM,C_EcM)
    N_ErMSOMp = 0.0!C_ErMSOMp*(N_ErM/C_ErM)
    N_ErMSOMa = 0.0!C_ErMSOMa*(N_ErM/C_ErM)
    N_ErMSOMc = 0.0!C_ErMSOMc*(N_ErM/C_ErM)

    N_AMSOMp = calc_parallel_Nrates(C_AMSOMp,N_AM,C_AM)
    N_AMSOMa = calc_parallel_Nrates(C_AMSOMa,N_AM,C_AM)
    N_AMSOMc = calc_parallel_Nrates(C_AMSOMc,N_AM,C_AM)

    !Dead saphrotroph biomass enters SOM pools
    N_SAPbSOMp = calc_parallel_Nrates(C_SAPbSOMp,N_SAPb,C_SAPb)
    N_SAPbSOMa = calc_parallel_Nrates(C_SAPbSOMa,N_SAPb,C_SAPb)
    N_SAPbSOMc = calc_parallel_Nrates(C_SAPbSOMc,N_SAPb,C_SAPb)
    N_SAPfSOMp = calc_parallel_Nrates(C_SAPfSOMp,N_SAPf,C_SAPf)
    N_SAPfSOMa = calc_parallel_Nrates(C_SAPfSOMa,N_SAPf,C_SAPf)
    N_SAPfSOMc = calc_parallel_Nrates(C_SAPfSOMc,N_SAPf,C_SAPf)

    !Desorption of SOMp to SOMa
    N_SOMpSOMa = calc_parallel_Nrates(C_SOMpSOMa,N_SOMp,C_SOMp)

    !Transport from SOMc to SOMa:
    N_SOMcSOMa = calc_parallel_Nrates(C_SOMcSOMa,N_SOMc,C_SOMc)

    !---------------------------------------------------------------------------------------------------------------------------------------
    !Calculate amount of inorganic N saprotrophs have access to: 
    N_for_sap  = (N_IN - ( N_INPlant + N_INEcM + N_INAM)*dt)*pctN_for_sap
    ! print*, "max_Nimmobilized: ", max_Nimmobilized, depth,N_INEcM
    ! print*, "Fluxmod: " ,k2*(N_NH4+Deposition*dt-nitrif_rate*dt- ( N_INPlant + N_INEcM + N_INAM)*dt)
    ! print*, "N for SAP", N_for_sap
    !print*, N_for_sap, max_Nimmobilized*dt
    !total C uptake (growth + respiration) of saprotrophs
    U_sb = C_LITmSAPb + C_LITsSAPb + C_SOMaSAPb  
    U_sf = C_LITmSAPf + C_LITsSAPf + C_SOMaSAPf
    
    !total N uptake by saprotrophs
    UN_sb = (N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb)*NUE 
    UN_sf = (N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf)*NUE

    !SAP demand for N:
    N_demand_SAPb =  CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
    N_demand_SAPf =  CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
    
    !How much N saprotrophs need from the inorganic pool
    N_INSAPb = N_demand_SAPb-UN_sb
    N_INSAPf = N_demand_SAPf-UN_sf
    
    !Determine exchange of N between inorganic pool and saprotrophs:
    if ( N_INSAPb >= 0. .and. N_INSAPf >= 0. ) then !immobilization
      if ( max_Nimmobilized < abs((N_INSAPb + N_INSAPf)*dt) ) then !Not enough mineral N to meet demand

        f_b = N_INSAPb/(N_INSAPb + N_INSAPf) ! Bac. and fungi want the same inorganic N. This fraction determines how much N is available to each pool.
        CUE_bacteria_vr(depth)=((f_b*max_Nimmobilized+UN_sb)*CN_ratio(3))/(U_sb)
        CUE_fungi_vr(depth) = (((1-f_b)*max_Nimmobilized+UN_sf)*CN_ratio(4))/(U_sf)

        !SAP demand for N:
        N_demand_SAPb =  CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
        N_demand_SAPf =  CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
        
        N_INSAPb = f_b*max_Nimmobilized 
        N_INSAPf = (1-f_b)*max_Nimmobilized 
        c1a=c1a+1
      else !Enough mineral N to meet demand
        c1b=c1b+1
        continue
      end if    
        
    elseif ( N_INSAPb < 0. .and. N_INSAPf < 0. ) then !mineralization
        c2=c2+1
        continue
      
    elseif ( N_INSAPb >= 0. .and. N_INSAPf < 0. ) then ! bacteria can use N mineralized by fungi
      max_Nimmobilized = max_Nimmobilized + N_INSAPf 
      if ( max_Nimmobilized < N_INSAPb ) then
        CUE_bacteria_vr(depth)=((max_Nimmobilized+UN_sb)*CN_ratio(3))/U_sb
        N_demand_SAPb =  CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
        N_INSAPb = N_demand_SAPb-UN_sb    
         c3a=c3a+1 
            
      else
         c3b=c3b+1 
      end if
    elseif ( N_INSAPb < 0. .and. N_INSAPf >= 0. ) then !fungi can use N mineralized by bacteria
      max_Nimmobilized = max_Nimmobilized + N_INSAPb       
      if ( max_Nimmobilized < N_INSAPf ) then
        CUE_fungi_vr(depth)=((max_Nimmobilized+UN_sf)*CN_ratio(4))/U_sf
        
        N_demand_SAPf =  CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
        N_INSAPf = N_demand_SAPf-UN_sf
         c4a=c4a+1 
        
      else
         c4b=c4b+1 
      end if
      
    end if
    
    !---------------------------------------------------------------------------------------------------------------         
    
    nullify( C_LITm,C_LITs,C_SOMp,C_SOMa,C_SOMc,C_EcM,C_ErM,C_AM, C_SAPb,C_SAPf)
    nullify( N_LITm,N_LITs,N_SOMp,N_SOMa,N_SOMc,N_EcM,N_ErM,N_AM, N_SAPb,N_SAPf,N_NH4,N_NO3)
  end subroutine calculate_fluxes

  subroutine vertical_diffusion(tot_diffusion_dummy,upper_diffusion_flux,lower_diffusion_flux,pool_matrix,vert,D) !This subroutine calculates the vertical transport of carbon through the soil layers.
    !IN 
    real(r8), intent(in)   :: pool_matrix(:,:)
    real(r8), intent(in)   :: D
    !OUT
    real(r8), intent(out)  :: upper_diffusion_flux, lower_diffusion_flux
    real(r8), intent(out)  :: tot_diffusion_dummy ![gC/h]
    real(r8), allocatable, intent(out)  :: vert(:,:)
    
    !Local
    integer                :: depth, pool !For iteration
    integer,dimension(1)   :: max_pool, max_depth !For iteration

    allocate (vert, mold = pool_matrix)

      !Get how many depth levels and pools we will loop over.
      max_depth=shape(pool_matrix(:,1)) !TODO: Easier way to do this?
      max_pool=shape(pool_matrix(1,:))
      !print*, max_depth,max_pool
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

  subroutine input_rates(layer_nr, LEAFC_TO_LIT,FROOTC_TO_LIT,LEAFN_TO_LIT,&
                        FROOTN_TO_LIT,&
                        N_CWD,C_CWD, &
                        C_inLITm,C_inLITs,&
                        N_inLITm,N_inLITs, &
                        C_inSOMp,C_inSOMa,C_inSOMc, &
                        N_inSOMp,N_inSOMa,N_inSOMc)
                        

    !NOTE: Which and how many layers that receives input from the "outside" (CLM history file) is hardcoded here. This may change in the future.
    !in:
    integer,  intent(in) :: layer_nr
    real(r8), intent(in) :: LEAFC_TO_LIT
    real(r8), intent(in) :: FROOTC_TO_LIT
    real(r8), intent(in) :: LEAFN_TO_LIT
    real(r8), intent(in) :: FROOTN_TO_LIT    
    real(r8), intent(in) :: N_CWD(:)
    real(r8), intent(in) :: C_CWD(:)
    
    !out:
    real(r8), intent(out) :: C_inLITm
    real(r8), intent(out) :: C_inLITs
    real(r8), intent(out) :: N_inLITm
    real(r8), intent(out) :: N_inLITs
    real(r8), intent(out) :: C_inSOMp
    real(r8), intent(out) :: C_inSOMa
    real(r8), intent(out) :: C_inSOMc
    real(r8), intent(out) :: N_inSOMp
    real(r8), intent(out) :: N_inSOMa
    real(r8), intent(out) :: N_inSOMc
    
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
    N_inLITm    = fMET*totN_LIT_input*(1-f_met_to_som)
    N_inLITs    = (1-fMET)*totN_LIT_input + N_CWD(layer_nr)
    N_inSOMp = fMET*totN_LIT_input*f_met_to_som*fPHYS(1)
    N_inSOMc = fMET*totN_LIT_input*f_met_to_som*fCHEM(1)
    N_inSOMa = fMET*totN_LIT_input*f_met_to_som*fAVAIL(1)        
    
  end subroutine input_rates
  
  function calc_parallel_Nrates(C_rate,N_pool,C_pool) result(N_rate)
    !in: 
    real(r8),intent(in) :: C_rate
    real(r8),intent(in) :: N_pool
    real(r8),intent(in) :: C_pool
    
    !out
    real(r8) :: N_rate
     
    if ( abs(C_pool) < 1E-18 ) then
      N_rate=0.0
    else
      N_rate=C_rate*(N_pool/C_pool)
    end if
  end function calc_parallel_Nrates


end module fluxMod
