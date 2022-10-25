module fluxMod
  use paramMod
  use shr_kind_mod   , only : r8 => shr_kind_r8
  
  use dispmodule, only: disp !External module to pretty print matrices (mainly for testing purposes)
  use initMod, only: nlevels
  implicit none
  PRIVATE
  real(r8) :: NH4_sol_final, NH4_sorp_final,NO3_final
  public :: calc_nitrification,calc_Leaching,set_N_dep,forward_MMK_flux,reverse_MMK_flux,input_rates,calculate_fluxes,vertical_diffusion,myc_to_plant, NH4_sol_final,NH4_sorp_final, NO3_final
contains 
  
  function calc_nitrification(nh4,t_scalar,w_scalar,soil_temp) result(f_nit)
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
  end function calc_nitrification
  
  function set_N_dep(CLMdep,const_dep) result(Dep)
    real(r8)           :: Dep
    real(r8), optional :: CLMdep
    real(r8), optional :: const_dep
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
    real(r8),intent(in)           :: drain       !mmH20/h = kgH20/m2 h
    real(r8),intent(in)           :: h2o_tot     !kgH20/m2
    real(r8),intent(in)           :: N_NO3       !gN/m3
    Leach = N_NO3*drain/h2o_tot
  end function calc_Leaching

  function forward_MMK_flux(C_SAP,C_SUBSTRATE,MMK_nr) result(flux)
    !Compute C flux from substrate pool to saprotroph pool by using Michaelis Menten Kinetics.
    !NOTE: On the way, a fraction 1-CUE is lost as respiration. This is handeled in the "decomp" subroutine.
    real(r8):: flux ![gC/(m3 hr)]
    real(r8), intent(in) :: C_SAP
    real(r8), intent(in) :: C_SUBSTRATE
    integer, intent (in) :: MMK_nr
    !TODO: this works, but should not depend on Vmax & Km from mycmim mod
    flux = C_SAP*Vmax(MMK_nr)*C_SUBSTRATE/(Km(MMK_nr)+C_SUBSTRATE)
  end function forward_MMK_flux

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
    !NOTE: V_max(T) in article, but not sure how this temperature dependence is?
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

    if ( .not. (C_substrate < epsilon(C_substrate)) ) then
      D_Cmine = K_MO*soil_depth*C_EcM*C_substrate*mining_mod
      D_Nmine = D_Cmine*N_substrate/C_substrate
    else 
      D_Cmine = 0.0_r8
      D_Nmine=0.0_r8
    end if
  end subroutine mining_rates_Baskaran

  subroutine myc_to_plant(CUE_EcM,CUE_AM,enzyme_prod,NAMPlant,NEcMPlant)
    
    !INOUT: 
    real(r8), intent(inout) :: CUE_EcM
    real(r8), intent(inout) :: CUE_AM
    real(r8), intent(inout) :: enzyme_prod
    
    !OUTPUT
    real(r8),intent(out)  ::  NAMPlant
    real(r8), intent(out) ::  NEcMPlant
    
    !LOCAL
    real(r8) ::     AM_N_demand
    real(r8) ::     AM_N_uptake
    real(r8) ::     EcM_N_demand
    real(r8) ::     EcM_N_uptake


    !All N the Mycorrhiza dont need for its own, it gives to the plant:
    AM_N_demand = CUE_AM*C_PlantAM/CN_ratio(6)
    AM_N_uptake = N_INAM     
    if ( AM_N_uptake >= AM_N_demand ) then   
      NAMPlant = AM_N_uptake - AM_N_demand
    else
      NAMPlant = (1-f_growth)*AM_N_uptake
      CUE_AM = f_growth*AM_N_uptake*CN_ratio(6)/(C_PlantAM)
    end if
    if ( abs(NAMPlant) < 1e-16 ) then
      save_N=save_N+NAMPlant
      NAMPlant=0.0
    end if
    !All N the Mycorrhiza dont need for its own, it gives to the plant:
    EcM_N_demand = (CUE_EcM*(1-enzyme_prod)*C_PlantEcM)/CN_ratio(5)
    EcM_N_uptake = N_INEcM + N_SOMpEcM + N_SOMcEcM 
    if ( EcM_N_uptake >= EcM_N_demand ) then   
        NEcMPlant=EcM_N_uptake-EcM_N_demand      
    else
        NEcMPlant = (1-f_growth)*EcM_N_uptake
        if ( use_ENZ ) then
          enzyme_prod = 1 - (f_growth*EcM_N_uptake*CN_ratio(5))/(CUE_EcM*C_PlantEcM)
        else
          CUE_EcM = (f_growth*EcM_N_uptake*CN_ratio(5))/((1-enzyme_prod)*C_PlantEcM)
        end if
    end if
    if ( abs(NEcMPlant) < 1e-16 ) then
      save_N=save_N+NEcMPlant
      
      NEcMPlant=0.0
    end if
  end subroutine myc_to_plant 

  subroutine calculate_fluxes(depth,Temp_Celsius,water_content,C_pool_matrix,N_pool_matrix,N_inorg_matrix,Deposition_rate, Leaching_rate, nitrification) !This subroutine calculates the fluxes in and out of the SOM pools.
    integer,intent(in)        :: depth !depth level    
    real(r8), intent(in)      :: Temp_Celsius
    real(r8),intent(in)       :: Deposition_rate
    real(r8),intent(in)       :: Leaching_rate
    real(r8),intent(in)       :: nitrification
    real(r8),intent(in)       :: water_content
    
    real(r8),target :: C_pool_matrix(nlevels, pool_types)
    real(r8),target :: N_pool_matrix(nlevels, pool_types_N)
    real(r8),target :: N_inorg_matrix(nlevels, inorg_N_pools)
    
    !LOCAL:
    real(r8)  :: NH4_sol_tmp
    real(r8)  :: nh4_sol_frac
    real(r8)  :: NH4_tot
    
    real(r8)  :: NO3_tmp
    real(r8)  :: N_IN
    real(r8)  :: f_b
    real(r8)  :: minedSOMp
    real(r8)  :: minedSOMc
    real(r8)  :: Temp_Kelvin
    real(r8)  :: U_sb, U_sf,UN_sb,UN_sf
    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: C_LITm, C_LITs, C_SOMp,C_SOMa,C_SOMc,C_EcM,C_AM, &
    C_SAPb, C_SAPf, N_LITm, N_LITs, N_SOMp,N_SOMa,N_SOMc,N_EcM,N_AM, N_SAPb, N_SAPf, &
    N_NH4_sol,N_NH4_sorp,N_NO3
    C_LITm => C_pool_matrix(depth, 1)
    C_LITs => C_pool_matrix(depth, 2)
    C_SAPb => C_pool_matrix(depth, 3)
    C_SAPf => C_pool_matrix(depth, 4)
    
    C_EcM =>  C_pool_matrix(depth, 5)
    C_AM =>   C_pool_matrix(depth, 6)
    C_SOMp => C_pool_matrix(depth, 7)
    C_SOMa => C_pool_matrix(depth, 8)
    C_SOMc => C_pool_matrix(depth, 9)
    
    N_LITm => N_pool_matrix(depth, 1)
    N_LITs => N_pool_matrix(depth, 2)
    N_SAPb => N_pool_matrix(depth, 3)
    N_SAPf => N_pool_matrix(depth, 4)
    N_EcM =>  N_pool_matrix(depth, 5)
    
    N_AM =>   N_pool_matrix(depth, 6)
    N_SOMp => N_pool_matrix(depth, 7)
    N_SOMa => N_pool_matrix(depth, 8)
    N_SOMc => N_pool_matrix(depth, 9)
    
    N_NH4_sol => N_inorg_matrix(depth, 1)
    N_NH4_sorp => N_inorg_matrix(depth, 2)
    N_NO3 => N_inorg_matrix(depth,3)
    
    Temp_Kelvin = Temp_Celsius+abs_zero
    
    !------------------CARBON FLUXES----------------------------:
    !Decomposition of LIT and SOMa by SAP:
    !On the way, a fraction 1-CUE is lost as respiration. This is handeled in the "decomp" subroutine.
    C_LITmSAPb=reverse_MMK_flux(C_SAPb,C_LITm,1) !C6
    C_LITsSAPb=reverse_MMK_flux(C_SAPb,C_LITs,2) !C7
    C_SOMaSAPb=reverse_MMK_flux(C_SAPb,C_SOMa,3) !C8
    C_LITmSAPf=reverse_MMK_flux(C_SAPf,C_LITm,4) !C9
    C_LITsSAPf=reverse_MMK_flux(C_SAPf,C_LITs,5) !C10
    C_SOMaSAPf=reverse_MMK_flux(C_SAPf,C_SOMa,6) !C11
    
    !Oxidation from SOMc to SOMa
    !From equations for decomposing structural litter in mimics,eq. A10
    !KO modifies Km which is used in the litter->SAP equations.
    C_SOMcSOMa    = ( C_SAPb * Vmax(2) * C_SOMc / (KO(1)*Km(2) + C_SAPb)) + &
    (C_SAPf * Vmax(5) * C_SOMc / (KO(2)*Km(5) + C_SAPf)) !C12
    
    !Desorbtion controls transport from physically protected to available SOM
    C_SOMpSOMa=C_SOMp*desorp !C13
    
    !Turnover from SAP to SOM. Based on the turnover equations used in mimics for flux from microbial pools to SOM pools (correspond to eq A4,A8 in Wieder 2015)
    C_SAPbSOMp=C_SAPb*k_sapsom(1)*fPHYS(1)   !gC/m3h !C14
    C_SAPbSOMc=C_SAPb*k_sapsom(1)*fCHEM(1) !C15
    C_SAPbSOMa=C_SAPb*k_sapsom(1)*fAVAIL(1)!C16
    
    C_SAPfSOMp=C_SAPf*k_sapsom(2)*fPHYS(2) !C17 
    C_SAPfSOMc=C_SAPf*k_sapsom(2)*fCHEM(2) !C18
    C_SAPfSOMa=C_SAPf*k_sapsom(2)*fAVAIL(2)!C19
    
    !Dead mycorrhizal biomass enters the SOM pools:  gC/m3h
    C_EcMSOMp=C_EcM*k_mycsom(1)*fEcMSOM(1)!somp !C20
    C_EcMSOMc=C_EcM*k_mycsom(1)*fEcMSOM(2)!somc !C21
    C_EcMSOMa=C_EcM*k_mycsom(1)*fEcMSOM(3)!soma !C22
    
    C_AMSOMp=C_AM*k_mycsom(2)*fAMSOM(1) !C23
    C_AMSOMc=C_AM*k_mycsom(2)*fAMSOM(2) !C24
    C_AMSOMa=C_AM*k_mycsom(2)*fAMSOM(3) !C25

    !Ectomycorrhizal mining options:
    if ( use_Sulman ) then
      call mining_rates_Sulman(C_EcM,C_SOMc,N_SOMc,r_moist(depth),Temp_Kelvin,input_mod, minedSOMc,N_SOMcEcM)
      call mining_rates_Sulman(C_EcM,C_SOMp,N_SOMp,r_moist(depth),Temp_Kelvin, input_mod, minedSOMp,N_SOMpEcM)
    else
      call mining_rates_Baskaran(C_EcM,C_SOMp,N_SOMp,input_mod,minedSOMp,N_SOMpEcM) !N26
      call mining_rates_Baskaran(C_EcM,C_SOMc,N_SOMc,input_mod,minedSOMc,N_SOMcEcM) !N27              
    end if
    
    C_EcMdecompSOMp = minedSOMp   ![gC/m3h] !C26 !NOTE Can drop minedSOMx and define C_EcMdecompSOMx directly
    C_EcMdecompSOMc = minedSOMc   ![gC/m3h] !C27
    
  
    !-----------------------------------NITROGEN FLUXES----------------------------:
    !Decomposition of LIT and SOMa by SAP
    N_LITmSAPb = calc_parallel_Nrates(C_LITmSAPb,N_LITm,C_LITm) !N6
    N_LITsSAPb = calc_parallel_Nrates(C_LITsSAPb,N_LITs,C_LITs) !N7
    N_SOMaSAPb = calc_parallel_Nrates(C_SOMaSAPb,N_SOMa,C_SOMa) !N8
    N_LITmSAPf = calc_parallel_Nrates(C_LITmSAPf,N_LITm,C_LITm) !N9
    N_LITsSAPf = calc_parallel_Nrates(C_LITsSAPf,N_LITs,C_LITs) !N10
    N_SOMaSAPf = calc_parallel_Nrates(C_SOMaSAPf,N_SOMa,C_SOMa) !N11
    !Transport from SOMc to SOMa:
    N_SOMcSOMa = calc_parallel_Nrates(C_SOMcSOMa,N_SOMc,C_SOMc) !N12
    !Desorption of SOMp to SOMa
    N_SOMpSOMa = calc_parallel_Nrates(C_SOMpSOMa,N_SOMp,C_SOMp) !N13
    !Dead saphrotroph biomass enters SOM pools
    N_SAPbSOMp = calc_parallel_Nrates(C_SAPbSOMp,N_SAPb,C_SAPb) !N14
    N_SAPbSOMc = calc_parallel_Nrates(C_SAPbSOMc,N_SAPb,C_SAPb) !N15
    N_SAPbSOMa = calc_parallel_Nrates(C_SAPbSOMa,N_SAPb,C_SAPb) !N16
    N_SAPfSOMp = calc_parallel_Nrates(C_SAPfSOMp,N_SAPf,C_SAPf) !N17
    N_SAPfSOMc = calc_parallel_Nrates(C_SAPfSOMc,N_SAPf,C_SAPf) !N18
    N_SAPfSOMa = calc_parallel_Nrates(C_SAPfSOMa,N_SAPf,C_SAPf) !N19
    
    !Dead mycorrhizal biomass enters SOM pools
    N_EcMSOMp = calc_parallel_Nrates(C_EcMSOMp,N_EcM,C_EcM) !N20
    N_EcMSOMa = calc_parallel_Nrates(C_EcMSOMa,N_EcM,C_EcM) !N21
    N_EcMSOMc = calc_parallel_Nrates(C_EcMSOMc,N_EcM,C_EcM) !N22
    
    N_AMSOMp = calc_parallel_Nrates(C_AMSOMp,N_AM,C_AM) !N23
    N_AMSOMa = calc_parallel_Nrates(C_AMSOMa,N_AM,C_AM) !N24
    N_AMSOMc = calc_parallel_Nrates(C_AMSOMc,N_AM,C_AM) !N25
    !*****************************************************************************
    !(1)Update inorganic pools to account for Leaching, deposition,  nitrification rate and gain from decomposition (1-NUE):
    NH4_sol_tmp = N_NH4_sol + (1-NUE)*(N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf+N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb)*dt + (Deposition_rate - nitrification)*dt
    NO3_tmp = N_NO3-Leaching_rate*dt + nitrification*dt
    call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
       
    !Inorganic N taken up directly by plant roots:
    N_InPlant = calc_plant_uptake(N_IN) !N34
    
    NH4_sol_tmp = NH4_sol_tmp - nh4_sol_frac*N_InPlant*dt
    NO3_tmp = max(NO3_tmp - (1-nh4_sol_frac)*N_InPlant*dt,0._r8)
    call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
          
    N_INEcM  = calc_myc_uptake(N_IN,C_EcM) !N28
    N_INAM   = calc_myc_uptake(N_IN,C_AM) !N29
    
    !(2)Update inorganic pools to account for uptake by plans and mycorrhizal fungi
    NH4_sol_tmp = NH4_sol_tmp - nh4_sol_frac*(N_INEcM+N_INAM)*dt
    NO3_tmp = max(NO3_tmp - (1-nh4_sol_frac)*(N_INEcM+N_INAM)*dt,0._r8)
    call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
  
    !total C uptake (growth + respiration) of saprotrophs
    U_sb = C_LITmSAPb + C_LITsSAPb + C_SOMaSAPb  
    U_sf = C_LITmSAPf + C_LITsSAPf + C_SOMaSAPf    
    ! N uptake by saprotrophs
    UN_sb = (N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb)*NUE 
    UN_sf = (N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf)*NUE
    !SAP demand for N:
    N_demand_SAPb =  CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
    N_demand_SAPf =  CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
    
    !How much N saprotrophs need from the inorganic pool
    N_INSAPb = N_demand_SAPb-UN_sb
    N_INSAPf = N_demand_SAPf-UN_sf
    
    !Determine exchange of N between inorganic pool and saprotrophs, N_INSAPb and N_INSAPf: !N36, !N37 is determined here.
    if ( N_INSAPb >= 0. .and. N_INSAPf >= 0. ) then !immobilization
      if ( N_IN < (N_INSAPb + N_INSAPf)*dt) then !Not enough mineral N to meet demand
        
        f_b = N_INSAPb/(N_INSAPb + N_INSAPf) ! Bac. and fungi want the same inorganic N. This fraction determines how much N is available to each pool.
        if ( U_sb ==0._r8 ) then !To avoid division by zero 
          N_INSAPb =0._r8
        else
          CUE_bacteria_vr(depth)=((f_b*N_IN+UN_sb*dt)*CN_ratio(3))/(U_sb*dt)
          N_demand_SAPb =  CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
          N_INSAPb = f_b*N_IN/dt
        end if
        if ( U_sf ==0._r8 ) then !To avoid division by zero
          N_INSAPf = 0._r8
        else
          CUE_fungi_vr(depth) = (((1-f_b)*N_IN+UN_sf*dt)*CN_ratio(4))/(U_sf*dt)
          N_demand_SAPf =  CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
          N_INSAPf = (1-f_b)*N_IN/dt
        end if
        
        !SAP demand for N:      
        c1a=c1a+1
      else !Enough mineral N to meet demand
        
        c1b=c1b+1
        continue
      end if    
      NO3_tmp = NO3_tmp - (1-nh4_sol_frac)*(N_INSAPb + N_INSAPf)
      NH4_sol_tmp = NH4_sol_tmp - nh4_sol_frac*(N_INSAPb + N_INSAPf)
      call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
        
    elseif ( N_INSAPb < 0. .and. N_INSAPf < 0. ) then !mineralization
      NH4_sol_tmp = NH4_sol_tmp - (N_INSAPb + N_INSAPf)
      call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)      
      c2=c2+1
      continue     
      
    elseif ( N_INSAPb >= 0. .and. N_INSAPf < 0. ) then ! bacteria can use N mineralized by fungi
      if ( (N_IN +abs(N_INSAPf)*dt) < N_INSAPb*dt ) then
        if ( U_sb ==0._r8 ) then !To avoid division by zero 
          N_INSAPb =0._r8
        else
          CUE_bacteria_vr(depth)=(((N_IN + abs(N_INSAPf)*dt)+UN_sb*dt)*CN_ratio(3))/(U_sb*dt)
          N_demand_SAPb =  CUE_bacteria_vr(depth)*U_sb/CN_ratio(3)
          N_INSAPb = N_demand_SAPb/dt-UN_sb   
        end if 
        c3a=c3a+1             
      else
        c3b=c3b+1 
      end if
      
      NO3_tmp = NO3_tmp - (1-nh4_sol_frac)*(N_INSAPb + N_INSAPf)
      NH4_sol_tmp = NH4_sol_tmp - nh4_sol_frac*(N_INSAPb + N_INSAPf)
      call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
      
    elseif ( N_INSAPb < 0. .and. N_INSAPf >= 0. ) then !fungi can use N mineralized by bacteria
      if ( (N_IN+ abs(N_INSAPb)*dt) < N_INSAPf*dt ) then
        if ( U_sf == 0._r8) then
          N_INSAPf = 0._r8
        else 
          CUE_fungi_vr(depth)=(( (N_IN + abs(N_INSAPb)*dt)+UN_sf*dt)*CN_ratio(4))/(U_sf*dt)
          N_demand_SAPf =  CUE_fungi_vr(depth)*U_sf/CN_ratio(4)
          N_INSAPf = N_demand_SAPf/dt-UN_sf
          N_INSAPf=N_IN/dt+abs(N_INSAPb)
          !nh4_sol_frac=calc_nh4_frac(NH4_sol_tmp+abs(N_INSAPb)*dt,NO3_tmp)
        end if
        
        c4a=c4a+1         
      else
        c4b=c4b+1 
      end if
      NO3_tmp = NO3_tmp - (1-nh4_sol_frac)*(N_INSAPb + N_INSAPf)
      NH4_sol_tmp = NH4_sol_tmp - nh4_sol_frac*(N_INSAPb + N_INSAPf)
      call update_inorganic_N(NO3_tmp,NH4_sol_tmp,N_IN,nh4_sol_frac)
                  
    else 
      print*, "No condition applies (this should not happen); ", C_LITmSAPf ,C_LITsSAPf , C_SOMaSAPf,C_LITmSAPb, C_LITsSAPb, C_SOMaSAPb, depth
      stop
      
    end if

    NH4_tot = NH4_sol_tmp + N_NH4_sorp    
    call calc_NH4_sol_sorp(NH4_tot,water_content,N_NH4_sorp,NH4_sorp_eq_vr(depth),NH4_sorp_final)
    NH4_sol_final = max(NH4_sol_tmp - (NH4_sorp_final-N_NH4_sorp),0._r8)
    NO3_final = max(NO3_tmp,0._r8)  
    
    N_IN = 0._r8 !reset values
    NH4_sol_tmp = 0._r8 !reset values
    NO3_tmp = 0._r8 !reset values
    nullify( C_LITm,C_LITs,C_SOMp,C_SOMa,C_SOMc,C_EcM,C_AM, C_SAPb,C_SAPf)
    nullify( N_LITm,N_LITs,N_SOMp,N_SOMa,N_SOMc,N_EcM,N_AM, N_SAPb,N_SAPf,N_NH4_sol,N_NH4_sorp,N_NO3)
  end subroutine calculate_fluxes

  subroutine calc_NH4_sol_sorp(NH4_tot,soil_water_frac,NH4_sorp_previous,NH4_sorp_eq,NH4_sorp)
    !IN:
    real(r8), intent(in)  :: NH4_tot   !g/m3, total NH4, both in soil solution and adsorbed
    real(r8), intent(in)  :: soil_water_frac   !m3water/m3soil (input from CLM data)
    real(r8), intent(in)  :: NH4_sorp_previous  !g/m3
    !Out:
    real(r8),intent(out)            :: NH4_sorp !g/m3, NH4 sorbed to particles
    real(r8),intent(out)            :: NH4_sorp_eq !g/m3, adsorbed NH4 at equilibrium
    
    real(r8), parameter :: BD_soil=1.6e6  !g/m3 (loam) soil from DOI: 10.3390/APP6100269 Table 1
    real(r8), parameter :: NH4_sorp_max = 0.09*BD_soil/mg_pr_g    !mg NH4 /g soil
    real(r8), parameter :: KL = 0.4      !L/mg
    real(r8)            :: KL_prime       !m3/g
    real(r8), parameter :: K_pseudo = 0.0167*mg_pr_g*60./BD_soil !m3/(g hour)

    !1) Calculate NH4_sorp_eq 
    KL_prime = KL*mg_pr_g*m3_pr_L/soil_water_frac 
    NH4_sorp_eq=(1+KL_prime*NH4_tot+NH4_sorp_max*KL_prime)/(2*KL_prime) - sqrt((1+KL_prime*NH4_tot+NH4_sorp_max*KL_prime)**2-4*KL_prime**2*NH4_sorp_max*NH4_tot)/(2*KL_prime)
    
    !2) Calculate NH4_sorp after adjusting towards equilibrium for 1 timestep
    if ( NH4_sorp_eq==NH4_sorp_previous ) then !Already at equilibrium
      NH4_sorp = NH4_sorp_previous
    elseif (NH4_sorp_eq > NH4_sorp_previous ) then !Adsorption
      NH4_sorp = NH4_sorp_eq - 1_r8/(1._r8/(NH4_sorp_eq-NH4_sorp_previous) + k_pseudo*dt)
    else !Desorption
      NH4_sorp = NH4_sorp_eq + 1_r8/(1._r8/(NH4_sorp_previous-NH4_sorp_eq) + k_pseudo*dt)
    end if
  end subroutine calc_NH4_sol_sorp

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
      !In a timestep, the fluxes between pools in the same layer is calculated before the vertical diffusion. Therefore, a loop over all the entries in
      !pool_matrix is used here to calculate the diffusion.
    do depth = 1,max_depth(1)
      do pool =1, max_pool(1)
        !eq. 6.18 and 6.20 from Soetaert & Herman, A practical guide to ecological modelling.
        if (depth == 1) then
          upper_diffusion_flux= 0._r8
          lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
        elseif (depth==max_depth(1)) then
          upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
          lower_diffusion_flux= 0._r8
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

    totC_LIT_input = FROOTC_TO_LIT*froot_prof(layer_nr) + LEAFC_TO_LIT*leaf_prof(layer_nr) !gC/m3h !NOTE: These do not include CWD
    totN_LIT_input = FROOTN_TO_LIT*froot_prof(layer_nr) + LEAFN_TO_LIT*leaf_prof(layer_nr)!gN/m3h 
    
    C_inLITm = fMET*totC_LIT_input*(1-f_met_to_som) !C1
    N_inLITm = fMET*totN_LIT_input*(1-f_met_to_som) !N1
    
    C_inLITs = ((1-fMET)*totC_LIT_input + C_CWD(layer_nr))*(1-f_struct_to_som) !C2
    N_inLITs = ((1-fMET)*totN_LIT_input + N_CWD(layer_nr))*(1-f_struct_to_som) !N2
    
    C_inSOMp = fMET*totC_LIT_input*f_met_to_som 
    C_inSOMc = ((1-fMET)*totC_LIT_input + C_CWD(layer_nr))*f_struct_to_som
    C_inSOMa = 0.0!fMET*totC_LIT_input*f_met_to_som*fAVAIL(1) 
    
    N_inSOMp = fMET*totN_LIT_input*f_met_to_som
    N_inSOMc = ((1-fMET)*totN_LIT_input + N_CWD(layer_nr))*f_struct_to_som
    N_inSOMa = 0.0!fMET*totN_LIT_input*f_met_to_som*fAVAIL(1)     
    
  end subroutine input_rates
  
  function calc_parallel_Nrates(C_rate,N_pool,C_pool) result(N_rate)
    !in: 
    real(r8),intent(in) :: C_rate
    real(r8),intent(in) :: N_pool
    real(r8),intent(in) :: C_pool
    
    !out
    real(r8) :: N_rate
     
    if ( C_pool <= epsilon(C_pool) ) then
      N_rate=0.0
    else
      N_rate=C_rate*(N_pool/C_pool)
    end if
  end function calc_parallel_Nrates

  subroutine update_inorganic_N(NO3,NH4_sol, &
    Ninorg_avail, ratio)
    real(r8), intent(in) :: NO3
    real(r8), intent(in) :: NH4_sol
    
    !out: 
    real(r8), intent(out) :: Ninorg_avail
    real(r8), intent(out) :: ratio
    
    Ninorg_avail = NO3+NH4_sol !The inorganic N available to microbes and plants (the rest, NH4_sorb_tmp is sorbed onto particles)
    if (NH4_sol+NO3 == 0._r8) Then 
      ratio = 0.5_8
    else
      ratio = NH4_sol/(NH4_sol+NO3)
    end if
  end subroutine update_inorganic_N
  
  function calc_plant_uptake(N_inorganic) result(N_INVeg)
    !IN 
    real(r8), INTENT(IN) :: N_inorganic
    
    real(r8)             :: N_INVeg
  
    N_INVeg = k_plant*N_inorganic
  end function calc_plant_uptake
  
  function calc_myc_uptake(N_inorganic,C_MYC) result(N_INMYC)
    !IN 
    real(r8), INTENT(IN) :: N_inorganic
    real(r8), INTENT(IN) :: C_MYC
    !OUT 
    real(r8) :: N_INMYC
    N_INMYC = V_max_myc*N_inorganic*(C_MYC/(C_MYC + Km_myc/soil_depth))*input_mod
    !NOTE: MMK parameters should maybe be specific to mycorrhizal type?

  end function calc_myc_uptake
  
end module fluxMod
