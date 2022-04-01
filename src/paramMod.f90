module paramMod
use shr_kind_mod   , only : r8 => shr_kind_r8
use initMod, only : nlevels
implicit none

!For conversion
integer, parameter :: sec_pr_hr = 60*60        
integer, parameter :: hr_pr_yr  = 365*24     
integer, parameter :: hr_pr_day = 24
integer, parameter :: days_in_year = 365
real(r8), parameter:: abs_zero=273.15 !Kelvin

integer, parameter, dimension(12)            :: days_in_month =(/31,28,31,30,31,30,31,31,30,31,30,31/)


!For calculating the Km parameter in Michaelis Menten kinetics (expressions based on mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and https://github.com/wwieder/MIMICS)
integer, parameter                           :: MM_eqs  = 6                     !Number of Michaelis-Menten parameters
real(kind=r8),dimension(MM_eqs),parameter    :: Kslope  = (/0.017, 0.027, 0.017, 0.017, 0.027, 0.017/) !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Vslope  = (/0.063, 0.063, 0.063, 0.063, 0.063, 0.063/) !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Kint    = 3.19      !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Vint    = 5.47      !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),parameter                      :: a_k     = 1e4 !Tuning parameter g/m3 (10 mg/cm3 from german et al 2012)
real(kind=r8),parameter                      :: a_v     = 8e-6 !Tuning parameter
real(kind=r8)                                :: pscalar 
real(kind=r8),dimension(MM_eqs)              :: Kmod   
real(kind=r8),dimension(MM_eqs)              :: Vmod    = (/10.0, 2.0, 10.0, 3.0,3.0, 2.0/) !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),parameter, dimension(2)        :: KO      =  4                 ![-]Increases Km (the half saturation constant for oxidation of chemically protected SOM, SOM_c) from mimics
real(kind=r8),dimension(MM_eqs)              :: Km                              ![mgC/cm3]*10e3=[gC/m3]
real(kind=r8),dimension(MM_eqs)              :: Vmax                            ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics.

!Pools: NOTE: This needs to be updated if pools are added to/removed from the system.
integer, parameter                           :: no_of_litter_pools = 2         !Metabolic and structural
integer, parameter                           :: no_of_sap_pools = 2            !SAP bacteria and SAP fungi
integer, parameter                           :: no_of_myc_pools = 3            !Mycorrhiza: Ecto, ericoid & arbuscular
integer, parameter                           :: no_of_som_pools = 3            !Physically protected, chemically protected, available carbon
integer, parameter                           :: pool_types = no_of_litter_pools + no_of_myc_pools + &
                                                no_of_sap_pools + no_of_som_pools
integer, parameter                           :: pool_types_N = pool_types + 2 !pool_types + NH4 + NO3

!For calculating turnover from SAP to SOM (expressions from mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and  https://github.com/wwieder/MIMICS)
real(r8),parameter                      :: fMET =0.6                       ![-] Fraction determining distribution of total litter production between LITm and LITs NOTE: Needs revision
real(r8), dimension(no_of_sap_pools)    :: fPHYS,fCHEM,fAVAIL              ![-]
real(r8), dimension(no_of_sap_pools)    :: k_sapsom  ![1/h] (k_sapsom in MIMICS)
real(r8)                                :: f_EcM! !fraction of present vegetation associated with EcM
real(r8),parameter                     :: pctN_for_sap=0.9 !NB: VERY ASSUMED Only this percentage of remaining inorganic N is avalable to SAPS

real(kind=r8)                          :: fCLAY                          ![-] fraction of clay in soil
real(kind=r8),dimension(3)             :: k_mycsom                        ![1/h] decay constants, MYC to SOM pools

real(r8), dimension(no_of_som_pools), parameter    :: fEcMSOM = (/0.4,0.4,0.2/) !somp,soma,somc. Fraction of flux from EcM to different SOM pools NOTE: assumed
real(r8), dimension(no_of_som_pools), parameter    :: fErMSOM = (/0.3,0.4,0.3/)
real(r8), dimension(no_of_som_pools), parameter    :: fAMSOM = (/0.3,0.3,0.4/)
real(r8)                                :: desorp ![1/h]From Mimics, used for the transport from physically protected SOM to available SOM pool

!Depth & vertical transport
real(r8)                             :: soil_depth           ![m] 
real(r8),dimension(25),parameter     :: node_z =  (/0.01,0.04,0.09,0.16,0.26,0.40,0.587,0.80,1.06,1.36,1.70,1.08,2.50,2.99,3.58,4.27,5.06,5.95,6.94,8.03,9.795,13.328,19.483,28.871,41.998/)!(/0.076,0.228, 0.380,0.532, 0.684,0.836,0.988,1.140,1.292,1.444/)!![m] Depth of center in each soil layer. Same as the first layers of default CLM5 with vertical resolution.
real(r8),dimension(25),parameter     :: delta_z = (/0.02, 0.04, 0.06, 0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.54,0.64,0.74,0.84,0.94,1.04,1.14,2.39,4.676,7.635,11.140,15.115/)!0.152![m] Thickness of each soil of the top layers in default clm5.
real(r8),parameter                   :: D_carbon = 0.0!1.14e-8![m2/h] Diffusivity. Based on Koven et al 2013, 1cm2/yr = 10e-4/(24*365)
real(r8),parameter                   :: D_nitrogen = 0.0!1.14e-8![m2/h] Diffusivity. Based on Koven et al 2013, 1cm2/yr = 10e-4/(24*365)

!counts: 
integer                             :: c1a
integer                             :: c1b
integer                             :: c2
integer                             :: c3a
integer                             :: c3b
integer                             :: c4a
integer                             :: c4b

real(r8)                            :: max_Nimmobilized
real(r8),PARAMETER                  :: k1 = 0.042 !hr-1 (day)
real(r8),PARAMETER                  :: k2 = 0.0014 !hr-1 (month)
real(r8),PARAMETER                  :: k3 = 0.00014 !hr-1 (year)

real(r8),dimension(:),allocatable    :: r_moist

real(r8),dimension(:),allocatable    :: CUE_bacteria_vr
real(r8),dimension(:),allocatable    :: CUE_fungi_vr
real(r8),dimension(:),allocatable    :: CUE_ecm_vr         !Growth efficiency of mycorrhiza 
real(r8),dimension(:),allocatable    :: CUE_am_vr         !Growth efficiency of mycorrhiza 
real(r8),dimension(:),allocatable    :: CUE_erm_vr        !Growth efficiency of mycorrhiza 
real(r8),parameter                   :: CUE_myc_0=0.25_r8 !Baskaran
real(r8),parameter                   :: NUE=0.7


real(r8),parameter                   :: CUE_0=0.5
real(r8),parameter                   :: CUE_slope=0.0!-0.016 !From German et al 2012

real(r8), parameter                  :: f_met_to_som=0.05_r8 ! fraction of metabolic litter flux that goes directly to SOM pools
real(r8)                             :: max_mining 
real(r8)                             :: input_mod 

real(r8),dimension(:),allocatable    :: f_enzprod 
real(r8),parameter                   :: f_enzprod_0=0.1_r8
real(r8), parameter                  :: f_use = 0.1_r8 !Fraction of C released during mining that is taken up by EcM

real(r8), parameter                  :: f_growth = 0.5_r8 !Fraction of mycorrhizal N uptake that must go to plant if there is too little N to support both giving 
                                                       !to plant and maximum growth. New CUE are calculated based on this. NB: VERY ASSUMED!!

real(r8), dimension(pool_types), parameter   :: CN_ratio = (/15,15,5,8,20,20,20,11,8,11/) !Fungi/bacteria: Tang, Riley, Maggi 2019 as in Mouginot et al. 2014
                                                                                          !NOTE: Wallander/Rousk may have data more suited for Boreal/Arctic conditions
                                                                                          !EcM: From Baskaran et al as in Wallander et al 2004
                                                                                          !SOM: From CLM documentation, table 21.3 (Mendeley version)
                                                                                          !LITm: MIMICS-CN manuscript
                                                                                          !LITs, ErM, AM: Guesses!

!From Baskaran et al 2016
real(r8), parameter :: Km_myc = 0.08            ![gNm-2] Half saturation constant of mycorrhizal uptake of inorganic N (called S_m in article) 
real(r8), parameter :: V_max_myc = 1.8/hr_pr_yr  ![g g-1 hr-1] Max mycorrhizal uptake of inorganic N (called K_mn in article) 



!Decomposition rates:
real(r8), parameter :: K_MO = 0.003_r8/hr_pr_yr ![m2gC-1hr-1] Mycorrhizal decay rate constant for oxidizable store     NOTE: vary from 0.0003 to 0.003 in article

!Moisture dependence (based on function used for MIMICS in the CASA-CNP testbed)
real(r8), parameter                          :: P = 44.247 !normalization of moisture function
real(r8)                                     :: gas_diffusion

real(r8),dimension(:),allocatable  :: ndep_prof
real(r8),dimension(:),allocatable  :: leaf_prof
real(r8),dimension(:),allocatable  :: froot_prof

integer(r8)     :: clock_rate,clock_start,clock_stop
integer(r8)     :: full_clock_rate,full_clock_start,full_clock_stop
!Fluxes etc:
real(r8) :: C_LITmSAPb, C_LITsSAPb, C_EcMSOMp, C_EcMSOMa, C_EcMSOMc, C_AMSOMp, &
            C_LITmSAPf, C_LITsSAPf, C_AMSOMa, C_AMSOMc, C_SOMaSAPb,C_SOMaSAPf, C_SOMpSOMa, C_SOMcSOMa, &
C_SAPbSOMa, C_SAPbSOMp, C_SAPbSOMc,C_SAPfSOMa, C_SAPfSOMp, C_SAPfSOMc, C_PlantSOMc,C_PlantSOMp,C_PlantSOMa, &
N_LITmSAPb, N_LITsSAPb, N_EcMSOMp, N_EcMSOMa, N_EcMSOMc,  N_AMSOMp, N_AMSOMa,&
N_AMSOMc, N_SOMaSAPb,N_SOMaSAPf, N_SOMpSOMa, N_SOMcSOMa, N_LITmSAPf, N_LITsSAPf, &
N_PlantLITs, N_PlantLITm, N_INPlant, N_INEcM,  N_INAM, N_EcMPlant, N_AMPlant, &
N_SAPbSOMa, N_SAPbSOMp, N_SAPbSOMc,N_SAPfSOMa, N_SAPfSOMp, N_SAPfSOMc,&
N_SOMcEcM,N_SOMpEcM, &
C_PlantEcM,  C_PlantAM, C_PlantLITm, C_PlantLITs, C_EcMdecompSOMp,C_EcMdecompSOMc, &
 Leaching, Deposition,nitrif_rate,f, U_sb, U_sf,UN_sb,UN_sf,N_demand_SAPf,N_demand_SAPb,N_INSAPb,N_INSAPf,&
 C_EcMdecompSOMa,N_PlantSOMp,N_PlantSOMa,N_PlantSOMc,C_SOMcEcM,C_SOMpEcM,C_EcMenz_prod,&
 C_ErMSOMp, C_ErMSOMa, C_ErMSOMc,C_PlantErM,N_INErM,N_ErMSOMc,N_ErMPlant, N_ErMSOMp, N_ErMSOMa
 
!For writing to file:
integer                                      :: ios = 0 !Changes if something goes wrong when opening a file
character (len=*),parameter                  :: output_path = './results/'

contains
  
  function calc_Kmod(clay_fraction) result(K_mod)
    !Input:
    real(r8), intent(in) :: clay_fraction
    !Output:
    real(r8), dimension(MM_eqs) :: K_mod
    !Local:
    real(r8)         :: p 
    
    p = 1.0/(2*exp(-2.0*sqrt(clay_fraction)))
    K_mod =  [real(r8) :: 0.125,0.5,0.25*p,0.5,0.25,0.167*p]
  end function calc_Kmod
  
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
  
  function ROI_function(N_aquired,C_myc, loss_rate) result(ROI) ! Based on Sulman et al 2019
    !INPUT
    real(r8),intent(in) :: N_aquired
    real(r8),intent(in) :: C_myc 
    real(r8),intent(in) :: loss_rate ![1/h]
    
    !OUTPUT
    real(r8) :: ROI
    
    !LOCAL
    real(r8), parameter :: eps = 0.5 !From Sulman et al supplementary: epsilon_mine, epsilon_scav
    real(r8) :: turnover ! [hour]
    
    turnover = 1/loss_rate 
    ROI=(N_aquired/C_myc)*turnover*eps      
  end function ROI_function 
  
  function r_input(C_input, max_input) result(mod) !Modifies N mining/scavegeing fluxes to avoid that mycorriza provides the plant with free N 
    !input
    real(r8) :: C_input
    real(r8) :: max_input
    !output
    real(r8) :: mod

    mod = C_input/(max_input)
  end function r_input 
    
  
  function calc_sap_to_som_fractions(clay_frac,met_frac) result(f_saptosom)
    !input
    real(r8),intent(in)      :: clay_frac
    real(r8),INTENT(IN)      :: met_frac
    !output
    real(r8),dimension(no_of_som_pools,no_of_sap_pools) :: f_saptosom
    
    f_saptosom(1,:) = [real(r8) :: 0.3*exp(1.3*clay_frac), 0.2*exp(0.8*clay_frac)]
    f_saptosom(2,:) = [real(r8) :: 0.1*exp(-3.0*met_frac), 0.3*exp(-3.0*met_frac)]
    f_saptosom(3,:) = 1 - (f_saptosom(1,:)+f_saptosom(2,:))
  end function calc_sap_to_som_fractions
  
  function calc_myc_mortality() result(myc_mortality)
    !NOTE: Is it better to call it turnover rate? Is there a difference?
    real(r8), dimension(no_of_myc_pools) :: myc_mortality
    myc_mortality=(/1.14_r8,1.14_r8,1.14_r8/)*1e-4  ![1/h]  1/yr  
  end function calc_myc_mortality

  function calc_sap_turnover_rate(met_frac,moist_modifier) result(turnover_rate)
    real(r8),INTENT(IN) :: met_frac
    real(r8),INTENT(IN) :: moist_modifier
    
    real(r8),dimension(no_of_sap_pools) :: turnover_rate

    turnover_rate = [real(r8) ::  5.2e-4*exp(0.3_r8*met_frac)*moist_modifier, 2.4e-4*exp(0.1_r8*met_frac)*moist_modifier]
  end function calc_sap_turnover_rate

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
  
  subroutine moisture_func(theta_l,theta_sat, theta_f,r_moist) !NOTE: Should maybe be placed somewhere else?
    real(r8), intent(out), dimension(nlevels) :: r_moist
    real(r8), intent(in), dimension(nlevels)  :: theta_l, theta_sat, theta_f
    real(r8), dimension(nlevels)  :: theta_frzn, theta_liq, air_filled_porosity
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
  
  function calc_Fmax(k,nh4) result(Fmax)
    !In:
    real(r8),intent(IN) :: k !loss rate [hr-1]
    real(r8),INTENT(IN) :: nh4 !nh4 consentration [gN/m3]
    !Out:
    real(r8)            :: Fmax !Maximum flux from NH4 to SAP in cases with too limited N (SAP immobilization)
    
    Fmax = k*nh4 
  end function calc_Fmax 


  
end module paramMod
