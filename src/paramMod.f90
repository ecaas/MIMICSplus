module paramMod
use shr_kind_mod   , only : r8 => shr_kind_r8

implicit none

!Define variables
real(kind=r8)                                :: TSOI =5.18                      ![degC]
real(kind=r8)                                :: GEP                             ![gC/(m2 h)] Gross ecosystem productivity
real(kind=r8),parameter                      :: fCLAY  = 0.30                   ![-] fraction of clay in soil
real(kind=r8),dimension(3)                   :: k_mycsom                        ![1/h] decay constants, MYC to SOM pools


!For calculating the Km parameter in Michaelis Menten kinetics (expressions based on mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and https://github.com/wwieder/MIMICS)
integer, parameter                           :: MM_eqs  = 6                     !Number of Michaelis-Menten parameters
real(kind=r8),dimension(MM_eqs),parameter    :: Kslope  = (/0.034, 0.034, 0.034, 0.034, 0.034, 0.034/) !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Vslope  = (/0.063, 0.063, 0.063, 0.063, 0.063, 0.063/) !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Kint    = (/3.19, 3.19, 3.19, 3.19, 3.19, 3.19/)       !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Vint    = (/5.47, 5.47, 5.47, 5.47, 5.47, 5.47/)       !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),parameter                      :: a_k     = 1e3 !Tuning parameter g/m3
real(kind=r8),parameter                      :: a_v     = 8e-6 !Tuning parameter
real(kind=r8),parameter                      :: pscalar = 1.0/(2*exp(-2.0*dsqrt(fCLAY)))
real(kind=r8),dimension(MM_eqs)              :: Kmod    = (/1d0, 1d0, 0.5d0*pscalar, 0.5d0, 0.5d0, 0.5d0*pscalar/)!LITm, LITs, SOMa entering SAPb, sapf
real(kind=r8),dimension(MM_eqs)              :: Vmod    = (/4.0,  4.0, 4.0, 3.0, 3.0, 3.0/)            !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),parameter, dimension(2)        :: KO      =  4                    ![-]Increases Km (the half saturation constant for oxidation of chemically protected SOM, SOM_c) from mimics
real(kind=r8),dimension(MM_eqs)              :: Km                              ![mgC/cm3]*10e3=[gC/m3]
real(kind=r8),dimension(MM_eqs)              :: Vmax                            ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics.

!Pools: NOTE: This needs to be updated if pools are added to/removed from the system.
integer, parameter                           :: no_of_litter_pools = 2         !Metabolic and structural
integer, parameter                           :: no_of_sap_pools = 2            !SAP bacteria and SAP fungi
integer, parameter                           :: no_of_myc_pools = 3            !Mycorrhiza: Ecto, ericoid & arbuscular
integer, parameter                           :: no_of_som_pools = 3            !Physically protected, chemically protected, available carbon
integer, parameter                           :: pool_types = no_of_litter_pools + no_of_myc_pools + &
                                                no_of_sap_pools + no_of_som_pools
integer, parameter                           :: pool_types_N = pool_types+1

!For calculating turnover from SAP to SOM (expressions from mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and  https://github.com/wwieder/MIMICS)
real(r8), dimension(no_of_sap_pools)    :: fPHYS,fCHEM,fAVAIL              ![-]
real(r8), dimension(no_of_sap_pools)    :: tau                             ![1/h]
real(r8)                                :: fMET =0.4                       ![-] Fraction determining distribution of total litter production between LITm and LITs NOTE: Needs revision

real(r8), dimension(no_of_som_pools), parameter    :: fEcMSOM = (/0.4,0.4,0.2/) !somp,soma,somc. Fraction of flux from EcM to different SOM pools NOTE: assumed
real(r8), dimension(no_of_som_pools), parameter    :: fErMSOM = (/0.3,0.4,0.3/)
real(r8), dimension(no_of_som_pools), parameter    :: fAMSOM = (/0.3,0.3,0.4/)
real(r8)                                :: desorb = 3e-4*exp(-4*(sqrt(fclay)))![1/h]From Mimics, used for the transport from physically protected SOM to available SOM pool

!Depth & vertical transport
real(r8),parameter                   :: soil_depth = 1.52            ![m] used if isVertical is False (sum(delta_z))
real(r8),dimension(10),parameter     :: node_z = (/0.01,0.04,0.09,0.16,0.26,0.40,0.587,0.80,1.06,1.36/) ![m] Depth of center in each soil layer. Same as the first layers of default CLM5 with vertical resolution.
real(r8),dimension(10)               :: delta_z = (/0.02, 0.04, 0.06, 0.08,0.12,0.16,0.20,0.24,0.28,0.32/)![m] Thickness of each soil of the top layers in default clm5.
real(r8),parameter                   :: D = 1.14e-8![m2/h] Diffusivity. Based on Koven et al 2013, 1cm2/yr = 10e-4/(24*365)


real(r8), dimension(pool_types), parameter   :: CN_ratio = (/15,15,5,8,20,20,20,11,8,11/) !Fungi/bacteria: Tang, Riley, Maggi 2019 as in Mouginot et al. 2014
                                                                                          !NOTE: Wallander/Rousk may have data more suited for Boreal/Arctic conditions
                                                                                          !EcM: From Baskaran et al as in Wallander et al 2004
                                                                                          !SOM: From CLM documentation, table 21.3 (Mendeley version)
                                                                                          !LITm: MIMICS-CN manuscript
                                                                                          !LITs, ErM, AM: Guesses!

!From Baskaran et al 2016
real(r8), parameter :: hr_pr_yr = 365*24        !For conversion
real(r8), parameter :: my_sap = 1/hr_pr_yr      ![1/hr]mortality rate sap
real(r8), parameter :: my_myc = 1/hr_pr_yr      ![1/hr] mortality rate myc
real(r8), parameter :: my_root = 0.15/hr_pr_yr  ![1/hr]mortality rate plant root
real(r8), parameter :: my_shoot = 0.15/hr_pr_yr ![1/hr]  mortality rate plant shoot
real(r8), parameter :: gamma_rs = 0.3           !Plant root:shoot ratio                      TODO: Vary with plant type/pft/myc?
real(r8)            :: delta=0.15               !Fraction of plant C allocated to mycorrhiza TODO: determine based on myc type?
real(r8), parameter :: a = 80/hr_pr_yr          ![gC(gN)⁻¹yr⁻¹], max plant N productivity
real(r8), parameter :: b = 0.09/hr_pr_yr        !Shading factor of plant productivity
real(r8), parameter :: Km_plant = 0.6           ![gNm-2] Half saturation constant of plant uptake of inorganic N (called S_p in article)
real(r8), parameter :: Km_myc = 0.08            ![gNm-2] Half saturation constant of mycorrhizal uptake of inorganic N (called S_m in article)
real(r8), parameter :: V_max_plant = 1.8/hr_pr_yr![g g-1 hr-1] Max plant root uptake of inorganic N (called K_pn in article)
real(r8), parameter :: V_max_myc = 1.8/hr_pr_yr  ![g g-1 hr-1] Max mycorrhizal uptake of inorganic N (called K_mn in article)
real(r8)  :: Leaching_rate = 3/hr_pr_yr          ![hr-1] Leaching rate
real(r8)  :: Deposition_rate =3/hr_pr_yr         ![gNm-2hr-1] Deposition rate  NOTE: varied from 0.3-3 in article
real(r8), parameter :: e_s = 0.25                !Growth efficiency of saprotrophs        TODO: Compare these to the efficiencies from Mimics
real(r8), parameter :: e_m = 0.25                !Growth efficiency of mycorrhiza NOTE: If efficiency is too high, SAPbIN will become negative bc.  e_s*U_sb/CN_ratio(3) will be too large. Problem??
!Decomposition rates:
real(r8), parameter :: K_SH = 0.006/hr_pr_yr ![m2gC-1hr-1] Saprotrophic decay rate constant for hydrolizable store. TODO: review these
                        !NOTE: K_SH is not used in fluxMod!!
real(r8), parameter :: K_MO = 0.0003/hr_pr_yr ![m2gC-1hr-1] Mycorrhizal decay rate constant for oxidizable store     NOTE: vary from 0.0003 to 0.003 in article


!Moisture dependence (based on function used for MIMICS in the CASA-CNP testbed)
real(r8), parameter                          :: P = 44.247 !normalization of moisture function
real(r8)                                     :: theta_l, theta_f, theta_sat ! liquid, frozen and saturation water content (Needs to come from some kind of forcing)
real(r8)                                     :: gas_diffusion
integer, parameter, dimension(12)            :: days_in_month =(/31,28,31,30,31,30,31,31,30,31,30,31/)
integer                                      :: current_month, previous_month

real(r8)                                     :: Loss_termN, Loss_termC, Loss_termNP, Loss_termCP, Plant_gainN,&
                                                Plant_GainC, Plant_lossN, Plant_lossC, a_NPlant, a_CPlant

!Fluxes etc:
real(r8) :: C_LITmSAPb, C_LITsSAPb, C_EcMSOMp, C_EcMSOMa, C_EcMSOMc, C_ErMSOMp, C_ErMSOMa, C_ErMSOMc, C_AMSOMp, &
C_LITmSAPf, C_LITsSAPf, C_AMSOMa, C_AMSOMc, C_SOMaSAPb,C_SOMaSAPf, C_SOMpSOMa, C_SOMcSOMa, &
C_SAPbSOMa, C_SAPbSOMp, C_SAPbSOMc,C_SAPfSOMa, C_SAPfSOMp, C_SAPfSOMc, &
N_LITmSAPb, N_LITsSAPb, N_EcMSOMp, N_EcMSOMa, N_EcMSOMc, N_ErMSOMp, N_ErMSOMa, N_ErMSOMc, N_AMSOMp, N_AMSOMa,&
N_AMSOMc, N_SOMaSAPb,N_SOMaSAPf, N_SOMpSOMa, N_SOMcSOMa, N_LITmSAPf, N_LITsSAPf, N_SOMaEcM, N_SOMaErM,N_SOMaAM,&
N_PlantLITs, N_PlantLITm, N_INPlant, N_INEcM, N_INErM, N_INAM, N_EcMPlant, N_ErMPlant, N_AMPlant, &
N_SAPbSOMa, N_SAPbSOMp, N_SAPbSOMc,N_SAPfSOMa, N_SAPfSOMp, N_SAPfSOMc, N_SAPfIN, N_SAPbIN,&
C_growth_rate, C_PlantEcM, C_PlantErM, C_PlantAM, C_PlantLITm, C_PlantLITs, Decomp_ecm, &
Decomp_erm, Decomp_am, Leaching, Deposition, C_PR, C_PS,N_PR, N_PS, Total_plant_mortality,f, U_sb, U_sf, CPlant,&
NPlant, P_N, Plant_CN, CPlant_tstep, NPlant_tstep, growth_rate_sum

character (len=*),parameter                  :: clm_data_file = &
'/home/ecaas/clm/cruncep_iso_hist/Dovre/clm50_clm50d001_1deg_CRUNCEPV7_iso_hist.clm2.h0.SOILLIQ_SOILICE_TSOI_W_SCALAR.185001-201412.nc_Dovre2014.nc'
!For writing to file:
character (len=*),parameter                  :: output_path = '/home/ecaas/decomposition_results/vertical/'
integer                                      :: ios = 0 !Changes if something goes wrong when opening a file
character (len=4), dimension(pool_types)     :: variables = &
(/  "LITm", "LITs", "SAPb","SAPf", "EcM ", "ErM ", "AM  ", "SOMp", "SOMa", "SOMc" /)
character (len=*), dimension(pool_types+1), parameter:: N_variables = &
(/  "N_LITm", "N_LITs", "N_SAPb","N_SAPf", "N_EcM ", "N_ErM ", "N_AM  ", "N_SOMp", "N_SOMa", "N_SOMc", "N_Inor"/)
character (len=10), dimension(pool_types):: change_variables = &
(/  "changeLITm", "changeLITs", "changeSAPb","changeSAPf", "changeEcM ", "changeErM ",&
    "changeAM  ", "changeSOMp", "changeSOMa", "changeSOMc" /)
character (len=*), dimension(pool_types),parameter:: an_variables = &
&(/  "anLITm", "anLITs", "anSAPb","anSAPf", "anEcM ", "anErM ", "anAM  ", "anSOMp", "anSOMa", "anSOMc" /)
character (len=*), dimension(*), parameter ::  C_name_fluxes = &
[character(len=11) ::"LITmSAPb","LITmSAPf","LITsSAPb","LITsSAPf", "SAPbSOMp","SAPfSOMp", "SAPbSOMa","SAPfSOMa", "SAPbSOMc","SAPfSOMc", &
  "EcMSAPb ", "EcMSAPf ","ErMSAPb ","ErMSAPf ", "AMSAPb  ","AMSAPf  ","EcMSOMp ", "EcMSOMa ","EcMSOMc ", "ErMSOMp ",&
  "ErMSOMa ","ErMSOMc ","AMSOMp  ","AMSOMa  ","AMSOMc  ","SOMaSAPb","SOMaSAPf","SOMpSOMa","SOMcSOMa","PlantLITm" &
  ,"PlantLITs","PlantEcM","PlantErM","PlantAM"]

character (len=*), dimension(*), parameter ::  N_name_fluxes = &
[character(len=11) ::"LITmSAPb","LITmSAPf","LITsSAPb","LITsSAPf", "SAPbSOMp","SAPfSOMp", "SAPbSOMa","SAPfSOMa", "SAPbSOMc","SAPfSOMc", &
  "EcMSAPb ", "EcMSAPf ","ErMSAPb ","ErMSAPf ", "AMSAPb  ","AMSAPf  ","EcMSOMp ", "EcMSOMa ","EcMSOMc ", "ErMSOMp ",&
  "ErMSOMa ","ErMSOMc ","AMSOMp  ","AMSOMa  ","AMSOMc  ","SOMaSAPb","SOMaSAPf","SOMpSOMa","SOMcSOMa","PlantLITm" &
  ,"PlantLITs","EcMPlant","ErMPlant","AMPlant", "Deposition", "Leaching", "INEcM", "INErM","INAM", &
  "SAPbIN", "SAPfIN"]

end module paramMod
