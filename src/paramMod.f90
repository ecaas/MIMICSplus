module paramMod
use shr_kind_mod   , only : r8 => shr_kind_r8

implicit none

!Define variables

real(r8), parameter :: sec_pr_hr = 60*60        !For conversion
real(r8), parameter :: hr_pr_yr = 365*24        !For conversion
real(kind=r8)                                :: T_soil_const =5.18              ![degC]
real(kind=r8)                                :: GEP                             ![gC/(m2 h)] Gross ecosystem productivity
real(kind=r8)                                :: fCLAY                          ![-] fraction of clay in soil
real(kind=r8),dimension(3)                   :: k_mycsom                        ![1/h] decay constants, MYC to SOM pools

!For calculating the Km parameter in Michaelis Menten kinetics (expressions based on mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and https://github.com/wwieder/MIMICS)
integer, parameter                           :: MM_eqs  = 6                     !Number of Michaelis-Menten parameters
real(kind=r8),dimension(MM_eqs),parameter    :: Kslope  = (/0.017, 0.027, 0.017, 0.017, 0.027, 0.017/) !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Vslope  = (/0.063, 0.063, 0.063, 0.063, 0.063, 0.063/) !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Kint    = 3.19      !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Vint    = 5.47      !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),parameter                      :: a_k     = 1e4 !Tuning parameter g/m3 (10 mg/cm3 from german et al 2012)
real(kind=r8),parameter                      :: a_v     = 8e-6 !Tuning parameter
real(kind=r8)                                :: pscalar != 1.0/(2*exp(-2.0*dsqrt(fCLAY))) !Value range:  0.5-3.7
real(kind=r8),dimension(MM_eqs)              :: Kmod    !LITm, LITs, SOMa entering SAPb, sapf
real(kind=r8),dimension(MM_eqs)              :: Vmod    = (/10.0, 2.0, 10.0, 3.0,3.0, 2.0/)            !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
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
real(r8),parameter                      :: fMET =0.6                       ![-] Fraction determining distribution of total litter production between LITm and LITs NOTE: Needs revision
real(r8), dimension(no_of_sap_pools)    :: fPHYS,fCHEM,fAVAIL              ![-]
real(r8), dimension(no_of_sap_pools)    :: tau  ![1/h]

real(r8), dimension(no_of_som_pools), parameter    :: fEcMSOM = (/0.4,0.4,0.2/) !somp,soma,somc. Fraction of flux from EcM to different SOM pools NOTE: assumed
real(r8), dimension(no_of_som_pools), parameter    :: fErMSOM = (/0.3,0.4,0.3/)
real(r8), dimension(no_of_som_pools), parameter    :: fAMSOM = (/0.3,0.3,0.4/)
real(r8)                                :: desorb ![1/h]From Mimics, used for the transport from physically protected SOM to available SOM pool

!Depth & vertical transport
real(r8)                             :: soil_depth           ![m] used if isVertical is False (sum(delta_z))
real(r8),dimension(10),parameter     :: node_z =  (/0.01,0.04,0.09,0.16,0.26,0.40,0.587,0.80,1.06,1.36/)!(/0.076,0.228, 0.380,0.532, 0.684,0.836,0.988,1.140,1.292,1.444/)!![m] Depth of center in each soil layer. Same as the first layers of default CLM5 with vertical resolution.
real(r8),dimension(10),parameter     :: delta_z = (/0.02, 0.04, 0.06, 0.08,0.12,0.16,0.20,0.24,0.28,0.32/)!0.152![m] Thickness of each soil of the top layers in default clm5.
real(r8),parameter                   :: D_carbon = 1.14e-8![m2/h] Diffusivity. Based on Koven et al 2013, 1cm2/yr = 10e-4/(24*365)
real(r8),parameter                   :: D_nitrogen = 1.14e-8![m2/h] Diffusivity. Based on Koven et al 2013, 1cm2/yr = 10e-4/(24*365)


real(r8),dimension(:),allocatable    :: CUE_bacteria_vr
real(r8),dimension(:),allocatable    :: CUE_fungi_vr

real(r8),parameter                   :: CUE_0=0.50
real(r8),parameter                   :: CUE_slope=0.0!-0.016 !From German et al 2012
real(r8)                             :: CUEmod_bacteria
real(r8)                             :: CUEmod_fungi

real(r8), parameter                  :: f_met_to_som=0.05 ! fraction of metabolic litter flux that goes directly to SOM pools
real(r8), parameter                  :: max_mining = 0.3
real(r8), parameter                  :: enzyme_pct = 0.1

real(r8), dimension(pool_types), parameter   :: CN_ratio = (/15,15,5,8,20,20,20,11,8,11/) !Fungi/bacteria: Tang, Riley, Maggi 2019 as in Mouginot et al. 2014
                                                                                          !NOTE: Wallander/Rousk may have data more suited for Boreal/Arctic conditions
                                                                                          !EcM: From Baskaran et al as in Wallander et al 2004
                                                                                          !SOM: From CLM documentation, table 21.3 (Mendeley version)
                                                                                          !LITm: MIMICS-CN manuscript
                                                                                          !LITs, ErM, AM: Guesses!

!From Baskaran et al 2016
real(r8), parameter :: Km_myc = 0.08            ![gNm-2] Half saturation constant of mycorrhizal uptake of inorganic N (called S_m in article) 
real(r8), parameter :: V_max_myc = 1.8/hr_pr_yr  ![g g-1 hr-1] Max mycorrhizal uptake of inorganic N (called K_mn in article) 
real(r8), parameter :: e_m = 0.25             !Growth efficiency of mycorrhiza 


!Decomposition rates:
real(r8), parameter :: K_MO = 0.003/hr_pr_yr ![m2gC-1hr-1] Mycorrhizal decay rate constant for oxidizable store     NOTE: vary from 0.0003 to 0.003 in article

!Moisture dependence (based on function used for MIMICS in the CASA-CNP testbed)
real(r8), parameter                          :: P = 44.247 !normalization of moisture function
real(r8)                                     :: gas_diffusion
integer, parameter, dimension(12)            :: days_in_month =(/31,28,31,30,31,30,31,31,30,31,30,31/)

real(r8),dimension(:),allocatable  :: ndep_prof
real(r8),dimension(:),allocatable  :: leaf_prof
real(r8),dimension(:),allocatable  :: froot_prof


!Fluxes etc:
real(r8) :: C_LITmSAPb, C_LITsSAPb, C_EcMSOMp, C_EcMSOMa, C_EcMSOMc, C_ErMSOMp, C_ErMSOMa, C_ErMSOMc, C_AMSOMp, &
C_LITmSAPf, C_LITsSAPf, C_AMSOMa, C_AMSOMc, C_SOMaSAPb,C_SOMaSAPf, C_SOMpSOMa, C_SOMcSOMa, &
C_SAPbSOMa, C_SAPbSOMp, C_SAPbSOMc,C_SAPfSOMa, C_SAPfSOMp, C_SAPfSOMc, C_PlantSOMc,C_PlantSOMp,C_PlantSOMa, &
N_LITmSAPb, N_LITsSAPb, N_EcMSOMp, N_EcMSOMa, N_EcMSOMc, N_ErMSOMp, N_ErMSOMa, N_ErMSOMc, N_AMSOMp, N_AMSOMa,&
N_AMSOMc, N_SOMaSAPb,N_SOMaSAPf, N_SOMpSOMa, N_SOMcSOMa, N_LITmSAPf, N_LITsSAPf, N_SOMaEcM, N_SOMaErM,N_SOMaAM,&
N_PlantLITs, N_PlantLITm, N_INPlant, N_INEcM, N_INErM, N_INAM, N_EcMPlant, N_ErMPlant, N_AMPlant, &
N_SAPbSOMa, N_SAPbSOMp, N_SAPbSOMc,N_SAPfSOMa, N_SAPfSOMp, N_SAPfSOMc, N_SAPfIN, N_SAPbIN,&
N_SOMcEcM,N_SOMpEcM, &
C_PlantEcM, C_PlantErM, C_PlantAM, C_PlantLITm, C_PlantLITs, C_EcMdecompSOMp,C_EcMdecompSOMc, &
 Leaching, Deposition,f, U_sb, U_sf

!For writing to file:
character (len=*),parameter                  :: output_path = '/home/ecaas/decomposition_results/sites/'
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
  "EcMSOMp ", "EcMSOMa ","EcMSOMc ", "ErMSOMp ",&
  "ErMSOMa ","ErMSOMc ","AMSOMp  ","AMSOMa  ","AMSOMc  ","SOMaSAPb","SOMaSAPf","SOMpSOMa","SOMcSOMa","PlantLITm" &
  ,"PlantLITs","PlantEcM","PlantErM","PlantAM","PlantSOMc  ","PlantSOMp  ","PlantSOMa  ", &
  "EcMdecoSOMp","EcMdecoSOMc"]

character (len=*), dimension(*), parameter ::  N_name_fluxes = &
[character(len=11) ::"LITmSAPb","LITmSAPf","LITsSAPb","LITsSAPf", "SAPbSOMp","SAPfSOMp", "SAPbSOMa","SAPfSOMa", "SAPbSOMc","SAPfSOMc" &
  ,"EcMSOMp ", "EcMSOMa ","EcMSOMc ", "ErMSOMp ",&
  "ErMSOMa ","ErMSOMc ","AMSOMp  ","AMSOMa  ","AMSOMc  ","SOMaSAPb","SOMaSAPf","SOMaEcM","SOMpSOMa","SOMcSOMa","PlantLITm" &
  ,"PlantLITs","EcMPlant","ErMPlant","AMPlant", "Deposition", "Leaching", "INEcM", "INErM","INAM", &
  "SAPbIN", "SAPfIN","SOMpEcM","SOMcEcM"]
  
character (len=*), dimension(38),parameter :: site_names = &
[character(len=19) :: &
 'NR31585_Flekkefjord','NR31581_Lyngdal    ','NR32361_Lyngdal    ','NR31682_Tysvar     ','32288_Sortland     ',&
 'NR32249_Vik        ','NR32182_Stryn      ','NR31881_Sande      ','NR31578_Kvinesdal  ',&
 '31463_Hurdal       ','31464_Hurdal       ','32087_Dovre        ','32379_Hemne        ','32441_Sel          ','32258_Maaselv      ',&
 '32032_VestreToten  ','31984_Namdalseid   ','31976_Namdalseid   ','31461_Nittedal     ','31513_Nes          ',&
 '31539_Modum        ','31652_Bygland      ','31767_Kongsvinger  ','31780_Vaaler       ','31941_Roeyrvik     ','31997_Verdal       ',&
 '32088_Lesja        ','32103_Halden       ','32139_Rennebu      ','32374_Saltdal      ','32404_Vinje        ','32409_Vang         ',&
 '32438_Porsanger    ','31519_Nissedal     ','31650_Valle        ','31714_Flaa         ','32085_Skjaak       ','NR31906_Voss       ' ]
! '32124_Engerdal     ', '32246_SoerVaranger '  ,'NR31927_Os         ',,'NR31590_Farsund    ', 'NR31908_Ullensvang ','NR32485_Stord      '
!'NR31579_Kvinesdal  ','NR31577_Kvinesdal  ',,'NR31677_Suldal     '
end module paramMod
