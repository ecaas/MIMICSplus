module paramMod
use shr_kind_mod   , only : r8 => shr_kind_r8
!use SoilDecompType, only : carbon_pool_type

implicit none
!Define variables


real(kind=r8)                                :: tsoi                            ![degC]
real(kind=r8)                                :: GEP                             ![gC/(m2 h)] Gross ecosystem productivity
real(kind=r8),parameter                      :: fCLAY  = 0.15                   ![-] fraction of clay in soil

real(r8)                                     :: f_som1=0.05, f_som2=0.05
real(kind=r8),dimension(3)                   :: k_mycsap                         ![1/h](EcM, ErM, AM) decay constants, MYC to SAP pool
real(kind=r8),dimension(3)                   :: k_mycsom                         ! [1/h] decay constants, MYC to SOM pools
real(kind=r8), parameter                     :: MYC_SAPb=0.4, MYC_SAPf=1-MYC_SAPb![-]Fraction of the flux from mycorrhizal pools to SAPf pool. The rest is going to SAPb

!For calculating the Km parameter in Michaelis Menten kinetics (expressions based on mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and https://github.com/wwieder/MIMICS)
integer, parameter                           :: MM_eqs  = 6                     !Number of Michaelis-Menten parameters
real(kind=r8),dimension(MM_eqs),parameter    :: Kslope  = (/0.034, 0.034, 0.034, 0.034, 0.034, 0.034/)!LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Vslope  = (/0.063, 0.063, 0.063, 0.063, 0.063, 0.063/) !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Kint    = (/3.19, 3.19, 3.19, 3.19, 3.19, 3.19/)       !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),dimension(MM_eqs),parameter    :: Vint    = (/5.47, 5.47, 5.47, 5.47, 5.47, 5.47/)       !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8), parameter                     :: a_k     = 1e3 !Tuning parameter g/m3
real(kind=r8), parameter                     :: a_v     = 8e-6 !Tuning parameter
real(r8),parameter                           :: pscalar = 1.0/(2*exp(-2.0*dsqrt(fclay)))
real(kind=r8),dimension(MM_eqs)              :: Kmod    = (/1d0, 1d0, 0.5d0*pscalar, 0.5d0, 0.5d0, 0.5d0*pscalar/)!LITm, LITs, SOMa entering SAPb, sapf
real(kind=r8),dimension(MM_eqs)              :: Vmod    = (/4.0,  4.0, 4.0, 3.0, 3.0, 3.0/)            !LITm, LITs, SOMa entering SAPb, LITm, LITs, SOMa entering SAPf
real(kind=r8),parameter, dimension(2)        :: KO      =  4                    ![-]Increases Km (the half saturation constant for oxidation of chemically protected SOM, SOM_c) from mimics
real(kind=r8),dimension(MM_eqs)              :: Km                              ![mgC/cm3]*10e3=[gC/m3]
real(kind=r8),dimension(MM_eqs)              :: Vmax                            ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics.

real(kind=r8),dimension(MM_eqs),save         :: MGE                             ![mg/mg] Microbial growth efficiency/Carbon Use efficiency. The fraction of the flux from litter pools that is used in microbial processes.
                                                                                !The rest is lost in respiration. Values should be determined carefully. pH/soil quality/N availability may be important
                                                                                !Maybe also vary with time/depth..? See DOI:10.1016/j.soilbio.2018.09.036 and DOI: 10.1016/J.SOILBIO.2019.03.008

!Pools: NOTE: This needs to be updated if pools are added to/removed from the system.
integer, parameter                           :: num_soilc = 1                  !Number of soil columns modeled. Only use one column here, but TODO code should be working for several columns as well (add loop)
integer, parameter                           :: no_of_litter_pools = 2         !Metabolic and structural
integer, parameter                           :: no_of_sap_pools = 2            !SAP bacteria and SAP fungi!
integer, parameter                           :: no_of_myc_pools = 3            !Mycorrhiza: Ecto, ericoid & arbuscular
integer, parameter                           :: no_of_som_pools = 3            !Physically protected, chemically protected, available carbon
integer, parameter                           :: pool_types = no_of_litter_pools + no_of_myc_pools + &
                                                no_of_sap_pools + no_of_som_pools

!For calculating turnover from SAP to SOM (expressions from mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and  https://github.com/wwieder/MIMICS)
real(r8)                                :: fMET                            ![-]
real(r8), dimension(no_of_sap_pools)    :: fPHYS,fCHEM,fAVAIL              ![-]
real(r8), dimension(no_of_sap_pools)    :: tau                             ![1/h]
real(r8), dimension(no_of_som_pools), parameter    :: fEcMSOM = (/0.4,0.4,0.2/) !somp,soma,somc. Fraction of flux from EcM to different SOM pools
real(r8), dimension(no_of_som_pools), parameter    :: fErMSOM = (/0.2,0.4,0.3/)
real(r8), dimension(no_of_som_pools), parameter    :: fAMSOM = (/0.3,0.3,0.4/)
real(r8)                                :: desorb = 3e-4 * exp(-4*(sqrt(fclay)))!1.5e-5 * exp(-1.5*(fclay)) ![1/h]From Mimics, used for the transport from physically protected SOM to available SOM pool

!Depth & vertical transport
real(r8),parameter                           :: depth = 0.56                    ![m] used if isVertical is False
integer,parameter                            :: nlevdecomp = 7                  ! number of vertical layers
real(r8),dimension(nlevdecomp),parameter     :: node_z = (/0.01,0.04,0.09,0.16,0.26,0.40,0.587/)!,0.80,1.06,1.36/) ![m] Depth of center in each soil layer. Same as the first 4 layers of default CLM5 with vertical resolution.
real(r8),dimension(nlevdecomp),parameter     :: delta_z = (/0.02, 0.04, 0.06, 0.08,0.12,0.16,0.20/)!,0.24,0.28,0.32/)![m] Thickness of each soil of the four top layers in default clm5.
real(r8),parameter                           :: D = 1.14e-8![m2/h] Diffusivity. Based on Koven et al 2013, 1cm2/yr = 10e-4/(24*365)
real(r8)                                     :: diffusive_source,diffusive_sink !For use in vertical_diffusion subroutine. Currently using alt_vertical diffusion subroutine
real(r8),dimension(nlevdecomp,pool_types)    :: net_diffusion
!real(r8), dimension(nlevdecomp), parameter   :: f_depth = 1!(/1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4/)
real(r8), parameter                          :: min_pool_value = 0.01

real(r8), dimension(nlevdecomp)          :: TSOIL
real(r8), dimension(nlevdecomp)          :: SOILLIQ
real(r8), dimension(nlevdecomp)          :: SOILICE
real(r8), dimension(nlevdecomp)          :: WATSAT
real(r8), dimension(nlevdecomp)          :: r_moist

!Moisture dependence (based on function used for MIMICS in the CASA-CNP testbed)
real(r8), parameter                          :: P = 44.247 !normalization of moisture function
real(r8)                                     :: theta_l, theta_f, theta_sat ! liquid, frozen and saturation water content (Needs to come from some kind of forcing)
real(r8)                                     :: gas_diffusion
integer, parameter, dimension(12)           :: days_in_month =(/31,28,31,30,31,30,31,31,30,31,30,31/)
integer                                      :: current_month, previous_month
!Fluxes between pools:
real(r8) :: LITmSAPb, LITsSAPb, EcMSAPb, ErMSAPb, AMSAPb, EcMSOMp, EcMSOMa, EcMSOMc, ErMSOMp, ErMSOMa, ErMSOMc, AMSOMp, AMSOMa, AMSOMc, SOMaSAPb,SOMaSAPf, SOMpSOMa, SOMcSOMa
real(r8) :: LITmSAPf, LITsSAPf, EcMSAPf, ErMSAPf, AMSAPf
real(r8) :: SAPbSOMa, SAPbSOMp, SAPbSOMc,SAPfSOMa, SAPfSOMp, SAPfSOMc
!End fluxes between pools
character (len=*),parameter                  :: clm_data_file = '/home/ecaas/clm/cruncep_iso_hist/Dovre/clm50_clm50d001_1deg_CRUNCEPV7_iso_hist.clm2.h0.SOILLIQ_SOILICE_TSOI.Dovre2014.nc'
!For writing to file:
character (len=*),parameter                  :: output_path = '/home/ecaas/decomposition_results/'
integer                                      :: ios = 0 !Changes if something goes wrong when opening a file
character (len=4), dimension(pool_types):: variables = (/  "LITm", "LITs", "SAPb","SAPf", "EcM ", "ErM ", "AM  ", "SOMp", "SOMa", "SOMc" /)
character (len=10), dimension(pool_types):: change_variables = (/  "changeLITm", "changeLITs", "changeSAPb","changeSAPf", "changeEcM ", "changeErM ", "changeAM  ", "changeSOMp", "changeSOMa", "changeSOMc" /)
character (len=*), dimension(pool_types),parameter:: an_variables = (/  "anLITm", "anLITs", "anSAPb","anSAPf", "anEcM ", "anErM ", "anAM  ", "anSOMp", "anSOMa", "anSOMc" /)
character (len=*), dimension(29), parameter ::  name_fluxes = (/"LITmSAPb","LITmSAPf","LITsSAPb","LITsSAPf", "SAPbSOMp","SAPfSOMp", "SAPbSOMa","SAPfSOMa", "SAPbSOMc","SAPfSOMc", &
      "EcMSAPb ", "EcMSAPf ","ErMSAPb ","ErMSAPf ", "AMSAPb  ","AMSAPf  ","EcMSOMp ", "EcMSOMa ","EcMSOMc ", "ErMSOMp ","ErMSOMa ","ErMSOMc ","AMSOMp  ","AMSOMa  " &
      ,"AMSOMc  ","SOMaSAPb","SOMaSAPf","SOMpSOMa","SOMcSOMa" /)

end module paramMod
