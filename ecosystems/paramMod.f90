module paramMod
use shr_kind_mod   , only : r8 => shr_kind_r8
!use SoilDecompType, only : carbon_pool_type

implicit none
!Define variables


real(kind=r8)                                :: tsoi                      ![degC]
real(kind=r8),parameter                      :: fCLAY  = 0.12                  ![-] fraction of clay in soil
real(kind=r8),dimension(6),save              ::  MGE![mg/mg] Microbial growth efficiency/Carbon Use efficiency. The fraction of the flux from litter pools that is used in microbial processes.
                                                                                !The rest is lost in respiration.
                                                                                !/Litm-sapr, litm->sapk, lits->sapr, lits->sapk/. These should be determined carefully. pH/soil quality/N availability may be important
                                                                                !Maybe also vary with time/depth..? See DOI:10.1016/j.soilbio.2018.09.036 and DOI: 10.1016/J.SOILBIO.2019.03.008
real(r8)                                     :: f_som1=0.05, f_som2=0.05
real(r8)                                     :: f_myc_levels=1, f_lit_1=1, f_lit_234=1
real(kind=r8),dimension(3)                   :: k                              ![1/day](EcM, ErM, AM) decay constants, MYC to SAP pools
real(kind=r8),dimension(3)                   :: k2                             ! [1/day] decay constants, MYC to SOM pools
!real(kind=r8), parameter                     :: Myc_SAPr=0.6, Myc_SAPk=1-Myc_SAPr![-]Fraction of the flux from mycorrhizal pools to SAPk pool. The rest is going to SAPr

!For calculating the Km parameter in Michaelis Menten kinetics (expressions based on mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and https://github.com/wwieder/MIMICS)
!TODO: Varying Kslope?
integer, parameter                           :: MM_eqs  = 3                  !Number of Michaelis-Menten parameters
real(kind=r8),dimension(MM_eqs),parameter    :: Kslope  = (/0.017, 0.027, 0.017/) !LITm, LITs, SOMa entering SAPr, LITm, LITs, SOMa entering SAPk
real(r8),parameter                           :: pscalar = 1.0/(2*exp(-2.0*dsqrt(fclay)))
real(kind=r8),dimension(MM_eqs)              :: Kmod!LITm, LITs, SOMa entering SAPr, LITm, LITs, SOMa entering sapk
real(kind=r8),dimension(MM_eqs)              :: Km   ![mgC/cm3]*10e3=[gC/m3]
real(kind=r8),parameter                      :: KO      =  10  !4                      ![-]Increases Km (the half saturation constant for oxidation of chemically protected SOM, SOM_c) from mimics
real(kind=r8),dimension(MM_eqs)              :: Vmod     !LITm, LITs, SOMa entering SAPr, LITm, LITs, SOMa entering sapk
real(kind=r8),dimension(MM_eqs)              :: Vmax ![mgC/((mgSAP)h)] For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?
!end MM parameters


!For calculating turnover from SAP to SOM (expressions from mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and  https://github.com/wwieder/MIMICS)
real(kind=r8)                                :: fMET                           ![-]
real(kind=r8)                                :: fPHYS,fCHEM,fAVAIL             ![-]
real(kind=r8)                                :: tau                            ![1/h]
real(kind=r8)                                :: desorb = 3e-4 * exp(-4*(sqrt(fclay)))!1.5e-5 * exp(-1.5*(fclay)) ![1/h]From Mimics, used for the transport from physically protected SOM to available SOM pool

!Pools: NOTE: This needs to be updated if pools are added to/removed from the system.
integer, parameter                           :: num_soilc = 1                  !Number of soil columns modeled. Only use one column here, but TODO code should be working for several columns as well (add loop)
integer, parameter                           :: no_of_litter_pools = 2         !Metabolic and structural
integer, parameter                           :: no_of_sap_pools = 1            !NOTE changed from two pools to one pool!
integer, parameter                           :: no_of_myc_pools = 3            !Mycorrhiza: Ecto, ericoid & arbuscular
integer, parameter                           :: no_of_som_pools = 3            !Physically protected, chemically protected, available carbon
integer, parameter                           :: pool_types = no_of_litter_pools + no_of_myc_pools + &
                                                no_of_sap_pools + no_of_som_pools

!Depth & vertical transport
real(r8),parameter                           :: depth = 0.56                   ![m] used if isVertical is False
real(r8),dimension(4)                        :: node_z = (/0.05,0.2,0.5,1.10/) ![m] Depth of center in each soil layer. Used if isVertical is True
real(r8),dimension(4)                        :: delta_z = (/0.1,0.2,0.4,0.8/)  ![m] Thickness of each soil layer NOTE: delta_z and node_z is used in the alt_vertical_diffusion routine
!real(r8),dimension(nlevdecomp)               :: node_z = (/0.01,0.04,0.09,0.16/) ![m] Depth of center in each soil layer. Same as the first 4 layers of default CLM5 with vertical resolution.
!real(r8),dimension(nlevdecomp)               :: delta_z = (/0.02, 0.04, 0.06, 0.08/)![m] Thickness of each soil of the four top layers in default clm5.
real(r8),parameter                           :: D = 1.14e-8![m2/h] Diffusivity. Based on Koven et al 2013, 1cm2/yr = 10e-4/(24*365)
real(r8)                                     :: diffusive_source,diffusive_sink !For use in vertical_diffusion subroutine. Currently using alt_vertical diffusion subroutine
real(r8)                                     :: net_diffusion(4,pool_types)
!end depth & vertical transport

!Fluxes between pools:
real(r8) :: LITmSAP, LITsSAP, EcMSAP, ErMSAP, AMSAP, EcMSOMp, EcMSOMa, EcMSOMc, ErMSOMp, ErMSOMa, ErMSOMc, AMSOMp, AMSOMa, AMSOMc, SOMaSAP, SOMpSOMa, SOMcSOMa
real(r8) :: SAPSOMa, SAPSOMp, SAPSOMc
!End fluxes between pools

integer                                      :: ios = 0 !Changes if something goes wrong when opening a file
character (len=4), dimension(pool_types):: variables = (/  "LITm", "LITs", "SAP ", "EcM ", "ErM ", "AM  ", "SOMp", "SOMa", "SOMc" /)
character (len=10), dimension(pool_types):: change_variables = (/  "changeLITm", "changeLITs", "changeSAP ", "changeEcM ", "changeErM ", "changeAM  ", "changeSOMp", "changeSOMa", "changeSOMc" /)
character (len=*), dimension(20), parameter ::  name_fluxes = (/"LITmSAP ","LITsSAP ", "SAPSOMp ", "SAPSOMa ", "SAPSOMc ","EcMSAP  ","ErMSAP  ", &
      "AMSAP   ","EcMSOMp ", "EcMSOMa ","EcMSOMc ", "ErMSOMp ","ErMSOMa ","ErMSOMc ","AMSOMp  ","AMSOMa  " &
      ,"AMSOMc  ","SOMaSAP ","SOMpSOMa","SOMcSOMa"/)

end module paramMod
