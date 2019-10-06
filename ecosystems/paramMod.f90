module paramMod
use shr_kind_mod   , only : r8 => shr_kind_r8
!use SoilDecompType, only : carbon_pool_type

implicit none
!Define variables

!Things that should be given as input to the decomposition model:
real(kind=r8),parameter                       :: tsoi   = 8.9                   ![degC]
real(kind=r8),parameter                       :: fclay  = 0.12                  ![-] fraction of clay in soil
real(kind=r8),dimension(4)                    :: MGE = (/ 0.59,0.50,0.50,0.65 /)![mg/mg] Microbial growth efficiency/Carbon Use efficiency. The fraction of the flux from litter pools that is used in microbial processes.
                                                 !The rest is lost in respiration.
                                                 !/Litm-sapr, litm->sapk, lits->sapr, lits->sapk/. These should be determined carefully. pH/soil quality/N availability may be important
                                                 !Maybe also vary with time/depth..? See DOI:10.1016/j.soilbio.2018.09.036 and DOI: 10.1016/J.SOILBIO.2019.03.008
real(r8), dimension(2)                        :: veg_input                      ![gC/(m2 day)] vegetational input to litter pools

!real(r8), dimension(2),parameter              :: IfracSOM = (/0.05,0.05/)     ![-]Fraction of vegetational input that goes directly to SOM pools
!real(r8), parameter                           :: I_trees=0.2, I_eric=0.5, I_gras=0.3 ![gC/(m2 day)] vegetational input to mycorrhizal pools. TODO: Should vary with time and somehow combine I_trees and veg_input

real(r8)                                      :: f_som1=0.05, f_som2=0.05
real(r8)                                      :: f_myc_levels=1, f_lit_1=1, f_lit_234=1

real(kind=r8),dimension(6),parameter          :: k = (/0.05,0.5,0.005,0.05,0.5,0.005/)![1/day](EcM_sapr, ErM_sapr, AM_sapr,EcM_sapr, ErM_sapr, AM_sapr) decay constants, MYC to SAP pools
real(kind=r8),dimension(3),parameter          :: k2 = (/0.0007,0.0007,0.00014/)! [1/day] decay constants, MYC to SOM pools
real, parameter                               :: Myc_SAPr=0.6, Myc_SAPk=1-Myc_SAPr    ![-]Fraction of the flux from mycorrhizal pools to SAPk pool. The rest is going to SAPr

!For calculating the Km parameter in Michaelis Menten kinetics (expressions based on mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and https://github.com/wwieder/MIMICS)
!TODO: Varying Kslope?
integer, parameter                            :: MM_eqs = 4       !Number of Michaelis-Menten parameters
real(kind=r8),dimension(MM_eqs),parameter     :: Kslope  = (/0.017, 0.027, 0.017, 0.027/) ![ln(mgC/cm3)/degC]
real(kind=r8),dimension(MM_eqs)               :: Km = exp(Kslope*tsoi+3.19)*10![mgC/cm3]*10e3=[gC/m3]
real(kind=r8),dimension(2),parameter          :: KO =(4)   ![-]Increases Km (the half saturation constant for oxidation of chemically protected SOM, SOM_c) from mimics
real(kind=r8),parameter,dimension(MM_eqs)     :: Vmax = (/ 0.012, 0.002, 0.012, 0.003/)*24/10![mgC/((mgSAP)h)]*24= [mgC/(mgSAP*day)]  For use in Michaelis menten kinetics. TODO: Is mgSAP only carbon?
!end MM parameters

real(kind=r8)                                 :: desorb = 3e-4 * exp(-4*(sqrt(fclay)))*24 ![1/h]*24=[1/day]From Mimics, used for the transport from physically protected SOM to available SOM pool
!For calculating turnover from SAP to SOM (expressions from mimics model: https://doi.org/10.5194/gmd-8-1789-2015 and  https://github.com/wwieder/MIMICS)

real(kind=r8)                                 :: fMET  !Assume half of the input goes to LITm
real(kind=r8),dimension(2)     :: fPHYS
real(kind=r8),dimension(2)     :: fCHEM
real(kind=r8),dimension(2)     :: fAVAIL
real(kind=r8),dimension(2)      :: tau ![1/day]


!Pools: NOTE: This needs to be updated if pools are added to/removed from the system.
integer, parameter                            :: num_soilc = 1  !Number of soil columns modeled. Only use one column here, but TODO code should be working for several columns as well (add loop)
integer, parameter                            :: no_of_litter_pools = 2     !Metabolic and structural
integer, parameter                            :: no_of_sap_pools = 2        !saprotrophs r and k
integer, parameter                            :: no_of_myc_pools = 3        !Mycorrhiza: Ecto, ericoid & arbuscular
integer, parameter                            :: no_of_som_pools = 3        !Physically protected, chemically protected, available carbon
integer, parameter                            :: pool_types = no_of_litter_pools + no_of_myc_pools + &
                                                 no_of_sap_pools + no_of_som_pools
real(r8),dimension(4)                :: node_z = (/0.05,0.2,0.5,1.10/) ![m] Depth of center in each soil layer
real(r8),dimension(4)                :: delta_z = (/0.1,0.2,0.4,0.8/)![m] Thickness of each soil layer NOTE: delta_z and node_z is used in the alt_vertical_diffusion routine
!real(r8),dimension(4)                :: node_z = (/0.1,0.3,0.5,0.7/) ![m] Depth of center in each soil layer
!real(r8),dimension(4)                :: delta_z = (/0.2,0.2,0.2,0.2/)![m] Thickness of each soil layer NOTE: delta_z and node_z is used in the alt_vertical_diffusion routine

!real(r8),dimension(nlevdecomp)               :: node_z = (/0.01,0.04,0.09,0.16/) ![m] Depth of center in each soil layer. Same as the first 4 layers of default CLM5 with vertical resolution.
!real(r8),dimension(nlevdecomp)               :: delta_z = (/0.02, 0.04, 0.06, 0.08/)![m] Thickness of each soil of the four top layers in default clm5.
real(r8),parameter    :: D = 0.0027*10e-4!(/0.0,0.0027,0.0027,0.0027,0.0/)*10e-4 ![cm2/day]*10e-4m2/cm2 = [m2/day] Diffusivity. Based on Koven et al 2013, 1cm2/yr = 0.0027cm2/day

real(r8)                                      :: diffusive_source,diffusive_sink !For use in vertical_diffusion subroutine. Currently using alt_vertical diffusion subroutine
real(r8)                                      :: net_diffusion(4,pool_types)
real(r8),parameter                            :: depth = 0.56![m]
!real(r8),target                               :: pool_matrix(1,pool_types)=0.0
!real(r8)                                      :: change_matrix(1,pool_types)=0.0
!real(r8)                                      :: a_noVert(1,pool_types)=0.0

!Fluxes between pools:
!TODO: change hardcoded dimensions
real(kind=r8),dimension(4)                    :: LITtoSAP ! (Lit_m->SAP_r, LIT_m->saP_k, LIT_s->SAP_r, LIT_s->SAP_k)
real(kind=r8),dimension(6)                    :: SAPtoSOM ! SAP_r to phys, avail, chem, SAP_k to phys, avail, chem
real(kind=r8),dimension(6)                    :: MYCtoSAP ! Ecm to SAP_r, SAPk, ErM to SAP_r, SAP_k, AM to SAP_r,SAP_k
real(kind=r8),dimension(9)                    :: MYCtoSOM ! Ecm to SOM_p, SOM_a, SOM_c, ErM to SOM_p, SOM_a, SOM_c, AM to SOM_p, SOM_a, SOM_c
real(kind=r8),dimension(2)                    :: SOMtoSAP ! SOM_a -> SAP_r, SAP_k
real(kind=r8),dimension(2)                    :: SOMtoSOM ! SOM_p->SOM_a, SOM_c->SOM_a

!integer                                       :: i,fc,t!, j        !For loops
integer                                       :: ios = 0 !Changes if something goes wrong when opening a file
end module paramMod
