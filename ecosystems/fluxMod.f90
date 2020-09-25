module fluxMod
  use paramMod
  use dispmodule !External module to pretty print matrices (mainly for testing purposes)
  implicit none

  contains

  subroutine N_fluxes(depth, Npool_matrix, Cpool_matrix)
    integer :: depth!depth level
    real(r8),target :: Npool_matrix(nlevdecomp, pool_types+1)
    real(r8),target :: Cpool_matrix(nlevdecomp, pool_types)
    real(r8), pointer :: SAPb, SAPf, EcM, ErM, AM


    SAPb => pool_matrix(depth, 3)
    SAPf => pool_matrix(depth, 4)

    !Inorganic Nitrogen
    N_SAPbDIN  =
    N_SAPfDIN  =
    N_DINPlant =
    N_DINEcM =
    N_DINErM =
    N_DINAM  =

    !SOMa nitrogen (Organic, available N)
    !Dead saprotrophs and mycorrhiza:
    N_SAPbSOMa =
    N_SAPfSOMa =
    N_EcMSOMa =
    N_ErMSOMa =
    N_AMSOMa =
    !Saprotrophs decompose SOMa:
    N_SOMaSAPb =
    N_SOMaSAPf =
    !Mycorrhizal mining (N only):
    N_SOMaEcM =
    N_SOMaErM =
    N_SOMaAM =

    !SOMp nitrogen:
    !Dead saprotrophs and mycorrhiza:
    N_SAPbSOMp =
    N_SAPfSOMp =
    N_EcMSOMp =
    N_ErMSOMp =
    N_AMSOMp =

    !SOMc nitrogen:
    !Dead saprotrophs and mycorrhiza:
    N_SAPbSOMc =
    N_SAPfSOMc =
    N_EcMSOMc =
    N_ErMSOMc =
    N_AMSOMc =

    !Mycorrhizal nitrogen:
    N_EcMPlant =
    N_ErMPlant =
    N_AMPlant =

    !N in saprotrophs: N/C*MMK eq. for Carbon
    N_LITmSAPb =
    N_LITmSAPf =
    N_LITsSAPb =
    N_LITsSAPf =

    !N in litter:
    N_PlantLITm =
    N_PlantLITs =


  end subroutine DIN_fluxes


  subroutine litter_fluxes(depth,pool_matrix) !This subroutine calculates the fluxes out of the litter pools. The input to litter pools comes from vegetation, and is not handeled in this subroutine.
    !This is the total flux from the litter pools. A fraction (1-MGE) is lost to respiration before it reaches the SAP pools. This is handeled in the "decomp" subroutine.

    integer :: depth!depth level
    real(r8),target :: pool_matrix(nlevdecomp, pool_types)
    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: SAPb, SAPf, LITm, LITs
    SAPb => pool_matrix(depth, 3)
    SAPf => pool_matrix(depth, 4)
    LITm => pool_matrix(depth, 1)
    LITs => pool_matrix(depth, 2)

    !From LIT to SAPb
    LITmSAPb=SAPb*Vmax(1)*LITm/(Km(1)+LITm)
    LITsSAPb=SAPb*Vmax(2)*LITs/(Km(2)+LITs)

    LITmSAPf=SAPf*Vmax(4)*LITm/(Km(4)+LITm)
    LITsSAPf=SAPf*Vmax(5)*LITs/(Km(5)+LITs)

    !NOTE: These correspond to eq. A1,A5,A2,A6 in Wieder 2015
    nullify(SAPb, SAPf, LITm,LITs)
  end subroutine litter_fluxes

  subroutine som_fluxes(depth,pool_matrix) !This subroutine calculates the fluxes in and out of the SOM pools.

    integer :: depth!depth level
    real(r8),target :: pool_matrix(nlevdecomp, pool_types)
    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: SOMp,SOMa,SOMc,EcM,ErM,AM, SAPb, SAPf
    SAPb => pool_matrix(depth, 3)
    SAPf => pool_matrix(depth, 4)
    EcM =>  pool_matrix(depth, 5)
    ErM =>  pool_matrix(depth, 6)
    AM =>   pool_matrix(depth, 7)
    SOMp => pool_matrix(depth, 8)
    SOMa => pool_matrix(depth, 9)
    SOMc => pool_matrix(depth, 10)
    !print*, 'EcM', EcM

    EcMSOMp=EcM*k_mycsom(1)*fEcMSOM(1)!somp
    EcMSOMa=EcM*k_mycsom(1)*fEcMSOM(2)!soma
    EcMSOMc=EcM*k_mycsom(1)*fEcMSOM(3)!somc

    ErMSOMp=ErM*k_mycsom(2)*fErMSOM(1)
    ErMSOMa=ErM*k_mycsom(2)*fErMSOM(2)
    ErMSOMc=ErM*k_mycsom(2)*fErMSOM(3)

    AMSOMp=AM*k_mycsom(3)*fAMSOM(1)
    AMSOMa=AM*k_mycsom(3)*fAMSOM(2)
    AMSOMc=AM*k_mycsom(3)*fAMSOM(3)

    !Turnover from SAP to SOM. Based on the turnover equations used in mimics for flux from microbial pools to SOM pools.
    !NOTE: correspond to eq A4,A8 in Wieder 2015
    SAPbSOMp=SAPb*tau(1)*fPHYS(1)
    SAPbSOMa=SAPb*tau(1)*fAVAIL(1)
    SAPbSOMc=SAPb*tau(1)*fCHEM(1) !No arrow on illustration by Haavard and Ella

    SAPfSOMp=SAPf*tau(2)*fPHYS(2)
    SAPfSOMa=SAPf*tau(2)*fAVAIL(2)
    SAPfSOMc=SAPf*tau(2)*fCHEM(2)

    !Based on the equations from SOMa to microbial pools in mimics. On the way, a fraction 1-MGE is lost as respiration. This is handeled in the "decomp" subroutine.
    SOMaSAPb=SAPb*Vmax(3)*SOMa/(Km(3)+SOMa)
    SOMaSAPf=SAPf*Vmax(6)*SOMa/(Km(6)+SOMa)

    !Between SOM pools
    !Desorption of SOMp to SOMa, from Mimics model, eq A9
    SOMpSOMa=SOMp*desorb

    !---Oxidation from SOMc to SOMa
    !From equations for decomposing structural litter in mimics,eq. A10
    !KO modifies Km which is used in the litter->SAP equations.
    SOMcSOMa    = ( SAPb * Vmax(2) * SOMc / (KO(1)*Km(2) + SOMc)) + &
                   (SAPf * Vmax(5) * SOMc / (KO(2)*Km(5) + SOMc))

    nullify( SOMp,SOMa,SOMc,EcM,ErM,AM, SAPb,SAPf)
  end subroutine som_fluxes

  !NOTE myc->sap Commented out for now, let all dead mycorrhiza go to SOMa. see Commit: 075363ec3c7b30471e717ca2a7c94585366cdf56
  ! subroutine microbial_fluxes(depth,pool_matrix,level_max) !This subroutine calculates the fluxes between the microbial pools (Mycorrhiza and SAPotrophs).
  !
  !   integer :: depth,level_max !depth level
  !   real(r8),target :: pool_matrix(level_max, pool_types)
  !
  !   !Creating these pointers improve readability of the flux equations.
  !   real(r8), pointer :: EcM,ErM,AM
  !   EcM => pool_matrix(depth, 5)
  !   ErM => pool_matrix(depth, 6)
  !   AM => pool_matrix(depth, 7)
  !
  !   From Mycorrhizal pools to SAProtroph pools
  !   Mycorrhizal pool*fraction to SAP*fraction to SAP_r*decay constant for mycorrhizal pool.
  !   EcMSAPf=EcM*MYC_SAPf*k_mycsap(1)
  !   EcMSAPf=EcM*MYC_SAPb*k_mycsap(1)
  !
  !   ErMSAPf=ErM*MYC_SAPf*k_mycsap(2)
  !   ErMSAPb=ErM*MYC_SAPb*k_mycsap(2)
  !
  !   AMSAPb=AM*MYC_SAPb*k_mycsap(3)
  !   AMSAPf=AM*MYC_SAPf*k_mycsap(3)
  !
  !   nullify(EcM,ErM,AM)
  ! end subroutine microbial_fluxes

  subroutine vertical_diffusion(tot_diffusion_dummy,upper_diffusion_flux,lower_diffusion_flux,pool_matrix,vert,t,counter,step_frac) !This subroutine calculates the vertical transport of carbon through the soil layers.

      integer               :: depth, pool,counter
      real(r8)              :: pool_matrix(nlevdecomp, pool_types)
      real(r8),intent(out)  :: upper_diffusion_flux, lower_diffusion_flux
      real(r8), intent(out) :: tot_diffusion_dummy ![gC/day]
      real(r8),intent(out)  :: vert(nlevdecomp, pool_types)
      real(r8)              :: t !t*dt in main routine
      real(r8)              :: sum_day=0.0
      real(r8)              :: step_frac

      !In a timestep, the fluxes between pools in the same layer is calculated before the vertical diffusion. Therefore, a loop over all the entries in
      !pool_matrix is used here to calculate the diffusion.
      do depth = 1,nlevdecomp
        do pool =1, pool_types

          !eq. 6.18 and 6.20 from Soetaert & Herman, A practical guide to ecological modelling.
          if (depth == 1) then
            upper_diffusion_flux= 0.0
            lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
          elseif (depth==nlevdecomp) then
            upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
            lower_diffusion_flux= 0.0
          else
            upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
            lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
          end if

          tot_diffusion_dummy=(upper_diffusion_flux-lower_diffusion_flux)/delta_z(depth)
          vert(depth,pool) = tot_diffusion_dummy
          sum_day=sum_day+tot_diffusion_dummy
        end do !pool
      end do !depth

  end subroutine vertical_diffusion

  subroutine moisture_func(theta_l,theta_sat, theta_f,r_moist)
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
end module fluxMod
