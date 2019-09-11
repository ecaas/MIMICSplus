module fluxMod
  use paramMod
  use dispmodule !External module to pretty print matrices (mainly for testing purposes)
  implicit none

  contains

  subroutine litter_fluxes(depth,pool_matrix,level_max) !This subroutine calculates the fluxes out of the litter pools. The input to litter pools comes from vegetation, and is not handeled in this subroutine.
    !This is the total flux from the litter pools. A fraction (1-MGE) is lost to respiration before it reaches the SAP pool. This is handeled in the "decomp" subroutine.

    integer :: depth,level_max !depth level
    real(r8),target :: pool_matrix(level_max, pool_types)
    !Creating these pointers improve readability of the flux equations. TODO: Do you really need the matrix?
    real(r8), pointer :: SAPr, SAPk, LITm, LITs
    SAPr => pool_matrix(depth, 3)
    SAPk => pool_matrix(depth, 4)
    LITm => pool_matrix(depth, 1)
    LITs => pool_matrix(depth, 2)

    !From LITm to sap
    LITtoSAP(1)=SAPr*Vmax(1)*LITm/(Km(1)+LITm)
    LITtoSAP(2)=SAPk*Vmax(2)*LITm/(Km(2)+LITm)
    !From LITs to SAP
    LITtoSAP(3)=SAPr*Vmax(3)*LITs/(Km(3)+LITs)
    LITtoSAP(4)=SAPk*Vmax(4)*LITs/(Km(4)+LITs)
    !NOTE: These correspond to eq. A1,A5,A2,A6 in Wieder 2015
    !TODO Maybe declare one name for each flux for better readability (compared to LITtoSAP etc arrays...)
    nullify(SAPr, SAPk, LITm,LITs)
  end subroutine litter_fluxes

  subroutine som_fluxes(depth,pool_matrix,level_max) !This subroutine calculates the fluxes in and out of the SOM pools.

    integer :: depth,level_max !depth level
    real(r8),target :: pool_matrix(level_max, pool_types)
    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: SOMp,SOMa,SOMc,EcM,ErM,AM, SAPr,SAPk
    SAPr => pool_matrix(depth, 3)
    SAPk => pool_matrix(depth, 4)
    EcM =>  pool_matrix(depth, 5)
    ErM =>  pool_matrix(depth, 6)
    AM =>   pool_matrix(depth, 7)
    SOMp => pool_matrix(depth, 8)
    SOMa => pool_matrix(depth, 9)
    SOMc => pool_matrix(depth, 10)
    !print*, 'EcM', EcM

    !TODO No idea if the Michalis Menten equations is the right thing to use here..
    !       !03.09.19: Prob. not bc. the mycorrhiza is not a substrate?
    ! !From Mycorrhiza pools to SOMp
    ! MYCtoSOM(1)=SOMp*Vmax(5)*EcM/(Km(5)+EcM)
    ! MYCtoSOM(2)=SOMp*Vmax(6)*ErM/(Km(6)+ErM)
    ! MYCtoSOM(3)=SOMp*Vmax(7)*AM/(Km(7)+AM)
    ! !From Mycorrhiza pools to SOMa
    ! MYCtoSOM(4)=SOMa*Vmax(8)*EcM/(Km(8)+EcM)
    ! !print*, MYCtoSOM(4)
    ! MYCtoSOM(5)=SOMa*Vmax(9)*ErM/(Km(9)+ErM)
    ! MYCtoSOM(6)=SOMa*Vmax(10)*AM/(Km(10)+AM)
    ! !From Mycorrhiza pools to SOMc
    ! MYCtoSOM(7)=SOMc*Vmax(11)*EcM/(Km(11)+EcM)
    ! MYCtoSOM(8)=SOMc*Vmax(12)*ErM/(Km(12)+ErM)
    ! MYCtoSOM(9)=SOMc*Vmax(13)*AM/(Km(13)+AM)
    MYCtoSOM(1)=EcM*k2(1)*0.25
    MYCtoSOM(2)=EcM*k2(1)*0.5
    MYCtoSOM(3)=EcM*k2(1)*0.25

    MYCtoSOM(4)=ErM*k2(2)*0.25
    MYCtoSOM(5)=ErM*k2(2)*0.65
    MYCtoSOM(6)=ErM*k2(2)*0.1

    MYCtoSOM(7)=AM*k2(3)*0.3
    MYCtoSOM(8)=AM*k2(3)*0.5
    MYCtoSOM(9)=AM*k2(3)*0.2
    !Turnover from SAP to SOM. Based on the turnover equations used in mimics for flux from microbial pools to SOM pools.
    !NOTE: correspond to eq A4,A8 in Wieder 2015
    SAPtoSOM(1)=SAPr*tau(1)*fPHYS(1)
    SAPtoSOM(2)=SAPr*tau(1)*fAVAIL(1)
    SAPtoSOM(3)=SAPr*tau(1)*fCHEM(1) !No arrow on illustration by Haavard and Ella

    SAPtoSOM(4)=SAPk*tau(2)*fPHYS(2)
    SAPtoSOM(5)=SAPk*tau(2)*fAVAIL(2)
    SAPtoSOM(6)=SAPk*tau(2)*fCHEM(2)

    !Based on the equations from SOMa to microbial pools in mimics. On the way, a fraction 1-MGE is lost as respiration. This is handeled in the "decomp" subroutine.
    SOMtoSAP(2)=SAPk*Vmax(4)*SOMa/(Km(3)+SOMa)
    SOMtoSAP(1)=SAPr*Vmax(3)*SOMa/(Km(3)+SOMa)

    !Between SOM pools
    !Desorption of SOMp to SOMa, from Mimics model, eq A9
    SOMtoSOM(1)=SOMp*desorb

    !---Oxidation from SOMc to SOMa
    !From equations for decomposing structural litter in mimics,eq. A10
    !KO modifies Km which is used in the litter->SAP equations.
    SOMtoSOM(2)    = (( SAPr * Vmax(2) * SOMc / (KO(1)*Km(2) + SOMc)) + &
                   (SAPk* Vmax(4) * SOMc / (KO(2)*Km(4) + SOMc)))

    nullify( SOMp,SOMa,SOMc,EcM,ErM,AM, SAPr,SAPk)
  end subroutine som_fluxes

  subroutine microbial_fluxes(depth,pool_matrix,level_max) !This subroutine calculates the fluxes between the microbial pools (Mycorrhiza and saprotrophs).

    integer :: depth,level_max !depth level
    real(r8),target :: pool_matrix(level_max, pool_types)

    !Creating these pointers improve readability of the flux equations. TODO: Do you really need the matrix?
    real(r8), pointer :: EcM,ErM,AM
    EcM => pool_matrix(depth, 5)
    ErM => pool_matrix(depth, 6)
    AM => pool_matrix(depth, 7)

    !From Mycorrhizal pools to saprotroph pools
    !Mycorrhizal pool*fraction to SAP*fraction to SAP_r*decay constant for mycorrhizal pool. TODO: Maybe reconsider these equations..
    MYCtoSAP(1)=EcM*0.3*Myc_SAPr*k(1)
    MYCtoSAP(2)=ErM*0.1*Myc_SAPr*k(2)
    MYCtoSAP(3)=AM*0.3*Myc_SAPr*k(3)

    MYCtoSAP(4)=EcM*0.3*Myc_SAPk*k(4)
    MYCtoSAP(5)=ErM*0.1*Myc_SAPk*k(5)
    MYCtoSAP(6)=AM*0.3*Myc_SAPk*k(6)

    nullify(EcM,ErM,AM)
  end subroutine microbial_fluxes

  subroutine alt_vertical_diffusion(depth, pool,tot_diffusion_dummy,upper_diffusion_flux,lower_diffusion_flux,pool_matrix,level_max) !This subroutine calculates the vertical transport of carbon through the soil layers.
      integer :: depth, pool,level_max
      real(r8) :: pool_matrix(level_max, pool_types)
      real(r8),intent(out) :: upper_diffusion_flux, lower_diffusion_flux
      real(r8), intent(out) :: tot_diffusion_dummy ![gC/day]
      !eq. 6.18 and 6.20 from Soetaert & Herman, A practical guide to ecological modelling. TODO: Check if the signs make sense here
      if (depth == 1) then
        upper_diffusion_flux= 0.0
        lower_diffusion_flux=-D(depth+1)*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
      elseif (depth==level_max) then
        upper_diffusion_flux=-D(depth)*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
        lower_diffusion_flux= 0.0
      else
        upper_diffusion_flux=-D(depth)*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
        lower_diffusion_flux=-D(depth+1)*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
      end if

      tot_diffusion_dummy=(-upper_diffusion_flux-lower_diffusion_flux)/delta_z(depth)
  end subroutine alt_vertical_diffusion

end module fluxMod
