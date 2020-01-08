module fluxMod
  use paramMod
  !use dispmodule !External module to pretty print matrices (mainly for testing purposes)
  implicit none

  contains

  subroutine litter_fluxes(depth,pool_matrix,level_max) !This subroutine calculates the fluxes out of the litter pools. The input to litter pools comes from vegetation, and is not handeled in this subroutine.
    !This is the total flux from the litter pools. A fraction (1-MGE) is lost to respiration before it reaches the SAP pool. This is handeled in the "decomp" subroutine.

    integer :: depth,level_max !depth level
    real(r8),target :: pool_matrix(level_max, pool_types)
    !Creating these pointers improve readability of the flux equations. TODO: Do you really need the matrix?
    real(r8), pointer :: SAP, LITm, LITs
    SAP => pool_matrix(depth, 3)
    !SAPk => pool_matrix(depth, 4)
    LITm => pool_matrix(depth, 1)
    LITs => pool_matrix(depth, 2)

    !From LIT to SAP
    LITmSAP=SAP*Vmax(1)*LITm/(Km(1)+LITm)
    LITsSAP=SAP*Vmax(2)*LITs/(Km(2)+LITs)

    !From LIT to sapk
    !LITtoSAP(2)=SAPk*Vmax(4)*LITm/(Km(4)+LITm)
    !LITtoSAP(4)=SAPk*Vmax(5)*LITs/(Km(5)+LITs)

    !NOTE: These correspond to eq. A1,A5,A2,A6 in Wieder 2015
    !TODO Maybe declare one name for each flux for better readability (compared to LITtoSAP etc arrays...)
    nullify(SAP, LITm,LITs)
  end subroutine litter_fluxes

  subroutine som_fluxes(depth,pool_matrix,level_max) !This subroutine calculates the fluxes in and out of the SOM pools.

    integer :: depth,level_max !depth level
    real(r8),target :: pool_matrix(level_max, pool_types)
    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: SOMp,SOMa,SOMc,EcM,ErM,AM, SAP
    SAP => pool_matrix(depth, 3)
    !SAPk => pool_matrix(depth, 4)
    EcM =>  pool_matrix(depth, 5)
    ErM =>  pool_matrix(depth, 6)
    AM =>   pool_matrix(depth, 7)
    SOMp => pool_matrix(depth, 8)
    SOMa => pool_matrix(depth, 9)
    SOMc => pool_matrix(depth, 10)
    !print*, 'EcM', EcM

    !TODO No idea if the Michalis Menten equations is the right thing to use here..
    !       !03.09.19: Probably. not bc. the mycorrhiza is not a substrate?
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
    EcMSOMp=EcM*k2(1)*0.40!somp
    EcMSOMa=EcM*k2(1)*0.40!soma
    EcMSOMc=EcM*k2(1)*0.2!somc

    ErMSOMp=ErM*k2(2)*0.2
    ErMSOMa=ErM*k2(2)*0.6
    ErMSOMc=ErM*k2(2)*0.2

    AMSOMp=AM*k2(3)*0.1
    AMSOMa=AM*k2(3)*0.8
    AMSOMc=AM*k2(3)*0.1
    !Turnover from SAP to SOM. Based on the turnover equations used in mimics for flux from microbial pools to SOM pools.
    !NOTE: correspond to eq A4,A8 in Wieder 2015

    SAPSOMp=SAP*tau(1)*fPHYS(1)
    SAPSOMa=SAP*tau(1)*fAVAIL(1)
    SAPSOMc=SAP*tau(1)*fCHEM(1) !No arrow on illustration by Haavard and Ella

    !SAPtoSOM(4)=SAPk*tau(2)*fPHYS(2)
    !SAPtoSOM(5)=SAPk*tau(2)*fAVAIL(2)
    !SAPtoSOM(6)=SAPk*tau(2)*fCHEM(2)

    !Based on the equations from SOMa to microbial pools in mimics. On the way, a fraction 1-MGE is lost as respiration. This is handeled in the "decomp" subroutine.
    SOMaSAP=SAP*Vmax(3)*SOMa/(Km(3)+SOMa)
    !SOMtoSAP(2)=SAPk*Vmax(6)*SOMa/(Km(6)+SOMa)


    !Between SOM pools
    !Desorption of SOMp to SOMa, from Mimics model, eq A9
    SOMpSOMa=SOMp*desorb

    !---Oxidation from SOMc to SOMa
    !From equations for decomposing structural litter in mimics,eq. A10
    !KO modifies Km which is used in the litter->SAP equations.
    SOMcSOMa    = (( SAP * Vmax(2) * SOMc / (KO(1)*Km(2) + SOMc)))! + &
                   !(SAPk* Vmax(5) * SOMc / (KO(2)*Km(5) + SOMc)))

    nullify( SOMp,SOMa,SOMc,EcM,ErM,AM, SAP)
  end subroutine som_fluxes

  subroutine microbial_fluxes(depth,pool_matrix,level_max) !This subroutine calculates the fluxes between the microbial pools (Mycorrhiza and SAPotrophs).

    integer :: depth,level_max !depth level
    real(r8),target :: pool_matrix(level_max, pool_types)

    !Creating these pointers improve readability of the flux equations. TODO: Do you really need the matrix?
    real(r8), pointer :: EcM,ErM,AM
    EcM => pool_matrix(depth, 5)
    ErM => pool_matrix(depth, 6)
    AM => pool_matrix(depth, 7)

    !From Mycorrhizal pools to SAPotroph pools
    !Mycorrhizal pool*fraction to SAP*fraction to SAP_r*decay constant for mycorrhizal pool. TODO: Maybe reconsider these equations..
    EcMSAP=EcM*Myc_SAP*k(1)
    ErMSAP=ErM*Myc_SAP*k(2)
    AMSAP=AM*Myc_SAP*k(3)

    !MYCtoSAP(4)=EcM*Myc_SAPk*k(1)
    !MYCtoSAP(5)=ErM*Myc_SAPk*k(2)
    !MYCtoSAP(6)=AM*Myc_SAPk*k(3)

    nullify(EcM,ErM,AM)
  end subroutine microbial_fluxes

  subroutine alt_vertical_diffusion(depth, pool,tot_diffusion_dummy,upper_diffusion_flux,lower_diffusion_flux,pool_matrix,level_max) !This subroutine calculates the vertical transport of carbon through the soil layers.
      integer               :: depth, pool,level_max
      real(r8)              :: pool_matrix(level_max, pool_types)
      real(r8),intent(out)  :: upper_diffusion_flux, lower_diffusion_flux
      real(r8), intent(out) :: tot_diffusion_dummy ![gC/day]

      !eq. 6.18 and 6.20 from Soetaert & Herman, A practical guide to ecological modelling.
      if (depth == 1) then
        upper_diffusion_flux= 0.0
        lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
      elseif (depth==level_max) then
        upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
        lower_diffusion_flux= 0.0
      else
        upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
        lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
      end if

      tot_diffusion_dummy=(upper_diffusion_flux-lower_diffusion_flux)/delta_z(depth)
  end subroutine alt_vertical_diffusion

  subroutine vertical_diffusion(tot_diffusion_dummy,upper_diffusion_flux,lower_diffusion_flux,pool_matrix,level_max,vert,t,counter,step_frac) !This subroutine calculates the vertical transport of carbon through the soil layers.

      integer               :: level_max, depth, pool,counter
      real(r8)              :: pool_matrix(level_max, pool_types)
      real(r8),intent(out)  :: upper_diffusion_flux, lower_diffusion_flux
      real(r8), intent(out) :: tot_diffusion_dummy ![gC/day]
      real(r8),intent(out)  :: vert(level_max, pool_types)
      real(r8)              :: t !t*dt in main routine
      real(r8)              :: sum_day=0.0
      real(r8)              :: step_frac
      !eq. 6.18 and 6.20 from Soetaert & Herman, A practical guide to ecological modelling.
      do depth = 1,level_max
        do pool =1, pool_types
          if (depth == 1) then
            upper_diffusion_flux= 0.0
            lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
          elseif (depth==level_max) then
            upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
            lower_diffusion_flux= 0.0
          else
            upper_diffusion_flux=-D*(pool_matrix(depth,pool)-pool_matrix(depth-1,pool))/(node_z(depth)-node_z(depth-1))
            lower_diffusion_flux=-D*(pool_matrix(depth+1,pool)-pool_matrix(depth,pool))/(node_z(depth+1)-node_z(depth))
          end if
          tot_diffusion_dummy=(upper_diffusion_flux-lower_diffusion_flux)/delta_z(depth)
          vert(depth,pool) = tot_diffusion_dummy
          sum_day=sum_day+tot_diffusion_dummy
          if (counter == step_frac) then
            write(unit=10,fmt='(F10.0,A2,I2,A2,I6,A2,F30.10,A2,F30.10,A2,F30.10)') &
            t,',',depth,',',pool,',',tot_diffusion_dummy,',', upper_diffusion_flux,',', lower_diffusion_flux
            sum_day=0.0
          end if !writing
        end do !pool
      end do !depth


  end subroutine vertical_diffusion


end module fluxMod
