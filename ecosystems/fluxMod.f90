module fluxMod
  use paramMod
  !use dispmodule !External module to pretty print matrices (mainly for testing purposes)
  implicit none

  contains

  subroutine litter_fluxes(depth,pool_matrix,level_max) !This subroutine calculates the fluxes out of the litter pools. The input to litter pools comes from vegetation, and is not handeled in this subroutine.
    !This is the total flux from the litter pools. A fraction (1-MGE) is lost to respiration before it reaches the SAP pools. This is handeled in the "decomp" subroutine.

    integer :: depth,level_max !depth level
    real(r8),target :: pool_matrix(level_max, pool_types)
    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: SAPb, SAPf, LITm, LITs
    SAPb => pool_matrix(depth, 3)
    SAPf => pool_matrix(depth, 4)
    LITm => pool_matrix(depth, 1)
    LITs => pool_matrix(depth, 2)

    !From LIT to SAPb
    LITmSAPb=f_depth(nlevdecomp +1 -depth)*SAPb*Vmax(1)*LITm/(Km(1)+LITm)
    LITsSAPb=SAPb*Vmax(2)*LITs/(Km(2)+LITs)*f_depth(nlevdecomp +1 -depth)
    !print*, "LITmSAPb",LITmSAPb, depth
    !print*, "LITsSAPb",LITsSAPb, depth
    !From LIT to sapf
    LITmSAPf=SAPf*Vmax(4)*LITm/(Km(4)+LITm)*f_depth(depth)
    LITsSAPf=SAPf*Vmax(5)*LITs/(Km(5)+LITs)*f_depth(depth)

    !NOTE: These correspond to eq. A1,A5,A2,A6 in Wieder 2015
    nullify(SAPb, SAPf, LITm,LITs)
  end subroutine litter_fluxes

  subroutine som_fluxes(depth,pool_matrix,level_max) !This subroutine calculates the fluxes in and out of the SOM pools.

    integer :: depth,level_max !depth level
    real(r8),target :: pool_matrix(level_max, pool_types)
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

    EcMSOMp=EcM*k_mycsom(1)*0.40!somp
    EcMSOMa=EcM*k_mycsom(1)*0.40!soma
    EcMSOMc=EcM*k_mycsom(1)*0.2!somc

    ErMSOMp=ErM*k_mycsom(2)*0.2
    ErMSOMa=ErM*k_mycsom(2)*0.4
    ErMSOMc=ErM*k_mycsom(2)*0.3

    AMSOMp=AM*k_mycsom(3)*0.3
    AMSOMa=AM*k_mycsom(3)*0.3
    AMSOMc=AM*k_mycsom(3)*0.4

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
                   (SAPf* Vmax(5) * SOMc / (KO(2)*Km(5) + SOMc))

    nullify( SOMp,SOMa,SOMc,EcM,ErM,AM, SAPb,SAPf)
  end subroutine som_fluxes

  subroutine microbial_fluxes(depth,pool_matrix,level_max) !This subroutine calculates the fluxes between the microbial pools (Mycorrhiza and SAPotrophs).

    integer :: depth,level_max !depth level
    real(r8),target :: pool_matrix(level_max, pool_types)

    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: EcM,ErM,AM
    EcM => pool_matrix(depth, 5)
    ErM => pool_matrix(depth, 6)
    AM => pool_matrix(depth, 7)

    !From Mycorrhizal pools to SAPotroph pools
    !Mycorrhizal pool*fraction to SAP*fraction to SAP_r*decay constant for mycorrhizal pool. TODO: Maybe reconsider these equations..
    EcMSAPb=EcM*MYC_SAPb*k_mycsap(1)
    EcMSAPf=EcM*MYC_SAPf*k_mycsap(1)

    ErMSAPf=ErM*MYC_SAPf*k_mycsap(2)
    ErMSAPb=ErM*MYC_SAPb*k_mycsap(2)

    AMSAPb=AM*MYC_SAPb*k_mycsap(3)
    AMSAPf=AM*MYC_SAPf*k_mycsap(3)


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
          ! if (counter == step_frac) then
          !   write(unit=10,fmt='(F10.0,A2,I2,A2,I6,A2,F30.10,A2,F30.10,A2,F30.10)') &
          !   t,',',depth,',',pool,',',tot_diffusion_dummy,',', upper_diffusion_flux,',', lower_diffusion_flux
          !   sum_day=0.0
          ! end if !writing
        end do !pool
      end do !depth

  end subroutine vertical_diffusion


  subroutine moisture_func(r_moist)
    real(r8), intent(out) :: r_moist

    r_moist = max(0.05, P*(theta_l/theta_sat)**3*(1-(theta_l/theta_sat)-(theta_f/theta_sat))**gas_diffusion)
  end subroutine moisture_func
end module fluxMod
