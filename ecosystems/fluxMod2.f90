module fluxMod2
  use paramMod
  use dispmodule !External module to pretty print matrices (mainly for testing purposes)
  implicit none

  contains
    !testchange
  subroutine calculate_fluxes(depth,C_pool_matrix,N_pool_matrix, C_plant, N_plant) !This subroutine calculates the fluxes in and out of the SOM pools.
    integer :: depth!depth level
    real(r8) :: C_plant, N_plant
    real(r8),target :: C_pool_matrix(nlevdecomp, pool_types)
    real(r8),target :: N_pool_matrix(nlevdecomp, pool_types+1)
    !Creating these pointers improve readability of the flux equations.
    real(r8), pointer :: C_LITm, C_LITs, C_SOMp,C_SOMa,C_SOMc,C_EcM,C_ErM,C_AM, &
    C_SAPb, C_SAPf, N_LITm, N_LITs, N_SOMp,N_SOMa,N_SOMc,N_EcM,N_ErM,N_AM, N_SAPb, N_SAPf, N_IN
    C_LITm => C_pool_matrix(depth, 1)
    C_LITs => C_pool_matrix(depth, 2)
    C_SAPb => C_pool_matrix(depth, 3)
    C_SAPf => C_pool_matrix(depth, 4)
    C_ErM =>  C_pool_matrix(depth, 6)
    C_EcM =>  C_pool_matrix(depth, 5)
    C_AM =>   C_pool_matrix(depth, 7)
    C_SOMp => C_pool_matrix(depth, 8)
    C_SOMa => C_pool_matrix(depth, 9)
    C_SOMc => C_pool_matrix(depth, 10)

    N_LITm => N_pool_matrix(depth, 1)
    N_LITs => N_pool_matrix(depth, 2)
    N_SAPb => N_pool_matrix(depth, 3)
    N_SAPf => N_pool_matrix(depth, 4)
    N_EcM =>  N_pool_matrix(depth, 5)
    N_ErM =>  N_pool_matrix(depth, 6)
    N_AM =>   N_pool_matrix(depth, 7)
    N_SOMp => N_pool_matrix(depth, 8)
    N_SOMa => N_pool_matrix(depth, 9)
    N_SOMc => N_pool_matrix(depth, 10)
    N_IN => N_pool_matrix(depth, 11)

    !------------------CARBON FLUXES----------------------------:
    C_PR = gamma_rs*C_Plant/(1+gamma_rs) !Carbon in plant roots
    C_PS = C_plant/(1+ gamma_rs)!Carbon in plant shoots

    N_PR = gamma_rs*N_Plant/(1+gamma_rs) !N in plant roots
    N_PS = N_plant/(1+ gamma_rs)!N in plant shoots
    !print*, "N_PS, N_PR, N_Plant", N_PS, N_PR, N_Plant
    P_N =  a-b*C_PS !Plant N productivity (as in Baskaran 2016, eq (12))
    if (P_N < 0) then
      print*, "P_N < 0: ", P_N
      !P_N = 0.0
    end if
    !Plant growth rate
    C_growth_rate = (1-delta)*P_N*N_PS                                 !NOTE Usikker paa enheter her
  !  print*, 'C_growth_rate', C_growth_rate
  !  print*, "C_PS", C_PS, "N_PS", N_PS

    !Plant Carbon to mycorrhiza: = delta*P_N*N_PS                     !TODO: differentiate between the different mycorrhizae (By using myc specific gamma_rs?)
    C_PlantEcM = delta*0.4*P_N*N_PS/delta_z(depth)                        !NOTE: Deler pa lagdybde for a fordele inputen fra planten likt over alle lagene
    C_PlantErM = delta*0.3*P_N*N_PS/delta_z(depth)                        !gC/m3h  (?)
    C_PlantAM = delta*0.3*P_N*N_PS/delta_z(depth)

    !Used to calculate litter production in flux subroutine:
    Total_plant_mortality = (my_shoot + gamma_rs*my_root)*(C_plant/(1+gamma_rs))!gC/m2h
    !Plant mortality/litter production:
    C_PlantLITm = fMET*Total_plant_mortality/delta_z(depth)                     !gC/m3h
    C_PlantLITs =(1-fMET)*Total_plant_mortality/delta_z(depth)                  !TODO: Blir det riktig a bruke fMET for a dele opp totalproduksjonen?
                                                                                !Tallene fra Baskaran

    !Decomposition of LIT by SAP:
    !On the way, a fraction 1-MGE is lost as respiration. This is handeled in the "decomp" subroutine.
    C_LITmSAPb=C_SAPb*Vmax(1)*C_LITm/(Km(1)+C_LITm)
    C_LITsSAPb=C_SAPb*Vmax(2)*C_LITs/(Km(2)+C_LITs)
    C_LITmSAPf=C_SAPf*Vmax(4)*C_LITm/(Km(4)+C_LITm)
    C_LITsSAPf=C_SAPf*Vmax(5)*C_LITs/(Km(5)+C_LITs)

    !Decomposition of SOMa by SAP. Based on the equations from SOMa to microbial pools in mimics.
    !On the way, a fraction 1-MGE is lost as respiration. This is handeled in the "decomp" subroutine.
    C_SOMaSAPb=C_SAPb*Vmax(3)*C_SOMa/(Km(3)+C_SOMa)
    C_SOMaSAPf=C_SAPf*Vmax(6)*C_SOMa/(Km(6)+C_SOMa)

    !Dead mycorrhizal biomass enters the SOM pools:                             gC/m3h
    C_EcMSOMp=C_EcM*k_mycsom(1)*fEcMSOM(1)!somp
    C_EcMSOMa=C_EcM*k_mycsom(1)*fEcMSOM(2)!soma
    C_EcMSOMc=C_EcM*k_mycsom(1)*fEcMSOM(3)!somc

    C_ErMSOMp=C_ErM*k_mycsom(2)*fErMSOM(1)
    C_ErMSOMa=C_ErM*k_mycsom(2)*fErMSOM(2)
    C_ErMSOMc=C_ErM*k_mycsom(2)*fErMSOM(3)

    C_AMSOMp=C_AM*k_mycsom(3)*fAMSOM(1)
    C_AMSOMa=C_AM*k_mycsom(3)*fAMSOM(2)
    C_AMSOMc=C_AM*k_mycsom(3)*fAMSOM(3)

    !Turnover from SAP to SOM. Based on the turnover equations used in mimics for flux from microbial pools to SOM pools (correspond to eq A4,A8 in Wieder 2015)
    C_SAPbSOMp=C_SAPb*tau(1)*fPHYS(1)                                           !gC/m3h
    C_SAPbSOMa=C_SAPb*tau(1)*fAVAIL(1)
    C_SAPbSOMc=C_SAPb*tau(1)*fCHEM(1)

    C_SAPfSOMp=C_SAPf*tau(2)*fPHYS(2)
    C_SAPfSOMa=C_SAPf*tau(2)*fAVAIL(2)
    C_SAPfSOMc=C_SAPf*tau(2)*fCHEM(2)

    !Desorbtion controls transport from physically protected to available SOM
    C_SOMpSOMa=C_SOMp*desorb

    !Oxidation from SOMc to SOMa
    !From equations for decomposing structural litter in mimics,eq. A10
    !KO modifies Km which is used in the litter->SAP equations.
    C_SOMcSOMa    = ( C_SAPb * Vmax(2) * C_SOMc / (KO(1)*Km(2) + C_SOMc)) + &
                   (C_SAPf * Vmax(5) * C_SOMc / (KO(2)*Km(5) + C_SOMc))

    !Baskaran et al: Rates of decomposition of available SOM mediated by mycorrhizal enzymes:
    Decomp_ecm = K_MO*delta_z(depth)*C_EcM*C_SOMa                               ![gC/m3h]
    Decomp_erm = K_MO*delta_z(depth)*C_ErM*C_SOMa                                !TODO: Somehow incorporate this carbon into SAP pools (?)
    Decomp_am  = K_MO*delta_z(depth)*C_AM*C_SOMa

    !-----------------------------------NITROGEN FLUXES----------------------------:
    !Nitrogen aquired bymycorrhiza via oxidation of SOMa                        gN/m3h
    N_SOMaEcM = Decomp_ecm*N_SOMa/C_SOMa!/CN_ratio(9)
    N_SOMaErM = Decomp_erm*N_SOMa/C_SOMa!/CN_ratio(9)
    N_SOMaAM  = Decomp_am*N_SOMa/C_SOMa!/CN_ratio(9)

    !Inorganic N taken up directly by plant roots                               !Usikker pa enheter
    N_InPlant = V_max_plant*N_in*(1-delta)*(C_PR/(C_PR + Km_plant/delta_z(depth)))
    !Deposition and leacing from the inorganic N pool
    Deposition = Deposition_rate/delta_z(depth)                                          !Usikker pa enheter
    Leaching = Leaching_rate*N_in/delta_z(depth)

    N_INEcM = V_max_myc*N_IN*(C_EcM/(C_EcM + Km_myc/delta_z(depth)))            !Bor MMK parametere vaere spesifikke til mycorrhiza type?
    N_INErM = V_max_myc*N_IN*(C_ErM/(C_ErM + Km_myc/delta_z(depth)))            !Usikker pa enheter
    N_INAM = V_max_myc*N_IN*(C_AM/(C_AM + Km_myc/delta_z(depth)))

    !Plant mortality
    N_PlantLITm = C_PlantLITm*N_Plant/C_Plant
    N_PlantLITs = C_PlantLITs*N_Plant/C_Plant

    !Decomposition of LIT and SOMa by SAP
    N_LITmSAPb = C_LITmSAPb*N_LITm/C_LITm!/CN_ratio(1)
    N_LITsSAPb = C_LITsSAPb*N_LITs/C_LITs!/CN_ratio(2)

    N_LITmSAPf = C_LITmSAPf*N_LITm/C_LITm!/CN_ratio(1)

    N_LITsSAPf = C_LITsSAPf*N_LITs/C_LITs!/CN_ratio(2)
    N_SOMaSAPb = C_SOMaSAPb*N_SOMa/C_SOMa!/CN_ratio(9)
    N_SOMaSAPf = C_SOMaSAPf*N_SOMa/C_SOMa!/CN_ratio(9)

    !Dead mycorrhizal biomass enters SOM pools
    N_EcMSOMp = C_EcMSOMp*(N_EcM/C_EcM)!/CN_ratio(5)
    N_EcMSOMa = C_EcMSOMa*(N_EcM/C_EcM)!/CN_ratio(5)
    N_EcMSOMc = C_EcMSOMc*(N_EcM/C_EcM)!/CN_ratio(5)
    N_ErMSOMp = C_ErMSOMp*(N_ErM/C_ErM)!/CN_ratio(6)
    N_ErMSOMa = C_ErMSOMa*(N_ErM/C_ErM)!/CN_ratio(6)
    N_ErMSOMc = C_ErMSOMc*(N_ErM/C_ErM)!/CN_ratio(6)
    N_AMSOMp = C_AMSOMp*(N_AM/C_AM)!/CN_ratio(7)
    N_AMSOMa = C_AMSOMa*(N_AM/C_AM)!/CN_ratio(7)
    N_AMSOMc = C_AMSOMc*(N_AM/C_AM)!/CN_ratio(7)

    !Dead saphrotroph biomass enters SOM pools
    N_SAPbSOMp = C_SAPbSOMp*N_SAPb/C_SAPb!/CN_ratio(3)
    N_SAPbSOMa = C_SAPbSOMa*N_SAPb/C_SAPb!CN_ratio(3)
    N_SAPbSOMc = C_SAPbSOMc*N_SAPb/C_SAPb!/CN_ratio(3)
    N_SAPfSOMp = C_SAPfSOMp*N_SAPf/C_SAPf!/CN_ratio(4)
    N_SAPfSOMa = C_SAPfSOMa*N_SAPf/C_SAPf!/CN_ratio(4)
    N_SAPfSOMc = C_SAPfSOMc*N_SAPf/C_SAPf!/CN_ratio(4)

    !Desorption of SOMp to SOMa
    N_SOMpSOMa = C_SOMpSOMa*N_SOMp/C_SOMp!/CN_ratio(8)

    !Transport from SOMc to SOMa:
    N_SOMcSOMa = C_SOMcSOMa*N_SOMc/C_SOMc!/CN_ratio(10)

    !"Leftover" N in saprotrophs. Given to inorganic pool to ensure constant C:N ratios:
    f = 0.5                                                                     !NOTE: A fraction, f, of the C made available by myc is decomposed by SAPb, the rest by SAPf
    U_sb = (C_LITmSAPb + C_LITsSAPb + C_SOMaSAPb+f*(Decomp_ecm + Decomp_erm + Decomp_am))    !The saprotrophs decompose the carbon that is made more available when the mycorrhiza take N from SOM.
    U_sf = (C_LITmSAPf + C_LITsSAPf + C_SOMaSAPf+ (1-f)*(Decomp_ecm + Decomp_erm + Decomp_am))
    !print*, "Saprotrophic uptake of C:", U_sb
    N_SAPbIN = N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb - e_s*U_sb/CN_ratio(3)
  !  print*, "NSAPbIN", N_SAPbIN, N_LITmSAPb, N_LITsSAPb, N_SOMaSAPb, e_s*U_sb/CN_ratio(3)
    N_SAPfIN = N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf - e_s*U_sf/CN_ratio(4)
    !If nothing is leftover, nothing is given:
    ! if (N_SAPbIN <= 0) then
    !   N_SAPbIN = 0
    ! end if
    ! if (N_SAPfIN <=0) then
    !   N_SAPfIN = 0
    ! end if

    !All N the Mycorrhiza dont need for its own, it gives to the plant:
    N_EcMPlant = N_INEcM + N_SOMaEcM - e_m*C_PlantEcM/CN_ratio(5) !gN/m3h
    N_ErMPlant = N_INErM + N_SOMaErM - e_m*C_PlantErM/CN_ratio(6)
    N_AMPlant = N_INAM + N_SOMaErM - e_m*C_PlantAM/CN_ratio(7)
    !If not enough N to cover the internal need, nothing is given to the plant
    ! if (N_EcMPlant <= 0) then
    !   N_EcMPlant = 0
    ! end if
    ! if (N_ErMPlant <= 0) then
    !   N_ErMPlant = 0
    ! end if
    ! if (N_AMPlant <= 0) then
    !   N_AMPlant = 0
    ! end if

    nullify( C_SOMp,C_SOMa,C_SOMc,C_EcM,C_ErM,C_AM, C_SAPb,C_SAPf)
  end subroutine calculate_fluxes


  subroutine vertical_diffusion(tot_diffusion_dummy,upper_diffusion_flux,lower_diffusion_flux,pool_matrix,vert) !This subroutine calculates the vertical transport of carbon through the soil layers.
      integer               :: depth, pool
      real(r8),intent(in)   :: pool_matrix(nlevdecomp, pool_types)
      real(r8),intent(out)  :: upper_diffusion_flux, lower_diffusion_flux
      real(r8), intent(out) :: tot_diffusion_dummy ![gC/day]
      real(r8),intent(out)  :: vert(nlevdecomp, pool_types)
      real(r8)              :: sum_day=0.0

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
end module fluxMod2
