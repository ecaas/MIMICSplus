!Change equations:

dDIN = N_SAPbDIN + N_SAPfDIN + Deposition - N_DINEcM - N_DINErM - N_DINAM - N_DINplant-Leaching

dPlant_N = N_EcMPlant + N_ErMPlant + N_AMPlant + N_DINplant - N_PlantLITm -N_PlantLITs
dPlant_C = C_growth_rate - C_PlantEcM - C_PlantErM - C_PlantAM  - C_PlantLITm -C_PlantLITs

dLITm_N = N_PlantLITm - N_LITmSAPb - N_LITmSAPf
dLITm_C = C_PlantLITm - C_LITmSAPb - C_LITmSAPf

dLITs_N = N_PlantLITs - N_LITsSAPb - N_LITsSAPf
dLITs_C = C_PlantLITs - C_LITsSAPb - C_LITsSAPf

dSAPb_N = N_LITmSAPb + N_LITsSAPb + N_SOMaSAPb - N_SAPbSOMp - N_SAPbSOMa - N_SAPbSOMc - N_SAPbDIN
dSAPb_C = C_LITmSAPb + C_LITsSAPb + C_SOMaSAPb - C_SAPbSOMp - C_SAPbSOMa - C_SAPbSOMc

dSAPf_N = N_LITmSAPf + N_LITsSAPf + N_SOMaSAPf - N_SAPfSOMp - N_SAPfSOMa - N_SAPfSOMc -N_SAPfDIN
dSAPf_C = C_LITmSAPf + C_LITsSAPf + C_SOMaSAPf - C_SAPfSOMp - C_SAPfSOMa - C_SAPfSOMc

dEcM_N = N_DINEcM + N_SOMaEcM - N_EcMPlant - N_EcMSOMa - N_EcMSOMc - N_EcMSOMp
dEcM_C = C_PlantEcM - C_EcMSOMp - C_EcMSOMa - C_EcMSOMc

dErM_N = N_DINErM + N_SOMaErM - N_ErMPlant - N_ErMSOMa - N_ErMSOMc - N_ErMSOMp
dErM_C = C_PlantErM - C_ErMSOMp - C_ErMSOMa - C_ErMSOMc

dAM_N = N_DINAM + N_SOMaAM -N_AMPlant - N_AMSOMa - N_AMSOMc - N_AMSOMp
dAM_C = C_PlantAM - C_AMSOMp - C_AMSOMa - C_AMSOMc

dSOMc_N = N_SAPbSOMc + N_SAPfSOMc + N_EcMSOMc + N_ErMSOMc + N_AMSOMc -N_SOMcSOMa
dSOMc_C = C_SAPbSOMc + C_SAPfSOMc + C_EcMSOMc + C_ErMSOMc + C_AMSOMc - C_SOMcSOMa

dSOMa_N = N_SAPbSOMa + N_SAPfSOMa + N_EcMSOMa + N_ErMSOMa + N_AMSOMa - N_SOMaSAPb - N_SOMaSAPf - N_SOMaEcM - N_SOMaErM - N_SOMaAM + N_SOMpSOMa + N_SOMcSOMa
dSOMa_C = C_SAPbSOMa + C_SAPfSOMa + C_EcMSOMa + C_ErMSOMa + C_AMSOMa +  C_SOMcSOMa + C_SOMpSOMa - C_SOMaSAPf - C_SOMaSAPb

dSOMp_N = N_SAPbSOMp + N_SAPfSOMp + N_EcMSOMp + N_ErMSOMp + N_AMSOMp - N_SOMpSOMa
dSOMp_C = C_SAPbSOMp + C_SAPfSOMp + C_EcMSOMp + C_ErMSOMp + C_AMSOMp - C_SOMpSOMa
