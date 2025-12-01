function [IM, sigma] = GMPM_calculation_main(M, R, Fhw, Vs30, Region, Fault_Type_name, gmpm_type, Dip, FVS30, fas, T)
    [rup] = NGA_GMPM_parameter_estimation_approximate_Rjb(M, R, Dip, Vs30, Region, Fault_Type_name, FVS30, fas, Fhw);
    
    if strcmp(gmpm_type,'BSSA_2014')
        [IM, sigma] = BSSA_2014_nga(rup.M, T, rup.Rjb, rup.Fault_Type, rup.Region, rup.Z1o0_CY, rup.Vs30);
    end
    if strcmp(gmpm_type,'CB_2014')
        [IM, sigma] = CB_2014_nga(rup.M, T, rup.Rrup, rup.Rjb, rup.Rx, rup.W, rup.Ztor, rup.Zbot, rup.Dip_degree, rup.lambda_degree, rup.Fhw, rup.Vs30, rup.Z25, rup.Zhyp, rup.Region);
    end
    if strcmp(gmpm_type,'CY_2014')
        [IM, sigma] = CY_2014_nga(rup.M, T, rup.Rrup, rup.Rjb, rup.Rx, rup.Ztor, rup.Dip_degree, rup.lambda_degree, rup.Z1o0_CY, rup.Vs30, rup.Fhw, rup.FVS30, rup.Region);
    end
    if strcmp(gmpm_type,'ASK_2014')
        [IM, sigma] = ASK_2014_nga(rup.M, T, rup.Rrup, rup.Rjb, rup.Rx, rup.Ry, rup.Ztor, rup.Dip_degree, rup.lambda_degree, rup.fas, rup.Fhw,  rup.W, rup.Z1o0_AS, rup.Vs30, rup.FVS30, rup.Region);
    end
    if strcmp(gmpm_type,'TestGMPM')

    end

end