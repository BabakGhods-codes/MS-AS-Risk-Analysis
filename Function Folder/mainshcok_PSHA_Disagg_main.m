function [f_IM_E, fx_IM, exceedance_rate_IM, Deag_M_R, LnIM, Deag_M_R_Exceed] = mainshcok_PSHA_Disagg_main(M_mid_vec, R_mid_vec, dIM, IM_vec_edge, IM_vec, P_M, landa, gmpm_type, Vs30, HW, Dip, FVS30, fas,  T, region, Fault_Type)
    for i = 1 : length(M_mid_vec)
        for j = 1 : length(R_mid_vec)
            % Estimating all rupture properties required for the GMPMs
            [rup] = NGA_GMPM_parameter_estimation_approximate_Rjb(M_mid_vec(i), R_mid_vec(j), Dip, Vs30, region, Fault_Type, FVS30, fas, HW);
            if strcmp(gmpm_type,'BSSA_2014')
                [IM(i,j), sigma(i,j)] = BSSA_2014_nga(rup.M, T, rup.Rjb, rup.Fault_Type, rup.Region, rup.Z1o0_CY, rup.Vs30);
            end
            if strcmp(gmpm_type,'CB_2014')
                [IM(i,j), sigma(i,j)] = CB_2014_nga(rup.M, T, rup.Rrup, rup.Rjb, rup.Rx, rup.W, rup.Ztor, rup.Zbot, rup.Dip_degree, rup.lambda_degree, rup.Fhw, rup.Vs30, rup.Z25, rup.Zhyp, rup.Region);
            end
            if strcmp(gmpm_type,'CY_2014')
                [IM(i,j), sigma(i,j)] = CY_2014_nga(rup.M, T, rup.Rrup, rup.Rjb, rup.Rx, rup.Ztor, rup.Dip_degree, rup.lambda_degree, rup.Z1o0_CY, rup.Vs30, rup.Fhw, rup.FVS30, rup.Region);
            end
            if strcmp(gmpm_type,'ASK_2014')
                [IM(i,j), sigma(i,j)] = ASK_2014_nga(rup.M, T, rup.Rrup, rup.Rjb, rup.Rx, rup.Ry, rup.Ztor, rup.Dip_degree, rup.lambda_degree, rup.fas, rup.Fhw,  rup.W, rup.Z1o0_AS, rup.Vs30, rup.FVS30, rup.Region);
            end
            if strcmp(gmpm_type,'TestGMPM')
               IM(i,j)  = -0.152+0.859*M_mid_vec(i)-1.803*log(R_mid_vec(j)+25);
               sigma(i,j) = 0.57;
               IM(i,j) = exp(IM(i,j));
            end
            if strcmp(gmpm_type,'idriss_2013')
                error('Babak: idriss_2013 code is old and should be checked before usage')
                [sa(:,m,jj), sigma(:,m,jj)] = I_2014_nga(M_mid_vec(i), rup.Vs30, knownPer, rup.Rjb_disaggregation(m), rup.F);
            end
        end
    end
    LnIM = log(IM);
    
    % Exccedance rate mainshock (Hazard curve)    
    exceedance_rate_IM = zeros(size(IM_vec_edge));
    f_IM_E = zeros(length(IM_vec),length(R_mid_vec),length(M_mid_vec)); 
    for s = 1 : length(IM_vec_edge)
        for i = 1 : length(M_mid_vec)
            for j = 1 : length(R_mid_vec)
                exceedance_rate_IM(s) = exceedance_rate_IM(s) + landa *  P_M(i) * 1/length(R_mid_vec) * (1 - normcdf( (log(IM_vec_edge(s))-LnIM(i,j))/sigma(i,j) ));
                if s < length(IM_vec_edge)
                    % used in EQ 18 of the Irevalio et al. 2020
                    f_IM_E(s,j,i) = f_IM_E(s,j,i) + (-normcdf( (log(IM_vec_edge(s))-LnIM(i,j))/sigma(i,j)) + normcdf( (log(IM_vec_edge(s+1))-LnIM(i,j))/sigma(i,j)) ); % This is the distribution of mainshocks IM conditioned on the mainshock M and R 
                end
            end
        end
    end
    fx_IM = diff(exceedance_rate_IM)*-1;
    
    for s = 1 : length(IM_vec)
        for i = 1 : length(M_mid_vec)
            for j = 1 : length(R_mid_vec)
                % Disaggregation of the hazard  (Occurence)
                  Deag_M_R(i,j,s) = landa * P_M(i) * 1/length(R_mid_vec) *(-normcdf( (log(IM_vec(s)-dIM/2)-LnIM(i,j))/sigma(i,j)) + normcdf( (log(IM_vec(s)+dIM/2)-LnIM(i,j))/sigma(i,j))) / fx_IM(s);
            end
        end
    end
    
    for s = 1 : length(IM_vec_edge)
        for i = 1 : length(M_mid_vec)
            for j = 1 : length(R_mid_vec)
                % Disaggregation of the hazard (Exceedance)
                Deag_M_R_Exceed(i,j,s) = landa * P_M(i) * 1/length(R_mid_vec) *(1-normcdf( (log(IM_vec_edge(s))-LnIM(i,j))/sigma(i,j))) / exceedance_rate_IM(s);
            end
        end
    end

    
end

