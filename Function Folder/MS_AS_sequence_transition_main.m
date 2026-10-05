function [PS_E, PS2_onlymain, PAS_MeRe, PA_ij, f_IM_A, EN_A] = MS_AS_sequence_transition_main(IM_vec, epicentral_location_info, M_mid_vec, M_vec, Mmin, dma, b, a, c, p, IM_vec_edge, Mu_mat, Beta_mat, DT,  f_IM_E, landa, P_M, gmpm_type, Vs30, HW, T, region, Fault_Type_name, dra, Dip, FVS30, fas, H, AS_Spatial_Distribution,S) 
    % Construction of f(IM_A|ME,RE), Eq13 of Irevalino et al. 2020
    warning('R_mid_vec_A is calculated by the formula proposed by Utsu, 1970 (Refer to Irevalino short note) ')
    f_IM_A = zeros(length(IM_vec),length(epicentral_location_info(1,:)),length(M_vec)-1); % used in EQ 13 of the Irevalino et al. 2020
    PA_Ex_IM = zeros(length(IM_vec_edge),length(epicentral_location_info(1,:)),length(M_vec)-1);
    for m = 1 : length(M_vec)-1
        for n = 1 : length(epicentral_location_info(1,:))
            
            ME = (M_vec(m+1)+M_vec(m))/2;
            % Uniform aftershock magnitude grid up to ME
            M_vec_A = Mmin + (0:floor((ME - Mmin)/dma)) * dma;
            tol     = 1e-12 * ME;
            if abs(M_vec_A(end) - ME) > tol
                M_vec_A = [M_vec_A, ME];                               % Append ME if not hit exactly
            end
    
            % Truncated GR (base-10) CDF over the grid
            % FMA(k) = P(M <= M_vec_A(k) | M <= ME, M >= Mmin)
            denom = (10.^(-b*Mmin) - 10.^(-b*ME));
            FMA   = (10.^(-b*Mmin) - 10.^(-b*M_vec_A)) ./ denom;

            % --- Midpoints for binning ---
            M_mid_vec_A = (M_vec_A(2:end) + M_vec_A(1:end-1)) / 2;

            % --- Probability per magnitude bin (finite difference of CDF) ---
            P_M_A = FMA(2:end) - FMA(1:end-1);

            % Loop over aftershocks M and R
            for s = 1 : length(IM_vec_edge)
                for i = 1 : length(M_mid_vec_A)
                    SA = 10^(M_mid_vec_A(i)-4.1);
                    
                    %%% Start of the Calculation of Distance %%%
                    if AS_Spatial_Distribution == 1 % Rectangular around the epicenter with 2 to 1 aspect ratio. The larger length is along the fault 
                        if strcmp(Fault_Type_name, 'Strike-Slip') == 1
                            ds = sqrt(SA/3); % ratio 1:3
                            xa_mesh = unique([-1.5*ds:dra:0, 0 , (-1.5*ds:dra:0)*-1]);
                            if length(xa_mesh) > 11
                                xa_mesh = unique([linspace(-1.5*ds,0,6), linspace(0,1.5*ds,6)]); % along the fault length direction
                            end
                            ya_mesh = unique([-ds/2:dra:0, 0 , (-ds/2:dra:0)*-1]);
                            if length(ya_mesh) > 11
                                ya_mesh = unique([linspace(-ds/2,0,6), linspace(0,ds/2,6)]); 
                            end
                        elseif strcmp(Fault_Type_name, 'Reverse') == 1
                            ds = sqrt(SA/2); % ratio 1:2
                            xa_mesh = unique([-ds:dra:0, 0 , (-ds:dra:0)*-1]);
                            if length(xa_mesh) > 11
                                xa_mesh = unique([linspace(-ds,0,6), linspace(0,ds,6)]); % along the fault length direction
                            end
                            ya_mesh = unique([-ds/2:dra:0, 0 , (-ds/2:dra:0)*-1]);
                            if length(ya_mesh) > 11
                                ya_mesh = unique([linspace(-ds/2,0,6), linspace(0,ds/2,6)]); 
                            end
                        end
                        R1 = S-(xa_mesh + epicentral_location_info(1,n));
                        R2 = H-ya_mesh ;
                        R_mat = (repmat(R1.^2,length(R2),1) + repmat((R2').^2,1,length(R1)) ).^0.5;
                        R_mid_vec_A = reshape(R_mat,1,[]);

                    elseif AS_Spatial_Distribution == 0  % Squared around the epicenter
                        ds = SA^0.5;
                        xa_mesh = unique([-ds/2:dra:0, 0 , (-ds/2:dra:0)*-1]);
                        if length(xa_mesh) > 11
                            xa_mesh = unique([linspace(-ds/2,0,6), linspace(0,ds/2,6)]); % along the fault length direction
                        end
                        ya_mesh = xa_mesh;
                        R1 = (xa_mesh + epicentral_location_info(1,n))-S ;
                        R2 = H-ya_mesh ;
                        R_mat = (repmat(R1.^2,length(R2),1) + repmat((R2').^2,1,length(R1)) ).^0.5;
                        R_mid_vec_A = reshape(R_mat,1,[]);
                    end
                    
                    % AS Hazard 
                    for j = 1 : length(R_mid_vec_A)      
                        [IM, sigma] =  GMPM_calculation_main(M_mid_vec_A(i), R_mid_vec_A(j), HW, Vs30, region, Fault_Type_name, gmpm_type, Dip, FVS30, fas, T);
                        lnIM = log(IM);
                        PA_Ex_IM(s,n,m) = PA_Ex_IM(s,n,m) + P_M_A(i) * 1/length(R_mid_vec_A) * (1 - normcdf( (log(IM_vec_edge(s))-lnIM)/sigma)); 
                    end
                end
            end
            f_IM_A(:,n,m) = diff(PA_Ex_IM(:,n,m))*-1; % This is the distribution of aftershocks IM conditioned on the mainshock M and R --- used in EQ 13 of the Irevalino et al. 2020
        end
    end

    %% P_A|Me,RE: Eqs 15 and 14 --- The transition matrix given the occurrence of one generic aftershock from a sequence of a mainshock {ME=mE, RE=rE}
    Nd = length(Mu_mat);
    PA_ij=repmat(zeros(Nd+1,Nd+1),1,1,length(epicentral_location_info(1,:)),length(M_mid_vec)); 
    for m = 1 : length(M_mid_vec)
        for n = 1 : length(epicentral_location_info(1,:))
            for i = 1 : Nd
                for j = i:Nd
                    PA_ij(i,j+1,n,m) = PDs_func(Mu_mat,Beta_mat,i,j,IM_vec,f_IM_A(:,n,m));
                end
                PA_ij(i,i,n,m) = 1 - sum(PA_ij(i,i+1:end,n,m));
            end  
            PA_ij(end,end,n,m) = 1;
        end
    end

    %% P_E|Me,RE: Eq 19 & 18 --- The transition matrix due to the occurrence of the {ME=mE, RE=rE} earthquake.
    PE_ij=eye(Nd+1,Nd+1);
    for m = 1 : length(M_mid_vec)
        for n = 1 : length(epicentral_location_info(1,:))
            for i = 1 : Nd
                for j = i:Nd
                    PE_ij(i,j+1,n,m) = PDs_func(Mu_mat,Beta_mat,i,j,IM_vec,f_IM_E(:,n,m));
                end
                PE_ij(i,i,n,m) = 1 - sum(PE_ij(i,i+1:end,n,m));
            end
            PE_ij(end,end,n,m) = 1;
        end
    end

    %% Ps|Me,RE Eq 20 without PE_ij (presented as PHI_E) or its equivalent presented by Shokrabadi borton (presented as PS_MeRe) --- This is transition matrix for several possible aftershocks
    PAS_MeRe = ones(size(PA_ij));
    EN_A = zeros(length(epicentral_location_info(1,:)),length(M_mid_vec));
    for m = 1 : length(M_mid_vec)
        for n = 1 : length(epicentral_location_info(1,:))
   
        if p ~= 1
            EN_A(n,m) = ((p-1)^-1) * (10^(a+b*(M_mid_vec(m)-Mmin))-10^a) * (c^(1-p) - (DT+c)^(1-p));
            PAS_MeRe(:,:,n,m) = PA_ij(:,:,n,m)^EN_A(n,m);
        else
            term1 = 10^(a + b * (M_mid_vec(m) - Mmin)) - 10^a;
            term2 = log((DT + c) / c);  % natural logarithm (ln)
            EN_A(n,m) = term1 * term2;
        end

        end

    end

    %% PS_E and Ps Eq 21 and 22 : I mixed the equation 21 and 22! and elimated the Eq 22. --- The unit time damage transition probability matrix for seismic sequences [PS]
    PS_E = zeros(Nd+1,Nd+1);
    PS2_onlymain = zeros(Nd+1,Nd+1);
    for m = 1 : length(M_mid_vec)
        for n = 1 : length(epicentral_location_info(1,:))
            PS_E =         PS_E         + landa* P_M(m) * 1/length(epicentral_location_info(1,:)) * PE_ij(:,:,n,m) * PAS_MeRe(:,:,n,m); % * PS_MeRe(:,:,n,m);
            PS2_onlymain = PS2_onlymain + landa* P_M(m) * 1/length(epicentral_location_info(1,:)) * PE_ij(:,:,n,m) ;
        end
    end
end




