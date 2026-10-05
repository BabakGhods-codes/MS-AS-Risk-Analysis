function RS = build_rupture_scaling( ...
    Fault_Type, M_mid_vec, ...
    variation_in_rupture_area, n_sampling)

%%% Rupture Properties
warning('Rupture length and downdip width are calculated directly from magnitude-dependent scaling relationships.')

if strcmp(Fault_Type,'Reverse') == 1

    % Rupture area
    a_RA = -3.99; b_RA = 0.98; sigm_RA = 0.26;

    % Downdip rupture width
    a_RW = -1.61; b_RW = 0.41; sigm_RW = 0.15;

    % Subsurface rupture length
    a_RL = -2.42; b_RL = 0.58; sigm_RL = 0.16;

elseif strcmp(Fault_Type,'Strike Slip') == 1

    % Rupture area
    a_RA = -3.42; b_RA = 0.90; sigm_RA = 0.22;

    % Downdip rupture width
    a_RW = -0.76; b_RW = 0.27; sigm_RW = 0.14;

    % Subsurface rupture length
    a_RL = -2.57; b_RL = 0.62; sigm_RL = 0.15;

else

    error('Unknown Fault_Type: %s', Fault_Type)

end

nM = length(M_mid_vec) + 1;

%%% Rupture Properties
for i = 1:nM-1

    if variation_in_rupture_area == 1

        error(['This methodology and the current code cannot consider ' ...
            'rupture-dimension variability - Refer to Monte Carlo code!'])

    elseif variation_in_rupture_area == 0

        % Median magnitude-dependent rupture dimensions
        w_vec_rup(i,:) = ...
            10^(a_RW + b_RW*M_mid_vec(i)) * ones(1,n_sampling);

        L_vec_rup(i,:) = ...
            10^(a_RL + b_RL*M_mid_vec(i)) * ones(1,n_sampling);

        % Area derived from the directly estimated L and W
        A_vec_rup(i,:) = L_vec_rup(i,:) .* w_vec_rup(i,:);

        P_rup_W(i,:) = ones(1,n_sampling);

    else

        error('variation_in_rupture_area must be 0 or 1.')

    end

end

%% Store outputs
RS = struct();

RS.Fault_Type = Fault_Type;

RS.variation_in_rupture_area = variation_in_rupture_area;
RS.n_sampling = n_sampling;

RS.a_RA = a_RA;
RS.b_RA = b_RA;
RS.sigm_RA = sigm_RA;

RS.a_RW = a_RW;
RS.b_RW = b_RW;
RS.sigm_RW = sigm_RW;

RS.a_RL = a_RL;
RS.b_RL = b_RL;
RS.sigm_RL = sigm_RL;

RS.M_mid_vec = M_mid_vec;
RS.nM = nM;

RS.A_vec_rup = A_vec_rup;
RS.w_vec_rup = w_vec_rup;
RS.L_vec_rup = L_vec_rup;
RS.P_rup_W = P_rup_W;

end