function [PA_Ex_IM, f_IM_A_deag, f_IM_A, EN_A] = aftershock_sequence_calculation_main( ...
    M_mid_vec, epicentral_location_info, DT, ...
    a, b, c, p, Mmin, IM_vec, IM_vec_edge, dma, ...
    gmpm_type, Vs30, HW, T, region, Fault_Type_name, ...
    dra, Dip, FVS30, fas, H, AS_Spatial_Distribution, S)

% AFTERSHOCK_SEQUENCE_CALCULATION_MAIN
%
% Outputs:
%   PA_Ex_IM       : Exceedance probability of IM caused by aftershocks (per IM bin edge, unique epicenter or R, MS magnitude bin)
%   f_IM_A_deag    : Aftershock IM PDF aggregated over deaggregation bins
%   f_IM_A         : Aftershock IM PDF conditioned on unique epicenter or R and MS magnitude 
%   EN_A           : Expected number of aftershocks (Omori parameters)
%
% Notes:
% - Uses truncated Gutenberg–Richter in base-10 for aftershock magnitudes.
% - Distance fields built over rectangular/square meshes depending on fault type or AS_Spatial_Distribution.
warning('R_mid_vec_A is calculated by the formula proposed by Utsu (1970). Refer to Irevalino short note.');

    %% ------------------------------ Initial unpacking ------------------------------
    epicentral_location_deag   = epicentral_location_info(1, :);
    weights_for_R_deag         = epicentral_location_info(2:end, :);

    epicentral_location_unique = unique(epicentral_location_deag);
    nIM        = numel(IM_vec);
    nIMedge    = numel(IM_vec_edge);
    nME        = numel(M_mid_vec);
    nR_unique  = numel(epicentral_location_unique);

    % Distribution of aftershock IM conditioned on MS 
    f_IM_A   = zeros(nIM,        nR_unique, nME);

    % Exceedance probability of IM caused by aftershocks conditioned on MS rupture
    PA_Ex_IM = zeros(nIMedge,    nR_unique, nME);

    %% ------------------------------ Main loop over MS magnitude and epicentral loc ------------------------------
    for m = 1:nME
        ME = M_mid_vec(m);                                         % Mid value for MS magnitude bin

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

        % Midpoints per magnitude bin & bin probabilities (finite difference of CDF)
        M_mid_vec_A = 0.5 * (M_vec_A(2:end) + M_vec_A(1:end-1));
        P_M_A       = FMA(2:end) - FMA(1:end-1);

        nMA = numel(M_mid_vec_A);

        for n = 1:nR_unique
            % Loop over aftershock magnitudes
            for i = 1:nMA
                % Fault (rupture) area scaling
                SA = 10^(M_mid_vec_A(i) - 4.1);

                % ------------------------------ Distance mesh ------------------------------
                if AS_Spatial_Distribution == 1
                    % Rectangular region around epicenter (aspect ratio depends on fault type)
                    if strcmp(Fault_Type_name, 'Strike-Slip')
                        % 1:3 aspect (width: length)
                        ds      = sqrt(SA/3);
                        xa_mesh = unique([ -1.5*ds:dra:0, 0, (-1.5*ds:dra:0)*-1 ]);
                        if numel(xa_mesh) > 11
                            xa_mesh = unique([ linspace(-1.5*ds,0,6), linspace(0,1.5*ds,6) ]);
                        end
                        ya_mesh = unique([ -ds/2:dra:0, 0, (-ds/2:dra:0)*-1 ]);
                        if numel(ya_mesh) > 11
                            ya_mesh = unique([ linspace(-ds/2,0,6), linspace(0,ds/2,6) ]);
                        end

                    elseif strcmp(Fault_Type_name, 'Reverse')
                        % 1:2 aspect
                        ds      = sqrt(SA/2);
                        xa_mesh = unique([ -ds:dra:0, 0, (-ds:dra:0)*-1 ]);
                        if numel(xa_mesh) > 11
                            xa_mesh = unique([ linspace(-ds,0,6), linspace(0,ds,6) ]);
                        end
                        ya_mesh = unique([ -ds/2:dra:0, 0, (-ds/2:dra:0)*-1 ]);
                        if numel(ya_mesh) > 11
                            ya_mesh = unique([ linspace(-ds/2,0,6), linspace(0,ds/2,6) ]);
                        end
                    else
                        % If other fault types appear, current logic implies no change.
                        % (Kept as-is; no math changes.)
                    end

                    R1     = S - (xa_mesh + epicentral_location_unique(1, n));
                    R2     = H - ya_mesh;
                    R_mat  = sqrt( repmat(R1.^2, numel(R2), 1) + repmat((R2').^2, 1, numel(R1)) );
                    R_mid_vec_A = reshape(R_mat, 1, []);

                elseif AS_Spatial_Distribution == 0
                    % Square region around epicenter (1:1 aspect)
                    ds      = sqrt(SA);
                    xa_mesh = unique([ -ds/2:dra:0, 0, (-ds/2:dra:0)*-1 ]);
                    if numel(xa_mesh) > 11
                        xa_mesh = unique([ linspace(-ds/2,0,6), linspace(0,ds/2,6) ]);
                    end
                    ya_mesh = xa_mesh;

                    R1     = S - (xa_mesh + epicentral_location_unique(1, n));
                    R2     = H - ya_mesh;
                    R_mat  = sqrt( repmat(R1.^2, numel(R2), 1) + repmat((R2').^2, 1, numel(R1)) );
                    R_mid_vec_A = reshape(R_mat, 1, []);
                else
                    % (Kept as-is; no alternative spatial distributions.)
                    R_mid_vec_A = []; 
                    error('Unsupported AS_Spatial_Distribution value: %g', AS_Spatial_Distribution);
                end
                % --------------------------------------------------------------------------

                % ------------------------------ Exceedance accumulation ------------------------------
                nR = numel(R_mid_vec_A);
                for s = 1:nIMedge
                    for j = 1:nR
                        [IM_val, sigma] = GMPM_calculation_main( ...
                            M_mid_vec_A(i), R_mid_vec_A(j), HW, Vs30, region, ...
                            Fault_Type_name, gmpm_type, Dip, FVS30, fas, T);

                        lnIM = log(IM_val);
                        z    = (log(IM_vec_edge(s)) - lnIM) / sigma;
                        PA_Ex_IM(s, n, m) = PA_Ex_IM(s, n, m) + P_M_A(i) * (1 / nR) * (1 - normcdf(z));
                    end
                end
                % --------------------------------------------------------------------------
            end

            % Convert exceedance curve over edges into PDF on IM bins
            f_IM_A(:, n, m) = -diff(PA_Ex_IM(:, n, m));
        end
    end

    %% ------------------------------ f_IM_A over deaggregation bins (Eq. 13) ------------------------------
    nDeagRows = size(weights_for_R_deag, 1);
    nDeagCols = size(weights_for_R_deag, 2);

    f_IM_A_deag = zeros(nIM, nDeagRows, nME);   % used in Eq. 13 (Iervolino et al., 2020)

    for m = 1:nME
        for n = 1:nDeagCols
            id_epicentral = (epicentral_location_unique == epicentral_location_deag(n));
            idr           = (weights_for_R_deag(:, n) == 1);
            if any(idr)
                sum_weights = sum(weights_for_R_deag(idr,:));  % sum across columns 
                f_IM_A_deag(:, idr, m) = f_IM_A_deag(:, idr, m) + f_IM_A(:, id_epicentral, m) ./ sum_weights';
            end
        end
    end

    %% ------------------------------ PS|ME,RE (Eq. 20) and EN_A ------------------------------
    EN_A    = zeros(size(M_mid_vec));

    for m = 1:nME
        % Expected number of aftershocks (Omori-type)
        EN_A(m) = ((p - 1)^-1) * (10^(a + b*(M_mid_vec(m) - Mmin)) - 10^a) * (c^(1 - p) - (DT + c)^(1 - p));
    end
end
