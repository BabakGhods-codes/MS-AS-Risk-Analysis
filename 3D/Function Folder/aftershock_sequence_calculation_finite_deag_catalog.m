function [PA_Ex_IM_deag, f_IM_A_deag, EN_A_deag, n_parent_deag] = ...
    aftershock_sequence_calculation_finite_deag_catalog( ...
    MS_rupture_deag_i, idm_active, M_vec_deagg_mid, ...
    DT, a, b, c, p, Mmin, ...
    IM_vec, IM_vec_edge, dma, gmpm_type, ...
    Vs30, Z1o0, Z2o5, Region, Fault_Type_name, ...
    FVS30, fas, HW, T, ...
    Dip_i, Lambda_i, FG_i, Fault_Hf_i, ...
    Site_Xs, Site_Ys, dra, ...
    variation_in_rupture_area, n_sampling_AS, rupture_placement_type, AS_rupture_domain, ...
    NA_parent_averaging_method)

% AFTERSHOCK_SEQUENCE_CALCULATION_FINITE_DEAG_CATALOG
%
% Finite-fault replacement for aftershock_sequence_calculation_main in the
% deaggregation-based MS-AS formulation. This version evaluates ONE GMPM
% per function call, preserving the original outer loop:
%
%   GMPM -> fault -> AS branch
%
% -------------------------------------------------------------------------
% CONDITIONING / OUTPUTS
% -------------------------------------------------------------------------
% For one fault, one AS branch, and one GMPM, the function calculates:
%
%   PA_Ex_IM_deag(s,r,mm)
%   f_IM_A_deag(s,r,mm)
%
% where:
%   s  = AS IM edge/midpoint index
%   r  = MS R-deaggregation bin
%   mm = active MS M-bin order, with idm = idm_active(mm)
%
% Thus there is ONE f_IM_A distribution for each M-R deaggregation bin for
% the selected fault + AS branch + GMPM.
%
% Several finite MS parent ruptures may be associated with one M-R bin.
% Their normalized INNER parent weights are stored in column 18 of
% MS_rupture_deag_i{idm,r}.
%
% Relevant parent-row columns:
%   col  1 = representative MS branch id (geometry lookup only)
%   col  5 = actual MS magnitude M_E
%   col  7 = L0
%   col  8 = L1
%   col  9 = W0
%   col 10 = W1
%   col 14 = Dip [rad]
%   col 15 = Lambda [rad]
%   col 18 = normalized inner parent averaging weight, pi_g
%
% Column 18 may represent either q-weighted averaging or the original
% uniform averaging, depending on the association-catalog option. The
% Figure-2 outer weight is NOT applied here.
%
% -------------------------------------------------------------------------
% ACTUAL PARENT MAINSHOCK MAGNITUDE AND N_A
% -------------------------------------------------------------------------
% Every retained parent uses its own actual mainshock magnitude:
%
%       ME = parent(5)
%
% Therefore N_A, the AS magnitude distribution, AS rupture scaling, and the
% parent-specific AS intensity distribution are all evaluated using that ME.
%
% IMPORTANT: f_IM_A and N_A may use DIFFERENT parent-averaging rules:
%
%   * f_IM_A uses column 18, which is built using averaging_method in
%     build_MS_rupture_deag_association_from_catalog3 ('uniform' or
%     'q_weighted').
%
%   * N_A uses NA_parent_averaging_method selected independently here:
%       'uniform'    -> multiplicity-based weights from column 17
%       'q_weighted' -> occurrence-mass weights from column 16
%
% Thus f_IM_A is NOT additionally weighted by N_A.  The two parent averages
% are formed independently, preserving the same averaging philosophy as the
% original finite-parent f_IM_A calculation.
%
% -------------------------------------------------------------------------
% FINITE AS MODEL
% -------------------------------------------------------------------------
% For each retained AS-unique MS parent geometry:
%   1) Build AS magnitude distribution on [Mmin, ME].
%   2) Build finite AS rupture scaling.
%   3) Build a uniform AS hypocenter grid inside [L0,L1] x [W0,W1].
%   4) Place finite AS ruptures inside the MS parent plane.
%   5) Calculate the selected GMPM for unique AS rupture-input rows.
%   6) Construct parent-specific PA_Ex_IM and f_IM_A using the same
%      probability structure as MS_AS_build_f_IM_A:
%
%        P(M_A) * P(hyp_A) * P(sample_A)
%
%   7) Average parent-specific f_IM_A using column 18.
%   8) Average parent-specific N_A independently using the selected
%      NA_parent_averaging_method.
%
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% Basic sizes and input checks
% -------------------------------------------------------------------------

nIM     = length(IM_vec);
nIMedge = length(IM_vec_edge);

n_M_deag_mid = size(MS_rupture_deag_i,1);
n_R_deag_mid = size(MS_rupture_deag_i,2);
nME_active   = length(idm_active);

if length(M_vec_deagg_mid) ~= n_M_deag_mid
    error(['M_vec_deagg_mid length must equal the number of M rows in ', ...
           'MS_rupture_deag_i.'])
end

if length(Dip_i) ~= numel(FG_i) || length(Lambda_i) ~= numel(FG_i)
    error('Dip_i, Lambda_i, and FG_i must have consistent MS-branch sizes.')
end

if ~iscell(FG_i)
    error('FG_i must be a cell array containing the MS-branch fault geometries.')
end

if ~(ischar(gmpm_type) || isstring(gmpm_type))
    error('gmpm_type must be a character vector or string.')
end

gmpm_type = char(gmpm_type);

if ~(ischar(NA_parent_averaging_method) || isstring(NA_parent_averaging_method))
    error('NA_parent_averaging_method must be ''uniform'' or ''q_weighted''.')
end

NA_parent_averaging_method = lower(strtrim(char(NA_parent_averaging_method)));

if ~ismember(NA_parent_averaging_method, {'uniform','q_weighted'})
    error('Unknown NA_parent_averaging_method. Use ''uniform'' or ''q_weighted''.')
end

% Output dimensions remain compatible with the deaggregation workflow.
PA_Ex_IM_deag = zeros(nIMedge, n_R_deag_mid, nME_active);
f_IM_A_deag   = zeros(nIM,     n_R_deag_mid, nME_active);

% Expected number of ASs is now M-R dependent because each retained parent
% uses its own actual MS magnitude. Orientation: [R-bin x active-M-bin].
EN_A_deag = zeros(n_R_deag_mid, nME_active);

% Number of retained AS-unique parents in each M-R bin.
n_parent_deag = zeros(n_R_deag_mid, nME_active);

log_IM_edges = log(IM_vec_edge(:));
weight_tol = 1e-10;

%% ------------------------------------------------------------------------
% Loop over active MS M-deaggregation bins
% -------------------------------------------------------------------------

for mm = 1:nME_active

    idm = idm_active(mm);

    %% --------------------------------------------------------------------
    % Loop over MS R-deaggregation bins
    % ---------------------------------------------------------------------

    for r = 1:n_R_deag_mid

        parents = MS_rupture_deag_i{idm,r};

        if isempty(parents)
            continue
        end

        n_parent = size(parents,1);
        n_parent_deag(r,mm) = n_parent;

        %% ----------------------------------------------------------------
        % Parent averaging weights for f_IM_A
        % -----------------------------------------------------------------
        % Column 18 is already normalized according to averaging_method used
        % in build_MS_rupture_deag_association_from_catalog3.
        parent_weights_fIM = parents(:,18);
        weight_sum_fIM = sum(parent_weights_fIM);

        if weight_sum_fIM <= 0
            % No meaningful parent mixture is available for f_IM_A.
            continue
        end

        if abs(weight_sum_fIM - 1) > weight_tol
            error(['f_IM_A parent weights do not sum to 1 for active ', ...
                   'M-R bin (idm=%d, r=%d). Sum = %.16g.'], ...
                   idm, r, weight_sum_fIM)
        end

        %% ----------------------------------------------------------------
        % Parent averaging weights for N_A
        % -----------------------------------------------------------------
        % These are selected independently so N_A can use either the same
        % uniform-parent philosophy or occurrence/q weighting.
        switch NA_parent_averaging_method

            case 'uniform'

                multiplicity = parents(:,17);
                mult_total = sum(multiplicity);

                if mult_total > 0
                    parent_weights_NA = multiplicity ./ mult_total;
                else
                    parent_weights_NA = zeros(n_parent,1);
                end

            case 'q_weighted'

                q_geometry = parents(:,16);
                q_total = sum(q_geometry);

                if q_total > 0
                    parent_weights_NA = q_geometry ./ q_total;
                else
                    parent_weights_NA = zeros(n_parent,1);
                end
        end

        weight_sum_NA = sum(parent_weights_NA);

        if weight_sum_NA > 0 && abs(weight_sum_NA - 1) > weight_tol
            error(['N_A parent weights do not sum to 1 for active ', ...
                   'M-R bin (idm=%d, r=%d). Sum = %.16g.'], ...
                   idm, r, weight_sum_NA)
        end

        % f_IM_A and N_A are averaged independently:
        %
        %   f_bin  = sum_g w_fIM,g * f_g
        %   EN_bin = sum_g w_NA,g  * EN_g
        %
        % Because differentiation is linear, PA can be averaged first and
        % f_bin obtained from -diff(PA_bin).
        EN_bin = 0;
        PA_mix = zeros(nIMedge,1);

        %% ----------------------------------------------------------------
        % Loop over AS-unique finite MS parents
        % Each parent uses its own ACTUAL MS magnitude ME = parent(5).
        % -----------------------------------------------------------------

        for pp = 1:n_parent

            w_fIM = parent_weights_fIM(pp);
            w_NA  = parent_weights_NA(pp);

            % Skip only if this parent contributes to neither average.
            if w_fIM == 0 && w_NA == 0
                continue
            end

            parent = parents(pp,:);

            % Representative branch is bookkeeping for the pre-built full
            % fault geometry. Branch id itself is not an AS-uniqueness
            % variable once M, Dip, Lambda and geometry are in the key.
            MS_branch_id = round(parent(1));

            % Actual MS magnitude retained in this parent row.
            ME = parent(5);

            if ME <= Mmin
                error('Parent MS magnitude %.6f must be larger than Mmin %.6f.', ...
                      ME, Mmin)
            end

            % Final finite MS parent-plane bounds.
            L0 = parent(7);
            L1 = parent(8);
            W0 = parent(9);
            W1 = parent(10);

            % Dip and rake are explicit AS-uniqueness variables.
            Dip_parent    = parent(14);
            Lambda_parent = parent(15);

            FG_parent = FG_i{MS_branch_id};

            %% ------------------------------------------------------------
            % Parent-specific expected number of aftershocks N_A(ME)
            % -------------------------------------------------------------

            EN_parent = local_expected_aftershock_count( ...
                ME, Mmin, a, b, c, p, DT);

            %% ------------------------------------------------------------
            % Parent-specific AS magnitude distribution and rupture scaling
            % -------------------------------------------------------------

            MAG_AS = local_build_AS_magnitude_distribution( ...
                Mmin, ME, b, dma);

            RS_AS = build_rupture_scaling( ...
                Fault_Type_name, ...
                MAG_AS.M_mid_vec, ...
                variation_in_rupture_area, ...
                n_sampling_AS);

            %% ------------------------------------------------------------
            % AS hypocenters inside this finite MS parent rupture plane
            % -------------------------------------------------------------

            HYP_AS = build_hypocenter_grid( ...
                L0, L1, ...
                W0, W1, ...
                dra, Site_Xs, Site_Ys, Fault_Hf_i, Dip_parent);

            %% ------------------------------------------------------------
            % AS rupture-placement domain
            % -------------------------------------------------------------

            HYP_AS_place = HYP_AS;

            switch lower(AS_rupture_domain)

                case 'ms_rupture'

                    % AS hypocenter and complete AS rupture are constrained
                    % to the finite MS rupture plane. No change is required.

                case 'full_fault'

                    % AS hypocenter remains inside the MS rupture, but the
                    % finite AS rupture may extend to the full host fault.
                    HYP_AS_place.L_start = 0;
                    HYP_AS_place.L_end   = FG_parent.L;
                    HYP_AS_place.W_start = 0;
                    HYP_AS_place.W_end   = FG_parent.W;

                otherwise

                    error(['Unknown AS_rupture_domain: %s. Use ', ...
                        '''MS_rupture'' or ''full_fault''.'], ...
                        AS_rupture_domain);

            end

            %% ------------------------------------------------------------
            % Place finite AS ruptures
            % -------------------------------------------------------------

            RUP_AS = place_ruptures_on_fault( ...
                FG_parent, MAG_AS, RS_AS, HYP_AS_place, ...
                rupture_placement_type);

            %% ------------------------------------------------------------
            % Calculate selected GMPM only for unique AS input rows
            % -------------------------------------------------------------

            [uniqueRows_AS, ~, ic_AS] = ...
                unique(RUP_AS.GMPM_input_mat, 'rows');

            [~, sigma_AS, LnIM_AS] = ...
                calculate_GMPM_for_unique_rows( ...
                uniqueRows_AS, 1, {gmpm_type}, ...
                Vs30, Z1o0, Z2o5, Region, Fault_Type_name, ...
                FVS30, fas, HW, T, Dip_parent, Lambda_parent);

            % Keep one-GMPM arrays unambiguously as column vectors.
            sigma_AS = sigma_AS(:);
            LnIM_AS  = LnIM_AS(:);

            %% ------------------------------------------------------------
            % Build P(IM_A > x | this finite parent)
            % Same probability structure as MS_AS_build_f_IM_A.
            % -------------------------------------------------------------

            PA_parent = zeros(nIMedge,1);

            nAS_M        = length(MAG_AS.M_mid_vec);
            nAS_hyp      = HYP_AS.n_hyp;
            nAS_sampling = RS_AS.n_sampling;

            expected_n_rows = nAS_M * nAS_hyp * nAS_sampling;

            if length(ic_AS) ~= expected_n_rows
                error(['The finite AS rupture catalog contains %d GMPM ', ...
                       'rows, but %d were expected from nM*nHyp*nSampling. ', ...
                       'The current f_IM_A indexing assumes every AS ', ...
                       'magnitude-hypocenter-sample combination is valid.'], ...
                       length(ic_AS), expected_n_rows)
            end

            for ma = 1:nAS_M

                for ha = 1:nAS_hyp

                    for ka = 1:nAS_sampling

                        id_local = ...
                            ka + ...
                            (ha-1)*nAS_sampling + ...
                            (ma-1)*nAS_sampling*nAS_hyp;

                        id_unique = ic_AS(id_local);

                        % Probability of this generic AS realization:
                        % P(M_A) * P(hyp_A) * P(sample_A)
                        weight_AS = ...
                            MAG_AS.P_M(ma) * ...
                            (1/nAS_hyp) * ...
                            RS_AS.P_rup_W(ma,ka);

                        z = (log_IM_edges - LnIM_AS(id_unique)) ./ ...
                            sigma_AS(id_unique);

                        PA_parent = PA_parent + ...
                            weight_AS .* (1 - normcdf(z));

                    end

                end

            end

            %% ------------------------------------------------------------
            % Independent parent averages
            % -------------------------------------------------------------

            EN_bin = EN_bin + w_NA * EN_parent;
            PA_mix = PA_mix + w_fIM .* PA_parent;

        end

        %% ----------------------------------------------------------------
        % Final equivalent N_A and generic-AS distribution for this M-R bin
        % -----------------------------------------------------------------

        EN_A_deag(r,mm) = EN_bin;

        PA_Ex_IM_deag(:,r,mm) = PA_mix;
        f_IM_A_deag(:,r,mm) = -diff(PA_mix);

    end

end

end

%% ========================================================================
% Local helper: expected number of aftershocks for one actual parent M_E
% ========================================================================
function EN = local_expected_aftershock_count(ME, Mmin, a, b, c, p, DT)

if p ~= 1
    EN = ((p - 1)^-1) * ...
        (10^(a + b*(ME - Mmin)) - 10^a) * ...
        (c^(1 - p) - (DT + c)^(1 - p));
else
    term1 = 10^(a + b*(ME - Mmin)) - 10^a;
    term2 = log((DT + c) / c);
    EN = term1 * term2;
end

end

%% ========================================================================
% Local helper: truncated AS magnitude distribution on [Mmin, ME]
% ========================================================================
function MAG_AS = local_build_AS_magnitude_distribution(Mmin, ME, b, dma)

if ME <= Mmin-0.000000001
    error('ME %.6f must be larger than Mmin %.6f.', ME, Mmin)
end

if dma <= 0
    error('dma must be positive.')
end

% Uniform increments from Mmin, with ME explicitly appended if dma does not
% land exactly on the actual parent mainshock magnitude.
M_vec_A = Mmin + (0:floor((ME - Mmin)/dma)) * dma;

tol = 1e-12 * max(1,abs(ME));

if abs(M_vec_A(end) - ME) > tol
    M_vec_A = [M_vec_A, ME];
end

% Truncated Gutenberg-Richter CDF in the same form used in the existing
% finite MS-AS / MS_AS_build_f_IM_A workflow.
denom = 10.^(-b*Mmin) - 10.^(-b*ME);

if denom <= 0
    error('Invalid AS magnitude-distribution denominator.')
end

FMA = (10.^(-b*Mmin) - 10.^(-b*M_vec_A)) ./ denom;

M_mid_vec_A = 0.5 * (M_vec_A(2:end) + M_vec_A(1:end-1));
P_M_A       = FMA(2:end) - FMA(1:end-1);

MAG_AS = struct();

MAG_AS.Mmin = Mmin;
MAG_AS.Mmax = ME;
MAG_AS.b_M  = b;
MAG_AS.dm   = dma;

MAG_AS.M_vec     = M_vec_A;
MAG_AS.FM        = FMA;
MAG_AS.M_mid_vec = M_mid_vec_A;
MAG_AS.P_M       = P_M_A;
MAG_AS.nM        = length(M_vec_A);

end
