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

% AS calculations for one fault, AS branch and GMPM.
% One distribution per M-R bin, averaged over its MS parents.
% PA_Ex_IM_deag(s,r,mm), f_IM_A_deag(s,r,mm):
% s = AS IM edge/midpoint, r = R bin, idm = idm_active(mm).
% EN_A_deag and n_parent_deag are [R bin x active M bin].
% Parent row cols used here:
% 1 = branch id for geometry lookup, 5 = actual MS magnitude
% 7:10 = L0,L1,W0,W1, 14 = Dip, 15 = Lambda (radians)
% 16 = occurrence mass, 17 = original row count
% 18 = normalized parent weight for f_IM_A
% f_IM_A weights comes from col 18, uniform or q-weighted.
% N_A weights are chosen separately: uniform uses col 17,
% q_weighted uses col 16. No extra N_A weighting of f_IM_A.
% Outer deaggregation weights are applied elsewhere.
% Each parent uses ME = parent(5) for AS counts, magnitudes and scaling.
% AS hypocenters stay in MS plane. Rupture placement can use
% MS plane or full fault, according to AS_rupture_domain.
% Within parent: P(M_A) * P(hyp_A) * P(sample_A).

%% sizes and input checks

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

% IM, R bin, active M bin
PA_Ex_IM_deag = zeros(nIMedge, n_R_deag_mid, nME_active);
f_IM_A_deag   = zeros(nIM,     n_R_deag_mid, nME_active);

% AS count varies with M-R bin; each parent use its own magnitude
% rows = R bins, cols = active M bins
EN_A_deag = zeros(n_R_deag_mid, nME_active);

% unique AS parents in each M-R bin
n_parent_deag = zeros(n_R_deag_mid, nME_active);

log_IM_edges = log(IM_vec_edge(:));
weight_tol = 1e-10;

%% active M bins

for mm = 1:nME_active

    idm = idm_active(mm);

    %% R bins

    for r = 1:n_R_deag_mid

        parents = MS_rupture_deag_i{idm,r};

        if isempty(parents)
            continue
        end

        n_parent = size(parents,1);
        n_parent_deag(r,mm) = n_parent;

        %% parent weights for f_IM_A
        % col 18 already normalized in catalog builder
        parent_weights_fIM = parents(:,18);
        weight_sum_fIM = sum(parent_weights_fIM);

        if weight_sum_fIM <= 0
            % no parent weight here
            continue
        end

        if abs(weight_sum_fIM - 1) > weight_tol
            error(['f_IM_A parent weights do not sum to 1 for active ', ...
                   'M-R bin (idm=%d, r=%d). Sum = %.16g.'], ...
                   idm, r, weight_sum_fIM)
        end

        %% parent weights for N_A
        % N_A weights chosen seperately: uniform or q
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

        % seperate averages for f_IM_A and N_A
        %   f_bin  = sum_g w_fIM,g * f_g
        %   EN_bin = sum_g w_NA,g  * EN_g
        % average PA first, then -diff gives f_bin
        EN_bin = 0;
        PA_mix = zeros(nIMedge,1);

        %% unique MS parents
        % ME comes from parent(5)

        for pp = 1:n_parent

            w_fIM = parent_weights_fIM(pp);
            w_NA  = parent_weights_NA(pp);

            % skip if both weights zero
            if w_fIM == 0 && w_NA == 0
                continue
            end

            parent = parents(pp,:);

            % branch id is for fault geometry lookup
            % uniqueness uses M, Dip, Lambda and rupture bounds
            MS_branch_id = round(parent(1));

            % actual parent M
            ME = parent(5);

            if ME <= Mmin
                error('Parent MS magnitude %.6f must be larger than Mmin %.6f.', ...
                      ME, Mmin)
            end

            % parent rupture limits
            L0 = parent(7);
            L1 = parent(8);
            W0 = parent(9);
            W1 = parent(10);

            % parent dip and rake
            Dip_parent    = parent(14);
            Lambda_parent = parent(15);

            FG_parent = FG_i{MS_branch_id};

            %% AS count for this parent

            EN_parent = local_expected_aftershock_count( ...
                ME, Mmin, a, b, c, p, DT);

            %% AS mags and rupture sizes

            MAG_AS = local_build_AS_magnitude_distribution( ...
                Mmin, ME, b, dma);

            RS_AS = build_rupture_scaling( ...
                Fault_Type_name, ...
                MAG_AS.M_mid_vec, ...
                variation_in_rupture_area, ...
                n_sampling_AS);

            %% AS hypocenters inside MS plane

            HYP_AS = build_hypocenter_grid( ...
                L0, L1, ...
                W0, W1, ...
                dra, Site_Xs, Site_Ys, Fault_Hf_i, Dip_parent);

            %% AS placement limits

            HYP_AS_place = HYP_AS;

            switch lower(AS_rupture_domain)

                case 'ms_rupture'

                    % keep hypocenters and whole AS rupture inside MS plane

                case 'full_fault'

                    % AS hypocenters stay inside MS; rupture may extend over host fault
                    HYP_AS_place.L_start = 0;
                    HYP_AS_place.L_end   = FG_parent.L;
                    HYP_AS_place.W_start = 0;
                    HYP_AS_place.W_end   = FG_parent.W;

                otherwise

                    error(['Unknown AS_rupture_domain: %s. Use ', ...
                        '''MS_rupture'' or ''full_fault''.'], ...
                        AS_rupture_domain);

            end

            %% place AS ruptures

            RUP_AS = place_ruptures_on_fault( ...
                FG_parent, MAG_AS, RS_AS, HYP_AS_place, ...
                rupture_placement_type);

            %% selected GMPM, unique AS inputs

            [uniqueRows_AS, ~, ic_AS] = ...
                unique(RUP_AS.GMPM_input_mat, 'rows');

            [~, sigma_AS, LnIM_AS] = ...
                calculate_GMPM_for_unique_rows( ...
                uniqueRows_AS, 1, {gmpm_type}, ...
                Vs30, Z1o0, Z2o5, Region, Fault_Type_name, ...
                FVS30, fas, HW, T, Dip_parent, Lambda_parent);

            % use column vectors here
            sigma_AS = sigma_AS(:);
            LnIM_AS  = LnIM_AS(:);

            %% P(IM_A > x) for this parent
            % same weights as MS_AS_build_f_IM_A

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

                        % AS realization weight
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

            %% add parent contributions

            EN_bin = EN_bin + w_NA * EN_parent;
            PA_mix = PA_mix + w_fIM .* PA_parent;

        end

        %% save M-R bin results

        EN_A_deag(r,mm) = EN_bin;

        PA_Ex_IM_deag(:,r,mm) = PA_mix;
        f_IM_A_deag(:,r,mm) = -diff(PA_mix);

    end

end

end

%% AS count from actual parent ME
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

%% AS mag distribution, Mmin to ME
function MAG_AS = local_build_AS_magnitude_distribution(Mmin, ME, b, dma)

if ME <= Mmin-0.000000001
    error('ME %.6f must be larger than Mmin %.6f.', ME, Mmin)
end

if dma <= 0
    error('dma must be positive.')
end

% step from Mmin; add ME at end if step misses it
M_vec_A = Mmin + (0:floor((ME - Mmin)/dma)) * dma;

tol = 1e-12 * max(1,abs(ME));

if abs(M_vec_A(end) - ME) > tol
    M_vec_A = [M_vec_A, ME];
end

% truncated Gutenberg-Richter CDF
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
