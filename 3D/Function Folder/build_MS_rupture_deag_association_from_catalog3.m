function [MS_rupture_deag_True, idm_cell, idr_cell, ...
          MS_deag_count, MS_unique_deag_count] = ...
    build_MS_rupture_deag_association_from_catalog3( ...
    RUP_MS, M_vec_deagg, R_vec_deagg, ...
    Dip, Lambda, Branch_weight_MS, activity, averaging_method)

% BUILD_MS_RUPTURE_DEAG_ASSOCIATION_FROM_CATALOG3
%
% Version 3: finite-fault / 3D parent-rupture association for the AS analysis.
%
% PURPOSE
% -------
% 1) Assign every original finite MS realization to an M-R deaggregation bin.
%
% 2) Within each fault and M-R bin, collapse original MS realizations that
%    produce the same AS-relevant parent state (actual M + geometry).
%
%    AS-uniqueness key:
%
%       [M_MS, Dip, Lambda, L0, L1, W0, W1]
%
%    The following are deliberately NOT uniqueness variables:
%       - MS_branch_id
%       - M_index
%       - hyp_index
%       - sample_index
%       - Rrup_MS
%
%    Reason:
%    * The actual MS magnitude M_MS is used later for the AS calculation,
%      therefore M_MS is part of the AS-uniqueness key.
%    * Different MS branches or hypocenter samples that lead to the same
%      [M_MS,Dip,Lambda,L0,L1,W0,W1] are AS-equivalent and can be collapsed.
%    * Rrup_MS has already served its role in assigning the realization to
%      the R deaggregation bin.
%
% 3) Provide one normalized inner-parent averaging weight for f_IM_A for
%    every AS-unique parent. Two averaging methods are available:
%
%    averaging_method = 'q_weighted'
%       General occurrence-weighted formulation:
%
%          q_r = w_MS(ii) * activity(ii) *
%                P_M * P_hyp * P_sample
%
%       For all original realizations r that collapse to geometry g:
%
%          q_geometry(g) = sum_{r in g} q_r
%
%       The normalized parent-mixture weight is:
%
%          pi_g = q_geometry(g) / sum_h q_geometry(h)
%
%    averaging_method = 'uniform'
%       Reproduces the original SIMPLE AVERAGE OVER ORIGINAL PARENT
%       REALIZATIONS in the bin.
%
%       Every original catalog realization receives equal weight 1/N.
%       If several original realizations collapse to the same geometry,
%       their multiplicity is retained:
%
%          pi_g = n_original_rows_collapsed(g) / N_original_rows_in_bin
%
%       Therefore collapsing duplicate AS calculations does NOT change the
%       original arithmetic average.
%
% IMPORTANT
% ---------
% Column 18 is the INNER finite-parent weight used to build f_IM_A for one
% fault / one M-R bin / one AS branch / one GMPM. It is NOT the outer
% Figure-2 W_step weight.
%
% Columns 16 and 17 are intentionally retained so the downstream AS helper
% can independently choose how to average N_A:
%   - col 17 -> uniform/multiplicity weighting
%   - col 16 -> q/occurrence weighting
%
% INPUTS
% ------
% RUP_MS            : cell(n_faults,n_branch_MS), output of
%                     place_ruptures_on_fault for the mainshocks.
% M_vec_deagg       : M deaggregation bin edges.
% R_vec_deagg       : Rrup deaggregation bin edges.
% Dip               : n_faults x n_branch_MS matrix [radians].
% Lambda            : n_faults x n_branch_MS rake matrix [radians].
% Branch_weight_MS  : n_faults x n_branch_MS MS logic-tree weights.
% activity          : n_faults x n_branch_MS MS activity/occurrence rates.
% averaging_method  : f_IM_A parent weighting: 'q_weighted' or 'uniform'.
%
% OUTPUT PARENT ROW FORMAT
% ------------------------
% MS_rupture_deag_True{i}{idm,idr} contains ONLY AS-unique geometries.
%
% col 1  = representative_MS_branch_id
%          (bookkeeping/geometry lookup only after branch collapse)
% col 2  = representative M_index
% col 3  = representative hyp_index
% col 4  = representative sample_index
% col 5  = representative original M_MS
% col 6  = representative Rrup_MS
% col 7  = L0
% col 8  = L1
% col 9  = W0
% col 10 = W1
% col 11 = rupture_length
% col 12 = rupture_width
% col 13 = sum_parent_weight_no_activity
%          = sum(P_M * P_hyp * P_sample) over collapsed original rows
% col 14 = Dip [rad]
% col 15 = Lambda [rad]
% col 16 = q_geometry
%          = sum(w_MS * activity * P_M * P_hyp * P_sample)
%            over collapsed original rows
% col 17 = n_original_rows_collapsed
% col 18 = f_IM_A parent_averaging_weight pi_g
%          * q-normalized weight for 'q_weighted'
%          * multiplicity/N_original for 'uniform'
%
% For N_A, the downstream AS helper may independently reconstruct either
% uniform weights from col 17 or q-weighted weights from col 16.
%
% COUNTS
% ------
% MS_deag_count{i}(idm,idr)
%     = number of ORIGINAL catalog realizations assigned to the bin.
%
% MS_unique_deag_count{i}(idm,idr)
%     = number of AS-unique [M_MS,Dip,Lambda,L0,L1,W0,W1] parents after collapse.
%
% NOTE ON COLUMN 13
% -----------------
% place_ruptures_on_fault defines the original catalog column 12 as:
%
%     parent_weight_no_activity = P_M * P_hyp * P_sample
%
% This function preserves that information in output column 13 by summing
% it over all original realizations represented by the retained geometry.
% The GENERAL q-weight including MS branch weight and activity is stored
% separately in column 16.

%% ------------------------------------------------------------------------
% Validate inputs
% -------------------------------------------------------------------------

n_faults    = size(RUP_MS,1);
n_branch_MS = size(RUP_MS,2);

% if ~isequal(size(Dip), [n_faults, n_branch_MS])
%     error('Dip must have size n_faults x n_branch_MS = %d x %d.', ...
%         n_faults, n_branch_MS)
% end
% 
% if ~isequal(size(Branch_weight_MS), [n_faults, n_branch_MS])
%     error(['Branch_weight_MS must have size n_faults x n_branch_MS = ', ...
%            '%d x %d.'], n_faults, n_branch_MS)
% end
% 
% if ~isequal(size(activity), [n_faults, n_branch_MS])
%     error('activity must have size n_faults x n_branch_MS = %d x %d.', ...
%         n_faults, n_branch_MS)
% end
% 
% if ~(ischar(averaging_method) || isstring(averaging_method))
%     error('averaging_method must be ''q_weighted'' or ''uniform''.')
% end

averaging_method = lower(strtrim(char(averaging_method)));

% if ~ismember(averaging_method, {'q_weighted','uniform'})
%     error('Unknown averaging_method. Use ''q_weighted'' or ''uniform''.')
% end
% 
% if any(Branch_weight_MS(:) < 0)
%     error('Branch_weight_MS cannot contain negative values.')
% end
% 
% if any(activity(:) < 0)
%     error('activity cannot contain negative values.')
% end

n_M_deag_mid = length(M_vec_deagg) - 1;
n_R_deag_mid = length(R_vec_deagg) - 1;

%% ------------------------------------------------------------------------
% Preallocate outputs
% -------------------------------------------------------------------------

MS_rupture_deag_True = cell(n_faults,1);
MS_deag_count        = cell(n_faults,1);
MS_unique_deag_count = cell(n_faults,1);

idm_cell = cell(n_faults,1);
idr_cell = cell(n_faults,1);

% Used only to make floating-point geometry comparisons robust.
round_scale = 1e8;

%% ------------------------------------------------------------------------
% Fault loop
% -------------------------------------------------------------------------

for i = 1:n_faults

    % Each temporary original row has 16 columns:
    %  1:13 = [MS_branch_id, rupture_catalog row]
    %  14   = Dip
    %  15   = Lambda
    %  16   = q_r = branch_weight * activity * catalog_parent_weight
    parents_all = cell(n_M_deag_mid, n_R_deag_mid);

    MS_rupture_deag_True{i} = cell(n_M_deag_mid, n_R_deag_mid);
    MS_deag_count{i}        = zeros(n_M_deag_mid, n_R_deag_mid);
    MS_unique_deag_count{i} = zeros(n_M_deag_mid, n_R_deag_mid);

    idm_cell{i} = [];
    idr_cell{i} = [];

    %% --------------------------------------------------------------------
    % Collect ALL original parent realizations from all MS branches
    % ---------------------------------------------------------------------

    for ii = 1:n_branch_MS

        if ~isfield(RUP_MS{i,ii}, 'rupture_catalog')
            error('RUP_MS{%d,%d} does not contain rupture_catalog.', i, ii)
        end

        cat = RUP_MS{i,ii}.rupture_catalog;

        if isempty(cat)
            continue
        end

        if size(cat,2) < 12
            error(['RUP_MS{%d,%d}.rupture_catalog must contain at least ', ...
                   '12 columns.'], i, ii)
        end

        for rr = 1:size(cat,1)

            M_current    = cat(rr,4);
            Rrup_current = cat(rr,5);

            idm = local_find_bin(M_vec_deagg, M_current, 'M');
            idr = local_find_bin(R_vec_deagg, Rrup_current, 'Rrup');

            % Original catalog weight from place_ruptures_on_fault:
            % P_M * P_hyp * P_sample
            parent_weight_no_activity = cat(rr,12);

            % General raw occurrence mass of this original MS realization.
            q_r = Branch_weight_MS(i,ii) * ...
                  activity(i,ii) * ...
                  parent_weight_no_activity;

            % [branch id, 12 catalog columns, Dip, Lambda, q_r]
            parent_original_row = [ ...
                    ii, ...             % 1
                    cat(rr,:), ...      % 2:13
                    Dip(i,ii), ...      % 14
                    Lambda(i,ii), ...   % 15
                    q_r];                % 16

            parents_all{idm,idr} = [ ...
                parents_all{idm,idr}; ...
                parent_original_row]; %#ok<AGROW>

            MS_deag_count{i}(idm,idr) = ...
                MS_deag_count{i}(idm,idr) + 1;

            idm_cell{i} = [idm_cell{i}, idm]; %#ok<AGROW>
            idr_cell{i} = [idr_cell{i}, idr]; %#ok<AGROW>

        end

    end

    idm_cell{i} = unique(idm_cell{i});
    idr_cell{i} = unique(idr_cell{i});

    %% --------------------------------------------------------------------
    % Collapse AS-equivalent parent geometries within each M-R bin
    % ---------------------------------------------------------------------

    for idm = 1:n_M_deag_mid

        for idr = 1:n_R_deag_mid

            parents = parents_all{idm,idr};

            if isempty(parents)
                continue
            end

            n_original_bin = size(parents,1);

            % -------------------------------------------------------------
            % AS-uniqueness key
            % -------------------------------------------------------------
            % parents(:,5)   = actual M_MS
            % parents(:,14)  = Dip
            % parents(:,15)  = Lambda
            % parents(:,7:10)= L0,L1,W0,W1 because:
            %   parent col 1 = branch id
            %   parent cols 2:13 = original rupture_catalog cols 1:12
            % Therefore M_MS is col 5 and L0..W1 are cols 7:10.
            %
            % Actual M_MS is included because each retained parent now uses
            % its own mainshock magnitude in the AS calculation.
            keys = [ ...
                parents(:,5), ...       % actual M_MS
                parents(:,14), ...      % Dip
                parents(:,15), ...      % Lambda
                parents(:,7:10)];       % L0 L1 W0 W1

            keys_round = round(keys * round_scale) / round_scale;

            % Stable: first occurrence is retained as representative.
            [~, ia, ic] = unique(keys_round, 'rows', 'stable');

            n_unique = length(ia);

            % Output has 18 columns, described in the header.
            unique_parents = zeros(n_unique,18);

            for uu = 1:n_unique

                members = find(ic == uu);
                rep     = ia(uu);

                % Representative original row for bookkeeping and geometry.
                unique_parents(uu,1:13) = parents(rep,1:13);

                % Dip is part of the AS-uniqueness key.
                unique_parents(uu,14) = parents(rep,14);

                % Lambda is part of the AS-uniqueness key.
                unique_parents(uu,15) = parents(rep,15);   

                % Sum raw no-activity parent masses represented by geometry.
                unique_parents(uu,13) = sum(parents(members,13));

                % q_geometry = sum of general occurrence masses represented
                % by this AS-unique geometry.
                unique_parents(uu,16) = sum(parents(members,16));

                % Number of original parent realizations collapsed here.
                unique_parents(uu,17) = numel(members);

            end

            % -------------------------------------------------------------
            % Select INNER parent-averaging method
            % -------------------------------------------------------------
            switch averaging_method

                case 'q_weighted'

                q_geometry = unique_parents(:,16);
                q_total = sum(q_geometry);
            
                if q_total > 0
                    unique_parents(:,18) = q_geometry ./ q_total;
                else
                    unique_parents(:,18) = 0;
                end
            
                case 'uniform'
                
                    unique_parents(:,18) = ...
                        unique_parents(:,17) ./ n_original_bin;
            end

            % Numerical check: selected weights should sum to 1 whenever
            % the bin has a meaningful positive weighting basis.
            selected_sum = sum(unique_parents(:,18));

            if strcmp(averaging_method,'uniform')
                if abs(selected_sum - 1) > 1e-12
                    error(['Uniform parent weights do not sum to 1 at ', ...
                           'fault=%d, M-bin=%d, R-bin=%d.'], i, idm, idr)
                end
            else
                if sum(unique_parents(:,16)) > 0 && ...
                        abs(selected_sum - 1) > 1e-12
                    error(['q-weighted parent weights do not sum to 1 at ', ...
                           'fault=%d, M-bin=%d, R-bin=%d.'], i, idm, idr)
                end
            end

            MS_rupture_deag_True{i}{idm,idr} = unique_parents;
            MS_unique_deag_count{i}(idm,idr) = n_unique;

        end

    end

end

%% ------------------------------------------------------------------------
% Summary
% -------------------------------------------------------------------------

fprintf('\nMS rupture association summary - Version 3\n')
fprintf('AS uniqueness key: [M_MS, Dip, Lambda, L0, L1, W0, W1]\n')
fprintf('f_IM_A parent averaging method: %s\n', averaging_method)

n_total_all  = 0;
n_unique_all = 0;

for i = 1:n_faults

    n_total_i  = sum(MS_deag_count{i}, 'all');
    n_unique_i = sum(MS_unique_deag_count{i}, 'all');

    n_total_all  = n_total_all  + n_total_i;
    n_unique_all = n_unique_all + n_unique_i;

    fprintf(['Fault %d: original parent realizations = %d, ', ...
             'AS-unique geometries = %d\n'], ...
             i, n_total_i, n_unique_i)

end

fprintf(['All faults: original parent realizations = %d, ', ...
         'AS-unique geometries = %d\n\n'], ...
         n_total_all, n_unique_all)

end

%% ========================================================================
% Local helper
% ========================================================================
function id = local_find_bin(edges, x, name)

% Match the bin convention used in MSAS_main_PSHA_loop_original:
%
%       id = sum(edges < x)
%
% Thus an interior-edge value belongs to the LOWER bin:
%
%       x = 10, edges = [... 7.5 10 12.5 ...]
%       --> bin [7.5,10]

tol = 1e-10 * max(1, max(abs(edges)));

if x < edges(1)-tol || x > edges(end)+tol
    error('%s value %.6f is outside deag edges [%.6f, %.6f].', ...
        name, x, edges(1), edges(end))
end

id = sum(edges < x);

% Protect the lower boundary
if id == 0 && abs(x-edges(1)) <= tol
    id = 1;
end

% Protect against tiny roundoff above the last edge
if id > length(edges)-1
    id = length(edges)-1;
end

end
