function [MS_rupture_deag_True, idm_cell, idr_cell, ...
          MS_deag_count, MS_unique_deag_count] = ...
    build_MS_rupture_deag_association_from_catalog3( ...
    RUP_MS, M_vec_deagg, R_vec_deagg, ...
    Dip, Lambda, Branch_weight_MS, activity, averaging_method)

% Group finite MS ruptures into M-R bins for AS calculations.
% Within each fault/bin merge parents with same
% [M_MS, Dip, Lambda, L0,L1,W0,W1], rounded before comparing.
% Branch and sample ids are just representative lookups after merging.
% Rrup is used for bin assignment, not for the uniqueness key.
% Parent weights for f_IM_A:
% q_weighted: sum(w_MS * activity * P_M * P_hyp * P_sample), normalize per bin
% uniform: original row multiplicity / total original rows in bin
% So uniform keeps the average over original realizations after merging.
% These are inner parent weights; outer deaggregation weights go elsewhere.
% N_A can use col 16 or 17 independently in the AS function.
% RUP_MS is cell(n_faults,n_branch_MS); M_vec_deagg and R_vec_deagg are edges.
% Dip, Lambda, Branch_weight_MS, activity: n_faults x n_branch_MS.
% Dip and Lambda in radians. averaging_method is uniform or q_weighted.
% MS_rupture_deag_True{i}{idm,idr}: one row per unique parent.
% cols 1:4 = representative branch, M index, hyp index, sample index
% 5 = actual M_MS, 6 = representative Rrup_MS
% 7:10 = L0,L1,W0,W1, 11 = rupture length, 12 = rupture width
% 13 = sum(P_M * P_hyp * P_sample), before branch weights and activity
% 14 = Dip, 15 = Lambda, 16 = summed occurrence mass q_geometry
% 17 = number of original rows merged, 18 = normalized f_IM_A parent weight
% col 13 comes from col 12 of the original rupture catalog.
% MS_deag_count{i}(idm,idr) counts original rows;
% MS_unique_deag_count{i}(idm,idr) counts parents after merging.

%% input sizes

n_faults    = size(RUP_MS,1);
n_branch_MS = size(RUP_MS,2);

% if ~isequal(size(Dip), [n_faults, n_branch_MS])
%     error('Dip must have size n_faults x n_branch_MS = %d x %d.', ...
%         n_faults, n_branch_MS)
% end
% if ~isequal(size(Branch_weight_MS), [n_faults, n_branch_MS])
%     error(['Branch_weight_MS must have size n_faults x n_branch_MS = ', ...
%            '%d x %d.'], n_faults, n_branch_MS)
% end
% if ~isequal(size(activity), [n_faults, n_branch_MS])
%     error('activity must have size n_faults x n_branch_MS = %d x %d.', ...
%         n_faults, n_branch_MS)
% end
% if ~(ischar(averaging_method) || isstring(averaging_method))
%     error('averaging_method must be ''q_weighted'' or ''uniform''.')
% end

averaging_method = lower(strtrim(char(averaging_method)));

% if ~ismember(averaging_method, {'q_weighted','uniform'})
%     error('Unknown averaging_method. Use ''q_weighted'' or ''uniform''.')
% end
% if any(Branch_weight_MS(:) < 0)
%     error('Branch_weight_MS cannot contain negative values.')
% end
% if any(activity(:) < 0)
%     error('activity cannot contain negative values.')
% end

n_M_deag_mid = length(M_vec_deagg) - 1;
n_R_deag_mid = length(R_vec_deagg) - 1;

%% output arrays

MS_rupture_deag_True = cell(n_faults,1);
MS_deag_count        = cell(n_faults,1);
MS_unique_deag_count = cell(n_faults,1);

idm_cell = cell(n_faults,1);
idr_cell = cell(n_faults,1);

% round geometry before comparing
round_scale = 1e8;

%% faults

for i = 1:n_faults

    % temporary rows, 16 cols
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

    %% collect parents from all MS branches

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

            % catalog weight before activity and branch weights
            % P_M * P_hyp * P_sample
            parent_weight_no_activity = cat(rr,12);

            % occurrence mass for this MS realization
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

    %% merge equivalent parents in each M-R bin

    for idm = 1:n_M_deag_mid

        for idr = 1:n_R_deag_mid

            parents = parents_all{idm,idr};

            if isempty(parents)
                continue
            end

            n_original_bin = size(parents,1);

            % key for equivalent AS parents
            % parents(:,5)   = actual M_MS
            % parents(:,14)  = Dip
            % parents(:,15)  = Lambda
            % parents(:,7:10) = L0,L1,W0,W1
            % include actual M, it controls the AS calculation
            keys = [ ...
                parents(:,5), ...       % actual M_MS
                parents(:,14), ...      % Dip
                parents(:,15), ...      % Lambda
                parents(:,7:10)];       % L0 L1 W0 W1

            keys_round = round(keys * round_scale) / round_scale;

            % keep first row as representative
            [~, ia, ic] = unique(keys_round, 'rows', 'stable');

            n_unique = length(ia);

            % 18 cols, see notes above
            unique_parents = zeros(n_unique,18);

            for uu = 1:n_unique

                members = find(ic == uu);
                rep     = ia(uu);

                % take representative row
                unique_parents(uu,1:13) = parents(rep,1:13);

                % dip
                unique_parents(uu,14) = parents(rep,14);

                % rake
                unique_parents(uu,15) = parents(rep,15);   

                % sum weights before activity
                unique_parents(uu,13) = sum(parents(members,13));

                % sum occurrence masses for this geometry
                unique_parents(uu,16) = sum(parents(members,16));

                % how many original rows merged here
                unique_parents(uu,17) = numel(members);

            end

            % parent averaging option
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

            % weights should sum to 1 when total weight is positive
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

%% print counts

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

%% bin lookup
function id = local_find_bin(edges, x, name)

% same bin convention as PSHA
%       id = sum(edges < x)
% interior edge goes in lower bin
%       x = 10, edges = [... 7.5 10 12.5 ...]
%       --> bin [7.5,10]

tol = 1e-10 * max(1, max(abs(edges)));

if x < edges(1)-tol || x > edges(end)+tol
    error('%s value %.6f is outside deag edges [%.6f, %.6f].', ...
        name, x, edges(1), edges(end))
end

id = sum(edges < x);

% include first edge
if id == 0 && abs(x-edges(1)) <= tol
    id = 1;
end

% roundoff near last edge
if id > length(edges)-1
    id = length(edges)-1;
end

end
