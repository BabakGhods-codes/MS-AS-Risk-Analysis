function RUP = place_ruptures_on_fault(FG, MAG, RS, HYP, rupture_placement_type)

%% Unpack geometry variables
W = FG.W;
L = FG.L;
Hf = FG.Hf;
Dip = FG.Dip;

points_fault_loc = FG.points_fault_loc;
points_loc = FG.points_loc;

Rrup_points = FG.Rrup_points;
Rjb_points = FG.Rjb_points;
Ztor_points = FG.Ztor_points;
Rx_points = FG.Rx_points;
Cases_points = FG.Cases_points;
Ry_points = FG.Ry_points;
alpha_points = FG.alpha_points;

%% Unpack magnitude variables
nM = MAG.nM;
M_mid_vec = MAG.M_mid_vec;

%% Unpack rupture variables
n_sampling = RS.n_sampling;
w_vec_rup = RS.w_vec_rup;
L_vec_rup = RS.L_vec_rup;

%% Unpack hypocenter variables
n_hyp = HYP.n_hyp;
hyp_points_loc = HYP.hyp_points_loc;

% These limits allow the same function to be used for:
% mainshock: full fault plane, [0,L] and [0,W]
% aftershock: mainshock rupture plane, [L_MS_0,L_MS_1] and [W_MS_0,W_MS_1]
L_domain_start = HYP.L_start;
L_domain_end   = HYP.L_end;
W_domain_start = HYP.W_start;
W_domain_end   = HYP.W_end;

L_domain = L_domain_end - L_domain_start;
W_domain = W_domain_end - W_domain_start;

%%% Ruptures Location on The Fault and the "Distance Paramters" or "R"s of
%%% each rupture
rup_lim_W = zeros(n_hyp, 2, n_sampling, nM-1);
rup_lim_L = zeros(n_hyp, 2, n_sampling, nM-1);

GMPM_input_mat = [];

Rrup = zeros(nM-1, n_hyp, n_sampling); 
Rjb = zeros(nM-1, n_hyp, n_sampling); 
Ztor = zeros(nM-1, n_hyp, n_sampling); 
Rx = zeros(nM-1, n_hyp, n_sampling);  
Cases = zeros(nM-1, n_hyp, n_sampling); 
Ry_test = zeros(nM-1, n_hyp, n_sampling); 
Ry = zeros(nM-1, n_hyp, n_sampling); 
Zhyp = zeros(nM-1, n_hyp, n_sampling); 
Alpha = zeros(nM-1, n_hyp, n_sampling);
error_Ry = zeros(nM-1, n_hyp, n_sampling);

valid_rupture_case = zeros(nM-1, n_hyp, n_sampling);
% Catalog of all valid ruptures:
% columns:
% 1  = M_index
% 2  = hyp_index
% 3  = sample_index
% 4  = M_mid
% 5  = Rrup
% 6  = L0
% 7  = L1
% 8  = W0
% 9  = W1
% 10 = rupture_length
% 11 = rupture_width
% 12 = parent_weight_no_activity = P_M * P_hyp * P_sample
max_catalog_rows = (nM-1) * n_hyp * n_sampling;
rupture_catalog = nan(max_catalog_rows, 12);
rupture_catalog_count = 0;
rupture_row_id = zeros(nM-1, n_hyp, n_sampling);


for i = 1 : nM-1

    for j = 1 : n_hyp

        for ii = 1 : n_sampling

            %% ============================================================
            %  Rupture placement option
            % =============================================================

            if strcmp(rupture_placement_type, 'original_uniform_hypocenter_clip')

                % ========================================================
                % YOUR ORIGINAL METHOD
                % Uniform hypocenter over the selected domain.
                % For MS, the selected domain is the full fault.
                % For AS, the selected domain can be the MS rupture plane.
                % If rupture exceeds the selected boundary, shift/clip it.
                % ========================================================

                % Rupture placement in width of the selected domain
                if w_vec_rup(i,ii) >= W_domain

                    rup_lim_W(j,:,ii,i) = [W_domain_start, W_domain_end];

                elseif hyp_points_loc(j,2)-w_vec_rup(i,ii)/2 >= W_domain_start  &&  ...
                       hyp_points_loc(j,2)+w_vec_rup(i,ii)/2 <= W_domain_end

                    rup_lim_W(j,:,ii,i) = [ ...
                        hyp_points_loc(j,2)-w_vec_rup(i,ii)/2, ...
                        hyp_points_loc(j,2)+w_vec_rup(i,ii)/2];

                elseif hyp_points_loc(j,2)-w_vec_rup(i,ii)/2 < W_domain_start

                    rup_lim_W(j,:,ii,i) = [ ...
                        W_domain_start, ...
                        W_domain_start + w_vec_rup(i,ii)];

                elseif hyp_points_loc(j,2)+w_vec_rup(i,ii)/2 > W_domain_end

                    rup_lim_W(j,:,ii,i) = [ ...
                        W_domain_end - w_vec_rup(i,ii), ...
                        W_domain_end];

                else

                    error('Something is wrong !!!')

                end

                % Rupture placement in length of the selected domain
                if L_vec_rup(i,ii) >= L_domain

                    rup_lim_L(j,:,ii,i) = [L_domain_start, L_domain_end];

                elseif hyp_points_loc(j,1)-L_vec_rup(i,ii)/2 >= L_domain_start  &&  ...
                       hyp_points_loc(j,1)+L_vec_rup(i,ii)/2 <= L_domain_end

                    rup_lim_L(j,:,ii,i) = [ ...
                        hyp_points_loc(j,1)-L_vec_rup(i,ii)/2, ...
                        hyp_points_loc(j,1)+L_vec_rup(i,ii)/2];

                elseif hyp_points_loc(j,1)-L_vec_rup(i,ii)/2 < L_domain_start

                    rup_lim_L(j,:,ii,i) = [ ...
                        L_domain_start, ...
                        L_domain_start + L_vec_rup(i,ii)];

                elseif hyp_points_loc(j,1)+L_vec_rup(i,ii)/2 > L_domain_end

                    rup_lim_L(j,:,ii,i) = [ ...
                        L_domain_end - L_vec_rup(i,ii), ...
                        L_domain_end];

                else

                    error('Something is wrong !!!')

                end

            elseif strcmp(rupture_placement_type, 'admissible_uniform_hypocenter')

                % ========================================================
                % Optional method
                % Uniform hypocenter, but only over the admissible region.
                % Invalid hypocenters are skipped.
                % ========================================================

                % Check width admissibility
                if w_vec_rup(i,ii) >= W_domain

                    rup_lim_W(j,:,ii,i) = [W_domain_start, W_domain_end];

                elseif hyp_points_loc(j,2)-w_vec_rup(i,ii)/2 >= W_domain_start && ...
                       hyp_points_loc(j,2)+w_vec_rup(i,ii)/2 <= W_domain_end

                    rup_lim_W(j,:,ii,i) = [ ...
                        hyp_points_loc(j,2)-w_vec_rup(i,ii)/2, ...
                        hyp_points_loc(j,2)+w_vec_rup(i,ii)/2];

                else

                    continue

                end

                % Check length admissibility
                if L_vec_rup(i,ii) >= L_domain

                    rup_lim_L(j,:,ii,i) = [L_domain_start, L_domain_end];

                elseif hyp_points_loc(j,1)-L_vec_rup(i,ii)/2 >= L_domain_start && ...
                       hyp_points_loc(j,1)+L_vec_rup(i,ii)/2 <= L_domain_end

                    rup_lim_L(j,:,ii,i) = [ ...
                        hyp_points_loc(j,1)-L_vec_rup(i,ii)/2, ...
                        hyp_points_loc(j,1)+L_vec_rup(i,ii)/2];

                else

                    continue

                end

            else

                error('Unknown rupture_placement_type.')

            end

            valid_rupture_case(i,j,ii) = 1;

            %% ============================================================
            %  The following part is your original distance/GMPM-input block
            % =============================================================

            % Identify the "grid points" associated with the rupture in the
            % "relative (Fault plane)" and "Global" coordinates
            id_rupture_surf = ...
                (points_fault_loc(:,1) <= rup_lim_L(j,2,ii,i) & ...
                 points_fault_loc(:,1) >= rup_lim_L(j,1,ii,i) & ...
                 points_fault_loc(:,2) <= rup_lim_W(j,2,ii,i) & ...
                 points_fault_loc(:,2) >= rup_lim_W(j,1,ii,i));

            Rrup_rupture = Rrup_points(id_rupture_surf);
            Rrup(i,j,ii) = min(Rrup_rupture);

            Rjb_rupture = Rjb_points(id_rupture_surf);
            Rjb(i,j,ii) = min(Rjb_rupture);

            Ztor_rupture = Ztor_points(id_rupture_surf);
            Ztor(i,j,ii) = min(abs(Ztor_rupture));

            Rx_rupture = Rx_points(id_rupture_surf);
            Cases_rupture = Cases_points(id_rupture_surf);

            % Ry_rupture = Ry_points(id_rupture_surf);
            % Ry(i,j,ii) = max(abs(Ry_rupture));

            L_site = FG.Site1_loc(1) - FG.Xs;

            L0 = rup_lim_L(j,1,ii,i);
            L1 = rup_lim_L(j,2,ii,i);

            if L_site < L0
                Ry(i,j,ii) = L0 - L_site;
            elseif L_site > L1
                Ry(i,j,ii) = L_site - L1;
            else
                Ry(i,j,ii) = 0;
            end

            Zhyp(i,j,ii) = hyp_points_loc(j,2)*sin(Dip) + Hf;

            % Identify the "line points" associated with the top edge of
            % rupture and calculate Rx, Cases, alpha, and Ry
            rupture_points_loc = points_loc(id_rupture_surf,:);

            z_rupture = abs(rupture_points_loc(:,3));
            top_z_point = min(z_rupture);

            id_top_z_rupture = (z_rupture == top_z_point);

            Rx_top_z_points = Rx_rupture(id_top_z_rupture);
            Cases_top_z_points = Cases_rupture(id_top_z_rupture);

            R_project_top_z = Rjb_rupture(id_top_z_rupture,:);

            [~, id_Rx] = min(R_project_top_z);

            Rx(i,j,ii) = Rx_top_z_points(id_Rx);
            Cases(i,j,ii) = Cases_top_z_points(id_Rx);
            % Alpha(i,j,ii) = alpha_points(id_Rx);

            Alpha_rupture = alpha_points(id_rupture_surf);
            Alpha_top_z_points = Alpha_rupture(id_top_z_rupture);

            Alpha(i,j,ii) = Alpha_top_z_points(id_Rx);


            if Alpha(i,j,ii) == pi/2 || Alpha(i,j,ii) == -pi/2

                Ry_test(i,j,ii) = 0;

            elseif Alpha(i,j,ii) == 0 || Alpha(i,j,ii) == pi || Alpha(i,j,ii) == -pi

                Ry_test(i,j,ii) = Rjb(i,j,ii);

            else

                Ry_test(i,j,ii) = abs(Rx(i,j,ii) * cot(Alpha(i,j,ii)));

            end

            error_Ry(i,j,ii) = abs((Ry_test(i,j,ii)-Ry(i,j,ii))/Ry(i,j,ii));
            
            GMPM_input_mat = [ ...
                GMPM_input_mat; ...
                M_mid_vec(i), ...
                Rrup(i,j,ii), ...
                Rjb(i,j,ii), ...
                Rx(i,j,ii), ...
                Ry(i,j,ii), ...
                Ztor(i,j,ii), ...
                Zhyp(i,j,ii), ...
                w_vec_rup(i,ii)];


            rupture_catalog_count = rupture_catalog_count + 1;

            rupture_row_id(i,j,ii) = rupture_catalog_count;

            if isfield(MAG, 'P_M')
                P_M_current = MAG.P_M(i);
            else
                P_M_current = NaN;
            end

            if isfield(RS, 'P_rup_W')
                P_sample_current = RS.P_rup_W(i,ii);
            else
                P_sample_current = NaN;
            end

            P_hyp_current = 1 / n_hyp;

            rupture_catalog(rupture_catalog_count,:) = [ ...
                i, ...                                  % 1  M_index
                j, ...                                  % 2  hyp_index
                ii, ...                                 % 3  sample_index
                M_mid_vec(i), ...                       % 4  M_mid
                Rrup(i,j,ii), ...                       % 5  Rrup
                rup_lim_L(j,1,ii,i), ...                % 6  L0
                rup_lim_L(j,2,ii,i), ...                % 7  L1
                rup_lim_W(j,1,ii,i), ...                % 8  W0
                rup_lim_W(j,2,ii,i), ...                % 9  W1
                L_vec_rup(i,ii), ...                    % 10 rupture length
                w_vec_rup(i,ii), ...                    % 11 rupture width
                P_M_current * P_hyp_current * P_sample_current]; % 12 parent weight

        end

    end

end

min_abs_Rx = min(min(abs(Rx)));
max_abs_Rx = max(max(abs(Rx)));
max_error_Ry = min(min(error_Ry));

Check_for_accurate_Rx_Plus = max(abs(Rx_points+points_loc(:,2)'));      
Check_for_accurate_Rx_Minus = max(abs(Rx_points-points_loc(:,2)'));

%% Store outputs
RUP = struct();

RUP.rup_lim_W = rup_lim_W;
RUP.rup_lim_L = rup_lim_L;

RUP.GMPM_input_mat = GMPM_input_mat;

RUP.Rrup = Rrup;
RUP.Rjb = Rjb;
RUP.Ztor = Ztor;
RUP.Rx = Rx;
RUP.Cases = Cases;
RUP.Ry_test = Ry_test;
RUP.Ry = Ry;
RUP.Zhyp = Zhyp;
RUP.Alpha = Alpha;
RUP.error_Ry = error_Ry;

RUP.valid_rupture_case = valid_rupture_case;

RUP.rupture_catalog = rupture_catalog(1:rupture_catalog_count,:);

RUP.rupture_catalog_colnames = { ...
    'M_index', ...
    'hyp_index', ...
    'sample_index', ...
    'M_mid', ...
    'Rrup', ...
    'L0', ...
    'L1', ...
    'W0', ...
    'W1', ...
    'rupture_length', ...
    'rupture_width', ...
    'parent_weight_no_activity'};

RUP.rupture_row_id = rupture_row_id;

RUP.min_abs_Rx = min_abs_Rx;
RUP.max_abs_Rx = max_abs_Rx;
RUP.max_error_Ry = max_error_Ry;

RUP.Check_for_accurate_Rx_Plus = Check_for_accurate_Rx_Plus;
RUP.Check_for_accurate_Rx_Minus = Check_for_accurate_Rx_Minus;

end