clc; clear; tic;

addpath('Function Folder');

warning('BSSA_2014 GMPM: Z1.0 from CY model is used (as in original).');
saving_folder = 'Outputs/Final8/LogicTree_4story_Site1_Dr2o5Dm0o1_IM4_Lin';
%% -----------------------------------------------------------------------
%  INPUTS
% ------------------------------------------------------------------------

Mmax = [6.6   6.6   7  7;   7.0 7.0 7.4 7.4];
% [6.5   6.5   7  7;   6.9 6.9 7.3 7.3]; %
b_M  = [1.15 1  1.15 1;    0.9 0.8 0.9 0.8];

% -----------------------------
% Fault constants
% -----------------------------
activity  = [0.4 0.4 0.4 0.4; 0.2 0.2 0.2 0.2];
% Fault_L = [60; 65];                % Fault lengths (km)
% Fault_W =  [15; 25];                % Fault lengths (km)
Fault_L = [60; 85];                % Fault lengths (km)
Fault_W =  [15; 30];                % Fault lengths (km)
Fault_Hf = [0; 0];
Region = 'California';
Fault_Type_name = {'Strike Slip'; 'Reverse'}; % 
Dip     = [pi/2 pi/2 pi/2 pi/2;   pi/3 pi/3 pi/3 pi/3];
Lambda =  [0 0 0 0; pi/2 pi/2 pi/2 pi/2];
% Global coordinates of the start of the upper edge of each fault
% Fault_X0 = [40; 70]; % [40; 40];%
% Fault_Y0 = [30; 60]; %[30; 30];%
Fault_X0 = [40; 70]; % [40; 40];%
Fault_Y0 = [30; 50]; %[30; 30];%

% Aftershock rupture spatial model
%
% 'MS_rupture' :
%     AS hypocenters are inside the MS rupture and the complete
%     AS rupture is also constrained to remain inside the MS rupture.
%
% 'full_fault' :
%     AS hypocenters are inside the MS rupture, but the finite AS
%     rupture may extend outside the MS rupture and is constrained
%     only by the full host fault boundaries.
AS_rupture_domain = 'full_fault';

% -----------------------------
% Aftershock Omori-type params
% -----------------------------
a = [-1.6 -1.4; -1.8 -1.6];   % [-1.6 -1.6 -1.4 -1.4;  -1.8 -1.8 -1.6 -1.6];
b_AS = [1.15 1.00;
        0.90 0.80]; 
c = [0.04 0.06; 0.02 0.03]; %[ 0.04 0.04 0.06 0.06;   0.02 0.02 0.03 0.03];
p = [1.03 1.15; 1.05 1.2];  %[ 1.03 1.03 1.15 1.15;   1.05 1.05 1.2  1.2];

% -----------------------------
% Site
% -----------------------------
FVS30 = 1;         % FVS30 = 0 for Vs30 is inferred from geology; FVS30 = 1 for measured  Vs30
fas_MS = 0;  % Flag for mainshocks  GMPMs
fas_AS = 1;  % Flag for aftershocks  GMPMs
Site_number = 1;   % Site 1 or 2 Properties    
if Site_number == 1
    HW     = [0; 0];                   % Hanging wall indicators

    % Global site coordinate
    Site_loc = [70,40,0];

    Vs30  = 540;   % m/s
    Z1o0  = 0.2;   % km
    Z2o5  = 1.0;   % km
    FVS30 = 1;

elseif Site_number == 2
    HW     = [0; 1];                   % Hanging wall indicators
    Site_loc = [110,80,0];

    Vs30  = 360;   % m/s
    Z1o0  = 0.2;   % km
    Z2o5  = 1.0;   % km
    FVS30 = 1;                     % m/s
   
end


% Spatial distribution flag for AS distribution flag around the MS epicenter 
% AS_Spatial_Distribution = [0 1 0 1; 0 1 0 1]; % 0 for squared and 1 for rectangulat (1:3 for Strike-Slip and 1:2 for Reverse faults)

% -----------------------------
% Logic-tree weights
% -----------------------------
Branch_weight_MS =  [0.2 0.2 0.3 0.3; 0.2 0.2 0.3 0.3];
Branch_weight_AS =  [0.6 0.4 ; 0.6 0.4]; % [0.15 0.35 0.15 0.35; 0.15 0.35 0.15 0.35];

% -----------------------------
% Structure
% -----------------------------
Structure_number = 2;   %  1 is the 8-story building and 2 is the 4 story building
if Structure_number == 1
    T     = 2.16;      % sec (8-story)
    % 8-story state-by-state fragility curves
    Mu_mat = [0.1320 0.3460 0.3910 0.4510 0.5030;
              0      0.3410 0.3900 0.4380 0.4920;
              0      0      0.3840 0.4240 0.4480;
              0      0      0      0.3440 0.4000;
              0      0      0      0      0.2270];   
    Beta_mat = [0.5070 0.4030 0.4030 0.4090 0.3870;
                0      0.4160 0.4160 0.4160 0.3990;
                0      0      0.4230 0.4230 0.4230;
                0      0      0      0.4390 0.4390;
                0      0      0      0      0.5030];

elseif Structure_number == 2
    T     = 1.57;      % sec (4-story)
    Mu_mat = [0.1490 0.3860 0.5180 0.5750 0.7320
              0      0.3860 0.5050 0.5720 0.7260
              0      0      0.4360 0.5120 0.6370
              0      0      0      0.3470 0.4980
              0      0      0      0      0.3200];   
    Beta_mat = [ 0.3310 0.3810 0.3640 0.3640 0.3640
                 0      0.4040 0.3670 0.3670 0.3670
                 0      0      0.4290 0.4360 0.4360
                 0      0      0      0.5830 0.5220
                 0      0      0      0      0.5690];
end
Nd = size(Mu_mat,1);  % number of damage states

% -----------------------------
% GMPMs (weights uniform)
% -----------------------------
GMPM_types  = {'CY_2014','BSSA_2014','CB_2014','ASK_2014','TestGMPM'};
n_GMPM      = 4;                         % using first 4 (as in your run)
GMPM_weight = ones(1,n_GMPM)/n_GMPM;

% -----------------------------
% PSHA / IM grid
% -----------------------------
dIM         = 0.1;
if Structure_number == 1
    IM_vec_edge = 0:dIM:3.5;       % 8-story as-is (0–3.5g) 
elseif Structure_number == 2
    IM_vec_edge = 0:dIM:4;         % 4-story as-is (0–4g) 
else
    IM_vec_edge = 0:dIM:4; 
end
IM_vec      = (IM_vec_edge(1:end-1)+IM_vec_edge(2:end))/2;
n_IM        = length(IM_vec);

% -----------------------------
% MS magnitude–distance mesh
% -----------------------------
dm        = 0.1;     % mag step
Mmin      = 5;       % min mag

% -----------------------------
% MS-AS Analysis Paramters and Methods
% -----------------------------
n_sampling = 1 ;
variation_in_rupture_area = 0;
rupture_placement_type = 'original_uniform_hypocenter_clip';

% -----------------------------
% Aftershock window and mesh
% -----------------------------
DT  = 60;            % days window
dma = dm;            % AS mag step (kept equal to MS)

% Derived sizes
n_faults     = 2 ;% size(Mmax,1);
n_branch_MS  = size(Branch_weight_MS,2);
n_branch_AS  = size(Branch_weight_AS,2);

% -----------------------------
% Aftershock Analysis Type
% -----------------------------
AS_Analysis = 'PARALLEL'; % SERIAL or PARALLEL
n_Core_Parallel = 16; % used only if AS_Analysis = 'PARALLEL'

dL_hyp  = 2.5;   % MS hypocenter spacing
dL_geom = 0.5;%0.5;   % resolution used to calculate rupture distances
dra     = 2.5;   % AS hypocenter spacing

dR_deagg = 2.5; %0.02;
dM_deagg = 0.1;
% Parent weighting used to average f_IM_A within each M-R bin:
parent_averaging_method = 'uniform'; % 'uniform' or 'q_weighted'

% Parent weighting used to average N_A within each M-R bin:
% This is intentionally independent of the f_IM_A option above.
NA_parent_averaging_method = 'uniform'; % 'uniform' or 'q_weighted'

%% -----------------------------------------------------------------------
%  MS: finite-fault rupture generation + PSHA + branch deaggregation
% ------------------------------------------------------------------------

% Containers
FG      = cell(n_faults, n_branch_MS);
MAG_MS  = cell(n_faults, n_branch_MS);
RS_MS   = cell(n_faults, n_branch_MS);
HYP_MS  = cell(n_faults, n_branch_MS);
RUP_MS  = cell(n_faults, n_branch_MS);

uniqueRows_MS = cell(n_faults, n_branch_MS);
ic_MS         = cell(n_faults, n_branch_MS);

IM_MS    = cell(n_faults, n_branch_MS);
sigma_MS = cell(n_faults, n_branch_MS);
LnIM_MS  = cell(n_faults, n_branch_MS);

PSHA = cell(n_faults, n_branch_MS);

M_vec       = cell(n_faults, n_branch_MS);
FM          = cell(n_faults, n_branch_MS);
P_M         = cell(n_faults, n_branch_MS);
M_mid_total = cell(n_faults, n_branch_MS);

f_IM_E              = cell(n_faults, n_branch_MS, n_GMPM);
fx_IM               = cell(n_faults, n_branch_MS, n_GMPM);
exceedance_rate_IM  = cell(n_faults, n_branch_MS, n_GMPM);
Deag_M_R_GMPM_fault = cell(n_faults, n_branch_MS, n_GMPM);
Deag_M_R_Exceed     = cell(n_faults, n_branch_MS, n_GMPM);

min_r_vec = inf(n_faults, n_branch_MS);
max_r_vec = zeros(n_faults, n_branch_MS);

%% -----------------------------------------------------------------------
% First pass: build MS finite ruptures and find global R limits
% ------------------------------------------------------------------------

for i = 1:n_faults


    for ii = 1:n_branch_MS


        % Fault geometry and site-distance grid
        FG{i,ii} = build_fault_geometry( ...
            Fault_Type_name{i}, ...
            Fault_W(i), Fault_L(i), Fault_Hf(i), Dip(i,ii), dL_geom, ...
            Fault_X0(i), Fault_Y0(i), Site_loc);

        % Mainshock magnitude distribution
        MAG_MS{i,ii} = build_magnitude_distribution( ...
            Mmin, Mmax(i,ii), b_M(i,ii), dm);

        % Mainshock rupture scaling
        RS_MS{i,ii} = build_rupture_scaling( ...
            Fault_Type_name{i}, ...
            MAG_MS{i,ii}.M_mid_vec, ...
            variation_in_rupture_area, ...
            n_sampling);

        % Mainshock hypocenter grid over the full fault plane
        HYP_MS{i,ii} = build_hypocenter_grid( ...
            0, Fault_L(i), ...
            0, Fault_W(i), ...
            dL_hyp, Fault_X0(i), Fault_Y0(i), Fault_Hf(i), Dip(i,ii));

        % Mainshock finite ruptures
        RUP_MS{i,ii} = place_ruptures_on_fault( ...
            FG{i,ii}, MAG_MS{i,ii}, RS_MS{i,ii}, HYP_MS{i,ii}, ...
            rupture_placement_type);

        % Keep previous-style variables
        M_vec{i,ii}       = MAG_MS{i,ii}.M_vec;
        FM{i,ii}          = MAG_MS{i,ii}.FM;
        P_M{i,ii}         = MAG_MS{i,ii}.P_M;
        M_mid_total{i,ii} = MAG_MS{i,ii}.M_mid_vec;

        % R limits are now from finite-rupture Rrup, not point-source distance
        min_r_vec(i,ii) = min(RUP_MS{i,ii}.Rrup(:));
        max_r_vec(i,ii) = max(RUP_MS{i,ii}.Rrup(:));

    end

end


%% -----------------------------------------------------------------------
%  MS Deaggregation mesh and PSHA Rrup bins
% ------------------------------------------------------------------------

R_limits = [min(min_r_vec(:)), max(max_r_vec(:))];
M_limits = [Mmin, max(Mmax(:))];

% Base mesh
R_base = 0:dR_deagg:(ceil(R_limits(2)/dR_deagg)*dR_deagg + dR_deagg);
% M_base = 0:dM_deagg:(ceil(M_limits(2)/dM_deagg)*dM_deagg + dM_deagg);

% Selected R mesh
id_R_min_deag = max(1, sum(R_base < R_limits(1)));
id_R_max_deag = min(length(R_base), sum(R_base <= R_limits(2)) + 1);

% Selected M mesh
% id_M_min_deag = max(1, sum(M_base < M_limits(1)));
% id_M_max_deag = min(length(M_base), sum(M_base <= M_limits(2)) + 1);

R_vec_deagg = R_base(id_R_min_deag:id_R_max_deag);

M_vec_deagg = Mmin:dM_deagg: ...
    (Mmin + ceil((max(Mmax(:))-Mmin)/dM_deagg)*dM_deagg);

if M_vec_deagg(end) <= max(Mmax(:))
    M_vec_deagg(end+1) = M_vec_deagg(end) + dM_deagg;
end

% M_vec_deagg = M_base(id_M_min_deag:id_M_max_deag);

R_vec_deagg_mid = 0.5 * (R_vec_deagg(1:end-1) + R_vec_deagg(2:end));
M_vec_deagg_mid = 0.5 * (M_vec_deagg(1:end-1) + M_vec_deagg(2:end));

n_R_deag     = length(R_vec_deagg);
n_R_deag_mid = length(R_vec_deagg_mid);

n_M_deag     = length(M_vec_deagg);
n_M_deag_mid = length(M_vec_deagg_mid);

% Use the same Rrup bins in PSHA and later deaggregation
Rrup_vec_PSHA = R_vec_deagg;
%% -----------------------------------------------------------------------
%  Mainshock GMPM assessment + PSHA
% ------------------------------------------------------------------------

for i = 1:n_faults


    for ii = 1:n_branch_MS


        %% Unique GMPM rows for this MS branch
        [uniqueRows_MS{i,ii}, ~, ic_MS{i,ii}] = ...
            unique(RUP_MS{i,ii}.GMPM_input_mat, 'rows');

        %% Calculate all GMPMs once for this MS branch
        [IM_MS{i,ii}, sigma_MS{i,ii}, LnIM_MS{i,ii}] = ...
            calculate_GMPM_for_unique_rows( ...
            uniqueRows_MS{i,ii}, n_GMPM, GMPM_types, ...
            Vs30, Z1o0, Z2o5, Region, Fault_Type_name{i}, ...
            FVS30, fas_MS, HW(i), T, Dip(i,ii), Lambda(i,ii));

        %% Mainshock PSHA
        PSHA{i,ii} = MSAS_main_PSHA_loop_original( ...
            LnIM_MS{i,ii}, sigma_MS{i,ii}, ic_MS{i,ii}, ...
            activity(i,ii), MAG_MS{i,ii}.P_M, RS_MS{i,ii}.P_rup_W, ...
            RUP_MS{i,ii}.Rrup, Rrup_vec_PSHA, ...
            IM_vec_edge, IM_vec, ...
            n_GMPM, MAG_MS{i,ii}.nM, HYP_MS{i,ii}.n_hyp, RS_MS{i,ii}.n_sampling);

        %% Store outputs in the old-code format
        for j = 1:n_GMPM

            f_IM_E{i,ii,j} = PSHA{i,ii}.f_IM_E(:,:,:,j);

            fx_IM{i,ii,j} = PSHA{i,ii}.fx_IM(:,j);

            exceedance_rate_IM{i,ii,j} = ...
                PSHA{i,ii}.exceedance_rate_IM(:,j);

            Deag_M_R_GMPM_fault{i,ii,j} = ...
                PSHA{i,ii}.Deag_M_R(:,:,:,j);

            Deag_M_R_Exceed{i,ii,j} = ...
                PSHA{i,ii}.Deag_M_R_Exceed(:,:,:,j);

        end

    end

end
%% -----------------------------------------------------------------------
%  MS Hazard Curve: combine faults, MS branches, and GMPMs
% ------------------------------------------------------------------------

exceedance_rate_IM_faults = zeros(n_faults, length(IM_vec_edge));

for i = 1:n_faults

    for ii = 1:n_branch_MS

        for j = 1:n_GMPM

            exceedance_rate_IM_faults(i,:) = exceedance_rate_IM_faults(i,:) +  Branch_weight_MS(i,ii) * GMPM_weight(j) * exceedance_rate_IM{i,ii,j}(:)';

        end

    end

end

exceedance_rate_IM_total = sum(exceedance_rate_IM_faults, 1);
fx_IM_total = -diff(exceedance_rate_IM_total);

fx_IM_fault = zeros(n_faults, length(IM_vec));

for i = 1:n_faults
    fx_IM_fault(i,:) = -diff(exceedance_rate_IM_faults(i,:));
end
 

%% -----------------------------------------------------------------------
%  Finite-fault MS Deaggregation
%  Output format is kept compatible with your previous code.
% ------------------------------------------------------------------------

idr_cell = cell(n_faults,1);
idm_cell = cell(n_faults,1);

Deag_M_R_GMPM_Exceed       = zeros(n_M_deag_mid, n_R_deag_mid, n_GMPM, length(IM_vec)+1, n_faults);
Total_Deag_M_R_GMPM_Exceed = zeros(n_M_deag_mid, n_R_deag_mid, n_GMPM, length(IM_vec)+1);

Deag_M_R_GMPM_Occur        = zeros(n_M_deag_mid, n_R_deag_mid, n_GMPM, length(IM_vec), n_faults);
Total_Deag_M_R_GMPM_Occur  = zeros(n_M_deag_mid, n_R_deag_mid, n_GMPM, length(IM_vec));

Total_Deag_M_R_Occure      = zeros(n_M_deag_mid, n_R_deag_mid, length(IM_vec));

for i = 1:n_faults

    idr_cell{i} = [];
    idm_cell{i} = [];

    for ii = 1:n_branch_MS

        for j = 1:n_GMPM

            w_ms_gmpm = Branch_weight_MS(i,ii) * GMPM_weight(j);

            D_occur_branch  = Deag_M_R_GMPM_fault{i,ii,j};
            D_exceed_branch = Deag_M_R_Exceed{i,ii,j};

            nM_branch = size(D_occur_branch,1);
            nR_branch = size(D_occur_branch,2);
            nIM_branch = min(length(IM_vec), size(D_occur_branch,3));

            for m = 1:nM_branch

                M_mid_current = MAG_MS{i,ii}.M_mid_vec(m);

                idm = sum(M_vec_deagg < M_mid_current);

                if idm < 1 || idm > n_M_deag_mid
                    error('M deaggregation bin index out of range.')
                end

                idm_cell{i} = [idm_cell{i}, idm];

                for r = 1:nR_branch

                    % Since Rrup_vec_PSHA = R_vec_deagg,
                    % the PSHA R bin r corresponds to R_vec_deagg_mid(r).
                    R_mid_current = R_vec_deagg_mid(r);

                    idr = sum(R_vec_deagg < R_mid_current);

                    if idr < 1 || idr > n_R_deag_mid
                        error('R deaggregation bin index out of range.')
                    end

                    idr_cell{i} = [idr_cell{i}, idr];

                    for s = 1:nIM_branch + 1

                        % ------------------------------------------------
                        % Occurrence-conditioned deaggregation
                        % -------------------------------------------------
                        if s < nIM_branch + 1
                            if fx_IM_total(s) > 0
    
                                mix_occur = ...
                                    w_ms_gmpm * fx_IM{i,ii,j}(s) / fx_IM_total(s);
    
                                Deag_M_R_GMPM_Occur(idm,idr,j,s,i) = ...
                                    Deag_M_R_GMPM_Occur(idm,idr,j,s,i) + ...
                                    mix_occur * D_occur_branch(m,r,s);
    
                                Total_Deag_M_R_GMPM_Occur(idm,idr,j,s) = ...
                                    Total_Deag_M_R_GMPM_Occur(idm,idr,j,s) + ...
                                    mix_occur * D_occur_branch(m,r,s);
    
                                Total_Deag_M_R_Occure(idm,idr,s) = ...
                                    Total_Deag_M_R_Occure(idm,idr,s) + ...
                                    mix_occur * D_occur_branch(m,r,s);
    
                            end
                        end

                        % ------------------------------------------------
                        % Exceedance-conditioned deaggregation
                        % -------------------------------------------------
            
                        if exceedance_rate_IM_total(s) > 0

                            mix_exceed = ...
                                w_ms_gmpm * ...
                                exceedance_rate_IM{i,ii,j}(s) / ...
                                exceedance_rate_IM_total(s);

                            Deag_M_R_GMPM_Exceed(idm,idr,j,s,i) = ...
                                Deag_M_R_GMPM_Exceed(idm,idr,j,s,i) + ...
                                mix_exceed * D_exceed_branch(m,r,s);

                            Total_Deag_M_R_GMPM_Exceed(idm,idr,j,s) = ...
                                Total_Deag_M_R_GMPM_Exceed(idm,idr,j,s) + ...
                                mix_exceed * D_exceed_branch(m,r,s);

                        end

                    end

                end

            end

        end

    end

    idr_cell{i} = unique(idr_cell{i});
    idm_cell{i} = unique(idm_cell{i});

end


%% -----------------------------------------------------------------------
%  Transition P(E|IMe) without AS (fragility-based, per IM)
%  Uses PDs_func (external) to fill upper-triangular transitions.
% ------------------------------------------------------------------------
PE_IMe_WoAS = eye(Nd+1,Nd+1);
for sa = 1:n_IM
    for i = 1:Nd
        for j = i:Nd
            PE_IMe_WoAS(i,j+1,sa) = PDs_func(Mu_mat, Beta_mat, i, j, IM_vec(sa), 1);
        end
        PE_IMe_WoAS(i,i,sa) = 1 - sum(PE_IMe_WoAS(i,i+1:end,sa));
    end
    PE_IMe_WoAS(end,end,sa) = 1;
end

% -----------------------------------------------------------------------
%  Associate finite MS rupture planes with each M-R deag bin
% ------------------------------------------------------------------------

[MS_rupture_deag_True, idm_cell, idr_cell, ...
    MS_deag_count, MS_unique_deag_count] = ...
    build_MS_rupture_deag_association_from_catalog3( ...
    RUP_MS, M_vec_deagg, R_vec_deagg, ...
    Dip, Lambda, Branch_weight_MS, activity, ...
    parent_averaging_method);
%% -----------------------------------------------------------------------
%  Aftershock analysis over (GMPM, Fault, AS-branch)
%  Finite 3D fault formulation
% ------------------------------------------------------------------------

PA_Ex_IM_deag = cell(n_branch_AS, n_faults, n_GMPM);
f_IM_A_deag   = cell(n_branch_AS, n_faults, n_GMPM);
EN_A_deag     = cell(n_branch_AS, n_faults, n_GMPM);
n_parent_deag = cell(n_branch_AS, n_faults, n_GMPM);

% AS rupture-scaling sampling
n_sampling_AS = n_sampling;


if strcmp(AS_Analysis, 'SERIAL')

    for j = 1:n_GMPM
        for i = 1:n_faults
            for jj = 1:n_branch_AS

                [PA_Ex_IM_deag{jj,i,j}, ...
                 f_IM_A_deag{jj,i,j}, ...
                 EN_A_deag{jj,i,j}, ...
                 n_parent_deag{jj,i,j}] = ...
                    aftershock_sequence_calculation_finite_deag_catalog( ...
                        MS_rupture_deag_True{i}, ...
                        idm_cell{i}, ...
                        M_vec_deagg_mid, ...
                        DT, ...
                        a(i,jj), b_AS(i,jj), c(i,jj), p(i,jj), ...
                        Mmin, ...
                        IM_vec, IM_vec_edge, ...
                        dma, ...
                        GMPM_types{j}, ...
                        Vs30, Z1o0, Z2o5, ...
                        Region, Fault_Type_name{i}, ...
                        FVS30, fas_AS, HW(i), T, ...
                        Dip(i,:), ...
                        Lambda(i,:), ...
                        FG(i,:), ...
                        Fault_Hf(i), ...
                        Fault_X0(i), Fault_Y0(i), ...
                        dra, ...
                        variation_in_rupture_area, ...
                        n_sampling_AS, ...
                        rupture_placement_type, AS_rupture_domain, ...
                        NA_parent_averaging_method);

            end
        end
    end


elseif strcmp(AS_Analysis, 'PARALLEL')

    pool = gcp('nocreate');

    if isempty(pool)
        parpool('local', n_Core_Parallel);
    end

    total_tasks = n_branch_AS * n_faults * n_GMPM;

    PA_flat = cell(total_tasks,1);
    fA_flat = cell(total_tasks,1);
    EN_flat = cell(total_tasks,1);
    NP_flat = cell(total_tasks,1);

    parfor task = 1:total_tasks

        [jj,i,j] = ind2sub( ...
            [n_branch_AS, n_faults, n_GMPM], task);

        [PA_flat{task}, ...
         fA_flat{task}, ...
         EN_flat{task}, ...
         NP_flat{task}] = ...
            aftershock_sequence_calculation_finite_deag_catalog( ...
                MS_rupture_deag_True{i}, ...
                idm_cell{i}, ...
                M_vec_deagg_mid, ...
                DT, ...
                a(i,jj), b_AS(i,jj), c(i,jj), p(i,jj), ...
                Mmin, ...
                IM_vec, IM_vec_edge, ...
                dma, ...
                GMPM_types{j}, ...
                Vs30, Z1o0, Z2o5, ...
                Region, Fault_Type_name{i}, ...
                FVS30, fas_AS, HW(i), T, ...
                Dip(i,:), ...
                Lambda(i,:), ...
                FG(i,:), ...
                Fault_Hf(i), ...
                Fault_X0(i), Fault_Y0(i), ...
                dra, ...
                variation_in_rupture_area, ...
                n_sampling_AS, ...
                rupture_placement_type, AS_rupture_domain, ...
                NA_parent_averaging_method);

    end

    PA_Ex_IM_deag = reshape(PA_flat, ...
        [n_branch_AS, n_faults, n_GMPM]);

    f_IM_A_deag = reshape(fA_flat, ...
        [n_branch_AS, n_faults, n_GMPM]);

    EN_A_deag = reshape(EN_flat, ...
        [n_branch_AS, n_faults, n_GMPM]);

    n_parent_deag = reshape(NP_flat, ...
        [n_branch_AS, n_faults, n_GMPM]);

else

    error('Unknown AS_Analysis option. Use SERIAL or PARALLEL.');

end


%% -----------------------------------------------------------------------
%  Extend parent-specific EN_A and f_IM_A_deag to all deaggregation M bins
%  EN_A is now M-R dependent because each parent uses its actual MS magnitude.
% ------------------------------------------------------------------------
EN_A_deag_ext     = cell(size(EN_A_deag));
f_IM_A_deag_ext   = cell(size(f_IM_A_deag));
EN_A_MR           = cell(size(EN_A_deag));

for j = 1:n_GMPM
    for i = 1:n_faults
        for jj = 1:n_branch_AS

            % Helper output orientation is [R-bin x active-M-bin].
            EN_local = EN_A_deag{jj,i,j};

            if size(EN_local,1) ~= n_R_deag_mid || ...
                    size(EN_local,2) ~= numel(idm_cell{i})
                error(['Unexpected EN_A_deag size for fault=%d, AS branch=%d, ', ...
                       'GMPM=%d. Expected [%d x %d], got [%d x %d].'], ...
                       i, jj, j, n_R_deag_mid, numel(idm_cell{i}), ...
                       size(EN_local,1), size(EN_local,2))
            end

            EN_A_deag_ext{jj,i,j} = zeros(n_M_deag_mid, n_R_deag_mid);
            EN_A_deag_ext{jj,i,j}(idm_cell{i},:) = EN_local.';

            % Direct M-R lookup used below in alpha_faults.
            EN_A_MR{jj,i,j} = EN_A_deag_ext{jj,i,j};

            f_IM_A_deag_ext{jj,i,j} = zeros(n_IM, n_R_deag_mid, n_M_deag_mid);
            f_IM_A_deag_ext{jj,i,j}(:,:,idm_cell{i}) = f_IM_A_deag{jj,i,j};

        end
    end
end


%% -----------------------------------------------------------------------
%  Weighted averaging across AS branches & faults (Linear )
% ------------------------------------------------------------------------
EN_A_ave_branch_faults_deag            = zeros(n_M_deag_mid, n_R_deag_mid, n_GMPM, n_IM); % Expected numebr of AS for each deag bin and IM, averaged on AS branches
f_IM_A_ave_branch_faults_deag_lin      = zeros(n_IM, n_R_deag_mid, n_M_deag_mid, n_GMPM); % f_IM_A each deag bin and IM, averaged on AS for each deag bin, averaged on AS branches, and faults
f_IM_A_ave_branch_faults_GMPM_deag_lin = zeros(n_IM, n_R_deag_mid, n_M_deag_mid);  % Expected numebr of AS for each deag bin and GMPM, averaged on AS branches, and faults and IM levels
EN_A_ave_branch_faults_GMPM_deag       = zeros(n_M_deag_mid, n_R_deag_mid, n_IM); % f_IM_A each deag bin, IM and GMPM, averaged on AS for each deag bin, averaged on AS branches, and faults
EN_A_mean                              = zeros(n_M_deag_mid, n_R_deag_mid); % Expected numebr of AS for each deag bin, averaged on AS branches, GMPMs, faults, and IM

for m = 1:n_M_deag_mid
    for n = 1:n_R_deag_mid
        for s = 1:n_IM
            for j = 1:n_GMPM
                for i = 1:n_faults
                    for jj = 1:n_branch_AS
                        if  Total_Deag_M_R_GMPM_Occur(m,n,j,s) ~= 0
                            alpha_faults(m,n,jj,i,j,s) = EN_A_MR{jj,i,j}(m,n) * ...
                                (Deag_M_R_GMPM_Occur(m,n,j,s,i)) ./ (Total_Deag_M_R_GMPM_Occur(m,n,j,s)); %  Total_Deag_M_R_GMPM_Occur(m,n,j,s) = sum(Deag_M_R_GMPM_Occur(m,n,j,s,:))
                        else
                            alpha_faults(m,n,jj,i,j,s) = 0;
                        end

                        fval = f_IM_A_deag_ext{jj,i,j}(s,n,m);
                        f_IM_A_ave_branch_faults_deag_lin(s,n,m,j)       = f_IM_A_ave_branch_faults_deag_lin(s,n,m,j)      + fval * Branch_weight_AS(i,jj) * alpha_faults(m,n,jj,i,j,s);
                        f_IM_A_ave_branch_faults_GMPM_deag_lin(s,n,m)    = f_IM_A_ave_branch_faults_GMPM_deag_lin(s,n,m)  + fval * Branch_weight_AS(i,jj) * GMPM_weight(j) * alpha_faults(m,n,jj,i,j,s);
                        
                        EN_A_ave_branch_faults_deag(m,n,j,s)    = EN_A_ave_branch_faults_deag(m,n,j,s)    + Branch_weight_AS(i,jj) * alpha_faults(m,n,jj,i,j,s);
                        EN_A_ave_branch_faults_GMPM_deag(m,n,s) = EN_A_ave_branch_faults_GMPM_deag(m,n,s) + Branch_weight_AS(i,jj) * GMPM_weight(j) * alpha_faults(m,n,jj,i,j,s);
                    end
                end
                if Total_Deag_M_R_GMPM_Occur(m,n,j,s) ~= 0
                    f_IM_A_ave_branch_faults_deag_lin(s,n,m,j) = f_IM_A_ave_branch_faults_deag_lin(s,n,m,j) ./ EN_A_ave_branch_faults_deag(m,n,j,s);
                end
            end
            if Total_Deag_M_R_Occure(m,n,s) ~= 0
                f_IM_A_ave_branch_faults_GMPM_deag_lin(s,n,m) = f_IM_A_ave_branch_faults_GMPM_deag_lin(s,n,m) ./ EN_A_ave_branch_faults_GMPM_deag(m,n,s);
            end
        end

        % Correct averaging of EN_A across GMPMs:
        % renormalize GMPM weights using only GMPMs that have support
        EN_s = [];
        
        for ss = 1:n_IM
        
            if Total_Deag_M_R_Occure(m,n,ss) <= 0
                continue
            end
        
            % GMPMs that actually contribute to this M-R-IM bin
            active_GMPM = squeeze( ...
                Total_Deag_M_R_GMPM_Occur(m,n,:,ss)) > 0;
        
            % Total prior weight of the active GMPMs
            w_active = sum(GMPM_weight(active_GMPM(:)'));
        
            if w_active <= 0
                continue
            end
        
            % Renormalize EN_A over the active GMPMs
            EN_s(end+1) = ...
                EN_A_ave_branch_faults_GMPM_deag(m,n,ss) / w_active;
        
        end
        
        if ~isempty(EN_s)
            EN_A_mean(m,n) = mean(EN_s);
        else
            EN_A_mean(m,n) = 0;
        end
     end
end


f_IM_A_ave_branch_faults_deag       = zeros(n_IM, n_R_deag_mid, n_M_deag_mid, n_GMPM);
f_IM_A_ave_branch_faults_GMPM_deag  = zeros(n_IM, n_R_deag_mid, n_M_deag_mid);

for m = 1:n_M_deag_mid
    for n = 1:n_R_deag_mid
        for s = 1:n_IM
            for j = 1:n_GMPM
                for i = 1:n_faults
                    for jj = 1:n_branch_AS
                        if Total_Deag_M_R_GMPM_Occur(m,n,j,s) ~= 0
                            alpha_faults(m,n,jj,i,j,s) = EN_A_MR{jj,i,j}(m,n) * Deag_M_R_GMPM_Occur(m,n,j,s,i) ./ Total_Deag_M_R_GMPM_Occur(m,n,j,s);
                        else
                            % Defensive checks as in your code
                            check1 = sum(Deag_M_R_GMPM_Occur(m,n,j,s,i));
                            check2 = sum(Total_Deag_M_R_GMPM_Occur(m,n,j,s));
                            if (check1==0 && check2~=0) || (check2==0 && check1~=0)
                                error('Inconsistent deaggregation mass at (m=%d,n=%d,j=%d,s=%d).',m,n,j,s);
                            end
                            alpha_faults(m,n,jj,i,j,s) = 0;
                        end
                        fval = f_IM_A_deag_ext{jj,i,j}(s,n,m);
                        if fval == 0, fval = 1e-50; end
                        f_IM_A_ave_branch_faults_deag(s,n,m,j)  = f_IM_A_ave_branch_faults_deag(s,n,m,j)  + log(fval) * Branch_weight_AS(i,jj)                   * alpha_faults(m,n,jj,i,j,s);
                        f_IM_A_ave_branch_faults_GMPM_deag(s,n,m) = f_IM_A_ave_branch_faults_GMPM_deag(s,n,m) + log(fval) * Branch_weight_AS(i,jj) * GMPM_weight(j) * alpha_faults(m,n,jj,i,j,s);
                    end
                end
                if Total_Deag_M_R_GMPM_Occur(m,n,j,s) ~= 0
                    f_IM_A_ave_branch_faults_deag(s,n,m,j) = exp(f_IM_A_ave_branch_faults_deag(s,n,m,j) ./ EN_A_ave_branch_faults_deag(m,n,j,s));
                end
            end
            if Total_Deag_M_R_Occure(m,n,s) ~= 0
                f_IM_A_ave_branch_faults_GMPM_deag(s,n,m) = exp(f_IM_A_ave_branch_faults_GMPM_deag(s,n,m) ./ EN_A_ave_branch_faults_GMPM_deag(m,n,s));
            end
        end
    end
end



%% -----------------------------------------------------------------------
disp('*****  Deag M & R for all IM (Linear weighted averaging) *****');

PS_Method_5_deag_Lin      = zeros(Nd+1,Nd+1);
PS_E_onlymain_5_deag_Lin  = zeros(Nd+1,Nd+1);
PS_MeRe_T_ave_branch_faults_GMPM_IM = zeros(Nd+1,Nd+1,n_IM);

for se = 1:n_IM
    sum_Mi_Rj = 0;
    for m = 1:n_M_deag_mid
        for n = 1:n_R_deag_mid
            if Total_Deag_M_R_Occure(m,n,se) ~= 0
                PA_ij_ave_branch_fault_GMPM(:,:,n,m) = P_ij_calculator(Mu_mat, Beta_mat, 1, 1, IM_vec, f_IM_A_ave_branch_faults_GMPM_deag_lin(:,n,m));
                sum_Mi_Rj = sum_Mi_Rj + PA_ij_ave_branch_fault_GMPM(:,:,n,m)^EN_A_mean(m,n) * Total_Deag_M_R_Occure(m,n,se);
            end
        end
    end
    PS_Method_5_deag_Lin     = PS_Method_5_deag_Lin     + PE_IMe_WoAS(:,:,se) * sum_Mi_Rj * fx_IM_total(se);
    PS_E_onlymain_5_deag_Lin = PS_E_onlymain_5_deag_Lin + PE_IMe_WoAS(:,:,se) * fx_IM_total(se);
end

P_E_5_deag_Lin            = PS_Method_5_deag_Lin + (1-PS_Method_5_deag_Lin(end,end))*eye(Nd+1,Nd+1)
P_E_5_OnlyMainshock_deag_Lin = PS_E_onlymain_5_deag_Lin + (1-PS_E_onlymain_5_deag_Lin(end,end))*eye(Nd+1,Nd+1)
vpa(P_E_5_deag_Lin,10)
vpa(P_E_5_OnlyMainshock_deag_Lin,10)


