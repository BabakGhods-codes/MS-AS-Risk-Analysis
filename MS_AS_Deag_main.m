clc; clear; tic;

%% =======================================================================
%  Multi-Branch MS–AS Risk Analysis with PSHA + Deaggregation
% ========================================================================

addpath('Function Folder');

warning('BSSA_2014 DMPM: Z1.0 from CY model is used (as in original).');

%% -----------------------------------------------------------------------
%  INPUTS
% ------------------------------------------------------------------------

% -----------------------------
% Faults MS logic-tree (epistemic)
% -----------------------------
Mmax = [8   8   7.6 7.6;   7.7 7.7 8.2 8.2];
b_M  = [1.15 1  1.15 1;    0.9 0.8 0.9 0.8];

% -----------------------------
% Fault constants
% -----------------------------
landa  = [0.4 0.4 0.4 0.4; 0.2 0.2 0.2 0.2];
Fault_L = [60; 40];                % Fault lengths (km)
region = 'California';
Fault_Type_name = {'Strike-Slip'; 'Reverse'};
Dip     = [pi/2 pi/2 pi/2 pi/2;   pi/3 pi/3 pi/3 pi/3];

% -----------------------------
% Site
% -----------------------------
FVS30 = 1;         % FVS30 = 0 for Vs30 is inferred from geology; FVS30 = 1 for measured  Vs30
fas   = 0;         % Flag for aftershocks in GMPMs
Site_number = 1;   % Site 1 or 2 Properties    
if Site_number == 1
    HW     = [1; 0];                   % Hanging wall indicators
    H      = [10,10];                  % Fault depth proxies (km)
    Site_alongLengthDist = [30 0];     % Site position along length (km)  
    Vs30  = 540;                       % m/s
elseif Site_number == 2
    HW     = [1; 1];                   % Hanging wall indicators
    H      = [52.43 ; 20];                  % Fault depth proxies (km)
    Site_alongLengthDist = [30 30];     % Site position along length (km)
    Vs30  = 360;                       % m/s
end
% -----------------------------
% Aftershock Omori-type params
% -----------------------------
a = [-1.6 -1.6 -1.4 -1.4;  -1.8 -1.8 -1.6 -1.6];
b = [ 1.15 1.15 1.0  1.0;    0.9  0.9  0.8  0.8];
c = [ 0.04 0.04 0.06 0.06;   0.02 0.02 0.03 0.03];
p = [ 0.85 0.85 0.95 0.95;   1.05 1.05 1.2  1.2];

% Spatial distribution flag for AS distribution flag around the MS epicenter 
AS_Spatial_Distribution = [0 1 0 1; 0 1 0 1]; % 0 for squared and 1 for rectangulat (1:3 for Strike-Slip and 1:2 for Reverse faults)

% -----------------------------
% Logic-tree weights
% -----------------------------
Branch_weight_MS = [0.2 0.2 0.3 0.3; 0.2 0.2 0.3 0.3];
Branch_weight_AS = [0.15 0.35 0.15 0.35; 0.15 0.35 0.15 0.35];

% -----------------------------
% Structure
% -----------------------------
Structure_number = 1;   %  1 is the 8-story building and 2 is the 4 story building
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
dm        = 0.05;     % mag step
Mmin      = 5;       % min mag
dL        = 2.5;     % along-strike step (km)

% -----------------------------
% Aftershock window and mesh
% -----------------------------
DT  = 60;            % days window
dma = dm;            % AS mag step (kept equal to MS)
dra = dL;            % AS distance step (kept equal to MS)

% -----------------------------
% Deaggregation mesh
% -----------------------------
dR_deagg = 20;
dM_deagg = 0.5;

% Derived sizes
n_faults     = size(Mmax,1);
n_branch_MS  = size(Branch_weight_MS,2);
n_branch_AS  = size(Branch_weight_AS,2);

% -----------------------------
% Aftershock Analysis Type
% -----------------------------
AS_Analysis = 'PARALLEL'; % SERIAL or PARALLEL
n_Core_Parallel = 8; % used only if AS_Analysis = 'PARALLEL'
%% -----------------------------------------------------------------------
%  MS: Loop over faults and MS branches -> PSHA + Deaggregation
% ------------------------------------------------------------------------

% Preallocations (cell because sizes differ by branch/fault)
R_vec            = cell(n_faults,1);
R_mid_vec_fault  = cell(n_faults,1);
R_mid_total      = cell(n_faults,n_branch_MS);

M_vec            = cell(n_faults,n_branch_MS);
FM               = cell(n_faults,n_branch_MS);
M_mid_total      = cell(n_faults,n_branch_MS);
P_M              = cell(n_faults,n_branch_MS);

fx_IM            = cell(n_faults,n_branch_MS,n_GMPM);
exceedance_rate_IM = cell(n_faults,n_branch_MS,n_GMPM);
Deag_M_R_GMPM_fault = cell(n_faults,n_branch_MS,n_GMPM);
Deag_M_R_Exceed  = cell(n_faults,n_branch_MS,n_GMPM);

min_r_vec        = zeros(1,n_faults);
max_r_vec        = zeros(1,n_faults);
% Main loop over faults and MS branches
for i = 1:n_faults
    % Along fault mesh and geometric distances to site
    R_vec{i} = dL/2 : dL : Fault_L(i);
    R_mid_vec_fault{i} = sqrt(H(i)^2 + (R_vec{i} - Site_alongLengthDist(i)).^2);

    for ii = 1:n_branch_MS
        % Magnitude distribution on [Mmin, Mmax(i,ii)]
        M_vec{i,ii}  = Mmin : dm : Mmax(i,ii);
        FM{i,ii}     = (1-10.^(-b_M(i,ii)*(M_vec{i,ii}-Mmin))) / (1-10^(-b_M(i,ii)*(Mmax(i,ii)-Mmin))); % CDF
        P_M{i,ii}    = diff(FM{i,ii});  % Prob masses
        M_mid_total{i,ii} = (M_vec{i,ii}(2:end)+M_vec{i,ii}(1:end-1))/2;

        for j = 1:n_GMPM
            [~, fx_IM{i,ii,j}, exceedance_rate_IM{i,ii,j}, ...
             Deag_M_R_GMPM_fault{i,ii,j}, ~, Deag_M_R_Exceed{i,ii,j}] = ...
                mainshcok_PSHA_Disagg_main( ...
                    M_mid_total{i,ii}, R_mid_vec_fault{i}, dIM, IM_vec_edge, IM_vec, ...
                    P_M{i,ii}, landa(i,ii), GMPM_types{j}, Vs30, HW(i), Dip(i,ii), ...
                    FVS30, fas, T, region, Fault_Type_name{i});
        end

        % Bookkeeping
        R_mid_total{i,ii}            = R_mid_vec_fault{i};
    end
    min_r_vec(i) = min(R_mid_vec_fault{i}); 
    max_r_vec(i) = max(R_mid_vec_fault{i});
end

%% -----------------------------------------------------------------------
%  MS Hazard Curve (Combine faults, branches, GMPMs)
% ------------------------------------------------------------------------
exceedance_rate_IM_faults = zeros(n_faults, n_IM+1);
fx_IM_fault               = zeros(n_faults, n_IM);
sum_fx_IM = 0;

for i = 1:n_faults
    for ii = 1:n_branch_MS
        for j = 1:n_GMPM
            w = Branch_weight_MS(i,ii) * GMPM_weight(j);
            exceedance_rate_IM_faults(i,:) = exceedance_rate_IM_faults(i,:) + w * exceedance_rate_IM{i,ii,j};
            sum_fx_IM = sum_fx_IM + w * sum(fx_IM{i,ii,j});
        end
    end
    fx_IM_fault(i,:) = -diff(exceedance_rate_IM_faults(i,:));
end

exceedance_rate_IM_total = sum(exceedance_rate_IM_faults,1);
fx_IM_total              = -diff(exceedance_rate_IM_total);
 

%% -----------------------------------------------------------------------
%  MS Deaggregation setup (M-R-GMPM-IM)
% ------------------------------------------------------------------------
R_limits = [min(min_r_vec), max(max_r_vec)];      % minimum and maximum limits of R from the PSHA inputs
M_limits = [Mmin, max(max(Mmax))];                % minimum and maximum limits of M from the PSHA inputs

% Based mesh for Deaggregation
R_base = 0:dR_deagg:(R_limits(2)+dR_deagg);          
M_base = 0:dM_deagg:(M_limits(2)+dM_deagg);

% Identifying the needed meshgird form the Base mesh
id_R_min_deag = sum(R_base < R_limits(1));
id_R_max_deag = sum(R_base <= R_limits(2)) + 1;

id_M_min_deag = sum(M_base < M_limits(1));
id_M_max_deag = sum(M_base <= M_limits(2)) + 1;

R_vec_deagg      = R_base(id_R_min_deag:id_R_max_deag);    % Selected R meshgird for deaggregation  
M_vec_deagg      = M_base(id_M_min_deag:id_M_max_deag);    % Selected M meshgird for deaggregation 
M_vec_deagg_mid  = 0.5*(M_vec_deagg(2:end)+M_vec_deagg(1:end-1)); % Selected mid-R meshgird for deaggregation 
R_vec_deagg_mid  = 0.5*(R_vec_deagg(2:end)+R_vec_deagg(1:end-1)); % Selected mid-M meshgird for deaggregation 


n_M_deag = length(M_vec_deagg);
n_M_deag_mid = length(M_vec_deagg_mid);
n_R_deag = length(R_vec_deagg);
n_R_deag_mid = length(R_vec_deagg_mid);

% Bins to place ruptures
idr_cell = cell(n_faults,1);
idm_cell = cell(n_faults,1);

% Containers for (M,R,GMPM,IM) deaggregation
Deag_M_R_GMPM_Exceed       = zeros(n_M_deag-1, n_R_deag-1, n_GMPM, n_IM, n_faults);
Total_Deag_M_R_GMPM_Exceed = zeros(n_M_deag-1, n_R_deag-1, n_GMPM, n_IM);
Total_Deag_M_R_GMPM_Occur  = zeros(n_M_deag-1, n_R_deag-1, n_GMPM, n_IM);
Total_Deag_M_R_Occure      = zeros(n_M_deag-1, n_R_deag-1, n_IM);
Deag_M_R_GMPM_Occur        = zeros(n_M_deag-1, n_R_deag-1, n_GMPM, n_IM, n_faults);

% Fill deaggregation cubes
for i = 1:n_faults
    for ii = 1:n_branch_MS
        for j = 1:n_GMPM
            for s = 1:n_IM
                for r = 1:length(R_mid_total{i,ii})
                    for m = 1:length(M_mid_total{i,ii})
                        idr = sum(R_vec_deagg < R_mid_total{i,ii}(r));
                        idm = sum(M_vec_deagg < M_mid_total{i,ii}(m));

                        idr_cell{i} = [idr_cell{i}, idr];
                        idm_cell{i} = [idm_cell{i}, idm];

                        w_ms_gmpm = Branch_weight_MS(i,ii) * GMPM_weight(j);

                        % Occurrence-conditioned contributions
                        Deag_M_R_GMPM_Occur(idm,idr,j,s,i) = Deag_M_R_GMPM_Occur(idm,idr,j,s,i) + ...
                            w_ms_gmpm * Deag_M_R_GMPM_fault{i,ii,j}(m,r,s) * fx_IM{i,ii,j}(s)/fx_IM_total(s);

                        Total_Deag_M_R_GMPM_Occur(idm,idr,j,s) = Total_Deag_M_R_GMPM_Occur(idm,idr,j,s) + ...
                            w_ms_gmpm * Deag_M_R_GMPM_fault{i,ii,j}(m,r,s) * fx_IM{i,ii,j}(s)/fx_IM_total(s);

                        Total_Deag_M_R_Occure(idm,idr,s) = Total_Deag_M_R_Occure(idm,idr,s) + ...
                            w_ms_gmpm * Deag_M_R_GMPM_fault{i,ii,j}(m,r,s) * fx_IM{i,ii,j}(s)/fx_IM_total(s);

                        % Exceedance-conditioned contributions
                        Deag_M_R_GMPM_Exceed(idm,idr,j,s,i) = Deag_M_R_GMPM_Exceed(idm,idr,j,s,i) + ...
                            w_ms_gmpm * Deag_M_R_Exceed{i,ii,j}(m,r,s) * exceedance_rate_IM{i,ii,j}(s)/exceedance_rate_IM_total(s);

                        Total_Deag_M_R_GMPM_Exceed(idm,idr,j,s) = Total_Deag_M_R_GMPM_Exceed(idm,idr,j,s) + ...
                            w_ms_gmpm * Deag_M_R_Exceed{i,ii,j}(m,r,s) * exceedance_rate_IM{i,ii,j}(s)/exceedance_rate_IM_total(s);
                    end
                end
            end
        end
        idr_cell{i} = unique(idr_cell{i});
        idm_cell{i} = unique(idm_cell{i});
    end
end

% Quick checks
for s = 1:n_IM
    checkOfDeagg_total(s) = sum(sum(sum(Total_Deag_M_R_GMPM_Exceed(:,:,:,s)))); %#ok<SAGROW>
    for i = 1:n_faults
        checkOfDeagg(s,i) = sum(sum(sum(Deag_M_R_GMPM_Exceed(:,:,:,s,i)))); %#ok<SAGROW>
    end
end
sum(checkOfDeagg,2);    
checkOfDeagg_total;    

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
%% -----------------------------------------------------------------------
%  Aftershock epicentral discretization (Method 2 in your code)
%  We mark the deaggregation R-bin for each MS rupture position.
% ------------------------------------------------------------------------
% --- Assessment of Epicenters associated with each Deaggregation bin --- %
% Given that M and R are assumed separate (No scaling assumption),
% epicenters are associated with each rupture bin of the deaggregation 
% ----------------------------------------------------------------------- %
for i = 1:n_faults
    epicentral_location_deag_True{i} = [];
    for ss = 1:length(R_mid_vec_fault{i})
         idr = sum(R_vec_deagg < R_mid_vec_fault{i}(ss));
         % Column vector with one-hot at (idr+1)
         col = [R_vec{i}(ss); zeros(n_R_deag-1, 1)];
         col(idr+1) = 1;
         epicentral_location_deag_True{i} = [epicentral_location_deag_True{i}, col]; 
    end
end
mmmmmmmm
%% -----------------------------------------------------------------------
%  Aftershock analysis over (GMPM, Fault, AS-branch)
% ------------------------------------------------------------------------
if  strcmp(AS_Analysis, 'SERIAL')
for j = 1:n_GMPM
    for i = 1:n_faults
        for jj = 1:n_branch_AS
            [~, f_IM_A_deag{jj,i,j}, ~, ...
             EN_A_deag{jj,i,j}] = ...
             aftershock_sequence_calculation_main( ...
                M_vec_deagg_mid(idm_cell{i}), epicentral_location_deag_True{i}, ...
                DT,  a(i,jj), b(i,jj), c(i,jj), p(i,jj), Mmin, ...
                IM_vec, IM_vec_edge, ...
                dma, GMPM_types{j}, Vs30, HW(i), T, region, ...
                Fault_Type_name{i}, dra, ...
                Dip(i,jj), FVS30, fas, H(i), AS_Spatial_Distribution(i,jj), ...
                Site_alongLengthDist(i));
        end
    end
end
% -----------------------------------------------------------------------
%  Aftershock analysis — PARALLEL version (kept, unchanged formulations)
% -----------------------------------------------------------------------
elseif strcmp(AS_Analysis, 'PARALLEL')
    
    if isempty(gcp('nocreate')); parpool(n_Core_Parallel); end
    
    total_tasks = n_GMPM * n_faults * n_branch_AS;
    
    PA_Ex_IM_deag_flat = cell(total_tasks,1);
    f_IM_A_deag_flat   = cell(total_tasks,1);
    f_IM_A_flat        = cell(total_tasks,1);
    EN_A_deag_flat     = cell(total_tasks,1);
    
    parfor task_idx = 1:total_tasks
        [jj, i, j] = ind2sub([n_branch_AS, n_faults, n_GMPM], task_idx);
        [PA_Ex_IM_deag_flat{task_idx}, ...
         f_IM_A_deag_flat{task_idx}, ...
         f_IM_A_flat{task_idx}, ...
         EN_A_deag_flat{task_idx}] = ...
         aftershock_sequence_calculation_main( ...
            M_vec_deagg_mid(idm_cell{i}), epicentral_location_deag_True{i}, ...
            DT, a(i,jj), b(i,jj), c(i,jj), p(i,jj), Mmin, ...
            IM_vec, IM_vec_edge, dma, ...
            GMPM_types{j}, Vs30, HW(i), T, region, Fault_Type_name{i}, ...
            dra,  Dip(i,jj), FVS30, fas, H(i), ...
            AS_Spatial_Distribution(i,jj), Site_alongLengthDist(i));
    end

    % Reshape back to (AS-branch, fault, GMPM)
    f_IM_A_deag   = reshape(f_IM_A_deag_flat,   [n_branch_AS, n_faults, n_GMPM]);
    EN_A_deag     = reshape(EN_A_deag_flat,     [n_branch_AS, n_faults, n_GMPM]);
end
%% -----------------------------------------------------------------------
%  Extend EN_A and f_IM_A_deag to total deaggregation bins
% (Only magnitude needs extension as the R bins are for total deaggregation)
% ------------------------------------------------------------------------
EN_A_deag_ext     = cell(size(EN_A_deag));
f_IM_A_deag_ext   = cell(size(f_IM_A_deag));

for j = 1:n_GMPM
    for i = 1:n_faults
        for jj = 1:n_branch_AS
            EN_A_deag_ext{jj,i,j} = zeros(size(M_vec_deagg_mid));
            EN_A_deag_ext{jj,i,j}(idm_cell{i}) = EN_A_deag{jj,i,j};

            f_IM_A_deag_ext{jj,i,j} = zeros(n_IM, n_R_deag_mid, n_M_deag_mid);
            f_IM_A_deag_ext{jj,i,j}(:,:,idm_cell{i}) = f_IM_A_deag{jj,i,j};
        end
    end
end

%% -----------------------------------------------------------------------
%  Weighted averaging across AS branches & faults (Linear)
% ------------------------------------------------------------------------
f_IM_A_ave_branch_faults_GMPM_deag_lin = zeros(n_IM, n_R_deag_mid, n_M_deag_mid);  % Expected numebr of AS for each deag bin and GMPM, averaged on AS branches, and faults and IM levels
EN_A_ave_branch_faults_GMPM_deag       = zeros(n_M_deag_mid, n_R_deag_mid, n_IM); % f_IM_A each deag bin, IM and GMPM, averaged on AS for each deag bin, averaged on AS branches, and faults
EN_A_mean                              = zeros(n_M_deag_mid, n_R_deag_mid); % Expected numebr of AS for each deag bin, averaged on AS branches, GMPMs, faults, and IM

for m = 1:n_M_deag_mid
    for n = 1:n_R_deag_mid
        for s = 1:n_IM
            for j = 1:n_GMPM
                for i = 1:n_faults
                    for jj = 1:n_branch_AS
                        EN_A_MR{jj,i,j} = repmat(EN_A_deag_ext{jj,i,j}',1,n_R_deag_mid); % Extending to all deaggregation bins
                        if  Total_Deag_M_R_GMPM_Occur(m,n,j,s) ~= 0
                            alpha_faults(m,n,jj,i,j,s) = EN_A_MR{jj,i,j}(m,n) * ...
                                (Deag_M_R_GMPM_Occur(m,n,j,s,i)) ./ (Total_Deag_M_R_GMPM_Occur(m,n,j,s)); %  Total_Deag_M_R_GMPM_Occur(m,n,j,s) = sum(Deag_M_R_GMPM_Occur(m,n,j,s,:))
                        else
                            alpha_faults(m,n,jj,i,j,s) = 0;
                        end

                        fval = f_IM_A_deag_ext{jj,i,j}(s,n,m);
                        f_IM_A_ave_branch_faults_GMPM_deag_lin(s,n,m)    = f_IM_A_ave_branch_faults_GMPM_deag_lin(s,n,m)  + fval * Branch_weight_AS(i,jj) * GMPM_weight(j) * alpha_faults(m,n,jj,i,j,s);
                        EN_A_ave_branch_faults_GMPM_deag(m,n,s) = EN_A_ave_branch_faults_GMPM_deag(m,n,s) + Branch_weight_AS(i,jj) * GMPM_weight(j) * alpha_faults(m,n,jj,i,j,s);
                    end
                end

            end
            if Total_Deag_M_R_Occure(m,n,s) ~= 0
                f_IM_A_ave_branch_faults_GMPM_deag_lin(s,n,m) = f_IM_A_ave_branch_faults_GMPM_deag_lin(s,n,m) ./ EN_A_ave_branch_faults_GMPM_deag(m,n,s);
            end
        end
        v = squeeze(EN_A_ave_branch_faults_GMPM_deag(m,n,:));
        EN_A_mean(m,n) = mean(v(v ~= 0));
    end
end

%% -----------------------------------------------------------------------
%  Linear average over AS branches, faults, & GMPMs
% ------------------------------------------------------------------------
disp('***** Deag M & R for all IM (Linear weighted averaging) *****');

PS_Method_deag_Lin      = zeros(Nd+1,Nd+1);
PS_E_onlymain_deag_Lin  = zeros(Nd+1,Nd+1);
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
    PS_Method_deag_Lin     = PS_Method_deag_Lin     + PE_IMe_WoAS(:,:,se) * sum_Mi_Rj * fx_IM_total(se);
    PS_E_onlymain_deag_Lin = PS_E_onlymain_deag_Lin + PE_IMe_WoAS(:,:,se) * fx_IM_total(se);
end

P_E_deag_Lin            = PS_Method_deag_Lin + (1-PS_Method_deag_Lin(end,end))*eye(Nd+1,Nd+1)
P_E_OnlyMainshock_deag_Lin = PS_E_onlymain_deag_Lin + (1-PS_E_onlymain_deag_Lin(end,end))*eye(Nd+1,Nd+1)
vpa(P_E_deag_Lin,10)
vpa(P_E_OnlyMainshock_deag_Lin,10)

toc;
