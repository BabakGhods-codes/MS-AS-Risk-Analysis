
clc; clear; tic;

%% =======================================================================
%  Multi-Branch MS–AS Risk Analysis with Iervolino et. al 2020 Paper Method 
% ========================================================================

addpath('Function Folder');

warning('BSSA_2014 GMPM: Z1.0 from CY model is used (as in original).');

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
dm        = 0.5;     % mag step
Mmin      = 5;       % min mag
dL        = 2.5;     % along-strike step (km)

% -----------------------------
% Aftershock window and mesh
% -----------------------------
DT  = 60;            % days window
dma = dm;            % AS mag step (kept equal to MS)
dra = dL;            % AS distance step (kept equal to MS)

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
% Main loop over MS logic-tree branches (per-fault PSHA & deaggregation)
% -----------------------------------------------------------------------
% Preallocations (cell because sizes differ by branch/fault)
R_vec               = cell(n_faults,1);
R_mid_vec_fault     = cell(n_faults,1);
R_mid_total         = cell(n_faults,n_branch_MS);
epicentral_location = cell(n_faults,1);
M_vec               = cell(n_faults,n_branch_MS);
FM                  = cell(n_faults,n_branch_MS);
M_mid_total         = cell(n_faults,n_branch_MS);
P_M                 = cell(n_faults,n_branch_MS);

fx_IM               = cell(n_faults,n_branch_MS,n_GMPM);
f_IM_E              = cell(n_faults,n_branch_MS,n_GMPM);
exceedance_rate_IM  = cell(n_faults,n_branch_MS,n_GMPM);
Deag_M_R_GMPM_fault = cell(n_faults,n_branch_MS,n_GMPM);
Deag_M_R_Exceed     = cell(n_faults,n_branch_MS,n_GMPM);

min_r_vec           = zeros(1,n_faults);
max_r_vec           = zeros(1,n_faults);
% Main loop over faults and MS branches

% Main loop over faults and MS branches
for i = 1:n_faults
    % Along fault mesh and geometric distances to site
    R_vec{i} = dL/2 : dL : Fault_L(i);
    R_mid_vec_fault{i} = sqrt(H(i)^2 + (R_vec{i} - Site_alongLengthDist(i)).^2);
    epicentral_location{i} = [R_vec{i}; H(i)*ones(size(R_vec{i}))];

    for ii = 1:n_branch_MS
        % Magnitude distribution on [Mmin, Mmax(i,ii)]
        M_vec{i,ii}  = Mmin : dm : Mmax(i,ii);
        FM{i,ii}     = (1-10.^(-b_M(i,ii)*(M_vec{i,ii}-Mmin))) / (1-10^(-b_M(i,ii)*(Mmax(i,ii)-Mmin))); % CDF
        P_M{i,ii}    = diff(FM{i,ii});  % Prob masses
        M_mid_total{i,ii} = (M_vec{i,ii}(2:end)+M_vec{i,ii}(1:end-1))/2;

        for j = 1:n_GMPM
            [f_IM_E{i,ii,j}, fx_IM{i,ii,j}, exceedance_rate_IM{i,ii,j}, ...
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
% Logic Tree Hazard (Total Hazard)
% -----------------------------------------------------------------------
exceedance_rate_IM_total = zeros(size(exceedance_rate_IM{i,1}));
for i = 1:n_faults
    for j = 1:n_GMPM
        for ii = 1:n_branch_MS
            exceedance_rate_IM_total = exceedance_rate_IM_total + ...
                Branch_weight_MS(i,ii) * GMPM_weight(j) * exceedance_rate_IM{i,ii,j};
        end
    end
end
fx_IM_total = -diff(exceedance_rate_IM_total);

%% Saving The Hazard Curve
saving_folder = 'Outputs/Final5/HazardCurves';
mkdir(saving_folder);

writematrix(exceedance_rate_IM_total, fullfile(saving_folder, 'exceedance_rate_IM_total_Site1_4Story.txt'));
writematrix(IM_vec_edge, fullfile(saving_folder, 'IM_vec_edge_Site1_4Story.txt'));
%% -----------------------------------------------------------------------
% AS Transition: Serial or Parallel workflow (Iervolino et al., 2020)
% -----------------------------------------------------------------------
if strcmp(AS_Analysis, 'SERIAL')
    % Preallocate cell containers (parfor-safe)
    PS            = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    PS2_onlymain  = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    PS_MeRe       = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    f_IM_A        = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    EN_A          = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);

    for i = 1:n_faults
        for j = 1:n_GMPM
            for ii = 1:n_branch_MS
                for jj = 1:n_branch_AS
                    [PS{jj,ii,i,j}, PS2_onlymain{jj,ii,i,j}, PS_MeRe{jj,ii,i,j}, ...
                     f_IM_A{jj,ii,i,j}, EN_A{jj,ii,i,j}] = ...
                        MS_AS_sequence_transition_main( ...
                            IM_vec, epicentral_location{i}, M_mid_total{i,ii}, M_vec{i,ii}, ...
                            Mmin, dma, b(i,jj), a(i,jj), c(i,jj), p(i,jj), IM_vec_edge, ...
                            Mu_mat, Beta_mat, DT,  ...
                            f_IM_E{i,ii,j}, landa(i,ii), P_M{i,ii}, GMPM_types{j}, ...
                            Vs30, HW(i), T, region, Fault_Type_name{i}, ...
                            dra, Dip(i,ii), FVS30, fas, H(i), ...
                            AS_Spatial_Distribution(i,jj), Site_alongLengthDist(i) );
                end
            end
        end
    end

elseif strcmp(AS_Analysis, 'PARALLEL')

    % % Build index combinations
    count = 0;
    for i = 1 : n_faults
        for j = 1 : n_GMPM
            for ii = 1 : n_branch_MS
                for jj = 1 : n_branch_AS
                    count = count + 1;
                    idx_map(count,:) = [i, j, ii, jj];
                end
            end
        end
    end
    
    % Create struct array to hold results
    results(count) = struct();
    
    % Start pool if not started
    if isempty(gcp('nocreate'))
        parpool('local', 8);  % Adjust core count
    end
    
    % Parallel loop
    parfor task = 1:count
        i  = idx_map(task,1);
        j  = idx_map(task,2);
        ii = idx_map(task,3);
        jj = idx_map(task,4);
        
        [ps, ps2, psmere, fima, ena] = MS_AS_sequence_transition_main(IM_vec, epicentral_location{i}, M_mid_total{i,ii}, M_vec{i,ii}, Mmin, dma, b(i,jj), a(i,jj), c(i,jj), p(i,jj), IM_vec_edge, Mu_mat, Beta_mat, DT, f_IM_E{i,ii,j}, landa(i,ii), P_M{i,ii}, GMPM_types{j}, Vs30, HW(i), T, region, Fault_Type_name{i}, dra, Dip(i,ii), FVS30, fas, H(i), AS_Spatial_Distribution(i,jj), Site_alongLengthDist(i));
    
        
        results(task).i = i;
        results(task).j = j;
        results(task).ii = ii;
        results(task).jj = jj;
        results(task).PS = ps;
        results(task).PS2 = ps2;
        results(task).PSA_MeRe = psmere;
        results(task).f_IM_A = fima;
        results(task).EN_A = ena;
    end
    
    % Preallocate output cells
    PS = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    PS2_onlymain = PS;
    PSA_MeRe = PS;
    f_IM_A = PS;
    EN_A = PS;
    
    for k = 1:length(results)
        i  = results(k).i;
        j  = results(k).j;
        ii = results(k).ii;
        jj = results(k).jj;
        
        PS{jj,ii,i,j} = results(k).PS;
        PS2_onlymain{jj,ii,i,j} = results(k).PS2;
        PSA_MeRe{jj,ii,i,j} = results(k).PSA_MeRe;
        f_IM_A{jj,ii,i,j} = results(k).f_IM_A;
        EN_A{jj,ii,i,j} = results(k).EN_A;
    end


end

%% -----------------------------------------------------------------------
% Weighted transition matrices (Iervolino et al., 2020)
% -----------------------------------------------------------------------
PS_total          = zeros(size(PS{1,1,1,1}));
PS_onlymain_total = zeros(size(PS{1,1,1,1}));

for i = 1:n_faults
    for j = 1:n_GMPM
        for ii = 1:n_branch_MS
            for jj = 1:n_branch_AS
                w = Branch_weight_MS(i,ii) * Branch_weight_AS(i,jj) * GMPM_weight(j);
                PS_total          = PS_total          + w * PS{jj,ii,i,j};
                PS_onlymain_total = PS_onlymain_total + w * PS2_onlymain{jj,ii,i,j};
            end
        end
    end
end

P_E = PS_total + (1 - PS_total(end, end)) * eye(Nd + 1, Nd + 1)
vpa(PS_onlymain_total)
vpa(P_E)
time = toc/3600

%% -----------------------------------------------------------------------
% Saving
% -----------------------------------------------------------------------
saving_folder = 'Outputs/Final5/Iervolino_LogicTree_Site1_8story_dl_2o5_dm_0o05_IM3o5_NewFaults';
mkdir(saving_folder);

writematrix(PS_onlymain_total, fullfile(saving_folder, 'PS_onlymain_total.txt'));
writematrix(PS_total,          fullfile(saving_folder, 'PS_total.txt'));

time = toc/3600;

save('f_IM_A.mat','f_IM_A');
save('EN_A.mat','EN_A');
save('PSA_MeRe.mat','PSA_MeRe');



