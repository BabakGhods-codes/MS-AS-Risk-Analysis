clc; clear; tic;

%% =======================================================================
%  Multi-Branch MS–AS Risk Analysis with Irevolino et. al 2020 Paper Method 
% ========================================================================

addpath('Function Folder');

warning('BSSA_2014 GMPM: Z1.0 from CY model is used (as in original).');
saving_folder = 'Outputs/Final8/Irevalino_LogicTree_Site1_4story_hyp2o5_dra2o5_dgeom0o5_dm_0o1_IM4_NewFaults';

%% -----------------------------------------------------------------------
%  INPUTS
% ------------------------------------------------------------------------

Mmax = [6.6   6.6   7  7;   7.0 7.0 7.4 7.4];
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
Lambda = [0 0 0 0; pi/2 pi/2 pi/2 pi/2];
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
for i = 1 : n_faults
i
    for ii = 1 : n_branch_MS
        ii
        % Magnitude distribution on [Mmin, Mmax(i,ii)]
        M_vec{i,ii}  = Mmin : dm : Mmax(i,ii);
        FM{i,ii}     = (1-10.^(-b_M(i,ii)*(M_vec{i,ii}-Mmin))) / (1-10^(-b_M(i,ii)*(Mmax(i,ii)-Mmin))); % CDF
        P_M{i,ii}    = diff(FM{i,ii});  % Prob masses
        M_mid_total{i,ii} = (M_vec{i,ii}(2:end)+M_vec{i,ii}(1:end-1))/2;
        
        % FG      % fault geometry and site-distance grid
        % MAG_MS  % magnitude distribution
        % RS_MS   % rupture scaling
        % HYP_MS  % hypocenter grid
        % RUP_MS  % placed ruptures and GMPM input matrix

        FG{i,ii} = build_fault_geometry( ...
            Fault_Type_name{i}, ...
            Fault_W(i), Fault_L(i), Fault_Hf(i), Dip(i,ii), dL_geom, ...
            Fault_X0(i), Fault_Y0(i), Site_loc);
        
        MAG_MS{i,ii} = build_magnitude_distribution( ...
            Mmin, Mmax(i,ii), b_M(i,ii), dm);
        
        RS_MS{i,ii} = build_rupture_scaling( ...
            Fault_Type_name{i}, ...
            MAG_MS{i,ii}.M_mid_vec, ...
            variation_in_rupture_area, n_sampling);
        
       
        HYP_MS{i,ii} = build_hypocenter_grid( ...
            0, Fault_L(i), ...
            0, Fault_W(i), ...
            dL_hyp, Fault_X0(i), Fault_Y0(i), Fault_Hf(i), Dip(i,ii));
        
        RUP_MS{i,ii} = place_ruptures_on_fault( ...
            FG{i,ii}, MAG_MS{i,ii}, RS_MS{i,ii}, HYP_MS{i,ii}, ...
            rupture_placement_type);
                

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%% Calculate GMPMs for unique rows %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        [uniqueRows, ~, ic] = unique(RUP_MS{i,ii}.GMPM_input_mat, 'rows');

        [IM, sigma, LnIM] = calculate_GMPM_for_unique_rows( ...
            uniqueRows, n_GMPM, GMPM_types, ...
            Vs30, Z1o0, Z2o5, Region, Fault_Type_name{i}, ...
            FVS30, fas_MS, HW(i), T, Dip(i,ii), Lambda(i,ii));
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% MAIN PSHA LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Rrup_vec = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200];%[0 min(min(Rrup))-dL/2: dL : max(max(Rrup))+dL];
        Rrup_mid_vec = (Rrup_vec(1:end-1) + Rrup_vec(2:end))/2;
        nR = length(Rrup_mid_vec);
        Deag_M_R = zeros(MAG_MS{i,ii}.nM-1, nR, n_IM-1, n_GMPM);
        % Deag_M_R_Exceed = zeros(FF{i,ii}.nM-1, nR, n_IM, n_GMPM); 
        % Rjb_vec = min(min(RUP_MS{i,ii}.Rjb))-dL/2: dL : max(max(RUP_MS{i,ii}.Rjb))+dL/2;
        % exceedance_rate_IM = zeros(n_IM, n_GMPM);
        % exceedance_rate_IM_midValues = zeros(n_IM-1, n_GMPM);
        % f_IM_E = zeros(n_IM-1, FF{i,ii}.n_hyp, FF{i,ii}.nM-1, n_GMPM);
        warning('The following Lines should be checked if "n_sampling" larger than 1 is used')

       
        PSHA{i,ii} = MSAS_main_PSHA_loop_original( ...
            LnIM, sigma, ic, ...
            activity(i,ii), MAG_MS{i,ii}.P_M, RS_MS{i,ii}.P_rup_W, ...
            RUP_MS{i,ii}.Rrup, Rrup_vec, ...
            IM_vec_edge, IM_vec, ...
            n_GMPM, MAG_MS{i,ii}.nM, HYP_MS{i,ii}.n_hyp, n_sampling);

        for j = 1 : n_GMPM
            f_IM_E{i,ii,j}               =  PSHA{i,ii}.f_IM_E(:,:,:,j);
            fx_IM{i,ii,j}                =  PSHA{i,ii}.fx_IM(:,j);   
            exceedance_rate_IM{i,ii,j}   =  PSHA{i,ii}.exceedance_rate_IM(:,j);   
            Deag_M_R_GMPM_fault{i,ii,j}  =  PSHA{i,ii}.Deag_M_R(:,:,:,j);
            Deag_M_R_Exceed{i,ii,j}      =  PSHA{i,ii}.Deag_M_R_Exceed(:,:,:,j);
        end

   
    
     end
    
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
saving_folder = 'Outputs/Final8/HazardCurves';
mkdir(saving_folder);

writematrix(exceedance_rate_IM_total, fullfile(saving_folder, 'exceedance_rate_IM_total_Site1_4Story.txt'));
writematrix(IM_vec_edge, fullfile(saving_folder, 'IM_vec_edge_Site1_4Story.txt'));

%% =======================================================================
% Build all AS GMPM input rows and calculate AS LnIM/sigma
% One calculation per fault i and MS branch ii
% ========================================================================

% All_AS_GMPM_input_mat = cell(n_faults, n_branch_MS);

AS_row_ids = cell(n_faults, n_branch_MS);
AS_MAG     = cell(n_faults, n_branch_MS);
AS_RS      = cell(n_faults, n_branch_MS);
AS_HYP     = cell(n_faults, n_branch_MS);
AS_RUP     = cell(n_faults, n_branch_MS);

uniqueRows_AS = cell(n_faults, n_branch_MS);
ic_AS_all     = cell(n_faults, n_branch_MS);

IM_AS    = cell(n_faults, n_branch_MS);
sigma_AS = cell(n_faults, n_branch_MS);
LnIM_AS  = cell(n_faults, n_branch_MS);

for i = 1:n_faults

    i

    for ii = 1:n_branch_MS

        ii

        %% ---------------------------------------------------------------
        % Local AS containers for this fault and this MS branch
        % ---------------------------------------------------------------

        All_AS_GMPM_input_mat_local = [];

        nMS_M        = MAG_MS{i,ii}.nM - 1;
        nMS_hyp      = HYP_MS{i,ii}.n_hyp;
        nMS_sampling = RS_MS{i,ii}.n_sampling;

        AS_row_ids_local = cell(nMS_M, nMS_hyp, nMS_sampling);
        AS_MAG_local     = cell(nMS_M, 1);
        AS_RS_local      = cell(nMS_M, 1);
        AS_HYP_local     = cell(nMS_M, nMS_hyp, nMS_sampling);
        AS_RUP_local     = cell(nMS_M, nMS_hyp, nMS_sampling);

        %% ---------------------------------------------------------------
        % Loop over mainshock magnitude bins
        % ---------------------------------------------------------------

        for m = 1:nMS_M
if m == nMS_M
    stop = 1
end
            ME = (MAG_MS{i,ii}.M_vec(m+1) + MAG_MS{i,ii}.M_vec(m)) / 2;

            %% -----------------------------------------------------------
            % AS magnitude distribution
            % Common AS magnitude grid; branch-specific probabilities are applied later
            % -----------------------------------------------------------
            M_vec_A = Mmin + (0:floor((ME - Mmin)/dma)) * dma;

            tol = 1e-12 * ME;

            if abs(M_vec_A(end) - ME) > tol
                M_vec_A = [M_vec_A, ME];
            end

            M_mid_vec_A = ...
                (M_vec_A(2:end) + M_vec_A(1:end-1)) / 2;

            MAG_AS_m = struct();

            MAG_AS_m.Mmin = Mmin;
            MAG_AS_m.Mmax = ME;
            MAG_AS_m.dm   = dma;

            MAG_AS_m.M_vec     = M_vec_A;
            MAG_AS_m.M_mid_vec = M_mid_vec_A;
            MAG_AS_m.nM        = length(M_vec_A);

            AS_MAG_local{m} = MAG_AS_m;
            
           

            %% -----------------------------------------------------------
            % AS rupture scaling
            % -----------------------------------------------------------

            RS_AS_m = build_rupture_scaling( ...
                Fault_Type_name{i}, ...
                MAG_AS_m.M_mid_vec, ...
                variation_in_rupture_area, ...
                n_sampling);

            AS_RS_local{m} = RS_AS_m;

            %% -----------------------------------------------------------
            % Loop over each MS rupture
            % -----------------------------------------------------------

            for n = 1:nMS_hyp

                for k = 1:nMS_sampling

                    %% Mainshock rupture limits
                    L_MS_0 = RUP_MS{i,ii}.rup_lim_L(n,1,k,m);
                    L_MS_1 = RUP_MS{i,ii}.rup_lim_L(n,2,k,m);

                    W_MS_0 = RUP_MS{i,ii}.rup_lim_W(n,1,k,m);
                    W_MS_1 = RUP_MS{i,ii}.rup_lim_W(n,2,k,m);

                    %% AS hypocenters inside this MS rupture plane
                    HYP_AS_mnk = build_hypocenter_grid( ...
                        L_MS_0, L_MS_1, ...
                        W_MS_0, W_MS_1, ...
                        dra, Fault_X0(i), Fault_Y0(i), Fault_Hf(i), Dip(i,ii));

                    AS_HYP_local{m,n,k} = HYP_AS_mnk;

                    if HYP_AS_mnk.n_hyp == 0
                        error('No AS hypocenter points generated. Reduce dra or check MS rupture size.')
                    end

                    %% AS rupture-placement domain

                    HYP_AS_place = HYP_AS_mnk;

                    switch lower(AS_rupture_domain)

                        case 'ms_rupture'

                            % Current model -- no modification

                        case 'full_fault'

                            HYP_AS_place.L_start = 0;
                            HYP_AS_place.L_end   = Fault_L(i);

                            HYP_AS_place.W_start = 0;
                            HYP_AS_place.W_end   = Fault_W(i);

                        otherwise

                            error(['Unknown AS_rupture_domain: %s. Use ', ...
                                '''MS_rupture'' or ''full_fault''.'], ...
                                AS_rupture_domain);

                    end

                    %% AS finite ruptures
                    RUP_AS_mnk = place_ruptures_on_fault( ...
                            FG{i,ii}, MAG_AS_m, RS_AS_m, HYP_AS_place, ...
                            rupture_placement_type);

                    % AS_RUP_local{m,n,k} = RUP_AS_mnk;

                    %% ---------------------------------------------------
                    % Store LOCAL row IDs
                    % These row IDs refer to All_AS_GMPM_input_mat_local.
                    % ---------------------------------------------------

                    row_start = size(All_AS_GMPM_input_mat_local,1) + 1;
                    row_end   = row_start + size(RUP_AS_mnk.GMPM_input_mat,1) - 1;

                    AS_row_ids_local{m,n,k} = row_start:row_end;

                    All_AS_GMPM_input_mat_local = [ ...
                        All_AS_GMPM_input_mat_local; ...
                        RUP_AS_mnk.GMPM_input_mat];

                end

            end

        end

        %% ---------------------------------------------------------------
        % Store local AS objects
        % ---------------------------------------------------------------

        % All_AS_GMPM_input_mat{i,ii} = All_AS_GMPM_input_mat_local;

        AS_row_ids{i,ii} = AS_row_ids_local;
        AS_MAG{i,ii}     = AS_MAG_local;
        AS_RS{i,ii}      = AS_RS_local;
        AS_HYP{i,ii}     = AS_HYP_local;
        % AS_RUP{i,ii}     = AS_RUP_local;

        %% ---------------------------------------------------------------
        % Unique AS GMPM rows for this fault/MS branch only
        % ---------------------------------------------------------------

        % [uniqueRows_AS{i,ii}, ~, ic_AS_all{i,ii}] = ...
        %     unique(All_AS_GMPM_input_mat{i,ii}, 'rows');

        [uniqueRows_AS{i,ii}, ~, ic_AS_all{i,ii}] = ...
            unique(All_AS_GMPM_input_mat_local, 'rows');

        %% ---------------------------------------------------------------
        % Calculate all GMPMs once for this fault/MS branch
        % No loop over j is needed here.
        % ---------------------------------------------------------------

        [IM_AS{i,ii}, sigma_AS{i,ii}, LnIM_AS{i,ii}] = ...
            calculate_GMPM_for_unique_rows( ...
            uniqueRows_AS{i,ii}, n_GMPM, GMPM_types, ...
            Vs30, Z1o0, Z2o5, Region, Fault_Type_name{i}, ...
            FVS30, fas_AS, HW(i), T, Dip(i,ii), Lambda(i,ii));

    end

end
%%

%% -----------------------------------------------------------------------
% AS Transition: Serial or Parallel workflow (Irevolino et al., 2020)
% -----------------------------------------------------------------------
if strcmp(AS_Analysis, 'SERIAL')
    AS_fIMA = cell(n_faults, n_branch_MS, n_branch_AS);

    for i = 1:n_faults

        for ii = 1:n_branch_MS
    
            for jj = 1:n_branch_AS
    
                % Copy common AS magnitude grids
                AS_MAG_current = AS_MAG{i,ii};
    
                % Apply the probability distribution of THIS AS branch
                for m = 1:length(AS_MAG_current)
    
                    M_vec_A = AS_MAG_current{m}.M_vec;
                    ME      = AS_MAG_current{m}.Mmax;
    
                    b_current = b_AS(i,jj);
    
                    denom = ...
                        10.^(-b_current*Mmin) - ...
                        10.^(-b_current*ME);
    
                    FMA = ...
                        (10.^(-b_current*Mmin) - ...
                         10.^(-b_current*M_vec_A)) ./ denom;
    
                    AS_MAG_current{m}.b_M = b_current;
                    AS_MAG_current{m}.FM  = FMA;
                    AS_MAG_current{m}.P_M = diff(FMA);
    
                end
    
                AS_fIMA{i,ii,jj} = MS_AS_build_f_IM_A( ...
                    MAG_MS{i,ii}, RS_MS{i,ii}, HYP_MS{i,ii}, ...
                    AS_MAG_current, AS_RS{i,ii}, AS_HYP{i,ii}, AS_row_ids{i,ii}, ...
                    ic_AS_all{i,ii}, LnIM_AS{i,ii}, sigma_AS{i,ii}, ...
                    IM_vec_edge, IM_vec);
    
            end
    
        end

    end



        PS            = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
        PS2_onlymain  = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
        PAS_MeRe      = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
        PA_ij         = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
        f_IM_A        = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
        EN_A          = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
        PE_ij         = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);

        for i = 1:n_faults

            for ii = 1:n_branch_MS

                for j = 1:n_GMPM

                    

                    for jj = 1:n_branch_AS
                        % f_IM_A for this GMPM
                        f_IM_A_current = AS_fIMA{i,ii,jj}.f_IM_A_3D_by_GMPM{j};

                        [PS{jj,ii,i,j}, ...
                            PS2_onlymain{jj,ii,i,j}, ...
                            PAS_MeRe{jj,ii,i,j}, ...
                            PA_ij{jj,ii,i,j}, ...
                            f_IM_A{jj,ii,i,j}, ...
                            EN_A{jj,ii,i,j}, ...
                            PE_ij{jj,ii,i,j}] = ...
                            MS_AS_sequence_transition_given_fIMA( ...
                            IM_vec, MAG_MS{i,ii}.M_mid_vec, Mmin, ...
                            a(i,jj), b_AS(i,jj), c(i,jj), p(i,jj), DT, ...
                            Mu_mat, Beta_mat, ...
                            f_IM_E{i,ii,j}, f_IM_A_current, ...
                            activity(i,ii), MAG_MS{i,ii}.P_M);

                    end

                end

            end

        end

    


elseif strcmp(AS_Analysis, 'PARALLEL')

            % -------------------------------------------------------------------
    % Start parallel pool
    % --------------------------------------------------------------------
    pool = gcp('nocreate');

    if isempty(pool)
        parpool('local', n_Core_Parallel);
    end


    % -------------------------------------------------------------------
    % Build f(IM_A | ME, rupture_E) for each
    % fault / MS branch / AS branch
    % --------------------------------------------------------------------
    AS_fIMA = cell(n_faults, n_branch_MS, n_branch_AS);

    n_fIMA_tasks = n_faults * n_branch_MS * n_branch_AS;
    fIMA_results = cell(n_fIMA_tasks,1);

    parfor task = 1:n_fIMA_tasks

        [i, ii, jj] = ind2sub( ...
            [n_faults, n_branch_MS, n_branch_AS], task);

        % Common AS magnitude/geometry information for this MS branch
        AS_MAG_current = AS_MAG{i,ii};

        % Apply the b-value of the current independent AS branch
        for m = 1:length(AS_MAG_current)

            M_vec_A = AS_MAG_current{m}.M_vec;
            ME      = AS_MAG_current{m}.Mmax;

            b_current = b_AS(i,jj);

            denom = ...
                10.^(-b_current*Mmin) - ...
                10.^(-b_current*ME);

            FMA = ...
                (10.^(-b_current*Mmin) - ...
                 10.^(-b_current*M_vec_A)) ./ denom;

            AS_MAG_current{m}.b_M = b_current;
            AS_MAG_current{m}.FM  = FMA;
            AS_MAG_current{m}.P_M = diff(FMA);

        end

        fIMA_results{task} = MS_AS_build_f_IM_A( ...
            MAG_MS{i,ii}, RS_MS{i,ii}, HYP_MS{i,ii}, ...
            AS_MAG_current, AS_RS{i,ii}, AS_HYP{i,ii}, AS_row_ids{i,ii}, ...
            ic_AS_all{i,ii}, LnIM_AS{i,ii}, sigma_AS{i,ii}, ...
            IM_vec_edge, IM_vec);

    end

    for task = 1:n_fIMA_tasks

        [i, ii, jj] = ind2sub( ...
            [n_faults, n_branch_MS, n_branch_AS], task);

        AS_fIMA{i,ii,jj} = fIMA_results{task};

    end


    % -------------------------------------------------------------------
    % Build all fault / MS branch / GMPM / AS branch combinations
    % --------------------------------------------------------------------
    count = 0;
    idx_map = zeros( ...
        n_faults*n_branch_MS*n_GMPM*n_branch_AS, 4);

    for i = 1:n_faults
        for ii = 1:n_branch_MS
            for j = 1:n_GMPM
                for jj = 1:n_branch_AS

                    count = count + 1;
                    idx_map(count,:) = [i, ii, j, jj];

                end
            end
        end
    end


    % -------------------------------------------------------------------
    % Parallel transition calculations
    % --------------------------------------------------------------------
    results = cell(count,1);

    parfor task = 1:count

        i  = idx_map(task,1);
        ii = idx_map(task,2);
        j  = idx_map(task,3);
        jj = idx_map(task,4);

        f_IM_A_current = ...
            AS_fIMA{i,ii,jj}.f_IM_A_3D_by_GMPM{j};

        [PS_current, ...
         PS2_current, ...
         PAS_current, ...
         PA_current, ...
         fIMA_current, ...
         ENA_current, ...
         PE_current] = ...
            MS_AS_sequence_transition_given_fIMA( ...
                IM_vec, MAG_MS{i,ii}.M_mid_vec, Mmin, ...
                a(i,jj), b_AS(i,jj), c(i,jj), p(i,jj), DT, ...
                Mu_mat, Beta_mat, ...
                f_IM_E{i,ii,j}, f_IM_A_current, ...
                activity(i,ii), MAG_MS{i,ii}.P_M);

        out = struct();

        out.i  = i;
        out.ii = ii;
        out.j  = j;
        out.jj = jj;

        out.PS           = PS_current;
        out.PS2_onlymain = PS2_current;
        out.PAS_MeRe     = PAS_current;
        out.PA_ij        = PA_current;
        out.f_IM_A       = fIMA_current;
        out.EN_A         = ENA_current;
        out.PE_ij        = PE_current;

        results{task} = out;

    end

    


    %% -------------------------------------------------------------------
    % Reassemble outputs in same format as SERIAL workflow
    % --------------------------------------------------------------------
    PS            = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    PS2_onlymain  = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    PAS_MeRe      = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    PA_ij         = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    f_IM_A        = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    EN_A          = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);
    PE_ij         = cell(n_branch_AS, n_branch_MS, n_faults, n_GMPM);

    for task = 1:count

        out = results{task};

        i  = out.i;
        ii = out.ii;
        j  = out.j;
        jj = out.jj;

        PS{jj,ii,i,j}           = out.PS;
        PS2_onlymain{jj,ii,i,j} = out.PS2_onlymain;
        PAS_MeRe{jj,ii,i,j}     = out.PAS_MeRe;
        PA_ij{jj,ii,i,j}        = out.PA_ij;
        f_IM_A{jj,ii,i,j}       = out.f_IM_A;
        EN_A{jj,ii,i,j}         = out.EN_A;
        PE_ij{jj,ii,i,j}        = out.PE_ij;

    end



end

%% -----------------------------------------------------------------------
% Weighted transition matrices (Irevolino et al., 2020)
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





