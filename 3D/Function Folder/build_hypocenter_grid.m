function HYP = build_hypocenter_grid( ...
    L_start, L_end, W_start, W_end, ...
    dL, Xs, Ys, Hf, Dip)

%%% Hypocenter Generation

L_domain = L_end - L_start;
W_domain = W_end - W_start;

if L_domain <= 0 || W_domain <= 0
    error('Hypocenter domain must have positive length and width.')
end

if dL <= 0
    error('Hypocenter spacing dL must be positive.')
end

% Number of cells required so that the entire domain is covered.
% The actual spacing is <= the requested dL.
nL = max(1, ceil(L_domain / dL));
nW = max(1, ceil(W_domain / dL));

% Cell edges.
% IMPORTANT: L_start/L_end and W_start/W_end are always included exactly.
L_vec_hyp = linspace(L_start, L_end, nL + 1);
W_vec_hyp = linspace(W_start, W_end, nW + 1);

% Uniform hypocenters at the centers of the cells.
Lhyp = 0.5 * (L_vec_hyp(1:end-1) + L_vec_hyp(2:end));
Whyp = 0.5 * (W_vec_hyp(1:end-1) + W_vec_hyp(2:end));

n_hyp = length(Whyp) * length(Lhyp);

Whyp_points_repmat = repmat(Whyp', 1, length(Lhyp));

hyp_points_loc = [ ...
    repmat(Lhyp', length(Whyp), 1), ...
    reshape(Whyp_points_repmat', [], 1)];

epi_points_loc = [ ...
    hyp_points_loc(:,1) + Xs, ...
    hyp_points_loc(:,2)*cos(Dip) + Ys];

hyp_points_loc_Global = [ ...
    hyp_points_loc(:,1) + Xs, ...
    hyp_points_loc(:,2)*cos(Dip) + Ys, ...
    -hyp_points_loc(:,2)*sin(Dip) - Hf];

%% Store outputs

HYP = struct();

HYP.L_start = L_start;
HYP.L_end   = L_end;
HYP.W_start = W_start;
HYP.W_end   = W_end;

HYP.dL  = dL;
HYP.Xs  = Xs;
HYP.Ys  = Ys;
HYP.Hf  = Hf;
HYP.Dip = Dip;

HYP.W_vec_hyp = W_vec_hyp;
HYP.L_vec_hyp = L_vec_hyp;

HYP.Whyp  = Whyp;
HYP.Lhyp  = Lhyp;
HYP.n_hyp = n_hyp;

HYP.hyp_points_loc        = hyp_points_loc;
HYP.epi_points_loc        = epi_points_loc;
HYP.hyp_points_loc_Global = hyp_points_loc_Global;

end

% 
% function HYP = build_hypocenter_grid( ...
%     L_start, L_end, W_start, W_end, ...
%     dL, Xs, Ys, Hf, Dip)
% 
% %%% Hypocenter Generation 
% W_vec_hyp = W_start : dL : W_end;
% L_vec_hyp = L_start : dL : L_end;
% 
% Whyp = (W_vec_hyp(1:end-1)+W_vec_hyp(2:end))/2; 
% Lhyp = (L_vec_hyp(1:end-1)+L_vec_hyp(2:end))/2;
% 
% n_hyp = length(Whyp)*length(Lhyp); 
% 
% Whyp_points_repmat = repmat(Whyp',1,length(Lhyp));
% 
% hyp_points_loc = [ ...
%     repmat(Lhyp',length(Whyp),1), ...
%     reshape(Whyp_points_repmat',[],1)];
% 
% epi_points_loc = [ ...
%     hyp_points_loc(:,1)+Xs , ...
%     hyp_points_loc(:,2)*cos(Dip) + Ys];
% 
% hyp_points_loc_Global = [ ...
%     hyp_points_loc(:,1)+Xs , ...
%     hyp_points_loc(:,2)*cos(Dip) + Ys, ...
%     -hyp_points_loc(:,2)*sin(Dip)-Hf];
% 
% %% Store outputs
% HYP = struct();
% 
% HYP.L_start = L_start;
% HYP.L_end = L_end;
% HYP.W_start = W_start;
% HYP.W_end = W_end;
% 
% HYP.dL = dL;
% HYP.Xs = Xs;
% HYP.Ys = Ys;
% HYP.Hf = Hf;
% HYP.Dip = Dip;
% 
% HYP.W_vec_hyp = W_vec_hyp;
% HYP.L_vec_hyp = L_vec_hyp;
% 
% HYP.Whyp = Whyp;
% HYP.Lhyp = Lhyp;
% HYP.n_hyp = n_hyp;
% 
% HYP.hyp_points_loc = hyp_points_loc;
% HYP.epi_points_loc = epi_points_loc;
% HYP.hyp_points_loc_Global = hyp_points_loc_Global;
% 
% end