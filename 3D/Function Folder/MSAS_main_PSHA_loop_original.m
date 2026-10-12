function PSHA = MSAS_main_PSHA_loop_original( ...
    LnIM, sigma, ic, ...
    activity, P_M, P_rup_W, ...
    Rrup, Rrup_vec, ...
    IM_vec_edge, IM_vec, ...
    n_GMPM, nM, n_hyp, n_sampling)
% PSHA using GMPMs calculated outside this function.
% LnIM, sigma: [unique input rows x n_GMPM], log median and log std.
% ic maps original input rows to unique rows.
% activity = annual source rate, P_M = mag bin probabilities
% P_rup_W: [nM-1 x n_sampling], rupture sample weights
% Rrup: [nM-1 x n_hyp x n_sampling]
% Rrup_vec = Rrup edges; IM_vec_edge = IM edges; IM_vec = midpoints
% nM is length(M_vec), n_hyp = hypocenter count, n_sampling = rupture samples.
% PSHA returns exceedance rates at edges and midpoints, fx_IM, f_IM_E,
% and occurrence/exceedance M-R deaggregation.

%%
%% sizes and arrays
%%

n_IM = length(IM_vec_edge);

Rrup_mid_vec = (Rrup_vec(1:end-1) + Rrup_vec(2:end))/2;
nR = length(Rrup_mid_vec);

Deag_M_R = zeros(nM-1, nR, n_IM-1, n_GMPM);
Deag_M_R_Exceed = zeros(nM-1, nR, n_IM, n_GMPM);

exceedance_rate_IM = zeros(n_IM, n_GMPM);
exceedance_rate_IM_midValues = zeros(n_IM-1, n_GMPM);

f_IM_E = zeros(n_IM-1, n_hyp, nM-1, n_GMPM);

warning('The following Lines should be checked if "n_sampling" larger than 1 is used')

%%
%% PSHA loop
%%

for iii = 1 : n_GMPM

    iii

    for s = 1 : n_IM

        for i = 1 : nM-1

            for j = 1 : n_hyp

                for ii = 1 : n_sampling

                    id_Rrup = sum((Rrup_vec < Rrup(i,j,ii)));

                    id_runned = ic(ii + (j-1)*n_sampling + (i-1)*n_sampling*n_hyp);

                    exceedance_rate_IM(s,iii) = exceedance_rate_IM(s,iii) + ...
                        activity * P_M(i) * 1/n_hyp * P_rup_W(i,ii) * ...
                        (1 - normcdf( ...
                        (log(IM_vec_edge(s))-LnIM(id_runned,iii)) / sigma(id_runned,iii)));

                    if s < n_IM

                        % Eq. 18, Iervolino et al. 2020
                        f_IM_E(s,j,i,iii) = f_IM_E(s,j,i,iii) + ...
                            ( -normcdf( ...
                            (log(IM_vec_edge(s))-LnIM(id_runned,iii)) / sigma(id_runned,iii)) ...
                            + normcdf( ...
                            (log(IM_vec_edge(s+1))-LnIM(id_runned,iii)) / sigma(id_runned,iii)) );

                        exceedance_rate_IM_midValues(s,iii) = exceedance_rate_IM_midValues(s,iii) + ...
                            activity * P_M(i) * 1/n_hyp * P_rup_W(i,ii) * ...
                            (1 - normcdf( ...
                            (log(IM_vec(s))-LnIM(id_runned,iii)) / sigma(id_runned,iii)));

                        % deaggregation for occurrence
                        Deag_M_R(i, id_Rrup, s, iii) = Deag_M_R(i, id_Rrup, s, iii) + ...
                            activity * P_M(i) * 1/n_hyp * P_rup_W(i,ii) * ...
                            ( -normcdf( ...
                            (log(IM_vec_edge(s))-LnIM(id_runned,iii)) / sigma(id_runned,iii)) ...
                            + normcdf( ...
                            (log(IM_vec_edge(s+1))-LnIM(id_runned,iii)) / sigma(id_runned,iii)) );

                    end

                    Deag_M_R_Exceed(i, id_Rrup, s, iii) = Deag_M_R_Exceed(i, id_Rrup, s, iii) + ...
                        activity * P_M(i) * 1/n_hyp * P_rup_W(i,ii) * ...
                        (1 - normcdf( ...
                        (log(IM_vec_edge(s))-LnIM(id_runned,iii)) / sigma(id_runned,iii)));

                end

            end

        end

    end

    fx_IM(:,iii) = diff(exceedance_rate_IM(:,iii)) * -1;

end

%%
%% normalize deaggregation
%%

for iii = 1 : n_GMPM

    iii

    for s = 1 : n_IM

        if s < n_IM

            if fx_IM(s,iii) > 0
                Deag_M_R(:, :, s, iii) = Deag_M_R(:, :, s, iii) / fx_IM(s,iii);
            end
        end

        if exceedance_rate_IM(s,iii) > 0
            Deag_M_R_Exceed(:,:,s,iii) = ...
                Deag_M_R_Exceed(:,:,s,iii) ./ exceedance_rate_IM(s,iii);
        else
            Deag_M_R_Exceed(:,:,s,iii) = 0;
        end

    end

end

%% save outputs

PSHA = struct();

PSHA.Rrup_vec = Rrup_vec;
PSHA.Rrup_mid_vec = Rrup_mid_vec;

PSHA.exceedance_rate_IM = exceedance_rate_IM;
PSHA.exceedance_rate_IM_midValues = exceedance_rate_IM_midValues;
PSHA.fx_IM = fx_IM;

PSHA.f_IM_E = f_IM_E;
PSHA.Deag_M_R = Deag_M_R;
PSHA.Deag_M_R_Exceed = Deag_M_R_Exceed;

end
