function AS_fIMA = MS_AS_build_f_IM_A( ...
    MAG_MS, RS_MS, HYP_MS, ...
    AS_MAG, AS_RS, AS_HYP, AS_row_ids, ...
    ic_AS_all, LnIM_AS, sigma_AS, ...
    IM_vec_edge, IM_vec)

% MS_AS_BUILD_F_IM_A
%
% This function does NOT calculate GMPMs.
% This function assumes that the following have already been calculated once:
%
%   [uniqueRows_AS, ~, ic_AS_all] = unique(All_AS_GMPM_input_mat, 'rows');
%   [IM_AS, sigma_AS, LnIM_AS] = calculate_GMPM_for_unique_rows(...);
%
% This function only builds:
%
%   PA_Ex_IM
%   f_IM_A
%
% by using:
%
%   ic_AS_all
%   LnIM_AS
%   sigma_AS
%
% Dimensions:
%
%   f_IM_A(s, n, m, k, g)
%
% where:
%   s = IM midpoint index
%   n = mainshock hypocenter / rupture-location index
%   m = mainshock magnitude-bin index
%   k = mainshock rupture-sampling index
%   g = GMPM index

%% ------------------------------------------------------------------------
% Basic sizes
% -------------------------------------------------------------------------

n_IM_edge = length(IM_vec_edge);
n_IM_mid  = length(IM_vec);

nMS_M        = MAG_MS.nM - 1;
nMS_hyp      = HYP_MS.n_hyp;
nMS_sampling = RS_MS.n_sampling;

n_GMPM_AS = size(LnIM_AS,2);

if n_GMPM_AS == 1
    LnIM_AS = LnIM_AS(:);
    sigma_AS = sigma_AS(:);
end

%% ------------------------------------------------------------------------
% Preallocation
% -------------------------------------------------------------------------

PA_Ex_IM = zeros(n_IM_edge, nMS_hyp, nMS_M, nMS_sampling, n_GMPM_AS);
f_IM_A   = zeros(n_IM_mid,  nMS_hyp, nMS_M, nMS_sampling, n_GMPM_AS);

%% ------------------------------------------------------------------------
% Build f(IM_A | M_E, rupture_E)
% -------------------------------------------------------------------------

for g = 1 : n_GMPM_AS

    for m = 1 : nMS_M

        MAG_AS = AS_MAG{m};
        RS_AS  = AS_RS{m};

        n_AS_M        = MAG_AS.nM - 1;
        n_AS_sampling = RS_AS.n_sampling;

        for n = 1 : nMS_hyp

            for k = 1 : nMS_sampling

                HYP_AS = AS_HYP{m,n,k};

                n_AS_hyp = HYP_AS.n_hyp;

                row_ids_local = AS_row_ids{m,n,k};
                ic_AS_local   = ic_AS_all(row_ids_local);

                for s = 1 : n_IM_edge

                    for iA = 1 : n_AS_M

                        for jA = 1 : n_AS_hyp

                            for iiA = 1 : n_AS_sampling

                                id_local = ...
                                    iiA + ...
                                    (jA-1)*n_AS_sampling + ...
                                    (iA-1)*n_AS_sampling*n_AS_hyp;

                                id_runned_AS = ic_AS_local(id_local);

                                PA_Ex_IM(s,n,m,k,g) = PA_Ex_IM(s,n,m,k,g) + ...
                                    MAG_AS.P_M(iA) * ...
                                    1/n_AS_hyp * ...
                                    RS_AS.P_rup_W(iA,iiA) * ...
                                    (1 - normcdf( ...
                                    (log(IM_vec_edge(s))-LnIM_AS(id_runned_AS,g)) / ...
                                    sigma_AS(id_runned_AS,g)));

                            end

                        end

                    end

                end

                f_IM_A(:,n,m,k,g) = diff(PA_Ex_IM(:,n,m,k,g)) * -1;

            end

        end

    end

end

%% ------------------------------------------------------------------------
% Store outputs
% -------------------------------------------------------------------------

AS_fIMA = struct();

AS_fIMA.PA_Ex_IM = PA_Ex_IM;
AS_fIMA.f_IM_A   = f_IM_A;

% Convenience output if mainshock n_sampling = 1
if nMS_sampling == 1

    AS_fIMA.PA_Ex_IM_4D = PA_Ex_IM(:,:,:,1,:);
    AS_fIMA.f_IM_A_4D   = f_IM_A(:,:,:,1,:);

    for g = 1 : n_GMPM_AS
        AS_fIMA.PA_Ex_IM_3D_by_GMPM{g} = PA_Ex_IM(:,:,:,1,g);
        AS_fIMA.f_IM_A_3D_by_GMPM{g}   = f_IM_A(:,:,:,1,g);
    end

end

AS_fIMA.n_GMPM_AS = n_GMPM_AS;

end