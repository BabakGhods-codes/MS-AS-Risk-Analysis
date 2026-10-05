function [PS_E, PS2_onlymain, PAS_MeRe, PA_ij, f_IM_A_out, EN_A, PE_ij] = ...
    MS_AS_sequence_transition_given_fIMA( ...
    IM_vec, M_mid_vec, Mmin, ...
    a, b, c, p, DT, ...
    Mu_mat, Beta_mat, ...
    f_IM_E, f_IM_A, ...
    activity_rate, P_M)

% MS_AS_SEQUENCE_TRANSITION_GIVEN_FIMA
%
% This function performs the MS-AS sequence transition calculations
% after f_IM_A has already been calculated.
%
% It does NOT:
%   - build AS magnitude distributions
%   - build AS spatial distributions
%   - calculate AS GMPMs
%   - calculate f_IM_A
%
% It only does:
%   1) PA_ij     = transition matrix for one generic aftershock
%   2) PE_ij     = transition matrix for the mainshock
%   3) EN_A      = expected number of aftershocks
%   4) PAS_MeRe  = aftershock-sequence transition matrix
%   5) PS_E      = total MS-AS sequence transition matrix
%   6) PS2_onlymain = mainshock-only transition matrix
%
% Dimensions expected:
%
%   f_IM_A : [n_IM x n_hyp x n_M]
%   f_IM_E : [n_IM x n_hyp x n_M]
%
% where:
%   n_IM  = length(IM_vec)
%   n_hyp = number of MS rupture locations
%   n_M   = number of MS magnitude bins

%% ------------------------------------------------------------------------
% Basic checks and sizes
% -------------------------------------------------------------------------

f_IM_A_out = f_IM_A;

Nd = size(Mu_mat,1);

n_IM = length(IM_vec);
n_hyp = size(f_IM_A,2);
n_M = size(f_IM_A,3);

if size(f_IM_A,1) ~= n_IM
    error('f_IM_A first dimension must be length(IM_vec).')
end

if size(f_IM_E,1) ~= n_IM
    error('f_IM_E first dimension must be length(IM_vec).')
end

if size(f_IM_E,2) ~= n_hyp
    error('f_IM_E and f_IM_A must have the same number of rupture-location bins.')
end

if size(f_IM_E,3) ~= n_M
    error('f_IM_E and f_IM_A must have the same number of magnitude bins.')
end

if length(M_mid_vec) ~= n_M
    error('length(M_mid_vec) must match size(f_IM_A,3).')
end

if length(P_M) ~= n_M
    error('length(P_M) must match size(f_IM_A,3).')
end

%% ========================================================================
% P_A | M_E, R_E
% Eqs. 15 and 14
% Transition matrix given one generic aftershock
% ========================================================================

PA_ij = zeros(Nd+1, Nd+1, n_hyp, n_M);

for m = 1:n_M

    for n = 1:n_hyp

        for i = 1:Nd

            for j = i:Nd

                PA_ij(i,j+1,n,m) = PDs_func( ...
                    Mu_mat, Beta_mat, ...
                    i, j, ...
                    IM_vec, f_IM_A(:,n,m));

            end

            PA_ij(i,i,n,m) = 1 - sum(PA_ij(i,i+1:end,n,m));

        end

        PA_ij(end,end,n,m) = 1;

    end

end

%% ========================================================================
% P_E | M_E, R_E
% Eqs. 19 and 18
% Transition matrix due to the mainshock event
% ========================================================================

PE_ij = zeros(Nd+1, Nd+1, n_hyp, n_M);

for m = 1:n_M

    for n = 1:n_hyp

        PE_ij(:,:,n,m) = eye(Nd+1, Nd+1);

        for i = 1:Nd

            for j = i:Nd

                PE_ij(i,j+1,n,m) = PDs_func( ...
                    Mu_mat, Beta_mat, ...
                    i, j, ...
                    IM_vec, f_IM_E(:,n,m));

            end

            PE_ij(i,i,n,m) = 1 - sum(PE_ij(i,i+1:end,n,m));

        end

        PE_ij(end,end,n,m) = 1;

    end

end

%% ========================================================================
% P_AS | M_E, R_E
% Eq. 20
% Transition matrix for several possible aftershocks
% ========================================================================

PAS_MeRe = zeros(size(PA_ij));
EN_A = zeros(n_hyp, n_M);

for m = 1:n_M

    for n = 1:n_hyp

        if p ~= 1

            EN_A(n,m) = ((p-1)^-1) * ...
                (10^(a + b*(M_mid_vec(m)-Mmin)) - 10^a) * ...
                (c^(1-p) - (DT+c)^(1-p));

        else

            term1 = 10^(a + b*(M_mid_vec(m)-Mmin)) - 10^a;
            term2 = log((DT + c) / c);

            EN_A(n,m) = term1 * term2;

        end

        PAS_MeRe(:,:,n,m) = PA_ij(:,:,n,m)^EN_A(n,m);

    end

end

%% ========================================================================
% PS_E and PS2_onlymain
% Eqs. 21 and 22
% Unit-time damage transition probability matrix
% ========================================================================

PS_E = zeros(Nd+1, Nd+1);
PS2_onlymain = zeros(Nd+1, Nd+1);

for m = 1:n_M

    for n = 1:n_hyp

        weight_MR = activity_rate * P_M(m) * 1/n_hyp;

        PS_E = PS_E + ...
            weight_MR * PE_ij(:,:,n,m) * PAS_MeRe(:,:,n,m);

        PS2_onlymain = PS2_onlymain + ...
            weight_MR * PE_ij(:,:,n,m);

    end

end

end