function [PS_E, PS2_onlymain, PAS_MeRe, PA_ij, f_IM_A_out, EN_A, PE_ij] = ...
    MS_AS_sequence_transition_given_fIMA( ...
    IM_vec, M_mid_vec, Mmin, ...
    a, b, c, p, DT, ...
    Mu_mat, Beta_mat, ...
    f_IM_E, f_IM_A, ...
    activity_rate, P_M)

% MS-AS transitions using f_IM_A already calculated outside.
% f_IM_A and f_IM_E: [n_IM x n_hyp x n_M]
% n_IM = length(IM_vec), n_hyp = MS rupture locations,
% n_M = number of MS magnitude bins.
% PA_ij: one AS, PE_ij: MS, EN_A: expected AS count
% PAS_MeRe: AS sequence, PS_E: full MS-AS sequence
% PS2_onlymain: MS only. f_IM_A_out returns input f_IM_A.
% GMPMs and AS distributions are inputs, no recalculation here.

%% sizes and checks

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

%% P_A given M_E, R_E
% Eqs. 15 and 14
% one generic AS transition

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

%% P_E given M_E, R_E
% Eqs. 19 and 18
% MS transition

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

%% P_AS given M_E, R_E
% Eq. 20
% AS sequence transition

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

%% combine MS and AS transitions
% Eqs. 21 and 22
% unit-time damage transitions

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
