function MAG = build_magnitude_distribution(Mmin, Mmax, b_M, dm)

%%% Mainshock Magnitudes or Aftershock Magnitudes
%%% Gutenberg-Truncated Relationship

M_vec = Mmin + (0:floor((Mmax-Mmin)/dm))*dm;

tol = 1e-12 * max(1,abs(Mmax));

if abs(M_vec(end)-Mmax) > tol
    M_vec = [M_vec, Mmax];
end

FM = (1-10.^(-b_M*(M_vec-Mmin))) / ...
    (1-10^(-b_M*(Mmax-Mmin))); % CDF of magnitude

M_mid_vec = (M_vec(2:end)+M_vec(1:end-1))/2;

P_M = FM(2:end)-FM(1:end-1);

nM = length(M_vec);

%% Store outputs
MAG = struct();

MAG.Mmin = Mmin;
MAG.Mmax = Mmax;
MAG.b_M = b_M;
MAG.dm = dm;

MAG.M_vec = M_vec;
MAG.FM = FM;
MAG.M_mid_vec = M_mid_vec;
MAG.P_M = P_M;
MAG.nM = nM;

end