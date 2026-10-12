function PDs_ij = PDs_func(Mu_mat,Beta_mat,i,j,IM_vec,f_IM_A)
    nDs = length(Mu_mat);

    if j == nDs
        Frag = normcdf(log(IM_vec/Mu_mat(i,j))/Beta_mat(i,j));
        PDs_ij = sum(Frag .* f_IM_A');
    else
        Frag1 = normcdf(log(IM_vec/Mu_mat(i,j))/Beta_mat(i,j));
        Frag2 = normcdf(log(IM_vec/Mu_mat(i,j+1))/Beta_mat(i,j+1));
        d_frag = Frag1-Frag2;        

 
        PDs_ij = sum(d_frag .* f_IM_A'); 
    end

end