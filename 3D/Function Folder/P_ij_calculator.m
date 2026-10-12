function PA_ij = P_ij_calculator(Mu_mat, Beta_mat, n_R, n_M, IM_vec, f_IM_A)
    Nd = length(Mu_mat);
    PA_ij=repmat(zeros(Nd+1,Nd+1),1,1,n_R,n_M); 
    for m = 1 : n_M
        for n = 1 : n_R
            for i = 1 : Nd
                for j = i:Nd
                    PA_ij(i,j+1,n,m) = PDs_func(Mu_mat,Beta_mat,i,j,IM_vec,f_IM_A(:,n,m));
                end
                PA_ij(i,i,n,m) = 1 - sum(PA_ij(i,i+1:end,n,m));
            end  
            PA_ij(end,end,n,m) = 1;
        end
    end

end
