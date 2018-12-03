
function [H_1_4_new,gamma_sqr_1_4_new]= ...
        PMM_eig_for_Maxwell(L_full)
   
    
    %L_full = L_HE*L_EH;
    [H_1_4, gamma_sqr_1_4] = eig(L_full);
    
    [gamma_sqr_1_4_diag,I] = sort(diag(gamma_sqr_1_4));
    [N,NN] = size(gamma_sqr_1_4_diag);

    H_1_4_new = zeros(N,N);
    for j=1:N
        H_1_4_new(:,j) = H_1_4(:,I(j));
    end
    
    gamma_sqr_1_4_new = diag(gamma_sqr_1_4_diag);
    