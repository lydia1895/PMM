
function [H_1_4,gamma_sqr_1_4]= ...
        PMM_eig_for_Maxwell(L_full)
   
    
    %L_full = L_HE*L_EH;
    [H_1_4, gamma_sqr_1_4] = eig(L_full);
    