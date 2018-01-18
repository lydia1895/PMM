function [H_perturb, gammasqr_perturb] = ...
        PMM_perturbation(Lfull_eig, H_eig, gammasqr_eig, Lfull_perturb)
    
    delta_L = Lfull_perturb - Lfull_eig;
    gammasqr_perturb = gammasqr_eig + diag(diag(transpose(H_eig)*delta_L*H_eig));
    
    
    
    