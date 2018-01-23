function [H_perturb, gammasqr_perturb] = ...
        eig_perturbation(Lfull_eig, H_eig, gammasqr_eig, Lfull_perturb, N_total_3)
    
    
    delta_L = Lfull_perturb - Lfull_eig;
    
    H_L_H = transpose(H_eig)*delta_L*H_eig;
    gammasqr_perturb = gammasqr_eig + diag(diag(H_L_H));
    
    lambda_row = transpose(diag(gammasqr_eig));
    lambda_matrix = repmat(lambda_row,2*N_total_3,1);
    diff_lambda = lambda_matrix - transpose(lambda_matrix);
    
    %%%small trick
    diff_lambda = diff_lambda + eye(2*N_total_3,2*N_total_3);
    H_L_H_divided = H_L_H./diff_lambda;
    H_L_H_divided = H_L_H_divided - diag(diag(H_L_H_divided));
    %%%
    
    H_perturb = H_eig + H_eig*H_L_H_divided;
    
    
                
                          