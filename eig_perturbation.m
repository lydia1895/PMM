function [H_perturb, gammasqr_perturb] = ...
        eig_perturbation(Lfull_eig, H_eig, gammasqr_eig, Lfull_perturb, N_total_3)
    
    delta_L = Lfull_perturb - Lfull_eig;
    gammasqr_perturb = gammasqr_eig + diag(diag(transpose(H_eig)*delta_L*H_eig));
    
    
    H_L_H = transpose(H_eig)*delta_L*H_eig;
    
    %H_L_H_divided = (H_L_H - diag(H_L_H))/
    
    H_perturb = H_eig;
    
    for i=1:2*N_total_3
        for j = 1:2*N_total_3
            if i~=j
                H_perturb(:,i) = H_perturb(:,i) + H_eig(:,j)*...
                    H_L_H(j,i)/(gammasqr_eig(i,i)-gammasqr_eig(j,j));
            end
        end
    end
    
                
                          