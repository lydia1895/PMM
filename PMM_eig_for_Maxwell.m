
function [W, pplus, pminus, gammaminus]= ...
        PMM_eig_for_Maxwell(L_HE,L_EH,h,N_total_3,verbose)
    if (verbose>5)
        title = 'enter PMM-eig'
    end
    
    L_full = L_HE*L_EH;
    [H_1_4, gamma_sqr_1_4] = eig(L_full);
    
    gamma_1_4 = diag(gamma_sqr_1_4.^0.5);
    n_plus_1_4 = 0;
    n_minus_1_4 = 0;
    for i=1:2*N_total_3
        if real(gamma_1_4(i))+imag(gamma_1_4(i))>0
            n_plus_1_4 = n_plus_1_4+1;
        else
            n_minus_1_4 = n_minus_1_4 + 1;
            gamma_1_4(i) = -gamma_1_4(i);
        end
    end
    
    E_1_4 = L_EH*H_1_4/diag(gamma_1_4);
    %H_1_4 = diag(gamma_1_4)\L_HE*E_1_4;
    W_up = cat(2,E_1_4,E_1_4);
    W_down = cat(2,H_1_4,-H_1_4);
    W = cat(1,W_up,W_down);
    EH = W;
    
    gammaplus = gamma_1_4;
    gammaminus = -gamma_1_4;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %test
    %{
    H_1_4_new = L_HE*E_1_4/diag(gamma_1_4);
    diff_H_1_4 = abs(H_1_4-H_1_4_new);
    max_diff_H_1_4 = max(diff_H_1_4(:))
    maxH = max(H_1_4(:))
    maxH_new = max(H_1_4_new(:))
    
    M_EH_minus_gamma_EH = abs(M*W-W*gammafull);
    max_M_EH_minus_gamma_EH = max(M_EH_minus_gamma_EH(:))
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NN = N_total_3;
    pplusv = zeros(2*NN,1);
    pminusv = zeros(2*NN,1);
    for m=1:(2*NN)
        pplusv(m) = exp(1i*gammaplus(m)*h);
        pminusv(m) = exp(-1i*gammaminus(m)*h);
    end
    
    pplus = diag(pplusv);
    pminus = diag(pminusv);

    if (verbose>5)
        title = 'escape PMM-gamma'
    end