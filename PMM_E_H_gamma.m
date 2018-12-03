
function [W, pplus, pminus, gammaminus]= ...
        PMM_E_H_gamma(H_1_4, gamma_sqr_1_4,L_EH,N_total_3,h)
    
    
    
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
    %%%%%%%%%%%%%%%this is new
    [gamma_1_4,I] = sort(gamma_1_4);
    H_1_4_new = zeros(2*N_total_3,2*N_total_3);
    for j=1:2*N_total_3;
        H_1_4_new(:,j) = H_1_4(:,I(j));
    end
    H_1_4 = H_1_4_new;
    %%%%%%%%%%%%%%%
    E_1_4 = L_EH*H_1_4/diag(gamma_1_4);
    %H_1_4 = diag(gamma_1_4)\L_HE*E_1_4;
    W_up = cat(2,E_1_4,E_1_4);
    W_down = cat(2,H_1_4,-H_1_4);
    W = cat(1,W_up,W_down);
    
    %gamma_1_4 = sort(gamma_1_4); %careful here!!!!!!!!!!!!!!!!!!!!!!
    gammaplus = gamma_1_4;
    gammaminus = -gamma_1_4;
    
    NN = N_total_3;
    pplusv = zeros(2*NN,1);
    pminusv = zeros(2*NN,1);
    for m=1:(2*NN)
        pplusv(m) = exp(1i*gammaplus(m)*h);
        pminusv(m) = exp(-1i*gammaminus(m)*h);
    end
    
    pplus = diag(pplusv);
    pminus = diag(pminusv);
