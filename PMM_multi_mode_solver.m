function [Rtotal] = PMM_multi_mode_solver(int_P1_Q1,int_P1_Q2, fx_coef, fy_coef,...
        Ex0, Ey0, alpha0,beta0,gamma0,k0, N_FMM, h, L, refIndices,...
        b_x1, b_x2, N_intervals_x, N_intervals_y, N_basis_x, N_basis_y,...
        verbose, H_1_4, gamma_sqr_1_4, L_EH)
    
    
    N_total_x = sum(N_basis_x);  %total number of basis functions
    N_total_y = sum(N_basis_y);  %total number of basis functions
    N_total_x3 = N_total_x - N_intervals_x;  %number of basis functions in "third" basis
    N_total_y3 = N_total_y - N_intervals_y;
    N_total_3 = N_total_x3*N_total_y3;
  
    n1 = refIndices(1);
    n2 = refIndices(2);

    gammaminus = zeros(2*N_total_3, L);    
    pplus = zeros(2*N_total_3, 2*N_total_3, L);
    pminus = zeros(2*N_total_3, 2*N_total_3, L);
    M = zeros(4*N_total_3, 4*N_total_3,L);
    W = zeros(4*N_total_3, 4*N_total_3, L);
    
    for nlayer=1:L 
        [Wt, pplust, pminust, gammaminus_t]= ...
            PMM_E_H_gamma(H_1_4(:,:,nlayer), gamma_sqr_1_4(:,:,nlayer),L_EH(:,:,nlayer),...
            N_total_3,h(nlayer));

            gammaminus(:,nlayer) = gammaminus_t;
            W(:,:,nlayer) = Wt;
            pplus(:,:,nlayer) = pplust;
            pminus(:,:,nlayer) = pminust;
    end
    
    %gammaminus(:,2)
    
    %S-matrix propagation
    
    Smin1 = eye(4*N_total_3,4*N_total_3);
    Stemp = Smin1;
    for i=1:(L-1)
        Si = new_recursion(Stemp, W(:,:,i), W(:,:,i+1), pplus(:,:,i), pminus(:,:,i+1));
        Stemp = Si;
    end
    Stotal = Stemp;
    Rtotal = Stotal\eye(4*N_total_3,4*N_total_3);
    %{
    WR = W;
    WR(:,:,1) = [W(:,2*N_total_3+1:4*N_total_3,1) W(:,1:2*N_total_3,1)];
    WR(:,:,L) = [W(:,2*N_total_3+1:4*N_total_3,L) W(:,1:2*N_total_3,L)];
    
    %R-matrix (reciprocal to S)
    Rmin1 = eye(4*N_total_3,4*N_total_3);
    Rtemp = Rmin1;
    for i=1:(L-1)
        Ri = new_recursion(Rtemp, WR(:,:,i), WR(:,:,i+1), pplus(:,:,i), pminus(:,:,i+1));
        Rtemp = Ri;
    end
    Rtotal = Rtemp;
    %}