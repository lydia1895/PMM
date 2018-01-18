
function [L_HE, L_EH]= ...
        PMM_Maxwell_matrix(alpha_ref, beta_ref, k0, alpha0, beta0,...
        N_intervals_x, N_intervals_y, N_basis_x, N_basis_y,...
        Dx, Dy, hx, hy, eps_total, mu_total)
    
    N_total_x = sum(N_basis_x);
    N_total_y = sum(N_basis_y);
    
    N_total_x3 = N_total_x - N_intervals_x; %in "third" basis
    N_total_y3 = N_total_y - N_intervals_y;
    
    N_total_3 = N_total_x3*N_total_y3;
    
    hxy = Kronecker_product(hx,hy);
    dx = hxy\Kronecker_product(1j*(alpha0 - alpha_ref)*hx + Dx, hy);
    dy = hxy\Kronecker_product(hx, 1j*(beta0  - beta_ref)*hy  + Dy);
    
    
    eps11 = eps_total(:,:,1);
    eps12 = eps_total(:,:,2);
    eps21 = eps_total(:,:,3);
    eps22 = eps_total(:,:,4);
    eps_inv33 = eps_total(:,:,5);
    
    mu11 = mu_total(:,:,1);
    mu12 = mu_total(:,:,2);
    mu21 = mu_total(:,:,3);
    mu22 = mu_total(:,:,4);
    mu_inv33 = mu_total(:,:,5);
    alpha = 1j*dx;
    beta = 1j*dy;
    
    M13 = k0*mu21 + (1/k0)*alpha*eps_inv33*beta;
    M14 = k0*mu22 - (1/k0)*alpha*eps_inv33*alpha;
    
    M23 = -k0*mu11 + (1/k0)*beta*eps_inv33*beta;
    M24 = -k0*mu12 - (1/k0)*beta*eps_inv33*alpha;
    
    M31 = -k0*eps21 - (1/k0)*alpha*mu_inv33*beta;
    M32 = -k0*eps22 + (1/k0)*alpha*mu_inv33*alpha;
    
    M41 = k0*eps11 - (1/k0)*beta*mu_inv33*beta;
    M42 = k0*eps12 + (1/k0)*beta*mu_inv33*alpha;
    
    %{
    
    mzero = zeros(N_total_3,N_total_3);
    
    M11 = mzero;
    M12 = mzero;
    M13 = k0*mu21 + (1/k0)*alpha*eps_inv33*beta;
    M14 = k0*mu22 - (1/k0)*alpha*eps_inv33*alpha;
    
    M21 = mzero;
    M22 = mzero;
    M23 = -k0*mu11 + (1/k0)*beta*eps_inv33*beta;
    M24 = -k0*mu12 - (1/k0)*beta*eps_inv33*alpha;
    
    M31 = -k0*eps21 - (1/k0)*alpha*mu_inv33*beta;
    M32 = -k0*eps22 + (1/k0)*alpha*mu_inv33*alpha;
    M33 = mzero;
    M34 = mzero;
    
    M41 = k0*eps11 - (1/k0)*beta*mu_inv33*beta;
    M42 = k0*eps12 + (1/k0)*beta*mu_inv33*alpha;
    M43 = mzero;
    M44 = mzero;
    
    M1 = cat(2, M11, M12, M13, M14);
    M2 = cat(2, M21, M22, M23, M24);
    M3 = cat(2, M31, M32, M33, M34);
    M4 = cat(2, M41, M42, M43, M44);
    
    M = cat(1, M1, M2, M3, M4)/k0;
   
    
    [EH, gamma_norm] = eig(M);
    
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    L_EH_up = cat(2,M13,M14);
    L_EH_down = cat(2,M23,M24);
    L_EH = cat(1,L_EH_up,L_EH_down);
    
    L_HE_up = cat(2,M31,M32);
    L_HE_down = cat(2,M41,M42);
    L_HE = cat(1,L_HE_up,L_HE_down);
    
    