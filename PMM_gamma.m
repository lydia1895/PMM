function [gamma_norm, EH, gamma_sorted, W, pplus, pminus, eps11,...
        M, M31, gammaminus, gamma_total, gamma_d1, gamma_u1] = ...
        PMM_gamma(alpha_ref, beta_ref, k0, alpha0, beta0, h,...
        N_intervals_x, N_intervals_y, N_basis_x, N_basis_y,...
        Dx, Dy, hx, hy, eps_total, mu_total,verbose)
    if (verbose>5)
        title = 'enter PMM-gamma'
    end
    N_total_x = sum(N_basis_x);
    N_total_y = sum(N_basis_y);
    
    N_total_x3 = N_total_x - N_intervals_x; %in "third" basis
    N_total_y3 = N_total_y - N_intervals_y;
    
    N_total_3 = N_total_x3*N_total_y3;
    
    hxy = Kronecker_product(hx,hy);
    dx = hxy\Kronecker_product(1j*(alpha0 - alpha_ref)*hx + Dx, hy);
    dy = hxy\Kronecker_product(hx, 1j*(beta0  - beta_ref)*hy  + Dy);
    
    mzero = zeros(N_total_3,N_total_3);
    
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
  
    gamma_v = diag(gamma_norm);
    gamma = gamma_norm*k0;
    
    NN = N_total_3;
    
    gammaplus = zeros(2*NN,1);
    gammaminus = zeros(2*NN,1);
    EHplus = zeros(4*NN,2*NN);
    EHminus = zeros(4*NN,2*NN);
    nplus = 0;
    nminus = 0;
    
    for i=1:4*NN
        
        if ( real(gamma(i,i))+imag(gamma(i,i)) )>0
            nplus = nplus + 1;
            gammaplus(nplus) = gamma(i,i);
            EHplus(:,nplus) = EH(:,i);
        else
            nminus = nminus + 1;
            gammaminus(nminus) = gamma(i,i);
            EHminus(:,nminus) = EH(:,i);
        end
    end
    
    W = cat(2, EHplus, EHminus);
    
    pplusv = zeros(2*NN,1);
    pminusv = zeros(2*NN,1);
    for m=1:(2*NN)
        pplusv(m) = exp(1i*gammaplus(m)*h);
        pminusv(m) = exp(-1i*gammaminus(m)*h);
    end
    
    pplus = diag(pplusv);
    pminus = diag(pminusv);
    gamma_sorted = cat(1,gammaplus,gammaminus);
    gamma_sorted = diag(gamma_sorted);
    gamma_norm = gamma_sorted/k0;
    
    
    gamma_total = zeros(N_total_3,4);
    gamma_d1 = zeros(N_total_3,1);
    gamma_u1 = zeros(N_total_3,1);
    if (verbose>5)
        title = 'escape PMM-gamma'
    end
