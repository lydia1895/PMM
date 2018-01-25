function [delta_lambda_new,delta_phi_new] =...
        eig_iteration(delta_K,K1,phi0,delta_phi,delta_phi0,lambda0,delta_lambda,N_total_3)
    
    %step 4
    phi1 = phi0 + delta_phi;
    delta_lambda_new = diag(diag(transpose(phi0)*delta_K*phi1)./(diag(transpose(phi0)*phi1)));
    lambda1 = lambda0 + delta_lambda_new;
    
    %step 5
    
    E = eye(2*N_total_3,2*N_total_3);
    Fnl = phi0*delta_lambda - delta_K*phi0;
    
    %%%test
    %{
    [row_phi0,col_phi0] = find(isnan(phi0));
    max_nan_phi0 = max(row_phi0)
    [row_delta_lambda,col_delta_lambda] = find(isnan(delta_lambda));
    max_nan_delta_lambda = max(row_delta_lambda)
    [row_delta_K,col_delta_K] = find(isnan(delta_K));
    max_nan_delta_K = max(row_delta_K)
    %}
    %[row_F,col_F]=find(isnan(Fnl));
    %max_nan_F = max(row_F)
    %delta_phi = zeros(2*N_total_3,2*N_total_3);
    
    
    V = zeros(2*N_total_3,2*N_total_3);
    for i=1:2*N_total_3
        D1 = K1 - lambda1(i,i)*E;
        [max_phi,n_max_phi] = max(abs(phi0(:,i)));
        j1 = n_max_phi;
        D1(j1,:) = zeros(1,2*N_total_3);
        D1(:,j1) = zeros(2*N_total_3,1);
        D1(j1,j1) = 1;
        Fnl(j1,i) = 0;
        V(:,i)=D1\Fnl(:,i);
        V(j1,i) = delta_phi0(j1,i);
    end
    
    %test
    %[row_V,col_V]=find(isnan(V))
    
    %step 6
    ksi = 2*(diag(transpose(phi0)*V))+(diag(transpose(V)*V));
    c = 1 - (1+ksi).^0.5;
    
    %step 7
    Erow=ones(2*N_total_3,1);
    delta_phi_new = phi0*diag(c./(1-c)) + V*diag(Erow./(1-c));
    %[row_phi,col_phi] = find(isnan(delta_phi_new))
    
