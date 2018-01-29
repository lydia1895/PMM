function [H_perturb, gammasqr_perturb] = ...
        eig_perturbation_iterative_High(Lfull_eig, H_eig, gammasqr_eig, Lfull_perturb, N_total_3)
    
    
    
    K0 = Lfull_eig;
    K1 = Lfull_perturb;
    delta_K = K1 - K0;
    phi0 = H_eig;
    lambda0 = gammasqr_eig;
    
    %step 1 - first-order eigenvalue perturbation
    D = transpose(phi0)*delta_K*phi0;
    delta_lambda = diag(diag(D));
    lambda1 = lambda0 + delta_lambda;
    
    %step 2 - first-order eigenvector perturbation
    lambda_row = transpose(diag(lambda0));
    lambda_matrix = repmat(lambda_row,2*N_total_3,1);
    diff_lambda = lambda_matrix - transpose(lambda_matrix);
    %%%small trick
    diff_lambda = diff_lambda + eye(2*N_total_3,2*N_total_3);
    D_divided = D./diff_lambda;
    X = D_divided - diag(diag(D_divided));
    delta_phi = phi0*X;
    
    %delete spurious solution as Eldred did?
    delta_phi = delta_phi + phi0;
    
    number_iterations = 0
    
    delta_lambda_previous = delta_lambda;
    delta_phi_previous = delta_phi;
    
    Y = delta_phi_previous;
    phi1 = phi0+Y;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    inv_K0 = eye(2*N_total_3,2*N_total_3)/K0;
    
    %%%test
    %{
    [row_K0,col_K0] = find(isnan(inv_K0));
    max_nan_inv__K0 = max(row_K0)
    %}
    while number_iterations<5
        number_iterations = number_iterations+1;
        
        %%%test
        %{
        [row_phi0,col_phi0] = find(isnan(delta_phi_previous));
        max_nan_delta_phi = max(row_phi0)
        [row_delta_lambda,col_delta_lambda] = find(isnan(delta_lambda_previous));
        max_nan_delta_lambda = max(row_delta_lambda)
        max_delta_lambda = max(delta_lambda_previous(:))
        %}
        
        %step 3 - new eigenvalue - check out again
        phi1 = phi0+delta_phi_previous;
        %delta_lambda_new = delta_lambda_previous;
        delta_lambda_new = diag( diag(transpose(phi0)*delta_K*phi1)./...
            (diag(transpose(phi0)*phi1)) );
        
        %step 4 - obtaining Y
        F = phi0*delta_lambda_new - delta_K*phi0;
        Y = inv_K0*(F + delta_phi_previous*lambda0);
        
        %step 5 - obtaining Z: Gram-Schmidt Ortogonalization
        
        phi1 = phi0+Y;
        R = -transpose(phi1)*phi1;
        R = triu(R,1) + eye(2*N_total_3,2*N_total_3);
        
        %delta_R = -transpose(delta_phi_previous)*phi1 -...
         %   transpose(phi1)*delta_phi_previous;
        delta_R = -transpose(Y)*phi1 -...
            transpose(phi1)*Y;
        delta_R = triu(delta_R,1);
        Z = Y*R + phi1*delta_R;
        
        %Z=Y;
        
        %Z = Y +phi0;
        
        %step 6 - mass normalization, Eldred 1992
        b = 2*transpose(Z)*phi0;
        a = transpose(Z)*Z;
        alpha = - diag(diag(b)./diag(2*a));
        
        %alpha = eye(2*N_total_3,2*N_total_3);
        
        delta_lambda_previous = delta_lambda_new;
        delta_phi_previous = Z*alpha;
        
        epsilon = zeros(2*N_total_3,1);
        for j=1:2*N_total_3
        epsilon(j) = (norm(K0*delta_phi_previous(:,j) - F(:,j) - delta_phi_previous(:,j)*lambda0(j,j)))^2;
        end
        epsilon_max(number_iterations) = max(epsilon);
    end
    eps_max = epsilon_max
    %delta_phi = delta_phi + phi0;
    %phi1 = phi0 + delta_phi;
 
    
   
    fin_number = number_iterations;
    H_perturb = phi0 + delta_phi_previous;
    gammasqr_perturb = lambda0 + delta_lambda_previous;
    
    
    
    
