
function [Rsum,Tsum] =...
        PMM_main_function(figure_shape, dispersion, lambda_full, theta_full, phi_full, delta,...
        h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
        b_x, b_y, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y, ellipse_parameters,...
        n_points, eta, f1, verbose,...
        Nlambda_eig, Nlambda_perturb, half_n_lambda, n_lambda_extra_perturb,...
        Ntheta_eig,  Ntheta_perturb,  half_n_theta,  n_theta_extra_perturb,...
        Nphi_eig,    Nphi_perturb,    half_n_phi,    n_phi_extra_perturb)
    
    %%%%%%%%%here starts the program
    
    [nx,Nx,N_total_x,N_total_x3] = PMM_number_of_basis_functions(N_intervals_x,N_basis_x);
    [ny,Ny,N_total_y,N_total_y3] = PMM_number_of_basis_functions(N_intervals_y,N_basis_y);
    N_total_3 = N_total_x3*N_total_y3;
    
    %for now
    
    b_x1 = b_x;
    b_x2 = b_y;
    
    eps_total = zeros(N_total_3,N_total_3,5,L);
    mu_total = zeros(N_total_3,N_total_3,5,L);
    
    %first compute coefficients [a] from boundary conditions
    %continuity and periodicity conditions are handled by sparse [a]
    
    if (verbose>5)
        title = 'enter boundary conditions'
    end
    ax = PMM_boundary_conditions(La, tau_x, N_intervals_x, N_basis_x, Nx, nx);
    ay = PMM_boundary_conditions(La, tau_y, N_intervals_y, N_basis_y, Ny, ny);
    if (verbose>5)
        title = 'escape boundary conditions'
    end
    %compute derivative matrices
    if (verbose>5)
        title = 'enter derivatives'
    end
    %for any Gegenbauer polynomial
    %{
[Dx, hx] = PMM_new_derivatives(La, N_intervals_x, N_basis_x, nx, Nx, ax, b_x1);
[Dy, hy] = PMM_new_derivatives(La, N_intervals_y, N_basis_y, ny, Ny, ay, b_x2);
    %}
    %for Legendre polynomial
    [Dx, hx] = PMM_Legendre_derivatives(La, N_intervals_x, N_basis_x, nx, Nx, ax, b_x1);
    [Dy, hy] = PMM_Legendre_derivatives(La, N_intervals_y, N_basis_y, ny, Ny, ay, b_x2);
    
    
    if (verbose>5)
        title = 'escape derivatives'
    end
    %%%%%%%%%%%%%%%for ellipse%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp (figure_shape,'ellipse')==1
        if (verbose>5)
            title = 'coordinates and derivatives for ellipse'
        end
        uni=1;
        
        %%%%%%%%%%last right
        [dx_x1,dx_x2,dy_x1,dy_x2] =...
            ellipse_coordinates_and_derivatives(ellipse_parameters,n_points,uni);
        
        if (verbose>5)
            title = 'integrals with metric tensor for eps and mu for ellipse'
        end
        
        
        %2D surface fit - working
        %{
    [int_g, int_for_ellipse] =...
        PMM_metric_integral_polyfit_matrices(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
        N_intervals_x,N_intervals_y,n_points,La,ax,ay,hx,hy,dx_x1,dx_x2,dy_x1,dy_x2,uni,b_x1,b_x2);
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%test
        %1D polyfit
        [int_g_1D, int_for_ellipse_1D] =...
            ellipse_coordinates_and_derivatives_for_1D_polyfit(ellipse_parameters,n_points,...
            La,Nx,nx,N_basis_x,ax,Ny,ny,N_basis_y,ay,N_total_3);
        
        int_g = int_g_1D;
        int_for_ellipse = int_for_ellipse_1D;
        
       
    end
    
    if strcmp (figure_shape,'ellipse')==1 && strcmp (dispersion,'no')==1
        for nlayer=1:L
            
            [eps_total_t, mu_total_t] =...
                PMM_epsilon_ellipse_matrices(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
                N_intervals_x,N_intervals_y,La,epsilon(nlayer,:),...
                int_g, int_for_ellipse);
            
            eps_total(:,:,:,nlayer) = eps_total_t;
            mu_total(:,:,:,nlayer) = mu_total_t;
            
        end
    end
    
    %%%%%%%%%%%%%%%for rectangle%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(figure_shape, 'rectangle')==1
        if (verbose>5)
            title = 'epsilon for rectangle'
        end
        %usual thing that works
        for i=1:L
            [eps_total_t, mu_total_t] =...
                PMM_epsilon_rectangle(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
                N_intervals_x,N_intervals_y,epsilon(:,:,i));
            eps_total(:,:,:,i) = eps_total_t;
            mu_total(:,:,:,i) = mu_total_t;
        end
        
    end
    
    if (verbose>5)
        title = 'incident integral for P0, Q0'
    end
    %precise solution that works
    [int_P1_Q1, int_P1_Q2] = incident_integral(La, b_x1, b_x2, alpha_ref, beta_ref,...
        N_basis_x, N_basis_y, N_intervals_x, N_intervals_y, Nx, Ny, nx, ny);
    
    if (verbose>5)
        title = 'derive matrix of transition from PMM to FMM'
    end
    N = N_FMM;
    NN = (2*N_FMM+1)*(2*N_FMM+1);
    
    %precise solution that works
    [fx_coef,fy_coef] = PMM_to_FMM_RT_La05_one_integral(N, NN, La, alpha_ref, beta_ref,...
        b_x1, b_x2, N_intervals_x, N_intervals_y, N_basis_x, N_basis_y, Nx, nx, Ny, ny,...
        ax, ay);
    
    
    
    %[Nll,Nlambda_perturb] = size(lambda_full);
    %[Ntt,Ntheta_perturb] = size(theta_full);
    [Npp,Nphi] = size(phi_full);
    
    Rsum = zeros(Nlambda_perturb, Ntheta_perturb);
    Tsum = zeros(Nlambda_perturb, Ntheta_perturb);
    gzero = zeros(Nlambda_perturb,Ntheta_perturb);
    gzero_norm = zeros(Nlambda_perturb,Ntheta_perturb);
    gamma_num = zeros(Nlambda_perturb,Ntheta_perturb,1);
    gamma00 = zeros(Ntheta_perturb,1);
    gammaminus = zeros(2*N_total_3,L);
    Lfull_eig = zeros(2*N_total_3,2*N_total_3,L);
    L_EH_eig = zeros(2*N_total_3,2*N_total_3,L);
    L_HE_eig = zeros(2*N_total_3,2*N_total_3,L);
    H_eig = zeros(2*N_total_3,2*N_total_3,L);
    gammasqr_eig = zeros(2*N_total_3,2*N_total_3,L);
    
    H_perturb = zeros(2*N_total_3,2*N_total_3,L);
    gammasqr_perturb = zeros(2*N_total_3,2*N_total_3,L);
    L_EH_perturb = zeros(2*N_total_3,2*N_total_3,L);
    %{
    M_eig_Au_layer=zeros(Nlambda_eig,Ntheta_eig,2*N_total_3,2*N_total_3);
    gammasqr_Au_layer=zeros(Nlambda_eig,Ntheta_eig,2*N_total_3);
    H_perturb_1_sorted_full=zeros(Nlambda_perturb,Ntheta_perturb,2*N_total_3,2*N_total_3);
    gammasqr_1_sorted_full=zeros(Nlambda_perturb,Ntheta_perturb,2*N_total_3);
    %}
    current_layer = 2;
    
    n1 = refIndices(1);
    for i_lambda_eig=1:Nlambda_eig
        for j_theta_eig=1:Ntheta_eig
            for k_phi_eig=1:Nphi_eig
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for points where we
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%calculate eig rigorously
                i_perturb = 0;
                i_lambda_perturb = half_n_lambda + 1 + i_perturb +...
                    (i_lambda_eig-1)*n_lambda_extra_perturb;
                j_perturb = 0;
                j_theta_perturb = half_n_theta + 1 + j_perturb +...
                    (j_theta_eig-1)*n_theta_extra_perturb;
                k_phi_perturb= k_phi_eig;
                %%%%%%%%%%%%%%%%%%%%%refIndices we calculate for different lambda
                if strcmp (figure_shape,'ellipse')==1 &&...
                        strcmp (dispersion,'yes')==1
                    refIndices_lambda = zeros(1,2);
                    refIndices_lambda(1) = refIndices(i_lambda_perturb,1);
                    refIndices_lambda(2) = refIndices(i_lambda_perturb,2);
                    n1 = refIndices_lambda(1);
                end
                if strcmp (dispersion,'no')==1
                    rrefIndices = refIndices;
                else
                    rrefIndices = refIndices_lambda;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%wave vector and E0 we calculate for different
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%lambda and theta and phi
                lambda = lambda_full(i_lambda_perturb)
                theta = theta_full(j_theta_perturb)
                phi = phi_full(k_phi_perturb);
                [alpha0,beta0,gamma0,k0,Ex0,Ey0] =...
                    PMM_incident_wave_vector_and_field(lambda,theta,phi,delta,n1);
                
                
                %%%%%%%%%%%%%%for each layer we calculate epsilon and matrix and then eig
                for nlayer=1:L
                    if strcmp (figure_shape,'ellipse')==1 &&...
                            strcmp (dispersion,'yes')==1
                        [eps_total, mu_total] =...
                            PMM_epsilon_ellipse_matrices(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
                            N_intervals_x,N_intervals_y,La,epsilon(nlayer,:,i_lambda_perturb),...
                            int_g, int_for_ellipse);
                        
                        
                        [L_HE_eig(:,:,nlayer), L_EH_eig(:,:,nlayer)]= ...
                            PMM_Maxwell_matrix(alpha_ref, beta_ref, k0, alpha0, beta0,...
                            N_intervals_x, N_intervals_y, N_basis_x, N_basis_y,...
                            Dx, Dy, hx, hy, eps_total, mu_total);
                        
                    else
                        
                        
                        [L_HE_eig(:,:,nlayer), L_EH_eig(:,:,nlayer)]= ...
                            PMM_Maxwell_matrix(alpha_ref, beta_ref, k0, alpha0, beta0,...
                            N_intervals_x, N_intervals_y, N_basis_x, N_basis_y,...
                            Dx, Dy, hx, hy, eps_total(:,:,:,nlayer), mu_total(:,:,:,nlayer));
                        
                    end
                    
                    Lfull_eig(:,:,nlayer) = L_HE_eig(:,:,nlayer)*L_EH_eig(:,:,nlayer);
                    
                    [H_eig(:,:,nlayer),gammasqr_eig(:,:,nlayer)] = PMM_eig_for_Maxwell(Lfull_eig(:,:,nlayer));
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%test!
                    %{
                    gammasqr_1 = diag(gammasqr_eig(:,:,current_layer));
                    H_perturb_1 = H_eig(:,:,current_layer);
                    [gammasqr_1_sorted, index] = sort(gammasqr_1);
                    H_perturb_1_sorted = zeros(2*N_total_3,2*N_total_3);
                    for i=1:2*N_total_3
                        H_perturb_1_sorted(:,i) = H_perturb(:,index(i),1);
                    end
                    H_perturb_1_sorted_full(i_lambda_perturb,j_theta_perturb,:,:) = H_perturb_1_sorted;
                    gammasqr_1_sorted_full(i_lambda_perturb,j_theta_perturb,:) = gammasqr_1;
                    %}
                end
                %M_eig_Au_layer(i_lambda_eig,j_theta_eig,:,:) = Lfull_eig(:,:,2);
                
                
                
                %%%%%%%%%%%%%here we start perturbation calculation
                i_current = 1209;
                for i_perturb = -half_n_lambda:half_n_lambda
                    for j_perturb = -half_n_theta:half_n_theta
                        title = 'perturbations'
                        jj=j_perturb
                        i_lambda_perturb = half_n_lambda + 1 + i_perturb +...
                            (i_lambda_eig-1)*n_lambda_extra_perturb
                        lambda = lambda_full(i_lambda_perturb)
                        
                        j_theta_perturb = half_n_theta + 1 + j_perturb +...
                            (j_theta_eig-1)*n_theta_extra_perturb
                        %j_theta_perturb = j_theta_eig
                        theta = theta_full(j_theta_perturb)
                        
                        phi = phi_full(k_phi_perturb);
                        
                        %%%%%%%%%%%%%we're in the same layer as we were
                        %%%%%%%%%%%%%%%%%%%%%epsilon we calculate for different lambda
                        if strcmp (figure_shape,'ellipse')==1 &&...
                                strcmp (dispersion,'yes')==1
                            refIndices_lambda = zeros(1,2);
                            refIndices_lambda(1) = refIndices(i_lambda_perturb,1);
                            refIndices_lambda(2) = refIndices(i_lambda_perturb,2);
                            n1 = refIndices_lambda(1);
                        end
                        if strcmp (dispersion,'no')==1
                            rrefIndices = refIndices;
                        else
                            rrefIndices = refIndices_lambda
                        end
                        
                        [alpha0,beta0,gamma0,k0,Ex0,Ey0] =...
                            PMM_incident_wave_vector_and_field(lambda,theta,phi,delta,n1);
                        
                        if (i_perturb==0)&&(j_perturb==0)
                            for nlayer=1:L
                                
                                L_EH_perturb(:,:,nlayer) = L_EH_eig(:,:,nlayer);
                                H_perturb(:,:,nlayer) = H_eig(:,:,nlayer);
                                gammasqr_perturb(:,:,nlayer) = gammasqr_eig(:,:,nlayer);
                            end
                            gammasqr_Au_layer(i_lambda_eig,j_theta_eig,:) = diag(gammasqr_perturb(:,:,2));
                        else
                            title = 'else'
                            for nlayer=1:L
                                if strcmp (figure_shape,'ellipse')==1 &&...
                                        strcmp (dispersion,'yes')==1 %&& (i_current~=i_perturb)
                                    [eps_total, mu_total] =...
                                        PMM_epsilon_ellipse_matrices(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
                                        N_intervals_x,N_intervals_y,La,epsilon(nlayer,:,i_lambda_perturb),...
                                        int_g, int_for_ellipse);
                                    i_current = i_perturb;
                                    
                                    [L_HE_perturb, L_EH_perturb(:,:,nlayer)]= ...
                                        PMM_Maxwell_matrix(alpha_ref, beta_ref, k0, alpha0, beta0,...
                                        N_intervals_x, N_intervals_y, N_basis_x, N_basis_y,...
                                        Dx, Dy, hx, hy, eps_total, mu_total);
                                else
                                    [L_HE_perturb, L_EH_perturb(:,:,nlayer)]= ...
                                        PMM_Maxwell_matrix(alpha_ref, beta_ref, k0, alpha0, beta0,...
                                        N_intervals_x, N_intervals_y, N_basis_x, N_basis_y,...
                                        Dx, Dy, hx, hy, eps_total(:,:,:,nlayer), mu_total(:,:,:,nlayer));
                                end
                                nnlayer = nlayer
                                
                                
                                Lfull_perturb = L_HE_perturb*L_EH_perturb(:,:,nlayer);
                                [H_perturb(:,:,nlayer), gammasqr_perturb(:,:,nlayer)] = ...
                                    eig_perturbation(Lfull_eig(:,:,nlayer), H_eig(:,:,nlayer),...
                                    gammasqr_eig(:,:,nlayer), Lfull_perturb, N_total_3);
                                
                                
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%test
                        %{
                        gammasqr_1 = diag(gammasqr_perturb(:,:,current_layer));
                        H_perturb_1 = H_perturb(:,:,current_layer);
                        [gammasqr_1_sorted, index] = sort(gammasqr_1);
                        H_perturb_1_sorted = zeros(2*N_total_3,2*N_total_3);
                        for i=1:2*N_total_3
                            H_perturb_1_sorted(:,i) = H_perturb_1(:,index(i),1);
                        end
                        H_perturb_1_sorted_full(i_lambda_perturb,j_theta_perturb,:,:) = H_perturb_1_sorted;
                        gammasqr_1_sorted_full(i_lambda_perturb,j_theta_perturb,:) = gammasqr_1;
                        %}
                        
                        [eta_R, eta_T, M,...
                            gzero_t,gzero_norm_t,gamma_num_t,gammaminus] =...
                            PMM_multi(int_P1_Q1,int_P1_Q2, fx_coef, fy_coef,...
                            Ex0, Ey0, alpha0,beta0,gamma0,k0, N_FMM, h, L, rrefIndices,...
                            b_x1, b_x2, N_intervals_x, N_intervals_y, N_basis_x, N_basis_y,...
                            verbose, H_perturb, gammasqr_perturb, L_EH_perturb);
                        
                        Rsum(i_lambda_perturb,j_theta_perturb) = sum(eta_R);
                        Tsum(i_lambda_perturb,j_theta_perturb) = sum(eta_T);
                        gzero(i_lambda_perturb,j_theta_perturb) = gzero_t;
                        gzero_norm(i_lambda_perturb,j_theta_perturb) = gzero_norm_t;
                        gamma00(j_theta_perturb)=gamma0;
                        gamma_num(i_lambda_perturb,j_theta_perturb)=gamma_num_t;
                    end
                end
            end
            
        end
    end
    
    
    
    %we can pack
    %alpha_ref,b_x1,N_intervals_x,N_basis_x,Dx,hx into one
    %x_object and analogous varibles into y_object
    %{
figure(2)
theta = theta_full*180/pi;

plot(theta,gamma00,'r',theta,gamma_num,'g',theta,gzero,'m','Linewidth', 2)
%plot(theta,gamma00,'r',theta,gzero,'m','Linewidth', 2)
ylabel('abs(min(kz0-gamma(i)))')
xlabel('theta')
hold off
    %}
    %{
                figure(3)
                pcolor(lambda_full,theta_full*180/pi,transpose(gzero))
                %pcolor(XI,YI*180/pi,ZI)
                xlabel('lambda for gzero');
                ylabel('theta');
                shading flat
                colorbar

                figure(2)
                pcolor(lambda_full,theta_full*180/pi,transpose(gzero_norm))
                %pcolor(XI,YI*180/pi,ZI)
                xlabel('lambda for gzero_norm');
                ylabel('theta');
                shading flat
                colorbar
    %}
    hold off
    if (verbose>5)
        title = 'escape eigenvalue solver and S-matrix'
    end
