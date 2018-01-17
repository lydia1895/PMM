
function [Rsum,Tsum, M, gammaminus,dq] =...
        PMM_main_function(figure_shape, dispersion, lambda_full, theta_full, phi_full, delta,...
        h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
        b_x, b_y, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y, ellipse_parameters,...
        n_points, eta, f1, verbose)
    
    %%%%%%%%%here starts the program
    
    nx = N_basis_x-1;   %because functions are p(0)...p(n(k)) on interval k
    Nx = zeros(N_intervals_x,1);
    for k=1:N_intervals_x
        Nx(k) = -nx(k);
        for p=1:k
            Nx(k) = Nx(k)+nx(p);
        end
    end
    
    ny = N_basis_y-1;
    Ny = zeros(N_intervals_y,1);
    for k=1:N_intervals_y
        Ny(k) = -ny(k);
        for p=1:k
            Ny(k) = Ny(k)+ny(p);
        end
    end
    
    N_total_x = sum(N_basis_x);  %total number of basis functions
    N_total_y = sum(N_basis_y);  %total number of basis functions
    N_total_x3 = N_total_x - N_intervals_x;  %number of basis functions in "third" basis
    N_total_y3 = N_total_y - N_intervals_y;
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
    [Dx, hx, P_dPx] = PMM_new_derivatives(La, N_intervals_x, N_basis_x, nx, Nx, ax, b_x1);
    [Dy, hy, P_dPy] = PMM_new_derivatives(La, N_intervals_y, N_basis_y, ny, Ny, ay, b_x2);
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
            
            
            [int_Ez_sqrt_g_full,int_Dz_unity_full,int_Dx_sqrt_g_full,int_Dy_sqrt_g_full,...
                int_Ex_g_down22_full,int_Ey_g_down12_full,int_Ex_g_down21_full,int_Ey_g_down11_full] =...
                PMM_metric_integral_polyfit_matrices(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
                N_intervals_x,N_intervals_y,n_points,La,ax,ay,hx,hy,dx_x1,dx_x2,dy_x1,dy_x2,uni,b_x1,b_x2);
            
    end
    
    if strcmp (figure_shape,'ellipse')==1 && strcmp (dispersion,'no')==1
        for nlayer=1:L
            
            [eps_total_t, mu_total_t] =...
                PMM_epsilon_ellipse_matrices(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
                N_intervals_x,N_intervals_y,La,epsilon(nlayer,:),...
                int_Ez_sqrt_g_full,int_Dz_unity_full,int_Dx_sqrt_g_full,int_Dy_sqrt_g_full,...
                int_Ex_g_down22_full,int_Ey_g_down12_full,int_Ex_g_down21_full,int_Ey_g_down11_full);
            
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
        
        %rectangle + ASR on x coordinate
        %{
        for i=1:L
            [eps_total_t, mu_total_t] =...
                PMM_metric_integral_polyfit_ASR(b_x,eta,f1,N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
                N_intervals_x,N_intervals_y,n_points,La,ax,ay,epsilon(:,:,i));
            eps_total(:,:,:,i) = eps_total_t;
            mu_total(:,:,:,i) = mu_total_t;
        end
        %}
        
    end
    if (verbose>5)
        title = 'incident integral for P0, Q0'
    end
    %precise solution that works
    [int_P1_Q1, int_P1_Q2] = incident_integral(La, b_x1, b_x2, alpha_ref, beta_ref,...
        N_basis_x, N_basis_y, N_intervals_x, N_intervals_y, Nx, Ny, nx, ny);
    %{
[int_P1_Q1, int_P1_Q2] = PMM_inc_coef(La, b_x1, b_x2, alpha_ref, beta_ref,...
    N_basis_x, N_basis_y, N_intervals_x, N_intervals_y, Nx, Ny, nx, ny);
        %}
        if (verbose>5)
            title = 'derive matrix of transition from PMM to FMM'
        end
        N = N_FMM;
        NN = (2*N_FMM+1)*(2*N_FMM+1);
        
        %precise solution that works
        [fx_coef,fy_coef] = PMM_to_FMM_RT_La05_one_integral(N, NN, La, alpha_ref, beta_ref,...
            b_x1, b_x2, N_intervals_x, N_intervals_y, N_basis_x, N_basis_y, Nx, nx, Ny, ny,...
            ax, ay);
        %{
[fx_coef,fy_coef] = PMM_to_FMM_RT_new(N, NN, La, alpha_ref, beta_ref,...
   b_x1, b_x2, N_intervals_x, N_intervals_y, N_basis_x, N_basis_y, Nx, nx, Ny, ny,...
   ax, ay);
            %}
            
            
            
            
            [Nll,Nlambda] = size(lambda_full);
            [Ntt,Ntheta] = size(theta_full);
            [Npp,Nphi] = size(phi_full);
            
            Rsum = zeros(Nlambda, Ntheta);
            Tsum = zeros(Nlambda, Ntheta);
            MM = zeros(4*N_total_3,4*N_total_3,Nlambda);
            u2d0FMM=zeros(NN,4,Ntheta);
            gzero = zeros(Nlambda,Ntheta);
            gzero_norm = zeros(Nlambda,Ntheta);
            gamma0 = zeros(Ntheta,1);
            gamma_num = zeros(Ntheta,1);
            gammaminus = zeros(2*N_total_3,L);
            
            n1 = refIndices(1);
            for i=1:Nlambda
                for j=1:Ntheta
                    for k=1:Nphi
                        if strcmp (figure_shape,'ellipse')==1 &&...
                                strcmp (dispersion,'yes')==1
                            refIndices_lambda = zeros(1,2);
                            for nlayer=1:L
                                
                                [eps_total_t, mu_total_t] =...
                                    PMM_epsilon_ellipse_matrices(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
                                    N_intervals_x,N_intervals_y,La,epsilon(nlayer,:,i),...
                                    int_Ez_sqrt_g_full,int_Dz_unity_full,int_Dx_sqrt_g_full,int_Dy_sqrt_g_full,...
                                    int_Ex_g_down22_full,int_Ey_g_down12_full,int_Ex_g_down21_full,int_Ey_g_down11_full);
                                
                                eps_total(:,:,:,nlayer) = eps_total_t;
                                mu_total(:,:,:,nlayer) = mu_total_t;
                                refIndices_lambda(1) = refIndices(i,1);
                                refIndices_lambda(2) = refIndices(i,2);
                                n1 = refIndices_lambda(1);
                            end
                        end
                        lambda = lambda_full(i)
                        theta = theta_full(j)
                        phi = phi_full(k);
                        
                        k0 = 2*pi/lambda;
                        
                        
                        alpha0 = k0*n1*sin(theta)*cos(phi);
                        beta0  = k0*n1*sin(theta)*sin(phi);
                        gamma0 = k0*n1*cos(theta);
                        
                        
                        %incident wave
                        
                        TETM = [0; cos(delta); sin(delta)];
                        TETMmatrix = [sin(theta)*cos(phi) cos(theta)*cos(phi) -sin(phi); ...
                            sin(theta)*sin(phi) cos(theta)*sin(phi) cos(phi);...
                            cos(phi) -sin(theta) 0];
                        E = TETMmatrix*TETM;
                        Ex = E(1);
                        Ey = E(2);
                        
                        k1 = k0*n1;
                        kz1v = gamma0;
                        A1_nul = ( k1^2 - alpha0^2)/(k0*kz1v);
                        B1_nul = ( k1^2 - beta0^2)/(k0*kz1v);
                        C1_nul = alpha0*beta0/(k0*kz1v);
                        
                        norm = A1_nul*abs(Ey)^2 + B1_nul*abs(Ex)^2 +...
                            C1_nul*( Ex*conj(Ey)+Ey*conj(Ex) );
                        
                        Ex0=Ex/sqrt(norm);
                        Ey0=Ey/sqrt(norm);
                        
                        %title = 'enter eigenvalue solver and S-matrix'
                        
                        if strcmp (dispersion,'no')==1
                            rrefIndices = refIndices;
                        else
                            rrefIndices = refIndices_lambda;
                        end
                        
                        [eta_R, eta_T, Stotal, ud_PMM, kz1v, kz2v,...
                            pplus, pminus, derx, eps, M, gamma_total,gamma_d1,...
                            gamma_norm, EH, gamma_sorted, W, u2d0_FMM_t,gzero_t,gzero_norm_t,...
                            gamma0_t,gamma_num_t,gammaminus,dq] =...
                            PMM_multi(int_P1_Q1,int_P1_Q2, fx_coef, fy_coef, Ex0, Ey0, lambda, theta, phi,...
                            N_FMM, h, L, rrefIndices, alpha_ref, beta_ref,...
                            b_x1, b_x2, N_intervals_x, N_intervals_y, N_basis_x, N_basis_y,...
                            Dx, Dy, hx, hy, eps_total, mu_total,verbose);
                        
                        %title = 'escape eigenvalue solver and S-matrix'
                        Wt(:,:,:,j)=W(:,:,L);
                        Rsum(i,j) = sum(eta_R);
                        Tsum(i,j) = sum(eta_T);
                        gzero(i,j) = gzero_t;
                        gzero_norm(i,j) = gzero_norm_t;
                        gamma00(j)=gamma0_t;
                        gamma_num(j)=gamma_num_t;
                        R00 = eta_R(N_FMM+1+N_FMM*(2*N_FMM+1));
                        T00 = eta_T(N_FMM+1+N_FMM*(2*N_FMM+1));
                        u2d0FMM(:,:,j)=u2d0_FMM_t;
                        
                    end
                end
            end
            
            
            figure(2)
            theta = theta_full*180/pi;
            plot(theta,gamma00,'r',theta,gamma_num,'g',theta,gzero,'m','Linewidth', 2)
            %plot(theta,gamma00,'r',theta,gzero,'m','Linewidth', 2)
            ylabel('abs(min(kz0-gamma(i)))')
            xlabel('theta')
            hold off
            
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
            %plot(theta_full*180/pi,gzero,'b','Linewidth', 2)
            %xlabel('theta')
            %ylabel('abs(min(kz0-gamma(i))/kz0)')
            hold off
            if (verbose>5)
                title = 'escape eigenvalue solver and S-matrix'
            end
