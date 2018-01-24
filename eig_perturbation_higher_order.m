        function [H_perturb, gammasqr_perturb] = ...
                eig_perturbation_higher_order(Lfull_eig, H_eig, gammasqr_eig, Lfull_perturb, N_total_3)



            K0 = Lfull_eig;
            K1 = Lfull_perturb;
            delta_K = K1 - K0;
            phi0 = H_eig;
            lambda0 = gammasqr_eig;

            %step 1
            delta_lambda = diag(diag(transpose(phi0)*delta_K*phi0));
            lambda1 = lambda0 + delta_lambda;

            %step 2
            E = eye(2*N_total_3,2*N_total_3);
            Fnl = phi0*delta_lambda - delta_K*phi0;
            delta_phi = zeros(2*N_total_3,2*N_total_3);
            for i=1:2*N_total_3
                D1 = K1 - lambda1(i,i)*E;
                delta_phi(:,i)=D1\Fnl(:,i);
            end

            %step 3
            delta_phi = delta_phi + phi0;
            phi1 = phi0 - delta_phi;

            
            number_iterations = 0
            
            fin = 1;
            while (fin>0.1)
                [delta_lambda_new,delta_phi_new] = eig_iteration(delta_K,K1,phi0,delta_phi,lambda0,delta_lambda,N_total_3);
                fin = max(diag(abs(delta_lambda_new-delta_lambda))./diag(delta_lambda))
                                
                delta_lambda = delta_lambda_new;
                delta_phi = delta_phi_new;
                
                number_iterations = number_iterations+1
                               
            end
            
            fin_number = number_iterations;
            H_perturb = phi0 + delta_phi;
            gammasqr_perturb = lambda0 + delta_lambda;
            



