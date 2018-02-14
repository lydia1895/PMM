
clc
clear all
% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);

verbose = 8;
tic;
format long
eta=0;
f1=0;
%%%%%%%%%
load MatParam_nk_Au_2_interpExportData.txt
load MatParam_nk_ZnSe_interpExportData.txt
load MatParam_nk_SiO2_interpExportData.txt
load MatParam_nk_Zep520A_interpExportData.txt
%%%%%%%%%
lamnkAu = MatParam_nk_Au_2_interpExportData;
lamnkZnSe = MatParam_nk_ZnSe_interpExportData;
lamnkSiO2 = MatParam_nk_SiO2_interpExportData;
lamnkZep = MatParam_nk_Zep520A_interpExportData;

Nlambda_eig = 1;
n_lambda_extra_perturb = 1;
Nlambda_perturb = n_lambda_extra_perturb * Nlambda_eig;
half_n_lambda = floor((n_lambda_extra_perturb-1)/2);

Ntheta_eig = 1;
n_theta_extra_perturb = 5;
Ntheta_perturb = n_theta_extra_perturb * Ntheta_eig;
half_n_theta = floor((n_theta_extra_perturb-1)/2);

Nphi_eig = 1;
n_phi_extra_perturb = 1;
Nphi_perturb = n_phi_extra_perturb * Nphi_eig;
half_n_phi = floor((n_phi_extra_perturb-1)/2);

lmin = 1200;
lmax = 1200;

%dlambda_eig = (lmax-lmin)/(Nlambda_eig-1);
%dlambda_perturb = (lmax-lmin)/n_lambda_extra_perturb;

lambda = linspace(lmin,lmax,Nlambda_perturb);


tmin = 50*pi/180;
tmax = 54*pi/180;

theta = linspace(tmin,tmax,Ntheta_perturb);
phi = 45*pi/180;

%{
%test
lambda_current = zeros(1,Nlambda_perturb);
for i =1:Nlambda_eig
    for i_perturb = -half_n_lambda:half_n_lambda
        i_lambda_perturb = 1 + half_n_lambda + i_perturb +...
            (i-1)*n_lambda_extra_perturb;
        lambda_current(i_lambda_perturb) = lambda(i_lambda_perturb);
    end
end
%}

Nl = Nlambda_perturb;
nAu = zeros(Nl,1);
nZnSe = zeros(Nl,1);
nSiO2 = zeros(Nl,1);
nZep = zeros(Nl,1);
for i=1:Nl
    num = floor(lambda(i)) - 500 + 1;
    nAu(i) = lamnkAu(num,2) + 1j*lamnkAu(num,3);
    nZnSe(i) = lamnkZnSe(num,2) + 1j*lamnkZnSe(num,3);
    nSiO2(i) = lamnkSiO2(num,2) + 1j*lamnkSiO2(num,3);
    nZep(i) = lamnkZep(num,2) + 1j*lamnkZep(num,3);
end
epsAu = nAu.^2;
epsZnSe = nZnSe.^2;
epsSiO2 = nSiO2.^2;
epsZep = nZep.^2;
epsAir = ones(Nl,1);
%%%%%%%%%conditions

figure_shape = 'ellipse';
dispersion = 'yes';

N_intervals_x = 3;
N_intervals_y = 3;
N_b = 6;
n_points = 150;
N_basis_x = N_b*ones(N_intervals_x,1);
N_basis_y = N_b*ones(N_intervals_y,1);



R1 = 140/2;
R2 = 175/2;
P1 = 200;
P2 = 200;
Q2 = R2/sqrt(2);
Q1 = (R1/R2)*sqrt(R2^2-Q2^2);

ellipse_parameters = [R1 R2 P1 P2 Q1 Q2];
b_x1 = zeros(N_intervals_x+1);
b_x2 = zeros(N_intervals_y+1);

x1_t_plus =  P1/2+Q1;
x1_t_minus = P1/2-Q1;
x2_t_plus =  P2/2+Q2;
x2_t_minus = P2/2-Q2;

b_x1 = [0 x1_t_minus x1_t_plus P1];
b_x2 = [0 x2_t_minus x2_t_plus P2];


[Nx1, NNxx1] = size(b_x1);
[Nx2, NNxx2] = size(b_x2);
periodx = b_x1(NNxx1)-b_x1(1);
periody = b_x2(NNxx2)-b_x2(1);

%delta is the angle between E and the incidence plane
%delta = 0 TM, delta = pi/2 TE

delta = pi/2;

refIndices = [nZnSe nSiO2];

L=5; %number of layers
%epsilon(iL,1,iNlambda) = eps outside the ellipse
%epsilon(iL,2,iNlambda) = eps inside the ellipse

epsilon = zeros(L,2,Nl);

epsilon(5,1,:) = epsZnSe;
epsilon(5,2,:) = epsZnSe;

epsilon(4,1,:) = epsAir;
epsilon(4,2,:) = epsAir;

epsilon(3,1,:) = epsZep;
epsilon(3,2,:) = epsZep;

epsilon(2,1,:) = epsZep;
epsilon(2,2,:) = epsAu;

epsilon(1,1,:) = epsSiO2;
epsilon(1,2,:) = epsSiO2;

h = zeros(L,1);
h(5) = 0.0;
h(4) = 25;
h(3) = 180;
h(2) = 20;
h(1) = 0.0;

alpha_ref = -sin(pi/6)/periodx;
beta_ref =  -sin(pi/6)/periody;

tau_x = exp(1j*alpha_ref*periodx);
tau_y = exp(1j*beta_ref*periody);

%La is lambda in Gegenbauer polynomials; La>-1/2

La = 0.5;

N_FMM = 1;
%{
save('ellipse_plasmonic_perturb_1order_input.mat','figure_shape', 'dispersion',...
    'lambda', 'theta', 'phi', 'delta',...
    'h', 'L', 'N_FMM', 'epsilon', 'refIndices', 'La', 'tau_x', 'tau_y',...
    'alpha_ref', 'beta_ref',...
    'b_x1', 'b_x2', 'N_basis_x', 'N_basis_y', 'N_intervals_x', 'N_intervals_y',...
    'ellipse_parameters', 'n_points', 'eta', 'f1', 'verbose',...
    'Nlambda_eig', 'Nlambda_perturb', 'half_n_lambda', 'n_lambda_extra_perturb',...
    'Ntheta_eig',  'Ntheta_perturb',  'half_n_theta',  'n_theta_extra_perturb',...
    'Nphi_eig',    'Nphi_perturb',    'half_n_phi',    'n_phi_extra_perturb')
%}
%{
load('ellipse_plasmonic.mat','figure_shape', 'dispersion', 'lambda', 'theta', 'phi', 'delta',...
    'h', 'L', 'N_FMM', 'epsilon', 'refIndices', 'La', 'tau_x', 'tau_y',...
    'alpha_ref', 'beta_ref',...
    'b_x1', 'b_x2', 'N_basis_x', 'N_basis_y', 'N_intervals_x', 'N_intervals_y',...
    'ellipse_parameters',...
    'n_points', 'eta', 'f1', 'verbose')
%}
%calculate reflection and transmission

[Rsum,Tsum, matrix_Au_layer, gamma_num,...
    H_sorted,gammasqr_sorted] = ...
    PMM_main_function(figure_shape, dispersion, lambda, theta, phi, delta,...
    h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
    b_x1, b_x2, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y,ellipse_parameters,...
    n_points, eta, f1, verbose,...
    Nlambda_eig, Nlambda_perturb, half_n_lambda, n_lambda_extra_perturb,...
    Ntheta_eig,  Ntheta_perturb,  half_n_theta,  n_theta_extra_perturb,...
    Nphi_eig,    Nphi_perturb,    half_n_phi,    n_phi_extra_perturb);

if (n_lambda_extra_perturb == 1)&&(n_theta_extra_perturb == 1)
    Rsum_nonperturbed=Rsum;
    gamma_num_nonperturbed = gamma_num;
    H_nonperturbed_sorted= H_sorted;
    gammasqr_nonperturbed_sorted = gammasqr_sorted;
else
    Rsum_perturbed=Rsum;
    gamma_num_perturbed = gamma_num;
    H_perturbed_sorted= H_sorted;
    gammasqr_perturbed_sorted = gammasqr_sorted;
end


%eig_max = max(eigenvalues_Au_layer);
load('eig_nonperturbed.mat', 'lambda','theta','gamma_num_nonperturbed',...
    'Rsum_nonperturbed','H_nonperturbed_sorted','gammasqr_nonperturbed_sorted');

[nx,Nx,N_total_x,N_total_x3] = PMM_number_of_basis_functions(N_intervals_x,N_basis_x);
[ny,Ny,N_total_y,N_total_y3] = PMM_number_of_basis_functions(N_intervals_y,N_basis_y);
N_total_3 = N_total_x3*N_total_y3;

%{
H_perturb_1_sorted_full(i_lambda_perturb,j_theta_perturb,:,:) = H_perturb_1_sorted;
gammasqr_1_sorted_full(i_lambda_perturb,j_theta_perturb,:) = gammasqr_1;
%}

difference_H = zeros(Ntheta_perturb,2*N_total_3);
difference_H_norm = zeros(Ntheta_perturb,2*N_total_3);
difference_eigenvalue = zeros(Ntheta_perturb,2*N_total_3);
difference_eigenvalue_norm = zeros(Ntheta_perturb,2*N_total_3);
eigenvector_nonperturb_for_max_diff_H = zeros(2*N_total_3,Ntheta_perturb);

for j_theta = 1:Ntheta_perturb
    for number_eig = 1:2*N_total_3
        difference(:,1)=H_perturbed_sorted(1,j_theta,:,number_eig)-...
            H_nonperturbed_sorted(1,j_theta,:,number_eig);
        H_nonperturbed_current(:,1) = H_nonperturbed_sorted(1,j_theta,:,number_eig);
        
        difference_H(j_theta,number_eig) = norm(difference);
        difference_H_norm(j_theta,number_eig) =...
            difference_H(j_theta,number_eig)/...
            norm(H_nonperturbed_current);
        
        difference_eigenvalue(j_theta,number_eig) =...
            abs(gammasqr_perturbed_sorted(1,j_theta,number_eig)-...
            gammasqr_nonperturbed_sorted(1,j_theta,number_eig));
        difference_eigenvalue_norm(j_theta,number_eig) =...
            difference_eigenvalue(j_theta,number_eig)/...
            abs(gammasqr_nonperturbed_sorted(1,j_theta,number_eig));
        
    end
    [difference_H_max(j_theta),num_diff_H(j_theta)] = max(difference_H(j_theta,:));
    eigenvalue_for_max_diff_H(j_theta) = gammasqr_nonperturbed_sorted(num_diff_H(j_theta));
    
    eigenvector_nonperturb_for_max_diff_H(:,j_theta) =...
        H_nonperturbed_sorted(1,j_theta,:,num_diff_H(j_theta));
    
    eigenvector_perturb_for_max_diff_H(:,j_theta) = H_perturbed_sorted(1,j_theta,:,num_diff_H(j_theta));
    
    [difference_H_max_norm(j_theta),num_diff_H_norm(j_theta)] = max(difference_H_norm(j_theta,:));
    eigenvalue_for_max_diff_H_norm(j_theta) = gammasqr_nonperturbed_sorted(num_diff_H_norm(j_theta));
    
    eigenvector_nonperturbed_for_max_diff_H_norm(:,j_theta) =...
        H_nonperturbed_sorted(1,j_theta,:,num_diff_H_norm(j_theta));
    
    [difference_eigenvalue_max(j_theta),num_diff_eigval(j_theta)] = max(difference_eigenvalue(j_theta,:));
    eigenvalue_for_max_diff_eigval(j_theta) =...
        gammasqr_nonperturbed_sorted(1,j_theta,num_diff_eigval(j_theta));
end

theta_deg = theta*180/pi;
for j_theta = 1:Ntheta_perturb
    for number_eig = 1:2*N_total_3
        Hc(:,1) = H_perturbed_sorted(1,j_theta,:,number_eig);
        H_perturbed_norm(number_eig,j_theta)=norm(Hc);
    end
    [diff_eig_minus_eig_with_max_H_diff(:,j_theta),numbers]=...
        sort(abs(gammasqr_nonperturbed_sorted(1,j_theta,:)-eigenvalue_for_max_diff_H_norm(j_theta)),'ascend');
    for i=1:2*N_total_3
        eig_closest_to_eig_with_max_H_diff(i,j_theta) = gammasqr_nonperturbed_sorted(1,j_theta,numbers(i));
    end
end
for i=1:2*N_total_3
    for j_theta=1:Ntheta_perturb
        gammasqr_nonperturbed_sorted_new(i,j_theta)=gammasqr_nonperturbed_sorted(1,j_theta,i);
        gammasqr_perturbed_sorted_new(i,j_theta)=gammasqr_perturbed_sorted(1,j_theta,i);
    end
end

for j_theta=1:Ntheta_perturb
    [H_perturbed_norm_sorted(:,j_theta),index] = sort(H_perturbed_norm(:,j_theta));
    for i=1:2*N_total_3
        eigenvalues_sorted_respect_to_H_norm(i,j_theta) = gammasqr_nonperturbed_sorted(1,j_theta,index(i));
    end
end



save('difference_for_theta.mat','theta_deg','difference_H','difference_H_norm',...
'difference_eigenvalue','difference_eigenvalue_norm',...
'gammasqr_perturbed_sorted_new','gammasqr_nonperturbed_sorted_new',...
'difference_H_max','eigenvalue_for_max_diff_H',...
'eigenvector_nonperturb_for_max_diff_H','eigenvector_perturb_for_max_diff_H',...
'difference_H_max_norm','eigenvalue_for_max_diff_H_norm',...
'difference_eigenvalue_max','eigenvalue_for_max_diff_eigval','H_perturbed_norm',...
'diff_eig_minus_eig_with_max_H_diff','eig_closest_to_eig_with_max_H_diff',...
'H_perturbed_norm_sorted','eigenvalues_sorted_respect_to_H_norm');

eigenvectors_perturbed_norms = H_perturbed_norm;
eigenvalues_nonperturbed_sorted = gammasqr_nonperturbed_sorted_new;
eigenvalues_perturbed_sorted = gammasqr_perturbed_sorted_new;

save('eigenvectors_norms.mat','eigenvectors_perturbed_norms',...
    'eigenvalues_nonperturbed_sorted','eigenvalues_perturbed_sorted','theta_deg')

%difference = abs(gamma_num_perturbed-gamma_num)./abs(gamma_num);

%{


difference_matrix_50_53_deg(:,:) = abs(matrix_Au_layer(1,1,:,:)-matrix_Au_layer(1,2,:,:));
relative_difference_matrix_50_53_deg = zeros(2*N_total_3,2*N_total_3);
for i=1:2*N_total_3
    for j=1:2*N_total_3
        if (matrix_Au_layer(1,1,i,j)~=0)
            relative_difference_matrix_50_53_deg(i,j) =...
                difference_matrix_50_53_deg(i,j)/abs(matrix_Au_layer(1,1,i,j));
        end
    end
end
matrix_50_deg_Au_layer = matrix_Au_layer(1,1,:,:);
matrix_53_deg_Au_layer = matrix_Au_layer(1,2,:,:);
max_matrix_50_deg = max(matrix_50_deg_Au_layer(:))
max_matrix_53_deg = max(matrix_53_deg_Au_layer(:))
max_diff_matrix_50_53_deg = max(difference_matrix_50_53_deg(:))
max_relative_diff_matrix_50_53_deg = max(relative_difference_matrix_50_53_deg(:))

eigenvalues_Au_layer_50_deg = (eigenvalues_Au_layer(1,1,:));
eigenvalues_Au_layer_53_deg = (eigenvalues_Au_layer(1,2,:));
diff_eigenvalues_Au_layer_50_53_deg = abs(eigenvalues_Au_layer_50_deg - eigenvalues_Au_layer_53_deg);
max_diff_eigenvalues_Au_layer_50_53_deg = max(diff_eigenvalues_Au_layer_50_53_deg(:))
relative_diff_eigenvalues_Au_layer_50_53_deg = diff_eigenvalues_Au_layer_50_53_deg./abs(eigenvalues_Au_layer_50_deg);
max_relative_diff_eigenvalues_Au_layer_50_53_deg = max(relative_diff_eigenvalues_Au_layer_50_53_deg(:))
%reflection_1200_nm = Rsum(1);
%reflection_1220_nm = Rsum(2);
save('matrix_difference_for_theta_50_53_deg.mat','max_matrix_50_deg','max_matrix_53_deg',...
    'eigenvalues_Au_layer_50_deg','eigenvalues_Au_layer_53_deg','max_diff_eigenvalues_Au_layer_50_53_deg');
    %}
%{
    difference_matrix_1200_1220_nm(:,:) = abs(matrix_Au_layer(1,1,:,:)-matrix_Au_layer(2,1,:,:));
relative_difference_matrix_1200_1220_nm = zeros(2*N_total_3,2*N_total_3);
for i=1:2*N_total_3
    for j=1:2*N_total_3
        if (matrix_Au_layer(1,1,i,j)~=0)
            relative_difference_matrix_1200_1220_nm(i,j) =...
                difference_matrix_1200_1220_nm(i,j)/abs(matrix_Au_layer(1,1,i,j));
        end
    end
end
matrix_1200_nm_Au_layer = matrix_Au_layer(1,1,:,:);
matrix_1220_nm_Au_layer = matrix_Au_layer(2,1,:,:);
max_matrix_1200_nm = max(matrix_1200_nm_Au_layer(:))
max_matrix_1220_nm = max(matrix_1220_nm_Au_layer(:))
max_diff_matrix_1200_1220_nm = max(difference_matrix_1200_1220_nm(:))
max_relative_diff_matrix_1200_1220_nm = max(relative_difference_matrix_1200_1220_nm(:))

eigenvalues_Au_layer_1200_nm = (eigenvalues_Au_layer(1,1,:));
eigenvalues_Au_layer_1220_nm = (eigenvalues_Au_layer(2,1,:));
diff_eigenvalues_Au_layer_1200_1220_nm = abs(eigenvalues_Au_layer_1200_nm - eigenvalues_Au_layer_1220_nm);
max_diff_eigenvalues_Au_layer_1200_1220_nm = max(diff_eigenvalues_Au_layer_1200_1220_nm(:))
relative_diff_eigenvalues_Au_layer_1200_1220_nm = diff_eigenvalues_Au_layer_1200_1220_nm./abs(eigenvalues_Au_layer_1200_nm);
max_relative_diff_eigenvalues_Au_layer_1200_1220_nm = max(relative_diff_eigenvalues_Au_layer_1200_1220_nm(:))
reflection_1200_nm = Rsum(1);
reflection_1220_nm = Rsum(2);
%}
    %{
save('matrix_difference_for_theta_50_53_deg.mat','max_matrix_50_deg','max_matrix_53_deg',...
    'eigenvalues_Au_layer_50_deg','eigenvalues_Au_layer_53_deg','max_diff_eigenvalues_Au_layer_50_53_deg');
%}
    %{
for i=1:L
    gammaminus(:,i)= sort(gammaminus(:,i));
end
%}
%load('ellipse_plasmonic_no_perturb_output.mat','Rsum','Tsum','theta');
%load('ellipse_plasmonic_perturb_1order_output.mat','Rsum_p','Tsum_p', 'gammaminus_p')
%diff_R = abs((Rsum-Rsum_p)./Rsum);
%diff_theta = abs(theta - theta_p);



%load('ellipse_plasmonic_output.mat','Rsum_ellipses','Tsum_ellipses',...
%    'M', 'gammaminus');
%difference_R_sum = abs(Rsum_1_4-Rsum);
%difference_T_sum = abs(Tsum_1_4-Tsum);
%{
[l, ll] = size(lambda)
[t,tt] = size(theta)
[r,rr] = size(Rsum_ellipses)
n_lambda_theta = 100;
tl = linspace (lmin,lmax,n_lambda_theta)/1000;
tt = linspace(tmin,tmax,n_lambda_theta);
[XI,YI] = meshgrid(tl,tt);
ZI = griddata(lambda/1000,theta,transpose(Rsum_ellipses),XI,YI);
figure(1);
pcolor(XI,YI*180/pi,ZI)

xlabel('lambda for R');
ylabel('theta');
shading flat
caxis([0 1])
colorbar
%}
%{
figure(4);
pcolor(lambda/1000,theta*180/pi,transpose(Rsum_p))

xlabel('lambda for R');
ylabel('theta');
shading flat
caxis([0 1])
colorbar
hold off
%}
toc

figure(1)
plot(lambda, Rsum, 'r', 'Linewidth', 2);
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off


