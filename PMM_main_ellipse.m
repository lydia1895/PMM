clc
clear all
% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);

%format long
%%%%%%%%%conditions

verbose = 7;

figure_shape = 'ellipse';
dispersion = 'no';
n_points = 150;
N_intervals_x = 3;
N_intervals_y = 3;
N_b = 6;
N_basis_x = N_b*ones(N_intervals_x,1);
N_basis_y = N_b*ones(N_intervals_y,1);
lambda = 2*pi;
%theta = linspace(26.5,26.9,10)*pi/180;
%theta = linspace(0,26.85,2)*pi/180;

Nlambda_eig = 1;
n_lambda_extra_perturb = 1;
Nlambda_perturb = n_lambda_extra_perturb * Nlambda_eig;
half_n_lambda = floor((n_lambda_extra_perturb-1)/2);

Ntheta_eig = 40;
n_theta_extra_perturb = 1;
Ntheta_perturb = n_theta_extra_perturb * Ntheta_eig;
half_n_theta = floor((n_theta_extra_perturb-1)/2);

Nphi_eig = 1;
n_phi_extra_perturb = 1;
Nphi_perturb = n_phi_extra_perturb * Nphi_eig;
half_n_phi = floor((n_phi_extra_perturb-1)/2);

theta = linspace(0,89,Ntheta_perturb)*pi/180;

phi = 0*pi/180;

%for ASR
eta=0;
f1=0;

R1 = 1;
R2 = 1;
P1 = 7;
P2 = 7;
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

delta = 0;


refIndices = [1.0 2.0];
epsilon(2,:)=[1.0 1.0];
epsilon(1,:)=[4.0 4.0];
%eps=epsilon;
%eps_out = epsilon(1);
%eps_in =  epsilon(2);

%L is number of layers including half-infinite medias
%numerate layers from lower one to the upper one

L=2; 
h(2) = 0.0;  %lower layer
h(1) = 0.0;

alpha_ref = -sin(pi/6)/periodx;
beta_ref =  -sin(pi/6)/periody;

tau_x = exp(1j*alpha_ref*periodx);
tau_y = exp(1j*beta_ref*periody);

%La is lambda in Gegenbauer polynomials; La>-1/2

La = 0.5;

%N_FMM is the number of Fourier harmonics in reflected or transmitted wave
%we go from Gegenbauer basis to Fourier basis
%to get diffraction in certain orders

N_FMM = 1;

%{
load('ellipse_2medias.mat','figure_shape', 'dispersion', 'lambda', 'theta', 'phi', 'delta',...
    'h', 'L', 'N_FMM', 'epsilon', 'refIndices', 'La', 'tau_x', 'tau_y',...
    'alpha_ref', 'beta_ref',...
    'b_x1', 'b_x2', 'N_basis_x', 'N_basis_y', 'N_intervals_x', 'N_intervals_y',...
    'ellipse_parameters',...
    'n_points', 'eta', 'f1', 'verbose')
%}
%{
load('ellipse_2medias_normal_incidence.mat','figure_shape', 'dispersion', 'lambda', 'theta', 'phi', 'delta',...
    'h', 'L', 'N_FMM', 'epsilon', 'refIndices', 'La', 'tau_x', 'tau_y',...
    'alpha_ref', 'beta_ref',...
    'b_x1', 'b_x2', 'N_basis_x', 'N_basis_y', 'N_intervals_x', 'N_intervals_y',...
    'ellipse_parameters',...
    'n_points', 'eta', 'f1', 'verbose')
%}
%{
load('ellipse_2medias.mat','figure_shape', 'dispersion', 'lambda', 'theta', 'phi', 'delta',...
    'h', 'L', 'N_FMM', 'epsilon', 'refIndices', 'La', 'tau_x', 'tau_y',...
    'alpha_ref', 'beta_ref',...
    'b_x1', 'b_x2', 'N_basis_x', 'N_basis_y', 'N_intervals_x', 'N_intervals_y',...
    'ellipse_parameters',...
    'n_points', 'eta', 'f1', 'verbose')
%}
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
[Rsum,Tsum,M,gammaminus,dq] = ...
    PMM_main_function(figure_shape, dispersion, lambda, theta, phi, delta,...
    h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
    b_x1, b_x2, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y,ellipse_parameters,...
    n_points, eta, f1, verbose);
save('ellipse_2medias_normal_incidence_output.mat','Rsum','Tsum',...
    'M', 'gammaminus','dq');
figure(1)
plot(theta*180/pi, Rsum, '-sr', theta*180/pi, Tsum, '-sg', 'Linewidth', 2);
xlabel('theta')
ylabel('R,T')
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[Rsum_1_4,Tsum_1_4,M_1_4,gammaminus_1_4] = ...
    PMM_main_function(figure_shape, dispersion, lambda, theta, phi, delta,...
    h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
    b_x1, b_x2, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y,ellipse_parameters,...
    n_points, eta, f1, verbose,...
     Nlambda_eig, Nlambda_perturb, half_n_lambda, n_lambda_extra_perturb,...
    Ntheta_eig,  Ntheta_perturb,  half_n_theta,  n_theta_extra_perturb,...
    Nphi_eig,    Nphi_perturb,    half_n_phi,    n_phi_extra_perturb);
for i=1:L
    gammaminus_1_4(:,i)= sort(gammaminus_1_4(:,i));
end
figure(1)
plot(theta*180/pi, Rsum_1_4, '-sr', theta*180/pi, Tsum_1_4, '-sg', 'Linewidth', 2);
xlabel('theta')
ylabel('R,T')
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off
%{
load('ellipse_2medias_normal_incidence_output.mat','Rsum','Tsum',...
    'M', 'gammaminus');
gamma_should_be_zero = abs(gammaminus_1_4 - gammaminus);
save('ellipse_2medias_gamma_difference','gammaminus','gammaminus_1_4','gamma_should_be_zero');
%}
