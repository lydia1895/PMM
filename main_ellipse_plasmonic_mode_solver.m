
clc
clear all
% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);

verbose = 2;

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


Nlambda_eig = 2;
n_lambda_extra_perturb = 1;
Nlambda_perturb = n_lambda_extra_perturb * Nlambda_eig;
half_n_lambda = floor((n_lambda_extra_perturb-1)/2);

Ntheta_eig = 1;
n_theta_extra_perturb = 1;
Ntheta_perturb = n_theta_extra_perturb * Ntheta_eig;
half_n_theta = floor((n_theta_extra_perturb-1)/2);

Nphi_eig = 1;
n_phi_extra_perturb = 1;
Nphi_perturb = n_phi_extra_perturb * Nphi_eig;
half_n_phi = floor((n_phi_extra_perturb-1)/2);


theta = 0;
phi = 0;


Nl = Nlambda_perturb;
nAu = zeros(Nl,1);
epsAir = ones(Nl,1);
nAir = ones(Nl,1);
%%%%%%%%%conditions

figure_shape = 'ellipse';
dispersion = 'yes';

N_intervals_x = 3;
N_intervals_y = 3;
N_b = 6;
n_points = 150;
N_basis_x = N_b*ones(N_intervals_x,1);
N_basis_y = N_b*ones(N_intervals_y,1);



R1 = 250*10^(-9);
R2 = 250*10^(-9);
P1 = 700*10^(-9);
P2 = 700*10^(-9);
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

refIndices = [nAir nAir];

L=3; %number of layers
%epsilon(iL,1,iNlambda) = eps outside the ellipse
%epsilon(iL,2,iNlambda) = eps inside the ellipse

epsilon = zeros(L,2,Nl);


epsilon(3,1,:) = epsAir;
epsilon(3,2,:) = epsAir;

epsilon(2,1,:) = epsAir;


epsilon(1,1,:) = epsAir;
epsilon(1,2,:) = epsAir;

h = zeros(L,1);
h(3) = 0.0;
h(2) = 30*10^(-9);
h(1) = 0.0;

alpha_ref = -sin(pi/6)/periodx;
beta_ref =  -sin(pi/6)/periody;

tau_x = exp(1j*alpha_ref*periodx);
tau_y = exp(1j*beta_ref*periody);

%La is lambda in Gegenbauer polynomials; La>-1/2

La = 0.5;

N_FMM = 3;

c = 3*10^8;
w = 332*10^12;
dw = 0.005*10^12;

for ii=1:1

    lambda1 = c/w;
    lambda2 = c/(w+dw);

    lambda = [lambda1 lambda2]

    for i=1:Nl
        num = floor(real(lambda(i)*10^9)) - 500 + 1;
        nAu(i) = lamnkAu(num,2) + 1j*lamnkAu(num,3);
    end
    epsAu = nAu.^2;
    epsilon(2,2,:) = epsAu;

    [R0, R1] = ...
        PMM_main_function_mode_solver(figure_shape, dispersion, lambda, theta, phi, delta,...
        h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
        b_x1, b_x2, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y,ellipse_parameters,...
        n_points, eta, f1, verbose,...
        Nlambda_eig, Nlambda_perturb, half_n_lambda, n_lambda_extra_perturb,...
        Ntheta_eig,  Ntheta_perturb,  half_n_theta,  n_theta_extra_perturb,...
        Nphi_eig,    Nphi_perturb,    half_n_phi,    n_phi_extra_perturb);

    dR = (R1-R0)/dw;
    MAXDR = max(dR(:))
    ddw = dw
    A=-dR\R0;
    [V, D] = eig(A);
    w_array = diag(D);

    dw0 = min(w_array);
    w0 = w + dw0;
    %[w0,dw0]= PMM_mode_solver(R0,R1,w,dw);
    dw00(ii) = dw0;
    w = w0;
end
%{
figure(1)
plot(lambda, Rsum, 'r', 'Linewidth', 2);
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off
%}

