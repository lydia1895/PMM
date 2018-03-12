
clc
clear all
% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);

verbose = 8;

format long
eta=0;
f1=0;
%%%%%%%%%
Si_dispersion = xlsread('silicon_cryst_500-1500nm.xlsx');
%%%%%%%%%


Nlambda_eig = 21;
n_lambda_extra_perturb = 1;
Nlambda_perturb = n_lambda_extra_perturb * Nlambda_eig;
half_n_lambda = floor((n_lambda_extra_perturb-1)/2);

Ntheta_eig = 26;
n_theta_extra_perturb = 1;
Ntheta_perturb = n_theta_extra_perturb * Ntheta_eig;
half_n_theta = floor((n_theta_extra_perturb-1)/2);

Nphi_eig = 1;
n_phi_extra_perturb = 1;
Nphi_perturb = n_phi_extra_perturb * Nphi_eig;
half_n_phi = floor((n_phi_extra_perturb-1)/2);

lmin = 1200;
lmax = 1500;
lambda = linspace(lmin,lmax,Nlambda_perturb);


tmin = 5*pi/180;
tmax = 80*pi/180;
theta = linspace(tmin,tmax,Ntheta_perturb);

phi = 0*pi/180;



Nl = Nlambda_perturb;
n_Si = zeros(Nl,1);
eps_Si = zeros(Nl,1);

n_media = 1.66*ones(Nl,1);
eps_media = n_media.^2;

n_prism = 3.5*ones(Nl,1);
eps_prism = n_prism.^2;

Si_lambda = Si_dispersion(:,1)*1000;

for i=1:Nl
    [ll,num] = min( abs (lambda(i)-Si_lambda(:) ) );
    llambda = lambda(i)
    si_llambda = Si_lambda(num)
    n_Si(i) = Si_dispersion(num,2) + 1j*Si_dispersion(num,3);
    eps_Si(i) = Si_dispersion(num,5) + 1j*Si_dispersion(num,6);
end
%%%%%%%%%conditions

figure_shape = 'ellipse';
dispersion = 'yes';

N_intervals_x = 3;
N_intervals_y = 3;
N_b = 6;
n_points = 150;
N_basis_x = N_b*ones(N_intervals_x,1);
N_basis_y = N_b*ones(N_intervals_y,1);



R1 = 242;
R2 = 242;
P1 = 666;
P2 = 666;
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

refIndices = [n_prism n_media];

L=3; %number of layers
%epsilon(iL,1,iNlambda) = eps outside the ellipse
%epsilon(iL,2,iNlambda) = eps inside the ellipse

epsilon = zeros(L,2,Nl);

epsilon(4,1,:) = eps_prism;
epsilon(4,2,:) = eps_prism;

epsilon(3,1,:) = eps_media;
epsilon(3,2,:) = eps_media;

epsilon(2,1,:) = eps_media;
epsilon(2,2,:) = eps_Si;

epsilon(1,1,:) = eps_media;
epsilon(1,2,:) = eps_media;

h = zeros(L,1);
h(4) = 0.0;
h(3) = 450;
h(2) = 220;
h(1) = 0.0;

alpha_ref = -sin(pi/6)/periodx;
beta_ref =  -sin(pi/6)/periody;

tau_x = exp(1j*alpha_ref*periodx);
tau_y = exp(1j*beta_ref*periody);

%La is lambda in Gegenbauer polynomials; La>-1/2

La = 0.5;

N_FMM = 1;

[Rsum,Tsum, matrix_Au_layer, eigenvalues_Au_layer] = ...
    PMM_main_function(figure_shape, dispersion, lambda, theta, phi, delta,...
    h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
    b_x1, b_x2, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y,ellipse_parameters,...
    n_points, eta, f1, verbose,...
    Nlambda_eig, Nlambda_perturb, half_n_lambda, n_lambda_extra_perturb,...
    Ntheta_eig,  Ntheta_perturb,  half_n_theta,  n_theta_extra_perturb,...
    Nphi_eig,    Nphi_perturb,    half_n_phi,    n_phi_extra_perturb);



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

figure(4);
pcolor(lambda/1000,theta*180/pi,transpose(Rsum))

xlabel('lambda, mkm');
ylabel('theta, deg');
colormap('jet');
colorbar;
set(gca,'fontsize', 16)
shading flat
caxis([0 1])
colorbar
hold off


figure(1)
plot(lambda, Rsum, 'r', 'Linewidth', 2);
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off

