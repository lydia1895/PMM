clc
clear all

%%%%%%%%%conditions

N_intervals_x = 2;
N_intervals_y = 2;
N_b = 9;
N_basis_x = N_b*ones(N_intervals_x,1);
N_basis_y = N_b*ones(N_intervals_y,1);

%boundaries in the periodic layer
%for ellipse here will be matched coordinates

%lambda = 0.01;
Nl = 1;
%lambda = 1*10^(-6);
lambda = 1;

b_x = zeros(N_intervals_x+1);
b_y = zeros(N_intervals_y+1);

b_x = [0.0 0.25*lambda 0.9*lambda];  
b_y = [0.0 0.25*lambda 0.9*lambda];
%b_x = [0 0.25 0.8 0.9]*lambda;
%b_y = [0 0.25 0.8 0.9]*lambda;

%{
b_x = [0.0 100*lambda 1000*lambda];    
b_y = [0.0 100*lambda 1000*lambda];
%}
%b_y = [0 0.9*lambda];

[Nxx, NNxx] = size(b_x);
[Nyy, NNyy] = size(b_y);
periodx = b_x(NNxx)-b_x(1);
periody = b_y(NNyy)-b_y(1);

%epsilon = zeros(N_intervals_x, N_intervals_y);

%{
refIndices = [1.0 2.0];
epsilon(:,:,2) = [1.0 1.0; 1.0 1.0];  %upper layer - wave comes from this media
epsilon(:,:,1) = 4.0*[1.0 1.0; 1.0 1.0];  %lower layer
%}
%{
refIndices = [1.0 1.0];
epsilon(:,:,1) = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]; 
epsilon(:,:,2) = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]; 
epsilon(:,:,3) = [1.0 1.0 1.0; 1.0 1.0 1.0; 1.0 1.0 1.0]; 
%}
refIndices = [1.0 1.0];
epsilon(:,:,1) = [1.0 1.0; 1.0 1.0];
epsilon(:,:,2) = [1.0 2.25; 2.25 2.25];
epsilon(:,:,3) = [1.0 1.0; 1.0 1.0];

%L is number of layers including half-infinite medias
%numerate layers from lower one to the upper one
%{
L=2; 
h(1) = 0.0;  %lower layer
h(2) = 0.0;
%}

L=3;
h(1) = 0.0;  %lower layer
h(2) = 0.25*lambda;
h(3) = 0.0;  %upper layer

%enter lambda in nm, theta, phi in radians

%{
lmin = 1400*10^(-9);
lmax = 1500*10^(-9);
lambda = linspace(lmin, lmax, 350);
[Nll,Nl] = size(lambda);
%}
theta = 0*pi/180;
phi = 0*pi/180;

%arbitrary boundary conditions tau

%alpha_ref = pi/(2*periodx); %1.0/periodx;
%beta_ref =  pi/(2*periodx); %1.0/periodx;

k0 = 2*pi/lambda;

alpha_ref = -k0*sin(pi/6);
beta_ref =  -k0*sin(pi/6);

tau_x = exp(1j*alpha_ref*periodx);
tau_y = exp(1j*beta_ref*periody);

%La is lambda in Gegenbauer polynomials; La>-1/2

La = 0.5;

%delta is the angle between E and the incidence plane
%delta = 0 TM, delta = pi/2 TE

delta = pi/2;

%N_FMM is the number of Fourier harmonics in reflected or transmitted wave
%we go from Gegenbauer basis to Fourier basis
%to get diffraction in certain orders

N_FMM = 12;

%calculate reflection and transmission

[R00, Rsum, T00, Tsum, eta_R, eta_T, gamma_norm, EH, gamma_sorted, W, Stotal,...
    ud_FMM, ud_PMM, kz1v, kz2v, hx, Dx, E_inc_PMM_to_FMM_1,E_inc_PMM_to_FMM_2,...
    E_PMM, nx, Nx, M, M31, pplus, pminus, derx, P_dPx, P_dPy, ax, ay, eps] = ...
    PMM_main_function(lambda, theta, phi, delta,...
    h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
    b_x, b_y, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y);

derxx = hx\Dx;

k0 = 2*pi/lambda;
%{
g1 = unique(diag(gamma(:,:,1)))/k0;
g2 = unique(diag(gamma(:,:,2)))/k0;
%g3 = unique(diag(gamma(:,:,3)))/k0;
%}
%kz1 = kz1v/k0;

N_total = sum(N_basis_x);
N_total_3=(N_total-N_intervals_x)^2;

sum = 0
[n,nn] = size(ud_PMM)
for i=1:n
    sum = sum+(ud_PMM(i))^2;
end

[row,col,v] = find(M(:,:,1)-M(:,:,2));
[row13,col13,v13] = find(M31);

nil_norm = M(:,:,1)*EH(:,:,1)-EH(:,:,1)*gamma_norm(:,:,1);
nil_sorted = M(:,:,1)*W(:,:,1)-W(:,:,1)*gamma_sorted(:,:,1)/k0;
g1 = sort(diag(gamma_norm(:,:,1)));
g2 = sort(diag(gamma_norm(:,:,2)));
%g3 = unique(diag(gamma_norm(:,:,3)));

