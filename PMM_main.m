clc
clear all

%%%%%%%%%conditions

figure_shape = 'rectangle';
dispersion = 'no';
ellipse_parameters = zeros(6,1);

N_intervals_x = 2;
N_intervals_y = 2;
N_b = 4;
N_basis_x = N_b*ones(N_intervals_x,1);
N_basis_y = N_b*ones(N_intervals_y,1);

%boundaries in the periodic layer
%for ellipse here will be matched coordinates

lambda = 2*pi;
theta = linspace(0,89,10)*pi/180;
phi = 0*pi/180;

%{
lambda = linspace(1400,1500,70);
theta = 4*pi/180;
phi = 0*pi/180;
%}
b_x = zeros(N_intervals_x+1);
b_y = zeros(N_intervals_y+1);

b_x = [0.0 150.0 1000.0];  
b_y = [0.0 150.0 1000.0];
%b_x = [0 0.25 0.9]*lambda;
%b_y = [0 0.25 0.9]*lambda;

[Nxx, NNxx] = size(b_x);
[Nyy, NNyy] = size(b_y);
periodx = b_x(NNxx)-b_x(1);
periody = b_y(NNyy)-b_y(1);

%delta is the angle between E and the incidence plane
%delta = 0 TM, delta = pi/2 TE

delta = pi/2;
%epsilon = zeros(N_intervals_x, N_intervals_y);
%{
nSi = 3.48;
nSiO2 = 1.44;
eSi = nSi^2;
eSiO2 = nSiO2^2;
refIndices = [1.0 nSi];
epsilon(:,:,5) = [1.0 1.0; 1.0 1.0];  %upper layer - wave comes from this media
epsilon(:,:,4) = eSiO2*[1.0 1.0; 1.0 1.0];
epsilon(:,:,3) = [eSi eSi; eSiO2 eSiO2];  
epsilon(:,:,2) = eSiO2*[1.0 1.0; 1.0 1.0];
epsilon(:,:,1) = eSi*[1.0 1.0; 1.0 1.0]; %lower layer
%}

refIndices = [1.0 2.0];
epsilon(:,:,2) = [1.0 1.0; 1.0 1.0];  %upper layer - wave comes from this media
epsilon(:,:,1) = 4.0*[1.0 1.0; 1.0 1.0];  %lower layer

%epsilon(:,:,1) = [1.0 1.0; 1.0 1.0];
%epsilon(:,:,2) = [1.0 2.25; 2.25 2.25];
%epsilon(:,:,3) = [1.0 1.0; 1.0 1.0];

%L is number of layers including half-infinite medias
%numerate layers from lower one to the upper one

L=2; 
h(1) = 0.0;  %lower layer
h(2) = 0.0;

%{
L=5;
h(5) = 0.0; %upper layer
h(4) = 3000;
h(3) = 350;
h(2) = 1000;
h(1) = 0.0;   %lower layer
%}

%arbitrary boundary conditions tau
  
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

%calculate reflection and transmission

[Rsum,Tsum] = ...
    PMM_main_function(figure_shape, dispersion, lambda, theta, phi, delta,...
    h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
    b_x, b_y, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y, ellipse_parameters);
%{
derxx = hx\Dx;

%k0 = 2*pi/lambda;
%{
g1 = unique(diag(gamma(:,:,1)))/k0;
g2 = unique(diag(gamma(:,:,2)))/k0;
%g3 = unique(diag(gamma(:,:,3)))/k0;
%}
%kz1 = kz1v/k0;

N_total = sum(N_basis_x);
N_total_3=(N_total-N_intervals_x)^2;


[row,col,v] = find(M(:,:,1)-M(:,:,2));
%[row13,col13,v13] = find(M31);

nil_norm = M(:,:,1)*EH(:,:,1)-EH(:,:,1)*gamma_norm(:,:,1);
%nil_sorted = M(:,:,1)*W(:,:,1)-W(:,:,1)*gamma_sorted(:,:,1)/k0;
g1 = sort(diag(gamma_norm(:,:,1)));
g2 = sort(diag(gamma_norm(:,:,2)));
%g3 = unique(diag(gamma_norm(:,:,3)));


kz1vsort = sort(kz1v)/k0;

NN=(2*N_FMM+1)^2;

u2_1_PMM = ud_PMM(1:N_total_3);
u2_2_PMM = ud_PMM(N_total_3+1   : 2*N_total_3);
d0_1_PMM = ud_PMM(2*N_total_3+1 : 3*N_total_3);
d0_2_PMM = ud_PMM(3*N_total_3+1 : 4*N_total_3);
ud_PMM_total = cat(2,u2_1_PMM,u2_2_PMM,d0_1_PMM,d0_2_PMM);

u2_1_FMM = u2d0_FMM(1:NN);
u2_2_FMM = u2d0_FMM(NN+1   : 2*NN);
d0_1_FMM = u2d0_FMM(2*NN+1 : 3*NN);
d0_2_FMM = u2d0_FMM(3*NN+1 : 4*NN);
ud_FMM_total = cat(2,u2_1_FMM,u2_2_FMM,d0_1_FMM,d0_2_FMM);
%}
figure(1)
plot(theta*180/pi, Rsum, 'r', theta*180/pi, Tsum, 'g', 'Linewidth', 2);
%[NNt, Nt] = size(theta)
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off

    
