clc
clear all

% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);

verbose = 5;

format long
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
theta = linspace(1,89,20)*pi/180;
%theta = linspace(40.1255,40.1257,2)*pi/180;  %for eta=0.001, b_x=[0 300 1000]
phi = 0*pi/180;


%ASR parameters
%eta = 0;
%eta = 0.0005;
eta = 0.001;
f1 = 0.5;
n_points=10;
%{
lambda = linspace(1400,1500,70);
theta = 4*pi/180;
phi = 0*pi/180;
%}
b_x = zeros(N_intervals_x+1);
b_y = zeros(N_intervals_y+1);

b_x = [0.0 300.0 1000.0];
b_y = [0.0 300.0 1000.0];
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
refIndices = [1.0 1.0];
epsilon(:,:,2) = [1.0 1.0; 1.0 1.0];  %upper layer - wave comes from this media
epsilon(:,:,1) = 1.0*[1.0 1.0; 1.0 1.0];  %lower layer

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
    b_x, b_y, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y, ellipse_parameters,...
    n_points, eta, f1, verbose);

figure(1)
plot(theta*180/pi, Rsum, '-sr', theta*180/pi, Tsum, '-sg', 'Linewidth', 2);
%[NNt, Nt] = size(theta)
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off


