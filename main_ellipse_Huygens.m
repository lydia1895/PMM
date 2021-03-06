
clc
clear all
% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);

verbose = 6;

format long
eta=0;
f1=0;
%%%%%%%%%conditions

figure_shape = 'ellipse';
dispersion = 'no';

N_intervals_x = 3;
N_intervals_y = 3;
N_b = 6;
n_points = 150;
N_basis_x = N_b*ones(N_intervals_x,1);
N_basis_y = N_b*ones(N_intervals_y,1);

lambda = linspace(1000,1700,100);
theta = 0*pi/180;
phi = 0*pi/180;

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

nSi = 3.5;
nSiO2 = 1.45;
epsSi = nSi^2;
epsSiO2 = nSiO2^2;
epsAir = 1.0;
refIndices = [1.0 nSi];

L=4;

epsilon = zeros(L,2);
epsilon(4,:) = epsAir*ones(1,2);
epsilon(3,:) = epsAir*ones(1,2);
epsilon(2,:) = epsSiO2*ones(1,2);
epsilon(1,:) = epsSi*ones(1,2);
epsilon(3,2) = epsSi; %ellipses

h = zeros(L,1); 
h(4) = 0.0;
h(3) = 220;
h(2) = 2000;
h(1) = 0.0;

alpha_ref = -sin(pi/6)/periodx;
beta_ref =  -sin(pi/6)/periody;

tau_x = exp(1j*alpha_ref*periodx);
tau_y = exp(1j*beta_ref*periody);

%La is lambda in Gegenbauer polynomials; La>-1/2

La = 0.5;

N_FMM = 1;

%calculate reflection and transmission
[Rsum_ellipses,Tsum_ellipses] = ...
    PMM_main_function(figure_shape, dispersion, lambda, theta, phi, delta,...
    h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
    b_x1, b_x2, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y,ellipse_parameters,...
    n_points, eta, f1, verbose);

L=3;

epsilon = zeros(L,2);
epsilon(3,:) = epsAir*ones(1,2);
epsilon(2,:) = epsSiO2*ones(1,2);
epsilon(1,:) = epsSi*ones(1,2);

h = zeros(L,1); 
h(3) = 0.0;
h(2) = 2000;
h(1) = 0.0;
[Rsum_etched,Tsum_etched] = ...
    PMM_main_function(figure_shape, dispersion, lambda, theta, phi, delta,...
    h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
    b_x1, b_x2, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y,ellipse_parameters,...
    n_points, eta, f1, verbose);

T_ref_T0 = Tsum_ellipses./Tsum_etched;
figure(1)
plot(lambda, transpose(T_ref_T0), 'r', 'Linewidth', 2);
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off

