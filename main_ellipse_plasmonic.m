
clc
clear all
% Find all windows of type figure, which have an empty FileName attribute.
allPlots = findall(0, 'Type', 'figure', 'FileName', []);
% Close.
delete(allPlots);

verbose = 5;

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

lmin = 600; %in nm
lmax = 1600;
lambda = linspace(lmin,lmax,50);
[Nll, Nl] = size(lambda);

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

%theta = 55*pi/180;
tmin = 0*pi/180;
tmax = 85*pi/180;
theta = linspace(tmin,tmax,50);
phi = 45*pi/180;

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

delta = 0;

refIndices = [nZnSe nSiO2];

L=5; %number of layers
%epsilon(iL,1,iNl) = eps outside the ellipse
%epsilon(iL,2,iNl) = eps inside the ellipse

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

%calculate reflection and transmission
[Rsum_ellipses,Tsum_ellipses] = ...
    PMM_main_function(figure_shape, dispersion, lambda, theta, phi, delta,...
    h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
    b_x1, b_x2, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y,ellipse_parameters,...
    n_points, eta, f1, verbose);

[l, ll] = size(lambda)
[t,tt] = size(theta)
[r,rr] = size(Rsum_ellipses)
tl = linspace (lmin,lmax,200)/1000;
tt = linspace(tmin,tmax,200);
[XI,YI] = meshgrid(tl,tt);
ZI = griddata(lambda/1000,theta,transpose(Rsum_ellipses),XI,YI);
figure(1);
pcolor(XI,YI*180/pi,ZI)
xlabel('lambda');
ylabel('theta');
shading flat
caxis([0 1])
hold off
%{
figure(1)
plot(lambda, Rsum_ellipses, 'r', 'Linewidth', 2);
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off
    %}

