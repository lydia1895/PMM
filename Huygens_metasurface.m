
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


<<<<<<< HEAD

Nlambda_eig = 50;

=======
Nlambda_eig = 15;
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
n_lambda_extra_perturb = 1;
Nlambda_perturb = n_lambda_extra_perturb * Nlambda_eig;
half_n_lambda = floor((n_lambda_extra_perturb-1)/2);

<<<<<<< HEAD

Ntheta_eig = 60;

=======
Ntheta_eig = 1;
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
n_theta_extra_perturb = 1;
Ntheta_perturb = n_theta_extra_perturb * Ntheta_eig;
half_n_theta = floor((n_theta_extra_perturb-1)/2);

Nphi_eig = 1;
n_phi_extra_perturb = 1;
Nphi_perturb = n_phi_extra_perturb * Nphi_eig;
half_n_phi = floor((n_phi_extra_perturb-1)/2);

<<<<<<< HEAD

lmin = 700;
lmax = 1200;
lambda = linspace(lmin,lmax,Nlambda_perturb);


tmin = 20*pi/180;
tmax = 80*pi/180;

=======
lmin = 1200;
lmax = 1500;
lambda = linspace(lmin,lmax,Nlambda_perturb);


tmin = 0*pi/180;
tmax = 0*pi/180;
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
theta = linspace(tmin,tmax,Ntheta_perturb);

phi = 0*pi/180;



Nl = Nlambda_perturb;
n_Si = zeros(Nl,1);
eps_Si = zeros(Nl,1);

<<<<<<< HEAD

n_media = 1.0*ones(Nl,1);
eps_media = n_media.^2;

n_prism = 2.5*ones(Nl,1);

=======
n_media = 1.66*ones(Nl,1);
eps_media = n_media.^2;

n_prism = 3.5*ones(Nl,1);
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
eps_prism = n_prism.^2;

Si_lambda = Si_dispersion(:,1)*1000;

for i=1:Nl
    [ll,num] = min( abs (lambda(i)-Si_lambda(:) ) );
<<<<<<< HEAD

    n_Si(i) = Si_dispersion(num,2) + 1j*Si_dispersion(num,3);
    eps_Si(i) = Si_dispersion(num,5) + 1j*Si_dispersion(num,6);

=======
    llambda = lambda(i)
    si_llambda = Si_lambda(num)
    n_Si(i) = Si_dispersion(num,2); %+ 1j*Si_dispersion(num,3);
    eps_Si(i) = Si_dispersion(num,5);% + 1j*Si_dispersion(num,6);
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
end

%n_Si = 4*ones(Nl,1);
%eps_Si = n_Si.^2;
%%%%%%%%%conditions

figure_shape = 'ellipse';
dispersion = 'yes';

N_intervals_x = 3;
N_intervals_y = 3;
N_b = 6;
n_points = 1000;
N_basis_x = N_b*ones(N_intervals_x,1);
N_basis_y = N_b*ones(N_intervals_y,1);

<<<<<<< HEAD
%{
a = 420 nm
R = 190 nm
h = 300 nm
n1 = n2 = 1

�������� �� 700 �� 1200 �� (50 �����) � ���� �� 20 �� 80 �������� (60 �����).
%}

R1 = 190;
R2 = 190;
P1 = 420;
P2 = 420;

=======


R1 = 200;
R2 = 200;
P1 = 666;
P2 = 666;
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
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

<<<<<<< HEAD

refIndices = [n_prism n_media];

L=4; %number of layers

=======
refIndices = [n_media n_media];

L=3; %number of layers
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
%epsilon(iL,1,iNlambda) = eps outside the ellipse
%epsilon(iL,2,iNlambda) = eps inside the ellipse

epsilon = zeros(L,2,Nl);

%epsilon(4,1,:) = eps_prism;
%epsilon(4,2,:) = eps_prism;

<<<<<<< HEAD
epsilon(4,1,:) = eps_prism;
epsilon(4,2,:) = eps_prism;

=======
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
epsilon(3,1,:) = eps_media;
epsilon(3,2,:) = eps_media;

epsilon(2,1,:) = eps_media;
epsilon(2,2,:) = eps_Si;

epsilon(1,1,:) = eps_media;
epsilon(1,2,:) = eps_media;

h = zeros(L,1);
<<<<<<< HEAD

h(4) = 0.0;
h(3) = 300;
h(2) = 300;
=======
%h(4) = 0.0;
h(3) = 0.0;
h(2) = 220;
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
h(1) = 0.0;

alpha_ref = -sin(pi/6)/periodx;
beta_ref =  -sin(pi/6)/periody;

tau_x = exp(1j*alpha_ref*periodx);
tau_y = exp(1j*beta_ref*periody);

%La is lambda in Gegenbauer polynomials; La>-1/2

La = 0.5;

N_FMM = 3;

<<<<<<< HEAD
[Rsum,Tsum] = ...

=======
[Rsum,Tsum, matrix_Au_layer, eigenvalues_Au_layer] = ...
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
    PMM_main_function(figure_shape, dispersion, lambda, theta, phi, delta,...
    h, L, N_FMM, epsilon, refIndices, La, tau_x, tau_y, alpha_ref, beta_ref,...
    b_x1, b_x2, N_basis_x, N_basis_y, N_intervals_x, N_intervals_y,ellipse_parameters,...
    n_points, eta, f1, verbose,...
    Nlambda_eig, Nlambda_perturb, half_n_lambda, n_lambda_extra_perturb,...
    Ntheta_eig,  Ntheta_perturb,  half_n_theta,  n_theta_extra_perturb,...
    Nphi_eig,    Nphi_perturb,    half_n_phi,    n_phi_extra_perturb);

<<<<<<< HEAD

figure(1);
=======
%save('Huygens_results_prism_TM.mat','lambda','theta','Rsum');

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
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c
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

<<<<<<< HEAD


c = physconst('LightSpeed');
h = 4.135666 * 10^(-15);
kx = (2*a*n_prism./lambda)'*sin(theta);
frequency = (c*10^(-3)./lambda');

figure(2);
pcolor(kx,frequency,Rsum)
xlabel('kx, \pi/a');
ylabel('frequency, THz');


=======
%}
figure(1)
plot(lambda, Rsum, 'g', lambda, Tsum, 'r', 'Linewidth', 2);
%plot(lambda, transpose(Rsum), 'r', 'Linewidth', 2);
hold off
>>>>>>> c1b42c226c43825d6c5ba0f1cf23c22d3e6f123c


