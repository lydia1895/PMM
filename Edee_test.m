clc
clear all

N = 1;                %number of Fourier orders
lambda = 1200*10^(-9);
Nl=1;

periodx = 0.9*lambda;
periody = 0.9*lambda;
dx = 0.25*lambda;
dy = 0.25*lambda;

L = 1;
h = zeros(L,1);
h(1) = 0.25*lambda;
M = 5001;
x = (1:1:M)*periodx/M;
epsilon = zeros(M, M, L);
epsnew = 2.25;
epsilon(:,:,1) = epsnew*ones(M,M);
for i=1:M
    for j=1:M
    if (x(i)<dx) && (x(j)<dy)
        epsilon(j,i,2) = 1.0;
    end
    end
end
refIndices = [1.0 1.0];
%{
periodx = 1000*(10^(-9));  %period of periodic layer
dx = 150*(10^(-9));     %ridge width
periody = 1000*(10^(-9));  %period of periodic layer
dy = 1000*(10^(-9));     
h = zeros(3,1);
h(3) = 3*10^(-6);
h(2) = 350*(10^(-9));       %thickness of periodic layer
h(1) = 1*10^(-6);
M = 101;               %number of modes for Fourier transform of epsilon
L = 3;                 %number of layers
x = (1:1:M)*periodx/M;
epsilon = zeros(M, M, L);
for i=1:M
    for j=1:M
    epsilon(i,j,1) = 1.44^2;
    if x(i)<dx
    %    epsilon(i,j,2) = 1.94^2;
        epsilon(j,i,2) = 3.48^2;
    else
        epsilon(j,i,2) = 1.44^2;
    end
    epsilon(i,j,3) = 1.44^2;
    end
end
refIndices = [1.0 3.48];
%}
%{
M = 301;
N = 2;                %number of Fourier orders
periodx = 1000*(10^(-9));  %period of periodic layer
dx = 150*(10^(-9));     %ridge width
periody = 1000*(10^(-9));  %period of periodic layer
dy = 1000*(10^(-9));     
epsilon = zeros(M,M,1);
L = 1;
h = 0.5*10^(-6);
for i=1:M
    for j = 1:M
        %{
        if ((i*periodx)/M<dx) && ((j*periody)/M<dy)
            epsilon(i,j,1) = 2.0;
        else
            epsilon(i,j,1) = 3.0;
    end
        %}
        epsilon(i,j,1) = 9.0;
    %epsilon(i,2) = 9.0;
    end
end

lambda = 1.064*(10^(-6));
refIndices = [3.0 1.0];     


Nl=1;
Nt=90;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);
theta=linspace(0,89,90)*pi/180;
[Ntt,Nt] = size(theta)
phi = 20*pi/180;
%}




%{
theta = linspace(1,15,15)*pi/180;
[Ntt, Nt] = size(theta)
%}
%{
lmin = 725*10^(-9);
lmax = 765*10^(-9);
lambda = linspace(lmin, lmax, 250);
[Nll,Nl] = size(lambda);

theta = [0 0.4]*pi/180;
Nt=2;
%}
theta = 0*pi/180;
Nt=1;
phi = 0*pi/180;
Np=1;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);

%lambda = 1400*10^(-9);
%Nl=1;
%theta = [0.4 2 4]*pi/180;
%Nt=3;
%{
phi = linspace(0,90,50)*pi/180;
[Npp,Np] = size(phi)
theta = linspace(0,89,50)*pi/180;
[Ntt, Nt] = size(theta)
lambda = 1400*10^(-9);
Nl = 1;
Rsum = zeros(Nt,Np);
Tsum = zeros(Nt,Np);
%}
%}
%{
Nl=1;
Nt=1;
Rsum=zeros(Nl,Nt);
Tsum=zeros(Nl,Nt);
lambda = 1400*10^(-9);
theta = 4*pi/180;
phi=0;
%}
P = 2*N+1;
Q = 2*N+1;
R = 1;

eps11=zeros(P*Q,P*Q,L);
eps22=zeros(P*Q,P*Q,L);
eps33=zeros(P*Q,P*Q,L);
for i=1:L
[eps11(:,:,i), eps22(:,:,i), eps33(:,:,i)] = FMM_eps123_new(epsilon(:,:,i),N,M);
end

for i=1:Nl
    for j=1:Nt
    for k=1:Np 
    [eta_R1, eta_T1, gamma] = FMM_2D_TE_RT_multi_test(eps11,eps22,eps33,periodx, periody, h, lambda(i), theta(j), phi(k), refIndices, N, M, L);
    %Rsum(j,k) = sum(eta_R1);
    %Tsum(j,k) = sum(eta_T1);
    Rsum(i,j) = sum(eta_R1);
    Tsum(i,j) = sum(eta_T1);
    end
    end
end
g = gamma
%{
figure(5)
hold on
%plot(lambda, Rsum(:,1), 'b', lambda, Rsum(:,2)+0.5, 'g', lambda, Rsum(:,3)+1, 'r', lambda, Rsum(:,4)+1.5, 'm', 'LineWidth', 2)
plot(lambda, Rsum(:,1), 'b',  lambda, Rsum(:,2)+1, 'r', 'LineWidth', 2)

%h5 = legend('theta=0','theta=1','theta=2','theta=3',4);

%title('W=',w*10^9,' nm, D=', period*10^9,' nm')
axis tight
axis([lmin lmax 0 2.5])
set(gca,'fontsize', 18)
%h5 = legend('theta=0','theta=0.2','theta=0.4','theta=0.6',4);
%set(h5,'Interpreter','none')
hold off
%}
%plot(lambda, Rsum)
%{
Collist = 'ybmgr'
col = randi([1 5])
figure(2)
hold on
plot(lambda, Rsum+0.5, Collist(col))
hold off
%}
%plot(lambda, Rsum)
%plot(lambda, Rsum(:,1,1), 'b', lambda, Rsum(:,2,1), 'g', lambda, Rsum(:,3,1), 'r');

%{
k0=2*pi/lambda;
kx = sin(theta')*cos(phi);
ky = sin(theta')*sin(phi);
figure;
pcolor(kx,ky,Rsum)
%shading interp
colormap gray;
%}
%pcolor(theta*180/pi,lambda,Rsum);
%imagesc(theta*180/pi,lambda,Rsum);
%set(gca,'Yscale','linear','Ydir','normal');


%plot([0:1:TH-1],Rsum,'b',[0:1:TH-1],Tsum,'r')

%num = zeros(2*N+1,1);
%for q=1:(2*N+1)
%   num(q) = q-N-1;
%end
%bar(num,eta,'stack')
%bar(num, eta_T1)