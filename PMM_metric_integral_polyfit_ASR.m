function [eps_total,mu_total] =...
    PMM_metric_integral_polyfit_ASR(b_x,eta,f1,N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
    N_intervals_x,N_intervals_y,n_points,La,ax,ay,epsilon)

w1 = b_x(2);
d1 = b_x(3);
%{
x1_1 = linspace(0,f1,n_points);
x1_2 = linspace(f1,1,n_points);
ksi = linspace(-1,1,n_points);
dx_dx1_1 = (w1/f1)*(1-eta*cos(2*pi*x1_1/f1));
dx_dx1_2 = ((d1-w1)/(1-f1))*( 1-eta*cos(2*pi*(x1_2-f1)/(1-f1)) );
%}
x1_1 = linspace(0,w1,n_points);
x1_2 = linspace(w1,d1,n_points);
ksi = linspace(-1,1,n_points);
x_1 = x1_1 - eta*(w1/(2*pi))*sin(2*pi*x1_1/w1);
x_2 = x1_2 - eta*((d1-w1)/(2*pi))*sin(2*pi*(x1_2-w1)/(d1-w1));

figure(7)
plot (x1_1, x_1,'m', 'Linewidth', 2)
xlabel('x');
ylabel('x1');
hold off

dx_dx1_1 = 1 - eta*cos(2*pi*x1_1/w1);
dx_dx1_2 = 1 - eta*cos(2*pi*(x1_2-w1)/(d1-w1));


g_sqrt_1 = dx_dx1_1;  %on the first interval  [0,f1]
g_sqrt_2 = dx_dx1_2;  %on the second interval [f1,1]
ksi = transpose(ksi);
g_sqrt_1 = transpose(g_sqrt_1);
g_sqrt_2 = transpose(g_sqrt_2);
g_sqrt_1_fit = fit(ksi,g_sqrt_1,'poly7')
g_sqrt_2_fit = fit(ksi,g_sqrt_2,'poly7');
g_sqrt_1_fit_num = g_sqrt_1_fit.p8+g_sqrt_1_fit.p7*ksi+g_sqrt_1_fit.p6*ksi.^2+...
    +g_sqrt_1_fit.p5*ksi.^3+g_sqrt_1_fit.p4*ksi.^4+g_sqrt_1_fit.p3*ksi.^5+...
    g_sqrt_1_fit.p2*ksi.^6+g_sqrt_1_fit.p1*ksi.^7;
%{
g_sqrt_1_fit = polyfit(ksi,g_sqrt_1,5);
g_sqrt_2_fit = polyfit(ksi,g_sqrt_2,5);
g_sqrt_1_fit_num = g_sqrt_1_fit(1)+g_sqrt_1_fit(2)*ksi+g_sqrt_1_fit(3)*ksi.^2+...
    +g_sqrt_1_fit(4)*ksi.^3+g_sqrt_1_fit(5)*ksi.^4+g_sqrt_1_fit(6)*ksi.^5;
g_sqrt_1_polyval = polyval(g_sqrt_1_fit,ksi);
g_sqrt_2_polyval = polyval(g_sqrt_2_fit,ksi);
    figure(11)
plot(ksi,g_sqrt_2,'r',ksi,g_sqrt_2_polyval,'g',...
    ksi,abs(g_sqrt_2-g_sqrt_2_polyval)./g_sqrt_2, 'm', 'Linewidth', 2)
hold off
%}
figure(10)
plot(ksi,g_sqrt_1,'r',ksi,g_sqrt_1_fit_num,'b',...
    ksi,abs(g_sqrt_1-g_sqrt_1_fit_num)./g_sqrt_1, 'm', 'Linewidth', 2)
xlabel('ksi');
ylabel('g sqrt 1, g sqrt 1 fit');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Tx, aTx, Ux] = PMM_T_matrices(Nx,nx,N_basis_x,ax,La,N_intervals_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_total_x = sum(N_basis_x);  %total number of basis functions
N_total_y = sum(N_basis_y);  %total number of basis functions
N_total_x3 = N_total_x - N_intervals_x;  %number of basis functions in "third" basis
N_total_y3 = N_total_y - N_intervals_y;
N_total_3 = N_total_x3*N_total_y3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ux1(:,:) = Ux(1,:,:);
Ux2(:,:) = Ux(2,:,:);

Tx1(:,:) = Tx(1,:,:);
Tx2(:,:) = Tx(2,:,:);

aTx1(:,:) = aTx(1,:,:);
aTx2(:,:) = aTx(2,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_P2_sqrt_g_1 = int_P2(Ux1,Tx1,aTx1,g_sqrt_1_fit);
int_P2_sqrt_g_2 = int_P2(Ux2,Tx2,aTx2,g_sqrt_2_fit);

int_P2_sqrt_g = int_P2_sqrt_g_1+int_P2_sqrt_g_2;

Nmax_x = max(N_basis_x);
Nmax_y = max(N_basis_y);
Nmax = max(Nmax_x,Nmax_y);
p = zeros(Nmax,1);     
norm = zeros(Nmax,1);
for i=0:(Nmax-1)
    p(i+1) = gamma(i+2*La)/(gamma(2*La)*gamma(i+1));
    %p(i)=Ci i-th Gegenbauer polynomial at 1
    norm(i+1) = pi^0.5*p(i+1)*gamma(La+0.5)/(gamma(La)*(i+La));
    %<Cn,Cm> = delta(n,m)*norm(n)
end
hx = zeros(N_total_x3,N_total_x3);
hy = zeros(N_total_y3,N_total_y3);
for k=1:N_intervals_x
    for i=(Nx(k)+1):(Nx(k)+nx(k))
    hx(i,i) = norm(i-Nx(k)); 
    end
end
for k=1:N_intervals_y
    for i=(Ny(k)+1):(Ny(k)+ny(k))
    hy(i,i) = norm(i-Ny(k)); 
    end
end
int_sqrt_g = Kr(int_P2_sqrt_g,hy);
%int_sqrt_g = Kr(hx,hy)\Kr(int_P2_sqrt_g,hy);
diag_g = diag(int_sqrt_g);

epsilon_inv = ones(N_intervals_x,N_intervals_y)./epsilon;
mu = ones(N_intervals_x,N_intervals_y);
mu_inv = ones(N_intervals_x,N_intervals_y);

eps_xx_inv = zeros(N_total_x3,N_total_x3);
eps_yy = zeros(N_total_x3,N_total_x3);
eps_zz = zeros(N_total_x3,N_total_x3);
mu_xx_inv = zeros(N_total_x3,N_total_x3);
mu_yy = zeros(N_total_x3,N_total_x3);
mu_zz = zeros(N_total_x3,N_total_x3);

for k1=1:N_intervals_x
for k2=1:N_intervals_y
    for m1=Nx(k1)+1:Nx(k1)+nx(k1)
    for m2=Ny(k2)+1:Ny(k2)+ny(k2)
        row = m2 + (m1-1)*N_total_y3;
        for i1=1:N_total_x3
        for i2=1:N_total_y3
        
        col = i2 + (i1-1)*N_total_y3;
        eps_xx_inv(row,col)=epsilon_inv(k1,k2)*int_sqrt_g(row,col);
        eps_yy(row,col)=epsilon(k1,k2)*int_sqrt_g(row,col);
        eps_zz(row,col)=epsilon(k1,k2)*int_sqrt_g(row,col);
        
        mu_xx_inv(row,col)=mu_inv(k1,k2)*int_sqrt_g(row,col);
        mu_yy(row,col)=mu(k1,k2)*int_sqrt_g(row,col);
        mu_zz(row,col)=mu(k1,k2)*int_sqrt_g(row,col);
        end
        end
    end
    end
end
end

Uxy=eye(N_total_3,N_total_3);

eps_xy = zeros(N_total_3,N_total_3);
eps_yx = zeros(N_total_3,N_total_3);
mu_xy = zeros(N_total_3,N_total_3);
mu_yx = zeros(N_total_3,N_total_3);

eps_total(:,:,1) = Uxy/eps_xx_inv;
eps_total(:,:,2) = eps_xy;
eps_total(:,:,3) = eps_yx;
eps_total(:,:,4) = eps_yy;
eps_total(:,:,5) = Uxy/eps_zz

mu_total(:,:,1) = Uxy/mu_xx_inv;
mu_total(:,:,2) = mu_xy;
mu_total(:,:,3) = mu_yx;
mu_total(:,:,4) = mu_yy;
mu_total(:,:,5) = Uxy/mu_zz;
