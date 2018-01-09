function [eps_total, mu_total] =...
    PMM_epsilon_rectangle(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
    N_intervals_x,N_intervals_y,epsilon)

N_total_x = sum(N_basis_x);
N_total_y = sum(N_basis_y);

N_total_x3 = N_total_x - N_intervals_x; %in "third" basis
N_total_y3 = N_total_y - N_intervals_y;

unity_x = zeros(N_total_x3, N_total_x3, N_intervals_x);
unity_y = zeros(N_total_y3, N_total_y3, N_intervals_y);

for j=1:N_intervals_x
    for i=(Nx(j)+1):(Nx(j)+nx(j))
    unity_x(i,i,j) = 1;
    end
end
for j=1:N_intervals_y
    for i=(Ny(j)+1):(Ny(j)+ny(j))
    unity_y(i,i,j) = 1;
    end
end

epsilon_inv = epsilon.^-1;
mmu = ones(N_intervals_x,N_intervals_y);
mmu_inv = mmu.^-1;

N_total3 = N_total_x3*N_total_y3;
epsm = zeros(N_total3, N_total3);
eps_inv = zeros(N_total3, N_total3);
mu  = zeros(N_total3, N_total3);
mu_inv = zeros(N_total3, N_total3);

for l = 1:N_intervals_x
    for q = 1:N_intervals_y
%eps = eps + epsilon(l,q)*Kronecker_product(hx\unity_x(:,:,l), hy\unity_y(:,:,q));
%mu  = mu  + Kronecker_product(hx\unity_x(:,:,l), hy\unity_y(:,:,q));
%!!!!!!!we will need the above expression for eps,mu when we introduce
%matched coordinates and sqrt(g)*(g)ii
Ilq = Kronecker_product(unity_x(:,:,l), unity_y(:,:,q));
epsm = epsm + epsilon(l,q)*Ilq;
eps_inv = eps_inv + epsilon_inv(l,q)*Ilq;
mu  = mu  + mmu(l,q)*Ilq;
mu_inv = mu_inv  + mmu_inv(l,q)*Ilq;
    end
end

%for now
g_g11 = 1;
g_g12 = 0;
g_g21 = 0;
g_g22 = 1;
g_g33 = 1;

eps11 = epsm*g_g11;
eps12 = epsm*g_g12;
eps21 = epsm*g_g21;
eps22 = epsm*g_g22;
eps_inv33 = eps_inv*g_g33;

eps_total(:,:,1) = eps11;
eps_total(:,:,2) = eps12;
eps_total(:,:,3) = eps21;
eps_total(:,:,4) = eps22;
eps_total(:,:,5) = eps_inv33;

mu11 = mu*g_g11;
mu12 = mu*g_g12;
mu21 = mu*g_g21;
mu22 = mu*g_g22;
mu_inv33 = mu_inv*g_g33;

mu_total(:,:,1)=mu11;
mu_total(:,:,2)=mu12;
mu_total(:,:,3)=mu21;
mu_total(:,:,4)=mu22;
mu_total(:,:,5)=mu_inv33;

%{
eps_inv11 = eps_inv*g_g11;
eps_inv12 = eps_inv*g_g12;
eps_inv21 = eps_inv*g_g21;
eps_inv22 = eps_inv*g_g22;

eps_inv_total(:,:,1) = eps_inv11;
eps_inv_total(:,:,2) = eps_inv12;
eps_inv_total(:,:,3) = eps_inv21;
eps_inv_total(:,:,4) = eps_inv22;
eps_inv_total(:,:,5) = eps_inv33;

mu_inv11 = mu_inv*g_g11;
mu_inv12 = mu_inv*g_g12;
mu_inv21 = mu_inv*g_g21;
mu_inv22 = mu_inv*g_g22;

mu_inv_total(:,:,1)=mu_inv11;
mu_inv_total(:,:,2)=mu_inv12;
mu_inv_total(:,:,3)=mu_inv21;
mu_inv_total(:,:,4)=mu_inv22;
mu_inv_total(:,:,5)=mu_inv33;
%}