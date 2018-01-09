function [eps_total, mu_total] =...
    PMM_epsilon_ellipse(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
    N_intervals_x,N_intervals_y,La,epsilon,...
    int_eps_xx_full, int_eps_xy_full, int_eps_yx_full, int_eps_yy_full)


N_total_x = sum(N_basis_x);  %total number of basis functions
N_total_y = sum(N_basis_y);  %total number of basis functions

N_total_x3 = N_total_x - N_intervals_x;
N_total_y3 = N_total_y - N_intervals_y;
N_total_3 = N_total_x3*N_total_y3;

Nmax_x = max(N_basis_x);   
Nmax_y = max(N_basis_y);   
Nmax = max(Nmax_x, Nmax_y);

p = zeros(Nmax,1);     
norm = zeros(Nmax,1);
for i=0:(Nmax-1)
    p(i+1) = gamma(i+2*La)/(gamma(2*La)*gamma(i+1));
    %p(i)=Ci i-th Gegenbauer polynomial at 1
    norm(i+1) = pi^0.5*p(i+1)*gamma(La+0.5)/(gamma(La)*(i+La));
    %<Cn,Cm> = delta(n,m)*norm(n)
end
hx = zeros(N_total_x3);
hy = zeros(N_total_y3);
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

eps_xx = zeros(N_total_3,N_total_3);
eps_xy = zeros(N_total_3,N_total_3);
eps_yx = zeros(N_total_3,N_total_3);
eps_yy = zeros(N_total_3,N_total_3);
mu_xx = zeros(N_total_3,N_total_3);
mu_xy = zeros(N_total_3,N_total_3);
mu_yx = zeros(N_total_3,N_total_3);
mu_yy = zeros(N_total_3,N_total_3);

epsilon_inv = epsilon.^-1;
mmu = ones(N_intervals_x,N_intervals_y);
mmu_inv = mmu.^-1;

for j1=1:N_intervals_x
for j2=1:N_intervals_y
    for m1=Nx(j1)+1:Nx(j1)+nx(j1)
    for m2=Ny(j2)+1:Ny(j2)+ny(j2)
        row = m2 + (m1-1)*N_total_y3;        
        for k2=1:N_intervals_y
            for i1 = Nx(j1)+1:Nx(j1)+nx(j1)
            for i2 = Ny(k2)+1:Ny(k2)+ny(k2)                
            col = i2 + (i1-1)*N_total_y3;
    eps_xx(row,col)=epsilon(j1,j2)*int_eps_xx_full(j1,j2,m1,m2,i1,i2)/(hx(m1,m1)*hy(m2,m2));
    eps_yx(row,col)=epsilon(j1,j2)*int_eps_yx_full(j1,j2,m1,m2,i1,i2)/(hx(m1,m1)*hy(m2,m2));
    mu_xx(row,col) =mmu(j1,j2)*int_eps_xx_full(j1,j2,m1,m2,i1,i2)/(hx(m1,m1)*hy(m2,m2));
    mu_yx(row,col) =mmu(j1,j2)*int_eps_yx_full(j1,j2,m1,m2,i1,i2)/(hx(m1,m1)*hy(m2,m2));
            end
            end
        end
        for k1=1:N_intervals_x
            for i1 = Nx(k1)+1:Nx(k1)+nx(k1)
            for i2 = Ny(j2)+1:Ny(j2)+ny(j2)
            col = i2 + (i1-1)*N_total_y3;
    eps_xy(row,col)=epsilon(j1,j2)*int_eps_xy_full(j1,j2,m1,m2,i1,i2)/(hx(m1,m1)*hy(m2,m2));
    eps_yy(row,col)=epsilon(j1,j2)*int_eps_yy_full(j1,j2,m1,m2,i1,i2)/(hx(m1,m1)*hy(m2,m2)); 
    mu_xy(row,col) =mmu(j1,j2)*int_eps_xy_full(j1,j2,m1,m2,i1,i2)/(hx(m1,m1)*hy(m2,m2));
    mu_yy(row,col) =mmu(j1,j2)*int_eps_yy_full(j1,j2,m1,m2,i1,i2)/(hx(m1,m1)*hy(m2,m2));
            end
            end
        end
    end
    end
end
end

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


N_total3 = N_total_x3*N_total_y3;
eps_inv = zeros(N_total3, N_total3);
mu_inv = zeros(N_total3, N_total3);

for l = 1:N_intervals_x
    for q = 1:N_intervals_y
Ilq = Kronecker_product(unity_x(:,:,l), unity_y(:,:,q));
%epsm = epsm + epsilon(l,q)*Ilq;
eps_inv = eps_inv + epsilon_inv(l,q)*Ilq;
%mu  = mu  + mmu(l,q)*Ilq;
mu_inv = mu_inv  + mmu_inv(l,q)*Ilq;
    end
end
    
eps_inv33 = eps_inv;
mu_inv33 = mu_inv;

eps_total(:,:,1) = eps_xx;
eps_total(:,:,2) = eps_xy;
eps_total(:,:,3) = eps_yx;
eps_total(:,:,4) = eps_yy;
eps_total(:,:,5) = eps_inv33;

mu_total(:,:,1)=mu_xx;
mu_total(:,:,2)=mu_xy;
mu_total(:,:,3)=mu_yx;
mu_total(:,:,4)=mu_yy;
mu_total(:,:,5)=mu_inv33;

