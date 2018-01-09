function [eps_total, mu_total] =...
    PMM_epsilon_ellipse_new(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
    N_intervals_x,N_intervals_y,La,epsilon,...
    int_epsmu_xx,int_epsmu_xy,int_epsmu_yx,int_epsmu_yy,int_epsmu_inv_zz)


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
eps_out = epsilon(1);
eps_in =  epsilon(2);
mu_in = 1.0;
mu_out = 1.0;

eps_xx=eps_out*int_epsmu_xx;
eps_xy=eps_out*int_epsmu_xy;
eps_yx=eps_out*int_epsmu_yx;
eps_yy=eps_out*int_epsmu_yy;

mu_xx =mu_out*int_epsmu_xx;
mu_xy =mu_out*int_epsmu_xy;
mu_yx =mu_out*int_epsmu_yx;    
mu_yy =mu_out*int_epsmu_yy;
    
eps_inv_zz = (1/eps_out)*int_epsmu_inv_zz;
mu_inv_zz = (1/mu_out)*int_epsmu_inv_zz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
miden = eye(N_total_3,N_total_3);
int_epsmu_inv_zz_2 = miden/int_epsmu_zz;
eps_inv_zz_2 = (1/eps_out)*int_epsmu_inv_zz_2;
mu_inv_zz_2 = (1/mu_out)*int_epsmu_inv_zz_2;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k1 = 2;
k2 = 2;

if (eps_in~=eps_out)
    
for m1=Nx(k1)+1:Nx(k1)+nx(k1)
for m2=Ny(k2)+1:Ny(k2)+ny(k2)
    row = m2 + (m1-1)*N_total_y3;
    for i1=1:N_total_x3
    for i2=1:N_total_y3
        %%%%%%%%%%??????????????????????
    %for i1=Nx(k1)+1:Nx(k1)+nx(k1)
    %for i2=Ny(k2)+1:Ny(k2)+ny(k2)
        col = i2 + (i1-1)*N_total_y3;
        eps_xx(row,col)=eps_in*int_epsmu_xx(row,col);
        eps_xy(row,col)=eps_in*int_epsmu_xy(row,col);
        eps_yx(row,col)=eps_in*int_epsmu_yx(row,col);
        eps_yy(row,col)=eps_in*int_epsmu_yy(row,col);
        
        mu_xx(row,col) =mu_in*int_epsmu_xx(row,col);
        mu_xy(row,col) =mu_in*int_epsmu_xy(row,col);
        mu_yx(row,col) =mu_in*int_epsmu_yx(row,col);
        mu_yy(row,col) =mu_in*int_epsmu_yy(row,col);
        
        eps_inv_zz(row,col) = (1/eps_in)*int_epsmu_inv_zz(row,col);
        mu_inv_zz(row,col) = (1/mu_in)*int_epsmu_inv_zz(row,col);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        eps_inv_zz_2(row,col) = (1/eps_in)*int_epsmu_inv_zz_2(row,col);
        mu_inv_zz_2(row,col) = (1/mu_in)*int_epsmu_inv_zz_2(row,col);
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    end
end
end

end

eps_total(:,:,1) = eps_xx;
eps_total(:,:,2) = eps_xy;
eps_total(:,:,3) = eps_yx;
eps_total(:,:,4) = eps_yy;
eps_total(:,:,5) = eps_inv_zz;

mu_total(:,:,1)=mu_xx;
mu_total(:,:,2)=mu_xy;
mu_total(:,:,3)=mu_yx;
mu_total(:,:,4)=mu_yy;
mu_total(:,:,5)=mu_inv_zz;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
eps_total(:,:,5) = eps_inv_zz_2;
mu_total(:,:,5) = mu_inv_zz_2;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

