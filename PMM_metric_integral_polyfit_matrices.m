
function [int_Ez_sqrt_g_full,int_Dz_unity_full,int_Dx_sqrt_g_full,int_Dy_sqrt_g_full,...
    int_Ex_g_down22_full,int_Ey_g_down12_full,int_Ex_g_down21_full,int_Ey_g_down11_full] =...
    PMM_metric_integral_polyfit_matrices(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
    N_intervals_x,N_intervals_y,n_points,La,ax,ay,hx,hy,dx_x1,dx_x2,dy_x1,dy_x2,uni,b_x1,b_x2)
    
g_down11 = dx_x1.^2 + dy_x1.^2;
g_down12 = dx_x1.*dx_x2 + dy_x1.*dy_x2;
g_down21 = g_down12;
g_down22 = dx_x2.^2 + dy_x2.^2;
g_det = g_down11.*g_down22 - g_down12.*g_down21;

g_sqrt = g_det.^0.5;

ksi1 = linspace(-uni,uni,n_points);
ksi2 = linspace(-uni,uni,n_points);

g_sqrt_12_raw(:,:) = g_sqrt(1,2,:,:);
g_sqrt_21_raw(:,:) = g_sqrt(2,1,:,:);
g_sqrt_22_raw(:,:) = g_sqrt(2,2,:,:);
g_sqrt_23_raw(:,:) = g_sqrt(2,3,:,:);
g_sqrt_32_raw(:,:) = g_sqrt(3,2,:,:);

g_down11_12_raw(:,:) = g_down11(1,2,:,:);
g_down11_21_raw(:,:) = g_down11(2,1,:,:);
g_down11_22_raw(:,:) = g_down11(2,2,:,:);
g_down11_23_raw(:,:) = g_down11(2,3,:,:);
g_down11_32_raw(:,:) = g_down11(3,2,:,:);

g_down22_12_raw(:,:) = g_down22(1,2,:,:);
g_down22_21_raw(:,:) = g_down22(2,1,:,:);
g_down22_22_raw(:,:) = g_down22(2,2,:,:);
g_down22_23_raw(:,:) = g_down22(2,3,:,:);
g_down22_32_raw(:,:) = g_down22(3,2,:,:);

g_down12_12_raw(:,:) = g_down12(1,2,:,:);
g_down12_21_raw(:,:) = g_down12(2,1,:,:);
g_down12_22_raw(:,:) = g_down12(2,2,:,:);
g_down12_23_raw(:,:) = g_down12(2,3,:,:);
g_down12_32_raw(:,:) = g_down12(3,2,:,:);

[g_sqrt_12_fit] = polyfit_procedure(ksi1,ksi2,g_sqrt_12_raw);
[g_sqrt_21_fit] = polyfit_procedure(ksi1,ksi2,g_sqrt_21_raw);
[g_sqrt_22_fit] = polyfit_procedure(ksi1,ksi2,g_sqrt_22_raw);
[g_sqrt_23_fit] = polyfit_procedure(ksi1,ksi2,g_sqrt_23_raw);
[g_sqrt_32_fit] = polyfit_procedure(ksi1,ksi2,g_sqrt_32_raw);

[g_down11_12_fit] = polyfit_procedure(ksi1,ksi2,g_down11_12_raw);
[g_down11_21_fit] = polyfit_procedure(ksi1,ksi2,g_down11_21_raw);
[g_down11_22_fit] = polyfit_procedure(ksi1,ksi2,g_down11_22_raw);
[g_down11_23_fit] = polyfit_procedure(ksi1,ksi2,g_down11_23_raw);
[g_down11_32_fit] = polyfit_procedure(ksi1,ksi2,g_down11_32_raw);

[g_down22_12_fit] = polyfit_procedure(ksi1,ksi2,g_down22_12_raw);
[g_down22_21_fit] = polyfit_procedure(ksi1,ksi2,g_down22_21_raw);
[g_down22_22_fit] = polyfit_procedure(ksi1,ksi2,g_down22_22_raw);
[g_down22_23_fit] = polyfit_procedure(ksi1,ksi2,g_down22_23_raw);
[g_down22_32_fit] = polyfit_procedure(ksi1,ksi2,g_down22_32_raw);

[g_down12_12_fit] = polyfit_procedure(ksi1,ksi2,g_down12_12_raw);
[g_down12_21_fit] = polyfit_procedure(ksi1,ksi2,g_down12_21_raw);
[g_down12_22_fit] = polyfit_procedure(ksi1,ksi2,g_down12_22_raw);
[g_down12_23_fit] = polyfit_procedure(ksi1,ksi2,g_down12_23_raw);
[g_down12_32_fit] = polyfit_procedure(ksi1,ksi2,g_down12_32_raw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Tx, aTx, Ux] = PMM_T_matrices(Nx,nx,N_basis_x,ax,La,N_intervals_x,b_x1);
[Ty, aTy, Uy] = PMM_T_matrices(Ny,ny,N_basis_y,ay,La,N_intervals_y,b_x2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_total_x = sum(N_basis_x);  %total number of basis functions
N_total_y = sum(N_basis_y);  %total number of basis functions
N_total_x3 = N_total_x - N_intervals_x;  %number of basis functions in "third" basis
N_total_y3 = N_total_y - N_intervals_y;
N_total_3 = N_total_x3*N_total_y3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ux1(:,:) = Ux(1,:,:);
Ux2(:,:) = Ux(2,:,:);
Ux3(:,:) = Ux(3,:,:);

Uy1(:,:) = Uy(1,:,:);
Uy2(:,:) = Uy(2,:,:);
Uy3(:,:) = Uy(3,:,:);

Tx1(:,:) = Tx(1,:,:);
Tx2(:,:) = Tx(2,:,:);
Tx3(:,:) = Tx(3,:,:);

Ty1(:,:) = Ty(1,:,:);
Ty2(:,:) = Ty(2,:,:);
Ty3(:,:) = Ty(3,:,:);

aTx1(:,:) = aTx(1,:,:);
aTx2(:,:) = aTx(2,:,:);
aTx3(:,:) = aTx(3,:,:);

aTy1(:,:) = aTy(1,:,:);
aTy2(:,:) = aTy(2,:,:);
aTy3(:,:) = aTy(3,:,:);


%%%%%%%%%%%%%%%%%%%%%%
int_Ez_sqrt_g = zeros(N_intervals_x,N_intervals_y,N_total_3,N_total_3);

int_Ez_sqrt_g(1,1,:,:) = Kr(Ux1,Uy1);
int_Ez_sqrt_g(1,3,:,:) = Kr(Ux1,Uy3);
int_Ez_sqrt_g(3,1,:,:) = Kr(Ux3,Uy1);
int_Ez_sqrt_g(3,3,:,:) = Kr(Ux3,Uy3);

int_Ez_sqrt_g(1,2,:,:) = int_P2_Q2(Ux1,Uy2,Tx1,aTx1,Ty2,aTy2,g_sqrt_12_fit);
int_Ez_sqrt_g(2,1,:,:) = int_P2_Q2(Ux2,Uy1,Tx2,aTx2,Ty1,aTy1,g_sqrt_21_fit);
int_Ez_sqrt_g(2,2,:,:) = int_P2_Q2(Ux2,Uy2,Tx2,aTx2,Ty2,aTy2,g_sqrt_22_fit);
int_Ez_sqrt_g(2,3,:,:) = int_P2_Q2(Ux2,Uy3,Tx2,aTx2,Ty3,aTy3,g_sqrt_23_fit);
int_Ez_sqrt_g(3,2,:,:) = int_P2_Q2(Ux3,Uy2,Tx3,aTx3,Ty2,aTy2,g_sqrt_32_fit);

int_Ez_sqrt_g_full = zeros(N_total_3,N_total_3);
for i=1:N_intervals_x
    for j=1:N_intervals_y
        int_Ez_sqrt_g_ij(:,:) = int_Ez_sqrt_g(i,j,:,:);
        int_Ez_sqrt_g_full(:,:) = int_Ez_sqrt_g_full(:,:)+int_Ez_sqrt_g_ij;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%{
hx = zeros(N_total_x3,N_total_x3);
hy = zeros(N_total_y3,N_total_y3);
lx=zeros(N_intervals_x,1);
suml = zeros(N_intervals_x,1);
for k=1:N_intervals_x
    lx(k) = (b_x(k+1) - b_x(k))/2;
    suml(k) = (b_x(k+1) + b_x(k))/2;
end
lx=zeros(N_intervals_x,1);
suml = zeros(N_intervals_x,1);
for k=1:N_intervals_x
    lx(k) = (b_x(k+1) - b_x(k))/2;
    suml(k) = (b_x(k+1) + b_x(k))/2;
end
for k=1:N_intervals_x
    for i=(Nx(k)+1):(Nx(k)+nx(k))
    hx(i,i) = norm(i-Nx(k))*lx(k); 
    end
end
for k=1:N_intervals_y
    for i=(Ny(k)+1):(Ny(k)+ny(k))
    hy(i,i) = norm(i-Ny(k)); 
    end
end
%}
int_Dz_unity_full=eye(N_total_3,N_total_3);%Kr(hx,hy); for PMM_T_matrices

%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_Dx_sqrt_g = zeros(N_intervals_x,N_intervals_y,N_total_3,N_total_3);
int_Dx_sqrt_g(1,1,:,:) = Kr(Ux1,Uy1);
int_Dx_sqrt_g(1,3,:,:) = Kr(Ux1,Uy3);
int_Dx_sqrt_g(3,1,:,:) = Kr(Ux3,Uy1);
int_Dx_sqrt_g(3,3,:,:) = Kr(Ux3,Uy3);

int_Dx_sqrt_g(1,2,:,:) = int_P2_Q3(Ux1,Uy2,Tx1,aTx1,Ty2,g_sqrt_12_fit);
int_Dx_sqrt_g(2,1,:,:) = int_P2_Q3(Ux2,Uy1,Tx2,aTx2,Ty1,g_sqrt_21_fit);
int_Dx_sqrt_g(2,2,:,:) = int_P2_Q3(Ux2,Uy2,Tx2,aTx2,Ty2,g_sqrt_22_fit);
int_Dx_sqrt_g(2,3,:,:) = int_P2_Q3(Ux2,Uy3,Tx2,aTx2,Ty3,g_sqrt_23_fit);
int_Dx_sqrt_g(3,2,:,:) = int_P2_Q3(Ux3,Uy2,Tx3,aTx3,Ty2,g_sqrt_32_fit);
int_Dx_sqrt_g_full = zeros(N_total_3,N_total_3);
for i=1:N_intervals_x
    for j=1:N_intervals_y
        int_Dx_sqrt_g_ij(:,:) = int_Dx_sqrt_g(i,j,:,:);
        int_Dx_sqrt_g_full(:,:) = int_Dx_sqrt_g_full+int_Dx_sqrt_g_ij;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_Dy_sqrt_g = zeros(N_intervals_x,N_intervals_y,N_total_3,N_total_3);
int_Dy_sqrt_g(1,1,:,:) = Kr(Ux1,Uy1);
int_Dy_sqrt_g(1,3,:,:) = Kr(Ux1,Uy3);
int_Dy_sqrt_g(3,1,:,:) = Kr(Ux3,Uy1);
int_Dy_sqrt_g(3,3,:,:) = Kr(Ux3,Uy3);

int_Dy_sqrt_g(1,2,:,:) = int_P3_Q2(Ux1,Uy2,Tx1,aTy2,Ty2,g_sqrt_12_fit);
int_Dy_sqrt_g(2,1,:,:) = int_P3_Q2(Ux2,Uy1,Tx2,aTy1,Ty1,g_sqrt_21_fit);
int_Dy_sqrt_g(2,2,:,:) = int_P3_Q2(Ux2,Uy2,Tx2,aTy2,Ty2,g_sqrt_22_fit);
int_Dy_sqrt_g(2,3,:,:) = int_P3_Q2(Ux2,Uy3,Tx2,aTy3,Ty3,g_sqrt_23_fit);
int_Dy_sqrt_g(3,2,:,:) = int_P3_Q2(Ux3,Uy2,Tx3,aTy2,Ty2,g_sqrt_32_fit);
int_Dy_sqrt_g_full = zeros(N_total_3,N_total_3);
for i=1:N_intervals_x
    for j=1:N_intervals_y
        int_Dy_sqrt_g_ij(:,:) = int_Dy_sqrt_g(i,j,:,:);
        int_Dy_sqrt_g_full(:,:) = int_Dy_sqrt_g_full+int_Dy_sqrt_g_ij;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_Ex_g_down22 = zeros(N_intervals_x,N_intervals_y,N_total_3,N_total_3);
int_Ex_g_down22(1,1,:,:) = Kr(Ux1,Uy1);
int_Ex_g_down22(1,3,:,:) = Kr(Ux1,Uy3);
int_Ex_g_down22(3,1,:,:) = Kr(Ux3,Uy1);
int_Ex_g_down22(3,3,:,:) = Kr(Ux3,Uy3);

int_Ex_g_down22(1,2,:,:) = int_P3_Q2(Ux1,Uy2,Tx1,aTy2,Ty2,g_down22_12_fit);
int_Ex_g_down22(2,1,:,:) = int_P3_Q2(Ux2,Uy1,Tx2,aTy1,Ty1,g_down22_21_fit);
int_Ex_g_down22(2,2,:,:) = int_P3_Q2(Ux2,Uy2,Tx2,aTy2,Ty2,g_down22_22_fit);
int_Ex_g_down22(2,3,:,:) = int_P3_Q2(Ux2,Uy3,Tx2,aTy3,Ty3,g_down22_23_fit);
int_Ex_g_down22(3,2,:,:) = int_P3_Q2(Ux3,Uy2,Tx3,aTy2,Ty2,g_down22_32_fit);
int_Ex_g_down22_full = zeros(N_total_3,N_total_3);
for i=1:N_intervals_x
    for j=1:N_intervals_y
        int_Ex_g_down22_ij(:,:) = int_Ex_g_down22(i,j,:,:);
        int_Ex_g_down22_full(:,:) = int_Ex_g_down22_full+int_Ex_g_down22_ij;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_Ey_g_down12 = zeros(N_intervals_x,N_intervals_y,N_total_3,N_total_3);
int_Ey_g_down12(1,1,:,:) = zeros(N_total_3,N_total_3);
int_Ey_g_down12(1,3,:,:) = zeros(N_total_3,N_total_3);
int_Ey_g_down12(3,1,:,:) = zeros(N_total_3,N_total_3);
int_Ey_g_down12(3,3,:,:) = zeros(N_total_3,N_total_3);

int_Ey_g_down12(1,2,:,:) = int_P2_Q3(Ux1,Uy2,Tx1,aTx1,Ty2,g_down12_12_fit);
int_Ey_g_down12(2,1,:,:) = int_P2_Q3(Ux2,Uy1,Tx2,aTx2,Ty1,g_down12_21_fit);
int_Ey_g_down12(2,2,:,:) = int_P2_Q3(Ux2,Uy2,Tx2,aTx2,Ty2,g_down12_22_fit);
int_Ey_g_down12(2,3,:,:) = int_P2_Q3(Ux2,Uy3,Tx2,aTx2,Ty3,g_down12_23_fit);
int_Ey_g_down12(3,2,:,:) = int_P2_Q3(Ux3,Uy2,Tx3,aTx3,Ty2,g_down12_32_fit);
int_Ey_g_down12_full = zeros(N_total_3,N_total_3);
for i=1:N_intervals_x
    for j=1:N_intervals_y
        int_Ey_g_down12_ij(:,:) = int_Ey_g_down12(i,j,:,:);
        int_Ey_g_down12_full(:,:) = int_Ey_g_down12_full+int_Ey_g_down12_ij;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_Ex_g_down21 = zeros(N_intervals_x,N_intervals_y,N_total_3,N_total_3);
int_Ex_g_down21(1,1,:,:) = zeros(N_total_3,N_total_3);
int_Ex_g_down21(1,3,:,:) = zeros(N_total_3,N_total_3);
int_Ex_g_down21(3,1,:,:) = zeros(N_total_3,N_total_3);
int_Ex_g_down21(3,3,:,:) = zeros(N_total_3,N_total_3);

int_Ex_g_down21(1,2,:,:) = int_P3_Q2(Ux1,Uy2,Tx1,aTy2,Ty2,g_down12_12_fit);
int_Ex_g_down21(2,1,:,:) = int_P3_Q2(Ux2,Uy1,Tx2,aTy1,Ty1,g_down12_21_fit);
int_Ex_g_down21(2,2,:,:) = int_P3_Q2(Ux2,Uy2,Tx2,aTy2,Ty2,g_down12_22_fit);
int_Ex_g_down21(2,3,:,:) = int_P3_Q2(Ux2,Uy3,Tx2,aTy3,Ty3,g_down12_23_fit);
int_Ex_g_down21(3,2,:,:) = int_P3_Q2(Ux3,Uy2,Tx3,aTy2,Ty2,g_down12_32_fit);
int_Ex_g_down21_full = zeros(N_total_3,N_total_3);
for i=1:N_intervals_x
    for j=1:N_intervals_y
        int_Ex_g_down21_ij(:,:) = int_Ex_g_down21(i,j,:,:);
        int_Ex_g_down21_full(:,:) = int_Ex_g_down21_full+int_Ex_g_down21_ij;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_Ey_g_down11 = zeros(N_intervals_x,N_intervals_y,N_total_3,N_total_3);
int_Ey_g_down11(1,1,:,:) = Kr(Ux1,Uy1);
int_Ey_g_down11(1,3,:,:) = Kr(Ux1,Uy3);
int_Ey_g_down11(3,1,:,:) = Kr(Ux3,Uy1);
int_Ey_g_down11(3,3,:,:) = Kr(Ux3,Uy3);

int_Ey_g_down11(1,2,:,:) = int_P2_Q3(Ux1,Uy2,Tx1,aTx1,Ty2,g_down11_12_fit);
int_Ey_g_down11(2,1,:,:) = int_P2_Q3(Ux2,Uy1,Tx2,aTx2,Ty1,g_down11_21_fit);
int_Ey_g_down11(2,2,:,:) = int_P2_Q3(Ux2,Uy2,Tx2,aTx2,Ty2,g_down11_22_fit);
int_Ey_g_down11(2,3,:,:) = int_P2_Q3(Ux2,Uy3,Tx2,aTx2,Ty3,g_down11_23_fit);
int_Ey_g_down11(3,2,:,:) = int_P2_Q3(Ux3,Uy2,Tx3,aTx3,Ty2,g_down11_32_fit);
int_Ey_g_down11_full = zeros(N_total_3,N_total_3);
for i=1:N_intervals_x
    for j=1:N_intervals_y
        int_Ey_g_down11_ij(:,:) = int_Ey_g_down11(i,j,:,:);
        int_Ey_g_down11_full(:,:) = int_Ey_g_down11_full+int_Ey_g_down11_ij;
    end
end
