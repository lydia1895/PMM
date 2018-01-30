function [int_Ex_g_down_12_1D,int_Ey_g_down_12_1D,int_Ey_g_down_11_1D,int_Ex_g_down_22_1D] =...
    ellipse_coordinates_and_derivatives_for_1D_polyfit(ellipse_parameters,n_points,...
    La,Nx,nx,N_basis_x,ax,Ny,ny,N_basis_y,ay,N_total_3)



R1 = ellipse_parameters(1);
R2 = ellipse_parameters(2);
P1 = ellipse_parameters(3);
P2 = ellipse_parameters(4);
Q1 = ellipse_parameters(5);
Q2 = ellipse_parameters(6);

N_intervals_x = 3;
N_intervals_y = 3;

x1_t_plus =  P1/2+Q1;
x1_t_minus = P1/2-Q1;
x2_t_plus =  P2/2+Q2;
x2_t_minus = P2/2-Q2;

b_x1 = [0 x1_t_minus x1_t_plus P1];
b_x2 = [0 x2_t_minus x2_t_plus P2];
[Tx, aTx, Ux] = PMM_T_matrices(Nx,nx,N_basis_x,ax,La,N_intervals_x,b_x1);
[Ty, aTy, Uy] = PMM_T_matrices(Ny,ny,N_basis_y,ay,La,N_intervals_y,b_x2);

ksi1 = linspace(-1,1,n_points);
ksi2 = linspace(-1,1,n_points);
ksi = linspace(-1,1,n_points);
x1 = zeros(N_intervals_x,n_points);
x2 = zeros(N_intervals_y,n_points);
for k1 = 1:N_intervals_x
    x1(k1,:) = ( (b_x1(k1+1)-b_x1(k1))*ksi1+(b_x1(k1+1)+b_x1(k1)) )/2;
end
for k2 = 1:N_intervals_y     
    x2(k2,:) = ( (b_x2(k2+1)-b_x2(k2))*ksi2+(b_x2(k2+1)+b_x2(k2)) )/2;
end

x1_plus = zeros(3,n_points);
x1_minus = zeros(3,n_points);
x2_plus = zeros(3,n_points);
x2_minus = zeros(3,n_points);


x1_plus(1,:)  = x1_t_plus*ones(n_points,1);
x1_minus(1,:) = x1_t_minus*ones(n_points,1);

x1_plus(2,:)  = P1/2 + (R1/R2)*(R2^2-(x2(2,:)-P2/2).^2).^0.5;
x1_minus(2,:) = P1/2 - (R1/R2)*(R2^2-(x2(2,:)-P2/2).^2).^0.5;

x1_plus(3,:)  = x1_t_plus*ones(n_points,1);
x1_minus(3,:) = x1_t_minus*ones(n_points,1);

x2_plus(1,:)  = x2_t_plus*ones(n_points,1);
x2_minus(1,:) = x2_t_minus*ones(n_points,1);

x2_plus(2,:)  = P2/2 + (R2/R1)*(R1^2-(x1(2,:)-P1/2).^2).^0.5;
x2_minus(2,:) = P2/2 - (R2/R1)*(R1^2-(x1(2,:)-P1/2).^2).^0.5;

x2_plus(3,:)  = x2_t_plus*ones(n_points,1);
x2_minus(3,:) = x2_t_minus*ones(n_points,1);

dx1_plus(1,:) = zeros(n_points,1);
dx1_plus(2,:) = -(R1/R2)*(x2(2,:)-P2/2)./(R2^2-(x2(2,:)-P2/2).^2).^0.5;
dx1_plus(3,:) = zeros(n_points,1);

dx2_plus(1,:) = zeros(n_points,1);
dx2_plus(2,:) = -(R2/R1)*(x1(2,:)-P1/2)./(R1^2-(x1(2,:)-P1/2).^2).^0.5;
dx2_plus(3,:) = zeros(n_points,1);

dx1_minus(1,:) = zeros(n_points,1);
dx1_minus(2,:) = (R1/R2)*(x2(2,:)-P2/2)./(R2^2-(x2(2,:)-P2/2).^2).^0.5;
dx1_minus(3,:) = zeros(n_points,1);

dx2_minus(1,:) = zeros(n_points,1);
dx2_minus(2,:) = (R2/R1)*(x1(2,:)-P1/2)./(R1^2-(x1(2,:)-P1/2).^2).^0.5;
dx2_minus(3,:) = zeros(n_points,1);

dx_x1 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
dx_x2 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
dy_x1 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
dy_x2 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);

gy_x_x1 = zeros(3,3,n_points);
gx_y_x2 = zeros(3,3,n_points);
vy_x_x2 = zeros(3,3,n_points);
g1y_x_x2 = zeros(3,3,n_points);
g2y_x_x2 = zeros(3,3,n_points);

ux_x_x2 = zeros(3,3,n_points);
f1x_x_x2 = zeros(3,3,n_points);
f2x_x_x2 = zeros(3,3,n_points);

vx_y_x1 = zeros(3,3,n_points);
g1x_y_x1 = zeros(3,3,n_points);
g2x_y_x1 = zeros(3,3,n_points);

uy_y_x1 = zeros(3,3,n_points);
f1y_y_x1 = zeros(3,3,n_points);
f2y_y_x1 = zeros(3,3,n_points);
%%%%%%%%%%%%%%%%%%%%%%%%dx/dx1
for ny=1:n_points
gy_x_x1(1,:,ny) = x1_minus(:,ny)/x1_t_minus;
gy_x_x1(2,:,ny) = (x1_plus(:,ny)-x1_minus(:,ny))/(x1_t_plus-x1_t_minus);
gy_x_x1(3,:,ny) = (x1_plus(:,ny)-P1)/(x1_t_plus-P1);
end
%%%%%%%%%%%%%%%%%%%%%%%%dy/dx2
for nx=1:n_points
gx_y_x2(:,1,nx) = x2_minus(:,nx)/x2_t_minus;
gx_y_x2(:,2,nx) = (x2_plus(:,nx)-x2_minus(:,nx))/(x2_t_plus-x2_t_minus);
gx_y_x2(:,3,nx) = (x2_plus(:,nx)-P2)/(x2_t_plus-P2);
end
%%%%%%%%%%%%%%%%%%%%%%%%dx/dx2
for ny=1:n_points
vy_x_x2(1,:,ny) = dx1_minus(:,ny);
vy_x_x2(3,:,ny) = dx1_plus(:,ny);

g1y_x_x2(2,:,ny) = dx1_plus(:,ny);
g2y_x_x2(2,:,ny) = dx1_minus(:,ny);
end
for nx=1:n_points
ux_x_x2(1,:,nx) = x1(1,nx)/x1_t_minus;
ux_x_x2(3,:,nx) = (x1(3,nx)-P1)/(x1_t_plus-P1);

f1x_x_x2(2,:,nx) = (x1(2,nx)-x1_t_minus)/(x1_t_plus-x1_t_minus);
f2x_x_x2(2,:,nx) = (x1(2,nx)-x1_t_plus)/(x1_t_minus-x1_t_plus);
end
%%%%%%%%%%%%%%%%%%%%%%%%dy/dx1
for nx=1:n_points
vx_y_x1(:,1,nx) = dx2_minus(:,nx);
vx_y_x1(:,3,nx) = dx2_plus(:,nx);

g1x_y_x1(:,2,nx) = dx2_plus(:,nx);
g2x_y_x1(:,2,nx) = dx2_minus(:,nx);
end
for ny=1:n_points
uy_y_x1(:,1,ny) = x2(1,ny)/x2_t_minus;
uy_y_x1(:,3,ny) = (x2(3,ny)-P2)/(x2_t_plus-P2);

f1y_y_x1(:,2,ny) = (x2(2,ny)-x2_t_minus)/(x2_t_plus-x2_t_minus);
f2y_y_x1(:,2,ny) = (x2(2,ny)-x2_t_plus)/(x2_t_minus-x2_t_plus);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comp_x = zeros(1,n_points);
comp_y = zeros(1,n_points);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%(dx/dx1)*(dx/dx2)
[matrix_dx_dx1_dx_dx2_for_Ex,ellipse_matrix_dx_dx1_dx_dx2_for_Ex,...
    matrix_dx_dx1_dx_dx2_for_Ey,ellipse_matrix_dx_dx1_dx_dx2_for_Ey]=...
    dx_dx1_dx_dx2(gy_x_x1,ux_x_x2,vy_x_x2, f1x_x_x2, f2x_x_x2, g1y_x_x2, g2y_x_x2,...
    N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi);            
%%%%%%%%%%%%%%(dy/dx1)*(dy/dx2)
[matrix_dy_dx1_dy_dx2_for_Ex,ellipse_matrix_dy_dx1_dy_dx2_for_Ex,...
    matrix_dy_dx1_dy_dx2_for_Ey,ellipse_matrix_dy_dx1_dy_dx2_for_Ey]=...
    dy_dx1_dy_dx2(gx_y_x2,uy_y_x1,vx_y_x1,g1x_y_x1,g2x_y_x1,f1y_y_x1,f2y_y_x1,...
    N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi);
%%%%%%
int_Ex_g_down_12_1D = matrix_dy_dx1_dy_dx2_for_Ex + matrix_dx_dx1_dx_dx2_for_Ex;
int_Ey_g_down_12_1D = matrix_dy_dx1_dy_dx2_for_Ey + matrix_dx_dx1_dx_dx2_for_Ey;

int_Ex_g_down_12_ellipse_1D = ellipse_matrix_dx_dx1_dx_dx2_for_Ex + ellipse_matrix_dy_dx1_dy_dx2_for_Ex;
int_Ey_g_down_12_ellipse_1D = ellipse_matrix_dx_dx1_dx_dx2_for_Ey + ellipse_matrix_dy_dx1_dy_dx2_for_Ey;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%(dx/dx1)^2
[matrix_dx_dx1_sqr_for_Ey,ellipse_matrix_dx_dx1_sqr_for_Ey]=...
    dx_dx1_sqr(gy_x_x1,N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi);
%%%%%%%%%%%%%%(dy/dx1)^2
[matrix_dy_dx1_sqr_for_Ey,ellipse_matrix_dy_dx1_sqr_for_Ey]=...
    dy_dx1_sqr(vx_y_x1,uy_y_x1,g1x_y_x1,g2x_y_x1,f1y_y_x1,f2y_y_x1,...
    N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi);
%%%%%%
int_Ey_g_down_11_1D = matrix_dx_dx1_sqr_for_Ey + matrix_dy_dx1_sqr_for_Ey;
int_Ey_g_down_11_ellipse_1D = ellipse_matrix_dx_dx1_sqr_for_Ey + ellipse_matrix_dy_dx1_sqr_for_Ey;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%(dy/dx2)^2
[matrix_dy_dx2_sqr_for_Ex,ellipse_matrix_dy_dx2_sqr_for_Ex]=...
    dy_dx2_sqr(gx_y_x2, N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi);
%%%%%%%%%%%%%%%%%%%(dx/dx2)^2
[matrix_dx_dx2_sqr_for_Ex,ellipse_matrix_dx_dx2_sqr_for_Ex]=...
    dx_dx2_sqr(ux_x_x2,vy_x_x2,f1x_x_x2,f2x_x_x2,g1y_x_x2, g2y_x_x2,...
    N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi);
%%%%%%
int_Ex_g_down_22_1D = matrix_dy_dx2_sqr_for_Ex + matrix_dx_dx2_sqr_for_Ex;
int_Ex_g_down_22_ellipse_1D = ellipse_matrix_dy_dx2_sqr_for_Ex + ellipse_matrix_dx_dx2_sqr_for_Ex;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[matrix_dx_dx1_sqr_for_Ez,ellipse_matrix_dx_dx1_sqr_for_Ez,...
    matrix_dx_dx1_sqr_for_Dx,ellipse_matrix_dx_dx1_sqr_for_Dx,...
    matrix_dx_dx1_sqr_for_Dy,ellipse_matrix_dx_dx1_sqr_for_Dy]=...
    dx_x1_sqr(gx_y_x2,gy_x_x1, N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi);
%%%%%%%%%%%%%%(dx/dx1)*(dy/dx2)
dx_dx1_dy_dx2_comp_x = gx_y_x2;
dx_dx1_dy_dx2_comp_y = gy_x_x1;
matrix_dx_dx1_dy_dx2_for_Ez = zeros(N_total_3,N_total_3);
for num_x = 1:3
    for num_y = 1:3
            %for Ez
            comp_x(1,:) = dx_dx1_dy_dx2_comp_x(num_x,num_y,:); 
            matrix_comp_x = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dx_dx1_dy_dx2_comp_y(num_x,num_y,:); 
            matrix_comp_y = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            
            matrix_x_y = kron(matrix_comp_x,matrix_comp_y);
            matrix_dx_dx1_dy_dx2_for_Ez = matrix_dx_dx1_dy_dx2_for_Ez + matrix_x_y;
            if (num_x==2)&&(num_y==2)
                ellipse_matrix_dx_dx1_dy_dx2_for_Ez = matrix_x_y;
            end
    end
end
%%%%%%%%%%%%%%(dx/dx2)*(dy/dx1)
dx_dx2_dy_dx1_comp_x = ux_x_x2.*vx_y_x1;
dx_dx2_dy_dx1_comp_y = vy_x_x2.*uy_y_x1;


dx_dx2_dy_dx1_on2_x1_comp_x_1 = vx_y_x1.*f1x_x_x2;
dx_dx2_dy_dx1_on2_x1_comp_x_2 = vx_y_x1.*f2x_x_x2;
dx_dx2_dy_dx1_on2_x1_comp_y_1 = uy_y_x1.*g1y_x_x2;
dx_dx2_dy_dx1_on2_x1_comp_y_2 = uy_y_x1.*g2y_x_x2;

dx_dx2_dy_dx1_on2_x2_comp_x_1 = ux_x_x2.*g1x_y_x1;
dx_dx2_dy_dx1_on2_x2_comp_x_2 = ux_x_x2.*g2x_y_x1;
dx_dx2_dy_dx1_on2_x2_comp_y_1 = vy_x_x2.*f1y_y_x1;
dx_dx2_dy_dx1_on2_x2_comp_y_2 = vy_x_x2.*f2y_y_x1;


dx_dx2_dy_dx1_on22_comp_x_1 = f1x_x_x2.*g1x_y_x1;
dx_dx2_dy_dx1_on22_comp_x_2 = f2x_x_x2.*g1x_y_x1;
dx_dx2_dy_dx1_on22_comp_x_3 = f1x_x_x2.*g2x_y_x1;
dx_dx2_dy_dx1_on22_comp_x_4 = f2x_x_x2.*g2x_y_x1;
dx_dx2_dy_dx1_on22_comp_y_1 = g1y_x_x2.*f1y_y_x1;
dx_dx2_dy_dx1_on22_comp_y_1 = g2y_x_x2.*f1y_y_x1;
dx_dx2_dy_dx1_on22_comp_y_1 = g1y_x_x2.*f2y_y_x1;
dx_dx2_dy_dx1_on22_comp_y_1 = g2y_x_x2.*f2y_y_x1;




