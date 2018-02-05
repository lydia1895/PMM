function [matrix_dx_dx2_dy_dx1,ellipse_matrix_dx_dx2_dy_dx1]=...
    dx_x2_dy_dx1(ux_x_x2,vx_y_x1,vy_x_x2,uy_y_x1,...
    f1x_x_x2,f2x_x_x2,g1x_y_x1,g2x_y_x1,...
    f1y_y_x1,f2y_y_x1,g1y_x_x2,g2y_x_x2,...
    N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi)

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
dx_dx2_dy_dx1_on22_comp_y_2 = g2y_x_x2.*f1y_y_x1;
dx_dx2_dy_dx1_on22_comp_y_3 = g1y_x_x2.*f2y_y_x1;
dx_dx2_dy_dx1_on22_comp_y_4 = g2y_x_x2.*f2y_y_x1;

matrix_dx_dx2_dy_dx1 = zeros(N_total_3,N_total_3);
for num_x = 1:3
    for num_y = 1:3
        %{
        if (num_x~=2)&&(num_y~=2)
            comp_x(1,:) = dx_dx2_dy_dx1_comp_x(num_x,num_y,:);
            matrix_comp_x = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dx_dx2_dy_dx1_comp_y(num_x,num_y,:);
            matrix_comp_y = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            
            matrix_x_y = kron(matrix_comp_x,matrix_comp_y);
            matrix_dx_dx2_dy_dx1 = matrix_dx_dx2_dy_dx1 + matrix_x_y;
        end
        %}
        
        if (num_x==2)&&(num_y~=2)
            
            comp_x(1,:) = dx_dx2_dy_dx1_on2_x1_comp_x_1(num_x,num_y,:);
            matrix_comp_x_1 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            comp_y(1,:) = dx_dx2_dy_dx1_on2_x1_comp_y_1(num_x,num_y,:);
            matrix_comp_y_1 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x1_y1 = kron(matrix_comp_x_1,matrix_comp_y_1);
            
            comp_x(1,:) = dx_dx2_dy_dx1_on2_x1_comp_x_2(num_x,num_y,:);
            matrix_comp_x_2 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            comp_y(1,:) = dx_dx2_dy_dx1_on2_x1_comp_y_2(num_x,num_y,:);
            matrix_comp_y_2 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x2_y2 = kron(matrix_comp_x_2,matrix_comp_y_2);
            matrix_dx_dx2_dy_dx1 = matrix_dx_dx2_dy_dx1 + matrix_x1_y1 + matrix_x2_y2;
        end
        if (num_x~=2)&&(num_y==2)
            
            comp_x(1,:) = dx_dx2_dy_dx1_on2_x2_comp_x_1(num_x,num_y,:);
            matrix_comp_x_1 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            comp_y(1,:) = dx_dx2_dy_dx1_on2_x2_comp_y_1(num_x,num_y,:);
            matrix_comp_y_1 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x1_y1 = kron(matrix_comp_x_1,matrix_comp_y_1);
            
            comp_x(1,:) = dx_dx2_dy_dx1_on2_x2_comp_x_2(num_x,num_y,:);
            matrix_comp_x_2 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            comp_y(1,:) = dx_dx2_dy_dx1_on2_x2_comp_y_2(num_x,num_y,:);
            matrix_comp_y_2 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x2_y2 = kron(matrix_comp_x_2,matrix_comp_y_2);
            matrix_dx_dx2_dy_dx1 = matrix_dx_dx2_dy_dx1 + matrix_x1_y1 + matrix_x2_y2;
        end
        if (num_x==2)&&(num_y==2)
          
            comp_x(1,:) = dx_dx2_dy_dx1_on22_comp_x_1(num_x,num_y,:);
            matrix_comp_x_1 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            comp_y(1,:) = dx_dx2_dy_dx1_on22_comp_y_1(num_x,num_y,:);
            matrix_comp_y_1 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x1_y1 = kron(matrix_comp_x_1,matrix_comp_y_1);
            
            comp_x(1,:) = dx_dx2_dy_dx1_on22_comp_x_2(num_x,num_y,:);
            matrix_comp_x_2 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            comp_y(1,:) = dx_dx2_dy_dx1_on22_comp_y_2(num_x,num_y,:);
            matrix_comp_y_2 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x2_y2 = kron(matrix_comp_x_2,matrix_comp_y_2);
            
            comp_x(1,:) = dx_dx2_dy_dx1_on22_comp_x_3(num_x,num_y,:);
            matrix_comp_x_3 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            comp_y(1,:) = dx_dx2_dy_dx1_on22_comp_y_3(num_x,num_y,:);
            matrix_comp_y_3 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x3_y3 = kron(matrix_comp_x_3,matrix_comp_y_3);
            
            comp_x(1,:) = dx_dx2_dy_dx1_on22_comp_x_4(num_x,num_y,:);
            matrix_comp_x_4 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            comp_y(1,:) = dx_dx2_dy_dx1_on22_comp_y_4(num_x,num_y,:);
            matrix_comp_y_4 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x4_y4 = kron(matrix_comp_x_4,matrix_comp_y_4);
            
            matrix_dx_dx2_dy_dx1 = matrix_dx_dx2_dy_dx1 +...
                matrix_x1_y1 + matrix_x2_y2 + + matrix_x3_y3 + matrix_x4_y4;
            ellipse_matrix_dx_dx2_dy_dx1 = ...
                matrix_x1_y1 + matrix_x2_y2 + + matrix_x3_y3 + matrix_x4_y4;
            
        end
    end
end