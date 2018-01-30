function [matrix_dy_dx1_sqr_for_Ey,ellipse_matrix_dy_dx1_sqr_for_Ey]=...
    dy_dx1_sqr(vx_y_x1,uy_y_x1,g1x_y_x1,g2x_y_x1,f1y_y_x1,f2y_y_x1,...
    N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi)
%%%%%%%%%%%%%%%%%%(dy/dx1)^2
dy_dx1_sqr_comp_x = vx_y_x1.^2;
dy_dx1_sqr_comp_y = uy_y_x1.^2;

dy_dx1_sqr_on2_x2_comp_x_1 = g1x_y_x1.^2;
dy_dx1_sqr_on2_x2_comp_x_2 = 2*g1x_y_x1.*g2x_y_x1;
dy_dx1_sqr_on2_x2_comp_x_3 = g2x_y_x1.^2;

dy_dx1_sqr_on2_x2_comp_y_1 = f1y_y_x1.^2;
dy_dx1_sqr_on2_x2_comp_y_2 = f1y_y_x1.*f2y_y_x1;
dy_dx1_sqr_on2_x2_comp_y_3 = f2y_y_x1.^2;
matrix_dy_dx1_sqr_for_Ey = zeros(N_total_3,N_total_3);
for num_x = 1:3
    for num_y = 1:3
        if num_y~=2
            %for Ey
            comp_x(1,:) = dy_dx1_sqr_comp_x(num_x,num_y,:); 
            matrix_comp_x = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dy_dx1_sqr_comp_y(num_x,num_y,:); 
            matrix_comp_y = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),Ty(num_y,:,:));
            
            matrix_x_y = kron(matrix_comp_x,matrix_comp_y);
            matrix_dy_dx1_sqr_for_Ey = matrix_dy_dx1_sqr_for_Ey + matrix_x_y;
        else
            comp_x(1,:) = dy_dx1_sqr_on2_x2_comp_x_1(num_x,num_y,:); 
            matrix_comp_x_1 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dy_dx1_sqr_on2_x2_comp_y_1(num_x,num_y,:); 
            matrix_comp_y_1 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),Ty(num_y,:,:));
            matrix_x1_y1 = kron(matrix_comp_x_1,matrix_comp_y_1);
            
            comp_x(1,:) = dy_dx1_sqr_on2_x2_comp_x_2(num_x,num_y,:); 
            matrix_comp_x_2 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dy_dx1_sqr_on2_x2_comp_y_2(num_x,num_y,:); 
            matrix_comp_y_2 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),Ty(num_y,:,:));
            matrix_x2_y2 = kron(matrix_comp_x_2,matrix_comp_y_2);
            
            comp_x(1,:) = dy_dx1_sqr_on2_x2_comp_x_3(num_x,num_y,:); 
            matrix_comp_x_3 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dy_dx1_sqr_on2_x2_comp_y_3(num_x,num_y,:); 
            matrix_comp_y_3 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),Ty(num_y,:,:));
            matrix_x3_y3 = kron(matrix_comp_x_3,matrix_comp_y_3);
            
            matrix_dy_dx1_sqr_for_Ey = matrix_dy_dx1_sqr_for_Ey +...
                matrix_x1_y1 + matrix_x2_y2 + matrix_x3_y3;
            if (num_x==2)&&(num_y==2)
                ellipse_matrix_dy_dx1_sqr_for_Ey = matrix_x1_y1 + matrix_x2_y2 + matrix_x3_y3;
            end
        end
    end
end