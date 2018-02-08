function [matrix_dx_dx1_dx_dx2,ellipse_matrix_dx_dx1_dx_dx2]=...
    dx_dx1_dx_dx2(gy_x_x1,ux_x_x2,vy_x_x2, f1x_x_x2, f2x_x_x2, g1y_x_x2, g2y_x_x2,...
    N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi)
dx_dx1_dx_dx2_comp_x = ux_x_x2;
dx_dx1_dx_dx2_comp_y = gy_x_x1.*vy_x_x2;
dx_dx1_dx_dx2_on2_x1_comp_x_1 = f1x_x_x2;
dx_dx1_dx_dx2_on2_x1_comp_x_2 = f2x_x_x2;
dx_dx1_dx_dx2_on2_x1_comp_y_1 = gy_x_x1.*g1y_x_x2;
dx_dx1_dx_dx2_on2_x1_comp_y_2 = gy_x_x1.*g2y_x_x2;
matrix_dx_dx1_dx_dx2 = zeros(N_total_3,N_total_3);

for num_x = 1:3
    for num_y = 1:3
        if num_x~=2
            
            comp_x(1,:) = dx_dx1_dx_dx2_comp_x(num_x,num_y,:); 
            matrix_comp_x = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dx_dx1_dx_dx2_comp_y(num_x,num_y,:); 
            matrix_comp_y = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            
            matrix_x_y = kron(matrix_comp_x,matrix_comp_y);
            matrix_dx_dx1_dx_dx2 = matrix_dx_dx1_dx_dx2 + matrix_x_y;
            
        else
            
            comp_x(1,:) = dx_dx1_dx_dx2_on2_x1_comp_x_1(num_x,num_y,:); 
            matrix_comp_x_1 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dx_dx1_dx_dx2_on2_x1_comp_y_1(num_x,num_y,:); 
            matrix_comp_y_1 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x1_y1 = kron(matrix_comp_x_1,matrix_comp_y_1);
            
            
            comp_x(1,:) = dx_dx1_dx_dx2_on2_x1_comp_x_2(num_x,num_y,:); 
            matrix_comp_x_2 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dx_dx1_dx_dx2_on2_x1_comp_y_2(num_x,num_y,:); 
            matrix_comp_y_2 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x2_y2 = kron(matrix_comp_x_2,matrix_comp_y_2);
            
            
            matrix_dx_dx1_dx_dx2 = matrix_dx_dx1_dx_dx2 + matrix_x1_y1 + matrix_x2_y2;
            if (num_x==2)&&(num_y==2)
                ellipse_matrix_dx_dx1_dx_dx2 = matrix_x1_y1 + matrix_x2_y2;
            end
            
        end
    end
end