function [matrix_dx_dx1_dy_dx2,ellipse_matrix_dx_dx1_dy_dx2]=...
    dx_x1_dy_dx2(gx_y_x2,gy_x_x1, N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi)
%%%%%%%%%%%%%%(dx/dx1)*(dy/dx2)
dx_dx1_dy_dx2_comp_x = gx_y_x2;
dx_dx1_dy_dx2_comp_y = gy_x_x1;
matrix_dx_dx1_dy_dx2 = zeros(N_total_3,N_total_3);
for num_x = 1:3
    for num_y = 1:3
            comp_x(1,:) = dx_dx1_dy_dx2_comp_x(num_x,num_y,:); 
            matrix_comp_x = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dx_dx1_dy_dx2_comp_y(num_x,num_y,:); 
            matrix_comp_y = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            if (num_x~=2)&&(num_y~=2)
                matrix_comp_x(:,:) = Ux(num_x,:,:);
                matrix_comp_y(:,:) = Uy(num_y,:,:);
            end
            matrix_x_y = kron(matrix_comp_x,matrix_comp_y);
            matrix_dx_dx1_dy_dx2 = matrix_dx_dx1_dy_dx2 + matrix_x_y;
            if (num_x==2)&&(num_y==2)
                ellipse_matrix_dx_dx1_dy_dx2 = matrix_x_y;
            end
    end
end