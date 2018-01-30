function [matrix_dx_dx1_sqr_for_Ey,ellipse_matrix_dx_dx1_sqr_for_Ey]=...
    dx_dx1_sqr(gy_x_x1,...
    N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi)
%%%%%%%%%%%%%%%%%%(dx/dx1)^2
dx_dx1_sqr_comp_y = gy_x_x1.^2;
matrix_dx_dx1_sqr_for_Ey = zeros(N_total_3,N_total_3);
for num_x = 1:3
    for num_y = 1:3
            %for Ey
            matrix_comp_x(:,:) = Ux(num_x,:,:);
            
            comp_y(1,:) = dx_dx1_sqr_comp_y(num_x,num_y,:); 
            matrix_comp_y = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),Ty(num_y,:,:));
            
            matrix_x_y = kron(matrix_comp_x,matrix_comp_y);
            matrix_dx_dx1_sqr_for_Ey = matrix_dx_dx1_sqr_for_Ey + matrix_x_y;
            if (num_x==2)&&(num_y==2)
                ellipse_matrix_dx_dx1_sqr_for_Ey = matrix_x_y;
            end
    end
end