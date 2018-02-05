function [matrix_dy_dx2_sqr_for_Ex,ellipse_matrix_dy_dx2_sqr_for_Ex]=...
    dy_dx2_sqr(gx_y_x2,...
    N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi)

dy_dx2_sqr_comp_x = gx_y_x2.^2;
matrix_dy_dx2_sqr_for_Ex = zeros(N_total_3,N_total_3);
for num_x = 1:3
    for num_y = 1:3
        %for Ex
        if (num_x~=2)&&(num_y~=2)
            matrix_comp_x(:,:) = Ux(num_x,:,:);
        else
            comp_x(1,:) = dy_dx2_sqr_comp_x(num_x,num_y,:);
            matrix_comp_x(:,:)= polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),Tx(num_x,:,:));
        end
        matrix_comp_y(:,:) = Uy(num_y,:,:);
        
        matrix_x_y = kron(matrix_comp_x,matrix_comp_y);
        matrix_dy_dx2_sqr_for_Ex = matrix_dy_dx2_sqr_for_Ex + matrix_x_y;
        if (num_x==2)&&(num_y==2)
            ellipse_matrix_dy_dx2_sqr_for_Ex = matrix_x_y;
        end
    end
end