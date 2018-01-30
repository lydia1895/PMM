function [matrix_dy_dx1_dy_dx2_for_Ex,ellipse_matrix_dy_dx1_dy_dx2_for_Ex,...
    matrix_dy_dx1_dy_dx2_for_Ey,ellipse_matrix_dy_dx1_dy_dx2_for_Ey]=...
    dy_dx1_dy_dx2(gx_y_x2,uy_y_x1,vx_y_x1,g1x_y_x1,g2x_y_x1,f1y_y_x1,f2y_y_x1,...
    N_total_3, Ux,Tx,aTx,Uy,Ty,aTy,ksi)

dy_dx1_dy_dx2_comp_x = gx_y_x2.*vx_y_x1;
dy_dx1_dy_dx2_comp_y = uy_y_x1;
dy_dx1_dy_dx2_on2_x2_comp_x_1 = gx_y_x2.*g1x_y_x1;
dy_dx1_dy_dx2_on2_x2_comp_x_2 = gx_y_x2.*g2x_y_x1;
dy_dx1_dy_dx2_on2_x2_comp_y_1 = f1y_y_x1;
dy_dx1_dy_dx2_on2_x2_comp_y_2 = f2y_y_x1;
matrix_dy_dx1_dy_dx2_for_Ex = zeros(N_total_3,N_total_3);
matrix_dy_dx1_dy_dx2_for_Ey = zeros(N_total_3,N_total_3);
for num_x = 1:3
    for num_y = 1:3
        if num_y~=2
            %for Ex
            comp_x(1,:) = dy_dx1_dy_dx2_comp_x(num_x,num_y,:); 
            matrix_comp_x = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),Tx(num_x,:,:));
            
            comp_y(1,:) = dy_dx1_dy_dx2_comp_y(num_x,num_y,:); 
            matrix_comp_y = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            
            matrix_x_y = kron(matrix_comp_x,matrix_comp_y);
            matrix_dy_dx1_dy_dx2_for_Ex = matrix_dy_dx1_dy_dx2_for_Ex + matrix_x_y;
            %for Ey
            comp_x(1,:) = dy_dx1_dy_dx2_comp_x(num_x,num_y,:); 
            matrix_comp_x = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dy_dx1_dy_dx2_comp_y(num_x,num_y,:); 
            matrix_comp_y = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),Ty(num_y,:,:));
            
            matrix_x_y = kron(matrix_comp_x,matrix_comp_y);
            matrix_dy_dx1_dy_dx2_for_Ey = matrix_dy_dx1_dy_dx2_for_Ey + matrix_x_y;
        else
            %for Ex
            comp_x(1,:) = dy_dx1_dy_dx2_on2_x2_comp_x_1(num_x,num_y,:); 
            matrix_comp_x_1 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),Tx(num_x,:,:));
            
            comp_y(1,:) = dy_dx1_dy_dx2_on2_x2_comp_y_1(num_x,num_y,:); 
            matrix_comp_y_1 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x1_y1 = kron(matrix_comp_x_1,matrix_comp_y_1);
            
            
            comp_x(1,:) = dy_dx1_dy_dx2_on2_x2_comp_x_2(num_x,num_y,:); 
            matrix_comp_x_2 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),Tx(num_x,:,:));
            
            comp_y(1,:) = dy_dx1_dy_dx2_on2_x2_comp_y_2(num_x,num_y,:); 
            matrix_comp_y_2 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),aTy(num_y,:,:));
            matrix_x2_y2 = kron(matrix_comp_x_2,matrix_comp_y_2);
            
            
            matrix_dy_dx1_dy_dx2_for_Ex = matrix_dy_dx1_dy_dx2_for_Ex + matrix_x1_y1 + matrix_x2_y2;
            if (num_x==2)&&(num_y==2)
                ellipse_matrix_dy_dx1_dy_dx2_for_Ex = matrix_x1_y1 + matrix_x2_y2;
            end
            %for Ey
            comp_x(1,:) = dy_dx1_dy_dx2_on2_x2_comp_x_1(num_x,num_y,:); 
            matrix_comp_x_1 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dy_dx1_dy_dx2_on2_x2_comp_y_1(num_x,num_y,:); 
            matrix_comp_y_1 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),Ty(num_y,:,:));
            matrix_x1_y1 = kron(matrix_comp_x_1,matrix_comp_y_1);
            
            
            comp_x(1,:) = dy_dx1_dy_dx2_on2_x2_comp_x_2(num_x,num_y,:); 
            matrix_comp_x_2 = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux(num_x,:,:),Tx(num_x,:,:),aTx(num_x,:,:));
            
            comp_y(1,:) = dy_dx1_dy_dx2_on2_x2_comp_y_2(num_x,num_y,:); 
            matrix_comp_y_2 = polyfit_1D_and_matrices(comp_y,ksi,...
                Uy(num_y,:,:),Ty(num_y,:,:),Ty(num_y,:,:));
            matrix_x2_y2 = kron(matrix_comp_x_2,matrix_comp_y_2);
            
            
            matrix_dy_dx1_dy_dx2_for_Ey = matrix_dy_dx1_dy_dx2_for_Ey + matrix_x1_y1 + matrix_x2_y2;
            if (num_x==2)&&(num_y==2)
                ellipse_matrix_dy_dx1_dy_dx2_for_Ey = matrix_x1_y1 + matrix_x2_y2;
            end
        end
    end
end