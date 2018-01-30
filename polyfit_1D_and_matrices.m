function [matrix_x] = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux_full,Tx_full,aTx_full)
            
 Ux(:,:) = Ux_full(1,:,:);
 Tx(:,:) = Tx_full(1,:,:);
 aTx(:,:) = aTx_full(1,:,:);
 
 fit_x = polyfit(ksi,comp_x,7);
 
 matrix_x = fit_x(8)*Ux + fit_x(7)*aTx + fit_x(6)*Tx*aTx+...
        + fit_x(5)*Tx^2*aTx + fit_x(4)*Tx^3*aTx + fit_x(3)*Tx^4*aTx +...
       fit_x(2)*Tx^5*aTx + fit_x(1)*Tx^6*aTx;