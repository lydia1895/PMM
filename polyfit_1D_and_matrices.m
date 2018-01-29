function [matrix_x] = polyfit_1D_and_matrices(comp_x,ksi,...
                Ux_full,Tx_full,aTx_full)
            
 Ux(:,:) = Ux_full(1,:,:);
 Tx(:,:) = Tx_full(1,:,:);
 aTx(:,:) = aTx_full(1,:,:);
 
 fit_x = fit(ksi,comp_x,'poly7');
 
 matrix_x = fit_x.p8*Ux + fit.p7*aTx + fit.p6*Tx*aTx+...
        + fit.p5*Tx^2*aTx + fit.p4*Tx^3*aTx + fit.p3*Tx^4*aTx +...
       fit.p2*Tx^5*aTx + fit.p1*Tx^6*aTx;