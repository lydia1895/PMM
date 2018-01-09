function [eps_xx_yx] = int_P3_Q2(Ux,Uy,Tx,aTy,Ty,g_fit)

eps_xx_yx = g_fit.p00 *Kr(Ux,Uy) + ...
    g_fit.p10 *Kr(Ux,aTy) + g_fit.p01 *Kr(Tx,Uy)+...
    g_fit.p20 *Kr(Ux,Ty*aTy) + g_fit.p11 *Kr(Tx,aTy) + g_fit.p02 *Kr(Tx^2,Uy)+...
    g_fit.p30 *Kr(Ux,Ty^2*aTy) + g_fit.p21 *Kr(Tx,Ty*aTy) + g_fit.p12 *Kr(Tx^2,aTy) +...
    g_fit.p03 *Kr(Tx^3,Uy)+...
    g_fit.p40 *Kr(Ux,Ty^3*aTy) + g_fit.p31 *Kr(Tx,Ty^2*aTy) + g_fit.p22 *Kr(Tx^2,Ty*aTy) +...
        g_fit.p13 *Kr(Tx^3,aTy) + g_fit.p04 *Kr(Tx^4,Uy)+...
    g_fit.p50 *Kr(Ux,Ty^4*aTy) + g_fit.p41 *Kr(Tx,Ty^3*aTy) + g_fit.p32 *Kr(Tx^2,Ty^2*aTy) +...
        g_fit.p23 *Kr(Tx^3,Ty*aTy) + g_fit.p14 *Kr(Tx^4,aTy) + g_fit.p05 *Kr(Tx^5,Uy);
    