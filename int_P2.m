function [int_P2_x] = int_P2(Ux,Tx,aTx,g_fit)

int_P2_x = g_fit.p8 *Ux + g_fit.p7 *aTx + g_fit.p6 *Tx*aTx +...
        g_fit.p5 *Tx^2*aTx + g_fit.p4 *Tx^3*aTx + g_fit.p3 *Tx^4*aTx+...
        g_fit.p2 *Tx^5*aTx+g_fit.p1 *Tx^6*aTx;
%{
int_P2_x = g_fit(1) *Ux + g_fit(2) *aTx + g_fit(3) *Tx*aTx +...
        g_fit(4) *Tx^2*aTx + g_fit(5) *Tx^3*aTx + g_fit(6) *Tx^4*aTx;
%}
    