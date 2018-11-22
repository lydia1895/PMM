function [R, R1] = PMM_reciprocal_S(w, dw)

c = 3*10^5;
lambda = c/w;
lambda1 = c/(w+dw);
R = Rmatrix(lambda)