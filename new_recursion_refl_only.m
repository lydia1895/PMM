function [Rud_p] = new_recursion_refl_only(Rud_p1, W_p, W_p1, pplus_p, pminus_p1)
%function [S_p] = new_recursion(S_p1, W_p, W_p1, pplus, pminus, N)

[N1, N2] = size(W_p);
NN = N1/2;

mzero = zeros(NN,NN);
W_p_plus = W_p(:,1:NN);
W_p_minus = W_p(:,(NN+1):2*NN);
W_p1_plus = W_p1(:,1:NN);
W_p1_minus = W_p1(:,(NN+1):2*NN);

%Tuu_p1 = S_p1( 1:NN, 1:NN );
%Rud_p1 = S_p1( 1:NN, (NN+1):2*NN );
%Rdu_p1 = S_p1( (NN+1):2*NN, 1:NN );
%Tdd_p1 = S_p1( (NN+1):2*NN, (NN+1):2*NN );

Z_plus = W_p1_plus;
Z_minus = -W_p_plus*pplus_p*Rud_p1 - W_p_minus;
%X1 = W_p_plus*pplus_p*Tuu_p1;
X2 = -W_p1_minus/pminus_p1;
Z = cat(2,Z_plus, Z_minus);

%ZX1 = Z\X1;
ZX2 = Z\X2;

%Tuu_p = ZX1(1:NN,:);
Rud_p = ZX2(1:NN,:);
%Rdu_p = Rdu_p1 + Tdd_p1 * ZX1((NN+1):(2*NN),:) ;
%Tdd_p = Tdd_p1 * ZX2((NN+1):(2*NN),:) ;

%S_p_up = cat(2, Tuu_p, Rud_p);
%S_p_down = cat(2, Rdu_p, Tdd_p);
%S_p = cat(1, S_p_up, S_p_down);
%S_p_up = cat(2, mzero, Rud_p);
%S_p_down = cat(2, mzero, mzero);
%S_p = cat(1, S_p_up, S_p_down);


