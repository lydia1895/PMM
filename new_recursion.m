function [S_p] = new_recursion(S_p1, W_p, W_p1, pplus_p, pminus_p1)
%function [S_p] = new_recursion(S_p1, W_p, W_p1, pplus, pminus, N)

[N1, N2] = size(W_p);
NN = N1/2;

W_p_plus = W_p(:,1:NN);
W_p_minus = W_p(:,(NN+1):2*NN);
W_p1_plus = W_p1(:,1:NN);
W_p1_minus = W_p1(:,(NN+1):2*NN);

Tuu_p1 = S_p1( 1:NN, 1:NN );
Rud_p1 = S_p1( 1:NN, (NN+1):2*NN );
Rdu_p1 = S_p1( (NN+1):2*NN, 1:NN );
Tdd_p1 = S_p1( (NN+1):2*NN, (NN+1):2*NN );

Z_plus = W_p1_plus;
Z_minus = -W_p_plus*pplus_p*Rud_p1 - W_p_minus;
X1 = W_p_plus*pplus_p*Tuu_p1;
X2 = -W_p1_minus*pminus_p1;
Z = cat(2,Z_plus, Z_minus);

unity = eye(2*NN,2*NN);
inv_Z=Z\unity;
ZX1 = inv_Z*X1;
ZX2 = inv_Z*X2;

%ZX1 = Z\X1;
%ZX2 = Z\X2;

Tuu_p = ZX1(1:NN,:);
%Tuu_p = zeros(NN,NN);
Rud_p = ZX2(1:NN,:);
%Rdu_p = zeros(NN,NN);
Rdu_p = Rdu_p1 + Tdd_p1 * ZX1((NN+1):(2*NN),:) ;
Tdd_p = Tdd_p1 * ZX2((NN+1):(2*NN),:) ;

S_p_up = cat(2, Tuu_p, Rud_p);
S_p_down = cat(2, Rdu_p, Tdd_p);
S_p = cat(1, S_p_up, S_p_down);

%{
NN = 2*(2*N+1)*(2*N+1);
W11_p = W_p( 1:NN, 1:NN );
W12_p = W_p( 1:NN, (NN+1):2*NN );
W21_p = W_p( (NN+1):2*NN, 1:NN );
W22_p = W_p( (NN+1):2*NN, (NN+1):2*NN );

W11_p1 = W_p1( 1:NN, 1:NN );
W12_p1 = W_p1( 1:NN, (NN+1):2*NN );
W21_p1 = W_p1( (NN+1):2*NN, 1:NN );
W22_p1 = W_p1( (NN+1):2*NN, (NN+1):2*NN );

Tuu_p1 = S_p1( 1:NN, 1:NN );
Rud_p1 = S_p1( 1:NN, (NN+1):2*NN );
Rdu_p1 = S_p1( (NN+1):2*NN, 1:NN );
Tdd_p1 = S_p1( (NN+1):2*NN, (NN+1):2*NN );

Rud_p1_w = pplus*Rud_p1/pminus;
Tdd_p1_w = Tdd_p1/pminus;
Tuu_p1_w = pplus*Tuu_p1;

Z11 = W11_p1;
Z12 = -W11_p*Rud_p1_w - W12_p;
Z21 = W21_p1;
Z22 = -W21_p*Rud_p1_w - W22_p;

Z_up = cat(2, Z11, Z12);
Z_down = cat(2, Z21, Z22);
Z = cat(1, Z_up, Z_down);


X11 = W11_p*Tuu_p1_w;
X12 = -W12_p1;
X21 = W21_p*Tuu_p1_w;
X22 = - W22_p1;

X1 = cat(1, X11, X21);
X2 = cat(1, X12, X22);

ZX1 = Z\X1;
ZX2 = Z\X2;

Rud_p = ZX2(1:NN,:);
Tdd_p = Tdd_p1_w * ZX2((NN+1):(2*NN),:) ;
Tuu_p = ZX1(1:NN,:);
Rdu_p = Rdu_p1 + Tdd_p1_w * ZX1((NN+1):(2*NN),:) ;

S_p_up = cat(2, Tuu_p, Rud_p);
S_p_down = cat(2, Rdu_p, Tdd_p);
S_p = cat(1, S_p_up, S_p_down);


%}
