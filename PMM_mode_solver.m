function [w0,dw0]= PMM_mode_solver(R0,R1,w,dw)

dR = (R1-R0)/dw;
A=-dR\R0;
[V, D] = eig(A);
w_array = diag(D);
%{
for i=1:size(w_array)
    if norm(w_array(i))<50 && norm(w_array(i))>2
        ii=i;
    end
end
dw0 = w_array(ii)
%}
dw0 = min(w_array);
w0 = w + dw0;