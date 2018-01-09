function [M1M2] = Kr(M1, M2)

[m,n] = size(M1);
[p,q] = size(M2);

M1M2 = zeros(m*p, n*q);

for r = 1:m
for s = 1:n
    for v = 1:p
    for w = 1:q
        row = p*(r-1)+v;
        col = q*(s-1)+w;
        M1M2(row, col) = M1(r,s)*M2(v,w);
    end
    end
end
end