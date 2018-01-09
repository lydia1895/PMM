function [a] = PMM_boundary_conditions(La, tau, N_intervals, N_basis, N, n)
%{
N_intervals = number of intervals
N_basis(i) = number of basis functions on interval i
%}

N_total = sum(N_basis);  %total number of basis functions
N_max = max(N_basis);
N_total_3 = N_total - N_intervals;

p = zeros(N_max,1);     %n-th Gegenbauer polynomial at 1
p_1 = zeros(N_max,1);   %n-th Gegenbauer polynomial at (-1)
for m = 0:(N_max-1)
    p(m+1) = gamma(m+2*La)/(gamma(2*La)*gamma(m+1));
    if mod(m,2)==0
        p_1(m+1) = p(m+1);
    else
        p_1(m+1) = -p(m+1);
    end
end

% M(p)*a = right(p) -> get a

M = zeros(N_intervals, N_total);
right = zeros(N_intervals, N_total_3);

for k = 1:(N_intervals-1)
    M(k,N(k+1)+n(k+1)+k +1) = p_1(n(k+1) +1);
    M(k,N(k)+n(k)+k-1 +1) = -p(n(k) +1);
    for j = 1:N_intervals
        for i=(N(j)):(N(j)+n(j)-1)
        if (j==k)
            right(k,i +1) = p(i-N(k) +1);
        end
        if (j==(k+1))
            right(k,i +1) = -p_1(i-N(k+1) +1);
        end
        end
    end
end

m = N_intervals;
M(m,N(m)+n(m)+m-1 +1) = p(n(m) +1);
M(m,n(1) +1) = -tau*p_1(n(1) +1);
for j=1:N_intervals
    for i=(N(j)):(N(j)+n(j)-1)
        if (j==m)
            right(m,i +1) = -p(i-N(m) +1);
        end
        if (j==1)
            right(m,i +1) = tau*p_1(i +1);
        end
    end
end
            
a = M\right;
for k = 1:N_intervals
    for i = (N(k)+1):(N(k)+n(k))
        a(i+k-1,i) = 1;
    end
end
    

%{
M = zeros(N_intervals, N_intervals);
right = zeros(N_intervals, N_total);

%conditions at internal boundaries of the unit cell

for k = 1:(N_intervals-1)
    for j = 1:N_intervals
        for i=(N(j)+1):(N(j)+n(j))
            M(k, k+1) = p_1(n(k+1)); %coef for a(N(k+1)+n(k+1)+k, i)
            M(k, k) = - p(n(k));     %coef for a(N(k)+n(k)+k-1, i)
            if j==k
                right(k,i) = p(i-N(k));
            end
            if j==(k+1)
                right(k,i) = -p_1(i-N(k+1));
            end
            if (j~=k)&&(j~=(k+1))
                right(k,i) = 0; 
            end
        end
    end
end

%conditions at boundaries of the unit cell

for j = 1:N_intervals
    for i=(N(j)+1):(N(j)+n(j))
        m = N_intervals;
        M(m, 1) = tau*p_1(n(1)); %a(N(m)+n(m)+m-1, i)
        M(m, m) = - p(n(m));     %a(n(1), i)
        if j==m
            right(m,i) = p(i-N(m));
        end
        if j==1
            right(m,i) = -tau*p_1(i);
        end
        if (j~=1)&&(j~=m)
            right(m,i) = 0; 
        end
            
    end
end

a_short = M\right;
a = zeros(N_total-N_intervals, N_total);
%a = zeros(N_total, N_total);
%a is the full matrix, a*Pn = ~Pn


for k = 1:N_intervals
    for i = (N(k)+1):(N(k)+n(k))
        a(i,i+k-1) = 1;
        for j = 1:N_intervals
            a(i,N(j)+n(j)+j) = a_short(j,i);
        end
    end
end 

a = transpose(a)
%}





