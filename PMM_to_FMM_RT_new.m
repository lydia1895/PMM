function [fx_coef,fy_coef] = PMM_to_FMM_RT_new(N, NN, La, alpha_ref, beta_ref,... 
   b_x1, b_x2, N_intervals_x, N_intervals_y, N_basis_x, N_basis_y, Nx, nx, Ny, ny,...
   ax, ay)

n_points_int = 200;
ksi1 = linspace(-1,1,n_points_int);
ksi2 = linspace(-1,1,n_points_int);
%x1 = zeros(n_points_int,N_intervals_x);
%x2 = zeros(n_points_int,N_intervals_y);
x1 = zeros(N_intervals_x,n_points_int);
x2 = zeros(N_intervals_y,n_points_int);
for k1 = 1:N_intervals_x
    x1(k1,:) = ( (b_x1(k1+1)-b_x1(k1))*ksi1+(b_x1(k1+1)+b_x1(k1)) )/2;
end
for k2 = 1:N_intervals_y     
    x2(k2,:) = ( (b_x2(k2+1)-b_x2(k2))*ksi2+(b_x2(k2+1)+b_x2(k2)) )/2;
end

%x(k1) = PMM_matched_coordinates_ellipse(x1(k1), x2(k2), P1, P2, R1, R2)
%y(k2) = PMM_matched_coordinates_ellipse(x1(k1), x2(k2), P1, P2, R1, R2)
x = x1;
y = x2;

Nmax_x = max(N_basis_x);   %maximum number of Gegenbauer polynomial on all intervals
Nmax_y = max(N_basis_y);   %maximum number of Gegenbauer polynomial on all intervals
Nmax = max(Nmax_x, Nmax_y);
C = zeros(Nmax, n_points_int);

for m=1:Nmax
    %for i = 1:n_points_int
    %    C(m,i) = mfun('G',m-1,La,ksi1(i));
    %end    
    C(m,:) = gegenbauerC(m-1,La,ksi1);
end


lx1=zeros(N_intervals_x,1);
for k=1:N_intervals_x
    lx1(k) = (b_x1(k+1) - b_x1(k))/2;

end
lx2=zeros(N_intervals_y,1);
for k=1:N_intervals_y
    lx2(k) = (b_x2(k+1) - b_x2(k))/2;
end


N_total_x = sum(N_basis_x);  %total number of basis functions
N_total_y = sum(N_basis_y);  %total number of basis functions

[Nxx, NNxx] = size(b_x1);
[Nyy, NNyy] = size(b_x2);
periodx = b_x1(NNxx)-b_x1(1);
periody = b_x2(NNyy)-b_x2(1);

alpha_mm = zeros(2*N+1,1);
beta_mm = zeros(2*N+1,1);

for m=1:(2*N+1)
    alpha_mm(m) = (m-N-1)*2*pi/periodx;
    beta_mm(m) =  (m-N-1)*2*pi/periody;
end

alpha_p = - alpha_ref - alpha_mm;
beta_p =  - beta_ref - beta_mm;

N_total_x3 = N_total_x - N_intervals_x;
N_total_y3 = N_total_y - N_intervals_y;

int_RT_x = zeros(2*N+1, 2*N+1, N_total_x3, N_total_y);
int_RT_y = zeros(2*N+1, 2*N+1, N_total_x, N_total_y3);

for p=1:(2*N+1)
for q=1:(2*N+1)
for i1 = 1:(n_points_int-1)   %points of integration on the interval
for i2 = 1:(n_points_int-1)
    %integral for Ex
    for k1 = 1:N_intervals_x   %for each interval
    for k2 = 1:N_intervals_y    
        for j1 = (Nx(k1)+1):(Nx(k1)+nx(k1))
        for j2 = (Ny(k2)+(k2-1)+1):(Ny(k2)+(k2-1)+ny(k2)+1)
            num1 = j1 - Nx(k1);      %true number of Gegenbauer polynomial
            num2 = j2 - Ny(k2) - (k2-1);  
            int_RT_x(p, q, j1, j2) = int_RT_x(p, q, j1, j2) + C(num1,i1)*C(num2,i2)*...
                exp(1j*alpha_p(p)*x(k1,i1))*exp(1j*beta_p(q)*y(k2,i2))*...
                lx1(k1)*lx2(k2)*(ksi1(i1+1)-ksi1(i1))*(ksi1(i2+1)-ksi1(i2));
                %(x(k1,i1+1)-x(k1,i1))*(y(k2,i2+1)-y(k2,i2));     
        end
        end
        %integral for Ey
        for jj1 = (Nx(k1)+(k1-1)+1):(Nx(k1)+(k1-1)+nx(k1)+1)
        for jj2 = (Ny(k2)+1):(Ny(k2)+ny(k2))
            nnum1 = jj1 - Nx(k1) - (k1-1);
            nnum2 = jj2 - Ny(k2);      %true number of Gegenbauer polynomial
            int_RT_y(p, q, jj1, jj2) = int_RT_y(p,q,jj1,jj2) + C(nnum1,i1)*C(nnum2,i2)*...
                exp(1j*alpha_p(p)*x(k1,i1))*exp(1j*beta_p(q)*y(k2,i2))*...
                (x(k1,i1+1)-x(k1,i1))*(y(k2,i2+1)-y(k2,i2));             
        end
        end 
    end
    end
end
end
end
end

N_total = N_total_x*N_total_y;
N_total_3 = (N_total_x3)*(N_total_y3);

fx = zeros(2*N+1, 2*N+1, N_total_x3, N_total_y3);
fy = zeros(2*N+1, 2*N+1, N_total_x3, N_total_y3);
fx_coef = zeros(NN, N_total_3);
fy_coef = zeros(NN, N_total_3);

%[na_x,nna_x] = size(ax);
%[na_y,nna_y] = size(ay);

for p = 1:(2*N+1)
for q = 1:(2*N+1)     
    for i = 1:N_total_x3
    for j = 1:N_total_y3 
        
    %coefficients for Ex    
    for sy = 1:N_total_y                                 
            fx(p,q,i,j) = fx(p,q,i,j) + int_RT_x(p,q,i,sy)*ay(sy,j);
    end
    %coefficients for Ey
    for sx = 1:N_total_x
            fy(p,q,i,j) = fy(p,q,i,j) + int_RT_y(p,q,sx,j)*ax(sx,i);
    end
    
    row_FMM = q + (p-1)*(2*N+1);
    col_PMM = j + (i-1)*(N_total_y3);
    fx_coef(row_FMM, col_PMM) = fx(p,q,i,j);
    fy_coef(row_FMM, col_PMM) = fy(p,q,i,j);
    end
    end
end
end

%fx_coef = Kronecker_delta

norm = periodx*periody;

fx_coef = fx_coef/norm;
fy_coef = fy_coef/norm;

%{
for p = 1:(2*N+1)
for q = 1:(2*N+1)
for u = 1:N_intervals_x
for v = 1:N_intervals_y 
    %coefficients for Ex
    for i = 1:N_total_x3
    for j = 1:N_total_y3  
    for s = 1:N_total_y                      
            if ((Nx(u)+1)<=i)&&(i<=(Nx(u)+nx(u)))&&...
               ((Ny(v)+(v-1)+1)<=s)&&(s<=(Ny(v)+(v-1)+ny(v)+1))&&(ay(s,j)~=0)            
            fx(p,q,i,j) = fx(p,q,i,j) + int_RT_x(p,q,i,s)*ay(s,j);
            end
    end
    row_FMM = q + (p-1)*(2*N+1);
    col_PMM = j + (i-1)*(N_total_y3);
    fx_coef(row_FMM, col_PMM) = fx(p,q,i,j);
    end 
    end    
 
    %coefficients for Ey
    for i = 1:N_total_x3
    for j = 1:N_total_y3
    for s = 1:N_total_x
            if ((Nx(u)+(u-1)+1)<=s)&&(s<=(Nx(u)+(u-1)+nx(u)+1))&&...
               ((Ny(v)+1)<=j)&&(j<=(Ny(v)+ny(v)))&&(ax(s,i)~=0)
            fy(p,q,i,j) = fy(p,q,i,j) + int_RT_y(p,q,s,j)*ax(s,i);
            end
    end
    row_FMM = q + (p-1)*(2*N+1);
    col_PMM = j + (i-1)*(N_total_y3);
    fy_coef(row_FMM, col_PMM) = fy(p,q,i,j);
    end
    end
end
end
end
end
%}
%{
u1_PMM = ud_PMM(1:N_total_3);        %u_last
u2_PMM = ud_PMM((N_total_3+1):2*N_total_3);
d1_PMM = ud_PMM((2*N_total_3+1):3*N_total_3);  %d0
d2_PMM = ud_PMM((3*N_total_3+1):4*N_total_3);

u1_FMM = fx_coef*u1_PMM;
u2_FMM = fy_coef*u2_PMM;
d1_FMM = fx_coef*d1_PMM;
d2_FMM = fy_coef*d2_PMM;

ud_FMM = cat(1, u1_FMM, u2_FMM, d1_FMM, d2_FMM);
%}
