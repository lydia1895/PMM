function [fx_coef,fy_coef] = PMM_to_FMM_RT_La05_one_integral(N, NN, La, alpha_ref, beta_ref,... 
   b_x1, b_x2, N_intervals_x, N_intervals_y, N_basis_x, N_basis_y, Nx, nx, Ny, ny,...
   ax, ay)
%{
n_points_int = 10000;
ksi = linspace(-1,1,n_points_int);
x1 = zeros(N_intervals_x,n_points_int);
x2 = zeros(N_intervals_y,n_points_int);
for k1 = 1:N_intervals_x
    x1(k1,:) = ( (b_x1(k1+1)-b_x1(k1))*ksi+(b_x1(k1+1)+b_x1(k1)) )/2;
end
for k2 = 1:N_intervals_y     
    x2(k2,:) = ( (b_x2(k2+1)-b_x2(k2))*ksi+(b_x2(k2+1)+b_x2(k2)) )/2;
end

%x(k1) = PMM_matched_coordinates_ellipse(x1(k1), x2(k2), P1, P2, R1, R2)
%y(k2) = PMM_matched_coordinates_ellipse(x1(k1), x2(k2), P1, P2, R1, R2)
x = x1;
y = x2;
C = zeros(Nmax, n_points_int);
for m=1:Nmax
    %for i = 1:n_points_int
    %    C(m,i) = mfun('G',m-1,La,ksi1(i));
    %end    
    C(m,:) = gegenbauerC(m-1,La,ksi);
end

%}
Nmax_x = max(N_basis_x);   %maximum number of Gegenbauer polynomial on all intervals
Nmax_y = max(N_basis_y);   %maximum number of Gegenbauer polynomial on all intervals
Nmax = max(Nmax_x, Nmax_y);


lx1=zeros(N_intervals_x,1);
sumx1=zeros(N_intervals_x,1);
for k=1:N_intervals_x
    lx1(k) = (b_x1(k+1) - b_x1(k))/2;
    sumx1(k) = (b_x1(k+1) + b_x1(k))/2;

end
lx2=zeros(N_intervals_y,1);
sumx2=zeros(N_intervals_y,1);
for k=1:N_intervals_y
    lx2(k) = (b_x2(k+1) - b_x2(k))/2;
    sumx2(k) = (b_x2(k+1) + b_x2(k))/2;
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

int_RT_Ex_P = zeros(N_total_x3, 2*N+1);
int_RT_Ex_Q = zeros(N_total_y,  2*N+1);
int_RT_Ey_P = zeros(N_total_x,  2*N+1);
int_RT_Ey_Q = zeros(N_total_y3, 2*N+1);



%for i = 1:(n_points_int-1)   %points of integration on the interval
for p=1:(2*N+1)

for k1 = 1:N_intervals_x   %for each interval  
    a = alpha_p(p)*lx1(k1);
    %integral for Ex, P
    for j1 = (Nx(k1)+1):(Nx(k1)+nx(k1))
        num1 = j1 - Nx(k1);      %true number of Gegenbauer polynomial
        integral_exp_Gegenbauer1 = int_exp_gegenbauer(a,num1-1,0.5);
        int_RT_Ex_P(j1, p) = exp(1j*alpha_p(p)*sumx1(k1))*...
            integral_exp_Gegenbauer1*lx1(k1); 
        %int_RT_Ex_P(j1, p) = int_RT_Ex_P(j1, p)+...
        %        C(num1,i)*exp(1j*alpha_p(p)*x(k1,i))*(x(k1,i+1)-x(k1,i));    
    end
    %integral for Ey, P
    for jj1 = (Nx(k1)+(k1-1)+1):(Nx(k1)+(k1-1)+nx(k1)+1)
         nnum1 = jj1 - Nx(k1) - (k1-1);
         integral_exp_Gegenbauer2 = int_exp_gegenbauer(a,nnum1-1,0.5);
         int_RT_Ey_P(jj1, p) = exp(1j*alpha_p(p)*sumx1(k1))*...
             integral_exp_Gegenbauer2*lx1(k1);
         %int_RT_Ey_P(jj1, p) = int_RT_Ey_P(jj1, p) +...
         %        C(nnum1,i)*exp(1j*alpha_p(p)*x(k1,i))*(x(k1,i+1)-x(k1,i)); 
    end
end

for k2 = 1:N_intervals_y  
    b = beta_p(p)*lx2(k2);
    %integral for Ex, Q      
    for j2 = (Ny(k2)+(k2-1)+1):(Ny(k2)+(k2-1)+ny(k2)+1)
         num2 = j2 - Ny(k2) - (k2-1);
         integral_exp_Gegenbauer3 = int_exp_gegenbauer(b,num2-1,0.5);
         int_RT_Ex_Q(j2, p) = exp(1j*beta_p(p)*sumx2(k2))*...
             integral_exp_Gegenbauer3*lx2(k2);  
         %int_RT_Ex_Q(j2, p) = int_RT_Ex_Q(j2, p)+...    
         %        C(num2,i)*exp(1j*beta_p(p)*y(k2,i))*(y(k2,i+1)-y(k2,i)); 
    end    
    %integral for Ey, Q
    for jj2 = (Ny(k2)+1):(Ny(k2)+ny(k2))
          nnum2 = jj2 - Ny(k2);      %true number of Gegenbauer polynomial
          integral_exp_Gegenbauer4 = int_exp_gegenbauer(b,nnum2-1,0.5);
          int_RT_Ey_Q(jj2, p) = exp(1j*beta_p(p)*sumx2(k2))*...
              integral_exp_Gegenbauer4*lx2(k2);  
          %int_RT_Ey_Q(jj2, p) = int_RT_Ey_Q(jj2, p) +...
          %       C(nnum2,i)*exp(1j*beta_p(p)*y(k2,i))*(y(k2,i+1)-y(k2,i));             
    end
end 

end
%end

N_total_3 = (N_total_x3)*(N_total_y3);

int_RT_Ex_Qnew = zeros(N_total_y3,2*N+1);
int_RT_Ey_Pnew = zeros(N_total_x3,2*N+1);


for p = 1:(2*N+1)          
    %coefficients for Ex, Qnew
    for i = 1:N_total_y3
    for sy = 1:N_total_y               
            int_RT_Ex_Qnew(i,p) = int_RT_Ex_Qnew(i,p) + int_RT_Ex_Q(sy,p)*ay(sy,i);
    end
    end
    
    %coefficients for Ey, Pnew
    for j = 1:N_total_x3 
    for sx = 1:N_total_x
            int_RT_Ey_Pnew(j,p) = int_RT_Ey_Pnew(j,p) + int_RT_Ey_P(sx,p)*ax(sx,j);
    end
    end
end

norm = periodx*periody;
int_RT_Ex_P = transpose(int_RT_Ex_P);
int_RT_Ex_Qnew = transpose(int_RT_Ex_Qnew);
int_RT_Ey_Pnew = transpose(int_RT_Ey_Pnew);
int_RT_Ey_Q = transpose(int_RT_Ey_Q);


fx_coef = (Kronecker_product(int_RT_Ex_P, int_RT_Ex_Qnew))/norm;
fy_coef = (Kronecker_product(int_RT_Ey_Pnew, int_RT_Ey_Q))/norm;

%MI = blkdiag(fx_coef,fy_coef,fx_coef,fy_coef);
