function [f_coef] = PMM_to_FMM_RT_inverse(N, NN, La, alpha_ref, beta_ref,...
    b_x1, b_x2, N_intervals_x, N_intervals_y, N_basis_x, N_basis_y, Nx, nx, Ny, ny)

n_points_int = 300;
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
%{
for m=1:Nmax
    for i = 1:n_points_int
        C(m,i) = mfun('G',m-1,La,ksi1(i));
    end    
    %C(m,:) = gegenbauerC(m-1,La,ksi1);
end
%}

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
%bb_x2 = b_x2
%llx2 = lx2

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

alpha_m =  alpha_mm + alpha_ref ;
beta_m =  beta_mm + beta_ref; 

N_total_x3 = N_total_x - N_intervals_x;
N_total_y3 = N_total_y - N_intervals_y;

%int_RT_x = zeros(2*N+1, 2*N+1, N_total_x3, N_total_y);
%int_RT_y = zeros(2*N+1, 2*N+1, N_total_x, N_total_y3);

int_x = zeros(2*N+1, N_total_x3);
int_y = zeros(2*N+1, N_total_y3);

for m=1:(2*N+1)
    for u = 1:N_intervals_x
        for k = (Nx(u)+1):(Nx(u)+nx(u))
            num1 = k - Nx(u);
            a = alpha_m(m)*lx1(u);
            integral_exp_Gegenbauer = int_exp_gegenbauer(a,num1-1,La);
            int_x(m,k) = exp(1j*alpha_m(m)*sumx1(u))*integral_exp_Gegenbauer;
        end
    end
    for v = 1:N_intervals_y
        for l = (Ny(v)+1):(Ny(v)+ny(v))
            num2 = l - Ny(v);
            b = beta_m(m)*lx2(v);
            integral_exp_Gegenbauer = int_exp_gegenbauer(b,num2-1,La);
            int_y(m,l) = exp(1j*beta_m(m)*sumx2(v))*integral_exp_Gegenbauer;
        end
    end
end       

%to define h=<Pi,Pj> we should start from here

p = zeros(Nmax,1);     
norm = zeros(Nmax,1);
for i=0:(Nmax-1)
    p(i+1) = gamma(i+2*La)/(gamma(2*La)*gamma(i+1));
    %p(i)=Ci i-th Gegenbauer polynomial at 1
    norm(i+1) = pi^0.5*p(i+1)*gamma(La+0.5)/(gamma(La)*(i+La));
    %<Cn,Cm> = delta(n,m)*norm(n)
end

%define h=<Pi,Pj>, integral over (1,-1), without lx, ly!
hx = zeros(N_total_x3, N_total_x3);
hy = zeros(N_total_y3, N_total_y3);
for k=1:N_intervals_x
    for i=(Nx(k)+1):(Nx(k)+nx(k))
    hx(i,i) = norm(i-Nx(k));   %unitless norm <Pi,Pi>
    end
end
for k=1:N_intervals_y
    for i=(Ny(k)+1):(Ny(k)+ny(k))
    hy(i,i) = norm(i-Ny(k));
    end
end


N_total_3 = (N_total_x3)*(N_total_y3);

int_xy = zeros(2*N+1, 2*N+1, N_total_x3, N_total_y3);
f_coef = zeros(N_total_3, NN);


for p=1:(2*N+1)
for q=1:(2*N+1)
    for i=1:N_total_x3
    for j=1:N_total_y3
        int_xy(p,q,i,j) = int_x(p,i)*int_y(q,j); 
        col_FMM = q + (p-1)*(2*N+1);
        row_PMM = j + (i-1)*(N_total_y3);
        f_coef(row_PMM, col_FMM) = int_xy(p,q,i,j)/(hx(i,i)*hy(j,j));      
    end
    end
end
end


%{
u1_PMM = ud_PMM(1:N_total_3);        %u_last
u2_PMM = ud_PMM((N_total_3+1):2*N_total_3);
d1_PMM = ud_PMM((2*N_total_3+1):3*N_total_3);  %d0
d2_PMM = ud_PMM((3*N_total_3+1):4*N_total_3);

hfull = Kronecker_product(hx,hy);
[nh,nnh] = size(hfull)
[nu,nnu] = size(u1_PMM)

u1_PMM_c = hfull*u1_PMM;
u2_PMM_c = hfull*u2_PMM;
d1_PMM_c = hfull*d1_PMM;
d2_PMM_c = hfull*d2_PMM;

u1_FMM = f_coef\u1_PMM_c;
u2_FMM = f_coef\u2_PMM_c;
d1_FMM = f_coef\d1_PMM_c;
d2_FMM = f_coef\d2_PMM_c;

ud_FMM = cat(1, u1_FMM, u2_FMM, d1_FMM, d2_FMM);
%}
