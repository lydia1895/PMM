function[int_P1_Q1, int_P1_Q2] = incident_integral(La, b_x1, b_x2, alpha_ref, beta_ref,...
    N_basis_x, N_basis_y, N_intervals_x, N_intervals_y, Nx, Ny, nx, ny)


n_points_int = 100;
ksi1 = linspace(-1,1,n_points_int);
ksi2 = linspace(-1,1,n_points_int);

Nmax_x = max(N_basis_x);   %maximum number of Gegenbauer polynomial on all intervals
Nmax_y = max(N_basis_y);   %maximum number of Gegenbauer polynomial on all intervals
Nmax = max(Nmax_x, Nmax_y);
C = zeros(Nmax, n_points_int);
%{
for m=1:Nmax
    for i = 1:n_points_int
        C(m,i) = mfun('G',m-1,La,ksi1(i));
        %C(m,i) = gegenbauerC(m-1,La,ksi1(i));
    end
end
%}

%that's where we submit matched coordinates,
%because x(k) and x(k+1) must be constant

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
N_total_x3 = N_total_x - N_intervals_x;
N_total_y3 = N_total_y - N_intervals_y;

int_inc = zeros(N_total_x3, N_total_y3);

int_x = zeros(N_total_x3,1);
int_y = zeros(N_total_y3,1);

for u = 1:N_intervals_x
        for k = (Nx(u)+1):(Nx(u)+nx(u))
            num1 = k - Nx(u);
            a = alpha_ref*lx1(u);
            integral_exp_Gegenbauer1 = int_exp_gegenbauer(a,num1-1,La);
            int_x(k) = exp(1j*alpha_ref*sumx1(u))*integral_exp_Gegenbauer1;
        end
end
for v = 1:N_intervals_y
        for l = (Ny(v)+1):(Ny(v)+ny(v))
            num2 = l - Ny(v);
            b = beta_ref*lx2(v);
            integral_exp_Gegenbauer2 = int_exp_gegenbauer(b,num2-1,La);
            int_y(l) = exp(1j*beta_ref*sumx2(v))*integral_exp_Gegenbauer2;
        end
end

p = zeros(Nmax,1);     
norm = zeros(Nmax,1);
for i=0:(Nmax-1)
    p(i+1) = gamma(i+2*La)/(gamma(2*La)*gamma(i+1));
    %p(i)=Ci i-th Gegenbauer polynomial at 1
    norm(i+1) = pi^0.5*p(i+1)*gamma(La+0.5)/(gamma(La)*(i+La));
    %<Cn,Cm> = delta(n,m)*norm(n)
end
%define h=<Pi,Pj>, integral over (1,-1), without lx, ly!
hx = zeros(N_total_x3);
hy = zeros(N_total_y3);
for k=1:N_intervals_x
    for i=(Nx(k)+1):(Nx(k)+nx(k))
    hx(i,i) = norm(i-Nx(k));
    end
end
for k=1:N_intervals_y
    for i=(Ny(k)+1):(Ny(k)+ny(k))
    hy(i,i) = norm(i-Ny(k));
    end
end
%define h=<Pi,Pj>, integral over (1,-1), without lx, ly!

int_P1_Q1 = int_x(1)*int_y(1)/(hx(1,1)*hy(1,1));
int_P1_Q2 = int_x(1)*int_y(2)/(hx(1,1)*hy(2,2));
