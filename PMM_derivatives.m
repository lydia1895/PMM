function [Dx, hx, P_dP] = PMM_derivatives(La, N_intervals, N_basis, n, N, ax, b_x)
                                  
%function to get derivative matrix for one dimension (x or y)


%N_total is the total number of basis functions

N_total = sum(N_basis); 
hx = zeros(N_total - N_intervals, N_total - N_intervals);

%Nmax is maximum number of Gegenbauer polynomial on all intervals 

Nmax = max(N_basis);   

%to define h=<Pi,Pj> we should start from here

p = zeros(Nmax,1);     
norm = zeros(Nmax,1);
for i=0:(Nmax-1)
    p(i+1) = gamma(i+2*La)/(gamma(2*La)*gamma(i+1));
    %p(i)=Ci i-th Gegenbauer polynomial at 1
    norm(i+1) = pi^0.5*p(i+1)*gamma(La+0.5)/(gamma(La)*(i+La));
    %<Cn,Cm> = delta(n,m)*norm(n)
end

 
%for l we submit matched coordinates x=x1 or x2,
%because x(k) and x(k+1) must be constant
    
lx=zeros(N_intervals,1);
for k=1:N_intervals
    lx(k) = (b_x(k+1) - b_x(k))/2;
end

%define h=<Pi,Pj>, last step

for k=1:N_intervals
    for i=(N(k)+1):(N(k)+n(k))
    hx(i,i) = norm(i-N(k))*lx(k);
    end
end

%to define d=<Pi,d(Pj)/d(ksi)> we should start from here
%{
npoints_int = 600;
ksi = linspace(-1,1,npoints_int);
C_La = zeros(Nmax, npoints_int);
C_La1 = zeros(Nmax, npoints_int);

for m=1:Nmax
    for i = 1:npoints_int       
        C_La(m,i) = mfun('G',m-1,La,ksi(i));   %C(m,lambda,ksi)
        C_La1(m,i) = mfun('G',m-1,La+1,ksi(i));%C(m,lambda+1,ksi)
    end
end

%define <Cm, d(Cn)/d(ksi)> for all m,n

core_int_C_dC = zeros(Nmax, Nmax, npoints_int);
int_C_dC = zeros(Nmax, Nmax);
for m=1:Nmax
for r=1:Nmax 
        for i = 1:(npoints_int-1)          
            if r>=2
            core_int_C_dC(m,r,i) = ( (1-ksi(i)^2)^(La-0.5) )*C_La(m,i)*(2*La*C_La1(r-1,i));
            end
            if r==1
            core_int_C_dC(m,r,i) = 0;
            end
        int_C_dC(m,r) = int_C_dC(m,r) + core_int_C_dC(m,r,i)*(ksi(i+1)-ksi(i));  
        end      
end
end
%}

%{
old derivatives
for m=1:Nmax
for r=1:Nmax 
        for i = 1:(npoints_int-1)          
            if r>1
                core_int_C_dC(m,r) = ((1-ksi(i)^2)^(La-1.5))*C(m,i)*...
                (-r*ksi(i)*C(r,i)+(r+2*La-1)*C(r-1,i));
            else
                core_int_C_dC(m,r) = ((1-ksi(i)^2)^(La-1.5))*C(m,i)*...
                (-r*ksi(i)*C(r,i)+(r+2*La-1)*0); 
            %C(0)=1
            end
        int_C_dC(m,r) = int_C_dC(m,r) + core_int_C_dC(m,r)*(ksi(i+1)-ksi(i));    
        end      
end
end
%}

%define P_dP=<Pi, d(Pj)/d(ksi)>, last step
%{
P_dP = zeros(N_total, N_total);
for j=1:N_intervals
    for i1=(N(j)+(j-1)+1):(N(j)+(j-1)+n(j)+1)
    for i2=(N(j)+(j-1)+1):(N(j)+(j-1)+n(j)+1)
    %for i1=(N(j1)+(j1-1)+1):(N(j1)+(j1-1)+n(j1)+1)
    %for i2=(N(j2)+(j2-1)+1):(N(j2)+(j2-1)+n(j2)+1)           
        %because in <P(i1),dP(i2)> both P(i1) and P(i2)
        %must be from the same interval
        ii1 = i1 - N(j) -(j-1); 
        ii2 = i2 - N(j) -(j-1); 
        P_dP(i1,i2) = int_C_dC(ii1,ii2);
        %check out the order of i1, i2 !!!!!!!! 
    end
    end
end
%}
P_dP(1,1)=0;
P_dP(1,2)=2.26462;
P_dP(1,3)=0;
P_dP(1,4)=2.26462;

P_dP(2,1)=0;
P_dP(2,2)=0;
P_dP(2,3)=2.71754;
P_dP(2,4)=0;

P_dP(3,1)=0;
P_dP(3,2)=0;
P_dP(3,3)=0;
P_dP(3,4)=2.98929;

P_dP(4,1)=0;
P_dP(4,2)=0;
P_dP(4,3)=0;
P_dP(4,4)=0;
for i=1:4
    for j=1:4
    P_dP(i+4,j+4) = P_dP(i,j);
    end
end
    
%from presentation: dx = derx*Iy, dy = Ix*dery
%we will compute dx, dy in main program
%derx = hx\Dx
%Dx = ((a)T)*P_dP*a
%P_dP = transpose(P_dP)

Dx = transpose(ax)*P_dP*ax;

