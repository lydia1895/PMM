function [T, aT, U] = PMM_T_matrices_new(Nx,nx,N_basis_x,ax,La,N_intervals_x)

N_total_x = sum(N_basis_x);
N_total_x3 = N_total_x - N_intervals_x;

Nmax = max(N_basis_x);   

p = zeros(Nmax,1);     
norm = zeros(Nmax,1);
for i=0:(Nmax-1)
    p(i+1) = gamma(i+2*La)/(gamma(2*La)*gamma(i+1));
    %p(i)=Ci i-th Gegenbauer polynomial at 1
    norm(i+1) = pi^0.5*p(i+1)*gamma(La+0.5)/(gamma(La)*(i+La));
    %<Cn,Cm> = delta(n,m)*norm(n)
end
hx = zeros(N_total_x3);

for k=1:N_intervals_x
    for i=(Nx(k)+1):(Nx(k)+nx(k))
    hx(i,i) = norm(i-Nx(k));
    end
end

U = zeros(N_intervals_x,N_total_x3,N_total_x3);
T = zeros(N_intervals_x,N_total_x3,N_total_x3);
Tnew = zeros(N_intervals_x,N_total_x3,N_total_x);
aT = zeros(N_intervals_x,N_total_x3,N_total_x3);

for j=1:N_intervals_x
    for i=(Nx(j)+1):(Nx(j)+nx(j))
    U(j,i,i) = hx(i,i);
    end
end

%integral for P 1st basis and P 3d basis
for k=1:N_intervals_x
    for m=Nx(k)+1:Nx(k)+nx(k)
        for i=Nx(k)+1:Nx(k)+nx(k)
            numm = m-Nx(k);
            numi = i-Nx(k);
            if (numi == numm-1)
                L = numm-1;
                T(k,m,i) = 2*L/((2*L-1)*(2*L+1));
            end
            if (numi == numm+1)
                L = numm-1;
                T(k,m,i) = 2*(L+1)/((2*L+1)*(2*L+3));
            end
        end
    end
end

%integral for P 3d basis and P 3d basis
for k=1:N_intervals_x
    for m=Nx(k)+1:Nx(k)+nx(k)
        for i=Nx(k)+(k-1)+1:Nx(k)+(k-1)+nx(k)+1
            numm = m-Nx(k);
            numi = i-Nx(k)-(k-1);
            if (numi == numm-1)
                L = numm-1;
                Tnew(k,m,i) = 2*L/((2*L-1)*(2*L+1));
            end
            if (numi == numm+1)
                L = numm-1;
                Tnew(k,m,i) = 2*(L+1)/((2*L+1)*(2*L+3));
            end
        end
    end
end

for i=1:N_intervals_x
    TTnew(:,:)=Tnew(i,:,:);
    aTT = TTnew*ax;
    aT(i,:,:) = aTT;
end
