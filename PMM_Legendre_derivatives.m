function [Dx, hx] = PMM_Legendre_derivatives(La, N_intervals, N_basis, n, N, ax, b_x)
    
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
    
    %d(m,nn) = <C(m)d(C(n))/dx>
    d = zeros(Nmax, Nmax);
    
    for mm=1:Nmax
        for nn=(mm+1):Nmax
            if mod((mm-nn),2)~=0
                d(mm,nn) = 2.0;
            end
        end
    end
    
    P_dP = zeros(N_total, N_total);
    for j=1:N_intervals
        for i1=(N(j)+(j-1)+1):(N(j)+(j-1)+n(j)+1)
            for i2=(N(j)+(j-1)+1):(N(j)+(j-1)+n(j)+1)
                %because in <P(i1),dP(i2)> both P(i1) and P(i2)
                %must be from the same interval
                ii1 = i1 - N(j) -(j-1);
                ii2 = i2 - N(j) -(j-1);
                P_dP(i1,i2) = d(ii1,ii2);
                %check out the order of i1, i2 !!!!!!!!
            end
        end
    end
    
    Dx = transpose(conj(ax))*P_dP*ax;
    
    
    
    
    
    
    
    
    
    
    
    
