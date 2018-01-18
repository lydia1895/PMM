function[nx,Nx,N_total_x,N_total_x3] = PMM_number_of_basis_functions(N_intervals_x,N_basis_x)
    
    nx = N_basis_x-1;   %because functions are p(0)...p(n(k)) on interval k
    Nx = zeros(N_intervals_x,1);
    for k=1:N_intervals_x
        Nx(k) = -nx(k);
        for p=1:k
            Nx(k) = Nx(k)+nx(p);
        end
    end
    N_total_x = sum(N_basis_x);  %total number of basis functions
    N_total_x3 = N_total_x - N_intervals_x;  %number of basis functions in "third" basis
