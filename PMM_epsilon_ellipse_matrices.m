function [eps_total, mu_total] =...
        PMM_epsilon_ellipse_matrices(N_basis_x,N_basis_y,Nx,nx,Ny,ny,...
        N_intervals_x,N_intervals_y,La,epsilon,...
        int_g, int_for_ellipse)
    
    int_Ez_sqrt_g_full = int_g(:,:,1);
    int_Dz_unity_full =  int_g(:,:,2);
    int_Dx_sqrt_g_full = int_g(:,:,3);
    int_Dy_sqrt_g_full = int_g(:,:,4);
    int_Ex_g_down22_full = int_g(:,:,5);
    int_Ey_g_down12_full = int_g(:,:,6);
    int_Ex_g_down21_full = int_g(:,:,7);
    int_Ey_g_down11_full = int_g(:,:,8);
    
    int_Ex_g_down22_ellipse = int_for_ellipse(:,:,1);
    int_Ey_g_down12_ellipse = int_for_ellipse(:,:,2);
    int_Ex_g_down21_ellipse = int_for_ellipse(:,:,3);
    int_Ey_g_down11_ellipse = int_for_ellipse(:,:,4);
    int_Ez_sqrt_g_ellipse = int_for_ellipse(:,:,5);
    
    N_total_x = sum(N_basis_x);  %total number of basis functions
    N_total_y = sum(N_basis_y);  %total number of basis functions
    
    N_total_x3 = N_total_x - N_intervals_x;
    N_total_y3 = N_total_y - N_intervals_y;
    N_total_3 = N_total_x3*N_total_y3;
    
    Nmax_x = max(N_basis_x);
    Nmax_y = max(N_basis_y);
    Nmax = max(Nmax_x, Nmax_y);
    
    p = zeros(Nmax,1);
    norm = zeros(Nmax,1);
    for i=0:(Nmax-1)
        p(i+1) = gamma(i+2*La)/(gamma(2*La)*gamma(i+1));
        %p(i)=Ci i-th Gegenbauer polynomial at 1
        norm(i+1) = pi^0.5*p(i+1)*gamma(La+0.5)/(gamma(La)*(i+La));
        %<Cn,Cm> = delta(n,m)*norm(n)
    end
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
    eps_out = epsilon(1);
    eps_in =  epsilon(2);
    mu_in = 1.0;
    mu_out = 1.0;
    
    eps_int_Ex_g_down22=eps_out*int_Ex_g_down22_full;
    eps_int_Ey_g_down12=eps_out*int_Ey_g_down12_full;
    eps_int_Ex_g_down21=eps_out*int_Ex_g_down21_full;
    eps_int_Ey_g_down11=eps_out*int_Ey_g_down11_full;
    
    mu_int_Ex_g_down22=mu_out*int_Ex_g_down22_full;
    mu_int_Ey_g_down12=mu_out*int_Ey_g_down12_full;
    mu_int_Ex_g_down21=mu_out*int_Ex_g_down21_full;
    mu_int_Ey_g_down11=mu_out*int_Ey_g_down11_full;
    
    eps_int_Ez_g_sqrt=eps_out*int_Ez_sqrt_g_full;
    mu_int_Ez_g_sqrt=mu_out*int_Ez_sqrt_g_full;
    
    
    if (eps_in~=eps_out)
        
        eps_int_Ex_g_down22=eps_int_Ex_g_down22+(eps_in-eps_out)*int_Ex_g_down22_ellipse;
        eps_int_Ey_g_down12=eps_int_Ey_g_down12+(eps_in-eps_out)*int_Ey_g_down12_ellipse;
        eps_int_Ex_g_down21=eps_int_Ex_g_down21+(eps_in-eps_out)*int_Ex_g_down21_ellipse;
        eps_int_Ey_g_down11=eps_int_Ey_g_down11+(eps_in-eps_out)*int_Ey_g_down11_ellipse;
        eps_int_Ez_g_sqrt=eps_int_Ez_g_sqrt+(eps_in-eps_out)*int_Ez_sqrt_g_ellipse;
    end
    if (mu_in~=mu_out)
        mu_int_Ex_g_down22=mu_int_Ex_g_down22+(mu_in-mu_out)*int_Ex_g_down22_ellipse;
        mu_int_Ey_g_down12=mu_int_Ey_g_down12+(mu_in-mu_out)*int_Ey_g_down12_ellipse;
        mu_int_Ex_g_down21=mu_int_Ex_g_down21+(mu_in-mu_out)*int_Ex_g_down21_ellipse;
        mu_int_Ey_g_down11=mu_int_Ey_g_down11+(mu_in-mu_out)*int_Ey_g_down11_ellipse;
        mu_int_Ez_g_sqrt=mu_int_Ez_g_sqrt+(mu_in-mu_out)*int_Ez_sqrt_g_ellipse;
    end
    
    %{

k1 = 2;
k2 = 2;

for m1=Nx(k1)+1:Nx(k1)+nx(k1)
for m2=Ny(k2)+1:Ny(k2)+ny(k2)
    row = m2 + (m1-1)*N_total_y3;
    for i1=1:N_total_x3
    for i2=1:N_total_y3
        %%%%%%%%%%??????????????????????
    %for i1=Nx(k1)+1:Nx(k1)+nx(k1)
    %for i2=Ny(k2)+1:Ny(k2)+ny(k2)
        col = i2 + (i1-1)*N_total_y3;
        eps_int_Ex_g_down22(row,col)=eps_in*int_Ex_g_down22_full(row,col);
        eps_int_Ey_g_down12(row,col)=eps_in*int_Ey_g_down12_full(row,col);
        eps_int_Ex_g_down21(row,col)=eps_in*int_Ex_g_down21_full(row,col);
        eps_int_Ey_g_down11(row,col)=eps_in*int_Ey_g_down11_full(row,col);
        
        mu_int_Ex_g_down22(row,col) =mu_in*int_Ex_g_down22_full(row,col);
        mu_int_Ey_g_down12(row,col) =mu_in*int_Ey_g_down12_full(row,col);
        mu_int_Ex_g_down21(row,col) =mu_in*int_Ex_g_down21_full(row,col);
        mu_int_Ey_g_down11(row,col) =mu_in*int_Ey_g_down11_full(row,col);
        
        eps_int_Ez_g_sqrt(row,col)=eps_in*int_Ez_sqrt_g_full(row,col);
        mu_int_Ez_g_sqrt(row,col)=mu_in*int_Ez_sqrt_g_full(row,col);
 
    end
    end
end
end

end
    %}
    
    unity = eye(N_total_3,N_total_3);
    
    inv_int_Dx_sqrt_g_full = unity/int_Dx_sqrt_g_full;
    inv_int_Dy_sqrt_g_full = unity/int_Dy_sqrt_g_full;
    inv_int_Dz_unity_full = unity/int_Dz_unity_full;
    
    
    eps_xx = inv_int_Dx_sqrt_g_full*eps_int_Ex_g_down22;
    eps_xy = -inv_int_Dx_sqrt_g_full*eps_int_Ey_g_down12;
    eps_yx = -inv_int_Dy_sqrt_g_full*eps_int_Ex_g_down21;
    eps_yy = inv_int_Dy_sqrt_g_full*eps_int_Ey_g_down11;
    eps_zz = inv_int_Dz_unity_full*eps_int_Ez_g_sqrt;
   
    
    mu_xx = inv_int_Dx_sqrt_g_full*mu_int_Ex_g_down22;
    mu_xy = -inv_int_Dx_sqrt_g_full*mu_int_Ey_g_down12;
    mu_yx = -inv_int_Dy_sqrt_g_full*mu_int_Ex_g_down21;
    mu_yy = inv_int_Dy_sqrt_g_full*mu_int_Ey_g_down11;
    mu_zz = inv_int_Dz_unity_full*mu_int_Ez_g_sqrt;
    %{
    eps_xx = int_Dx_sqrt_g_full\eps_int_Ex_g_down22;
    eps_xy = -int_Dx_sqrt_g_full\eps_int_Ey_g_down12;
    eps_yx = -int_Dy_sqrt_g_full\eps_int_Ex_g_down21;
    eps_yy = int_Dy_sqrt_g_full\eps_int_Ey_g_down11;
    eps_zz = int_Dz_unity_full\eps_int_Ez_g_sqrt;
   
    
    mu_xx = int_Dx_sqrt_g_full\mu_int_Ex_g_down22;
    mu_xy = -int_Dx_sqrt_g_full\mu_int_Ey_g_down12;
    mu_yx = -int_Dy_sqrt_g_full\mu_int_Ex_g_down21;
    mu_yy = int_Dy_sqrt_g_full\mu_int_Ey_g_down11;
    mu_zz = int_Dz_unity_full\mu_int_Ez_g_sqrt;
   %}
    
    eps_total(:,:,1) = eps_xx;
    eps_total(:,:,2) = eps_xy;
    eps_total(:,:,3) = eps_yx;
    eps_total(:,:,4) = eps_yy;
    eps_total(:,:,5) = unity/eps_zz;
    
    mu_total(:,:,1) = mu_xx;
    mu_total(:,:,2) = mu_xy;
    mu_total(:,:,3) = mu_yx;
    mu_total(:,:,4) = mu_yy;
    mu_total(:,:,5) = unity/mu_zz;
