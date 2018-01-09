function [gamma] = PMM_gamma(La, tau_x, tau_y, x, y, N_intervals_x, N_intervals_y, N_basis_x, N_basis_y)

%first compute a from boundary conditions
%continuity and periodicity conditions are handled by sparse [a]

ax = PMM_boundary_conditions(La, tau_x, N_intervals_x, N_basis_x);
ay = PMM_boundary_conditions(La, tau_y, N_intervals_y, N_basis_y);

%compute derivative matrices

derx = PMM_derivatives(La, N_intervals_x, N_basis_x, ax, x);
dery = PMM_derivatives(La, N_intervals_y, N_basis_y, ay, y);

N_total_x = sum(N_basis_x);
N_total_y = sum(N_basis_y);

Ix = eye(N_total_x, N_total_x);
Iy = eye(N_total_y, N_total_y);

dx = Kronecker_product(derx, Iy);
dy = Kronecker_product(Ix, dery);

N_xy = N_total_x*N_total*y;
miden = eye(N_xy, N_xy);

L_EH_11 = -dx/eps*dy;
L_EH_12 = dx/eps*dx + k^2*miden;
L_EH_21 = -dy/eps*dy - k^2*miden;
L_EH_22 = dy/eps*dx;
L_EH_up = cat(2, L_EH_11, L_EH_12);
L_EH_down = cat(2, L_EH_21, L_EH_22);
L_EH = cat(1, L_EH_up, L_EH_down);
L_EH = L_EH/(k^2);


L_HE_11 = -dx*dy;
L_HE_12 = dx*dx + (k^2)*eps;
L_HE_21 = -dy*dy - (k^2)*eps;
L_HE_22 = dy*dx;
L_HE_up = cat(2, L_HE_11, L_HE_12);
L_HE_down = cat(2, L_HE_21, L_HE_22);
L_HE = cat(1, L_HE_up, L_HE_down);
L_HE = L_HE/(k^2);

M = L_EH*L_HE;

[E, gamma2] = eig(M);
gamma = sqrt(-gamma2);
coef = (k/(w*eps0))*gamma;

H = -coef\L_HE*E;

