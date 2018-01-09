function [g_sqrt, g_down11, g_down22, g_down12] = ellipse_metric_integral_polyfit(ellipse_parameters,n_points)

%{
x = linspace(1,2,20);
y = 2*x.^2+1;
polyfit(x,y,2)

y = linspace(1,2,20);
z = x.^2+3*y;
fit([x, y],z,'poly23')
%}
R1 = ellipse_parameters(1);
R2 = ellipse_parameters(2);
P1 = ellipse_parameters(3);
P2 = ellipse_parameters(4);
Q1 = ellipse_parameters(5);
Q2 = ellipse_parameters(6);

N_intervals_x = 3;
N_intervals_y = 3;

x1_t_plus =  P1/2+Q1;
x1_t_minus = P1/2-Q1;
x2_t_plus =  P2/2+Q2;
x2_t_minus = P2/2-Q2;

b_x1 = [0 x1_t_minus x1_t_plus P1];
b_x2 = [0 x2_t_minus x2_t_plus P2];


ksi1 = linspace(-1,1,n_points);
ksi2 = linspace(-1,1,n_points);
x1 = zeros(N_intervals_x,n_points);
x2 = zeros(N_intervals_y,n_points);
for k1 = 1:N_intervals_x
    x1(k1,:) = ( (b_x1(k1+1)-b_x1(k1))*ksi1+(b_x1(k1+1)+b_x1(k1)) )/2;
end
for k2 = 1:N_intervals_y     
    x2(k2,:) = ( (b_x2(k2+1)-b_x2(k2))*ksi2+(b_x2(k2+1)+b_x2(k2)) )/2;
end

g_sqrt = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
g_sqrt(1,1,:,:) = ones(n_points, n_points);
g_sqrt(1,3,:,:) = ones(n_points, n_points);
g_sqrt(3,1,:,:) = ones(n_points, n_points);
g_sqrt(3,3,:,:) = ones(n_points, n_points);

g_down11 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
g_down11(1,1,:,:) = ones(n_points, n_points);
g_down11(1,3,:,:) = ones(n_points, n_points);
g_down11(3,1,:,:) = ones(n_points, n_points);
g_down11(3,3,:,:) = ones(n_points, n_points);

g_down22 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
g_down22(1,1,:,:) = ones(n_points, n_points);
g_down22(1,3,:,:) = ones(n_points, n_points);
g_down22(3,1,:,:) = ones(n_points, n_points);
g_down22(3,3,:,:) = ones(n_points, n_points);

g_down12 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
g_down12(1,1,:,:) = zeros(n_points, n_points);
g_down12(1,3,:,:) = zeros(n_points, n_points);
g_down12(3,1,:,:) = zeros(n_points, n_points);
g_down12(3,3,:,:) = zeros(n_points, n_points);

for nx = 1:n_points
    for ny=1:n_points
        
g_sqrt(1,2,nx,ny) = P1/2-(R1/R2)*(R2^2-(x2(2,ny)-P2/2)^2)^0.5/x1_t_minus;
g_sqrt(2,1,nx,ny) = P2/2-(R2/R1)*(R1^2-(x1(2,nx)-P1/2)^2)^0.5/x2_t_minus;

g_sqrt(2,2,nx,ny) = ( 1/(x1_t_plus-x1_t_minus)*(x2_t_plus-x2_t_minus) )*...
    ( 4*(R1^2-(x1(2,nx)-P1/2)^2)^0.5 * (R2^2-(x2(2,ny)-P2/2)^2)^0.5 - ...
    (2*x1(2,nx)-x1_t_plus-x1_t_minus)*(2*x2(2,ny)-x2_t_plus-x2_t_minus)*...
    (x1(2,nx)-P1/2)*(x2(2,ny)-P2/2)/((R1^2-(x1(2,nx)-P1/2)^2)^0.5*(R2^2-(x2(2,ny)-P2/2)^2)^0.5) );

g_sqrt(2,3,nx,ny) = ( (R2/R1)*(R1^2-(x1(2,nx)-P1/2)^2)^0.5 - P2/2 )/(x2_t_plus - P2);
g_sqrt(3,2,nx,ny) = ( (R1/R2)*(R2^2-(x2(2,ny)-P2/2)^2)^0.5 - P1/2 )/(x1_t_plus - P1);

g_down11(1,2,nx,ny) = (P1/2 - (R1/R2)*(R2^2-(x2(2,ny)-P2/2)^2)^0.5)^2/x1_t_minus^2;
g_down11(2,1,nx,ny) = 1 + ( (R2/R1)*(x1(2,nx)-P1/2)/(R1^2-(x1(2,nx)-P1/2)^2)^0.5)^2*...
    (x2(1,ny)/x2_t_minus)^2;

g_down11(2,2,nx,ny) = (1/(x1_t_plus-x1_t_minus)^2)*4*(R1/R2)^2*(R2^2-(x2(2,ny)-P2/2)^2) + ...
    ( (2*x2(2,ny)-x2_t_plus-x2_t_minus)/(x2_t_plus-x2_t_minus)*(R2/R1)*(x1(2,nx)-P1/2) )^2/...
    (R1^2-(x1(2,nx)-P1/2)^2);

g_down11(2,3,nx,ny) = 1 + ( (R2/R1)*(x2(3,ny)-P2)/(x2_t_plus-P2)*(x1(2,nx)-P1/2) )^2/...
    (R1^2-(x1(2,nx)-P1/2)^2);
g_down11(3,2,nx,ny) = ( ((R1/R2)*(R2^2-(x2(2,ny)-P2/2)^2)^0.5 - P1/2)/(x1_t_plus-P1) )^2;

g_down22(1,2,nx,ny) = ( R1/R2*(x2(2,ny)-P2/2)*x1(1,nx)/x1_t_plus )^2/(R2^2-(x2(2,ny)-P2/2)^2)+1;
g_down22(2,1,nx,ny) = (P2/2 - (R2/R1)*(R1^2-(x1(2,nx)-P1/2)^2)^0.5)^2/x2_t_minus^2;

g_down22(2,2,nx,ny) = ((2*x1(2,nx)-x1_t_plus-x1_t_minus)/(x1_t_plus-x1_t_minus))^2*...
    ((R1/R2)*(x2(2,ny)-P2/2))^2/(R2^2-(x2(2,ny)-P2/2)^2);

g_down22(2,3,nx,ny) = ((R2/R1)*(R1^2-(x1(2,nx)-P1/2)^2)^0.5 - P2/2)^2/(x2_t_plus-P2)^2;
g_down22(3,2,nx,ny) = ( ((x1(3,nx)-P1)/(x1_t_plus-P1))*( R1/R2*(x2(2,ny)-P2/2)) )^2/...
    (R2^2-(x2(2,ny)-P2/2)^2);
    
g_down12(1,2,nx,ny) = (P1/2 - (R1/R2)*(R2^2-(x2(2,ny)-P2/2)^2)^0.5)/x1_t_minus *...
    (R1/R2)*(x2(2,ny)-P2/2)/(R2^2-(x2(2,ny)-P2/2)^2)^0.5*x1(1,nx)/x1_t_minus;
g_down12(2,1,nx,ny) = (x2(1,ny)/x2_t_minus)*(R2/R1)*(x1(2,nx)-P1/2)/(R1^2-(x1(2,nx)-P1/2)^2)^0.5*...
    (P2/2-(R2/R1)*(R1^2-(x1(2,nx)-P1/2)^2)^0.5)/x2_t_minus;

g_down12(2,2,nx,ny) = 2*(R1/R2)^2/(x1_t_plus-x1_t_minus)^2 * ...
    (2*x1(2,nx)-x1_t_plus-x1_t_minus)*(-(x2(2,ny)-P2/2)) + ...
    (2*x2(2,ny)-x2_t_plus-x2_t_minus)/(x2_t_plus-x2_t_minus)^2 * ...
    2*(R2/R1)^2*(-(x1(2,nx)-P1/2));

g_down12(2,3,nx,ny) = (x2(3,ny)-P2)/(x2_t_plus-P2)^2*(R2/R1)*(-(x1(2,nx)-P1/2))/...
    (R1^2-(x1(2,nx)-P1/2)^2)^0.5*((R2/R1)*(R1^2-(x1(2,nx)-P1/2)^2)^0.5-P2/2);
g_down12(3,2,nx,ny) = ((R1/R2)*(R2^2-(x2(2,ny)-P2/2)^2)^0.5 - P1/2)*(x1(3,nx)-P1)/...
    (x1_t_plus-P1)^2*(R1/R2)*(-(x2-P2/2))/(R2^2-(x2(2,ny)-P2/2)^2)^0.5;


    end
end

g_sqrt = abs(g_sqrt);


