function [dx_x1,dx_x2,dy_x1,dy_x2] =...
    ellipse_coordinates_and_derivatives(ellipse_parameters,n_points,uni)

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


ksi1 = linspace(-uni,uni,n_points);
ksi2 = linspace(-uni,uni,n_points);
x1 = zeros(N_intervals_x,n_points);
x2 = zeros(N_intervals_y,n_points);
for k1 = 1:N_intervals_x
    x1(k1,:) = ( (b_x1(k1+1)-b_x1(k1))*ksi1+(b_x1(k1+1)+b_x1(k1)) )/2;
end
for k2 = 1:N_intervals_y     
    x2(k2,:) = ( (b_x2(k2+1)-b_x2(k2))*ksi2+(b_x2(k2+1)+b_x2(k2)) )/2;
end

x1_plus = zeros(3,n_points);
x1_minus = zeros(3,n_points);
x2_plus = zeros(3,n_points);
x2_minus = zeros(3,n_points);


x1_plus(1,:)  = x1_t_plus*ones(n_points,1);
x1_minus(1,:) = x1_t_minus*ones(n_points,1);

x1_plus(2,:)  = P1/2 + (R1/R2)*(R2^2-(x2(2,:)-P2/2).^2).^0.5;
x1_minus(2,:) = P1/2 - (R1/R2)*(R2^2-(x2(2,:)-P2/2).^2).^0.5;

x1_plus(3,:)  = x1_t_plus*ones(n_points,1);
x1_minus(3,:) = x1_t_minus*ones(n_points,1);

x2_plus(1,:)  = x2_t_plus*ones(n_points,1);
x2_minus(1,:) = x2_t_minus*ones(n_points,1);

x2_plus(2,:)  = P2/2 + (R2/R1)*(R1^2-(x1(2,:)-P1/2).^2).^0.5;
x2_minus(2,:) = P2/2 - (R2/R1)*(R1^2-(x1(2,:)-P1/2).^2).^0.5;

x2_plus(3,:)  = x2_t_plus*ones(n_points,1);
x2_minus(3,:) = x2_t_minus*ones(n_points,1);

dx1_plus(1,:) = zeros(n_points,1);
dx1_plus(2,:) = -(R1/R2)*(x2(2,:)-P2/2)./(R2^2-(x2(2,:)-P2/2).^2).^0.5;
dx1_plus(3,:) = zeros(n_points,1);

dx2_plus(1,:) = zeros(n_points,1);
dx2_plus(2,:) = -(R2/R1)*(x1(2,:)-P1/2)./(R1^2-(x1(2,:)-P1/2).^2).^0.5;
dx2_plus(3,:) = zeros(n_points,1);

dx1_minus(1,:) = zeros(n_points,1);
dx1_minus(2,:) = (R1/R2)*(x2(2,:)-P2/2)./(R2^2-(x2(2,:)-P2/2).^2).^0.5;
dx1_minus(3,:) = zeros(n_points,1);

dx2_minus(1,:) = zeros(n_points,1);
dx2_minus(2,:) = (R2/R1)*(x1(2,:)-P1/2)./(R1^2-(x1(2,:)-P1/2).^2).^0.5;
dx2_minus(3,:) = zeros(n_points,1);
%{
x = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
y = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
%}
dx_x1 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
dx_x2 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
dy_x1 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);
dy_x2 = zeros(N_intervals_x,N_intervals_y,n_points,n_points);

for nx=1:n_points
for ny=1:n_points
    
    %{
%%%%%%%%%%%%%%%%%%%%%% coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %x(x1,x2)
x(1,1,nx,ny) = x1(1,nx)*x1_plus(1,ny);
x(1,2,nx,ny) = x1(1,nx)*x1_plus(2,ny);
x(1,3,nx,ny) = x1(1,nx)*x1_plus(3,ny);

x(2,1,nx,ny) = (x1(2,nx)-x1_t_minus)*x1_plus(1,ny)/(x1_t_plus-x1_t_minus)+...
               (x1(2,nx)-x1_t_plus)*x1_minus(1,ny)/(x1_t_minus-x1_t_plus);
x(2,2,nx,ny) = (x1(2,nx)-x1_t_minus)*x1_plus(2,ny)/(x1_t_plus-x1_t_minus)+...
               (x1(2,nx)-x1_t_plus)*x1_minus(2,ny)/(x1_t_minus-x1_t_plus);
x(2,3,nx,ny) = (x1(2,nx)-x1_t_minus)*x1_plus(3,ny)/(x1_t_plus-x1_t_minus)+...
               (x1(2,nx)-x1_t_plus)*x1_minus(3,ny)/(x1_t_minus-x1_t_plus);
           
x(3,1,nx,ny) = (x1(3,nx)-P1)*x1_plus(1,ny)/(x1_t_plus-P1)+...
               (x1(3,nx)-x1_t_plus)*P1/(P1-x1_t_plus);
x(3,2,nx,ny) = (x1(3,nx)-P1)*x1_plus(2,ny)/(x1_t_plus-P1)+...
               (x1(3,nx)-x1_t_plus)*P1/(P1-x1_t_plus);
x(3,3,nx,ny) = (x1(3,nx)-P1)*x1_plus(3,ny)/(x1_t_plus-P1)+...
               (x1(3,nx)-x1_t_plus)*P1/(P1-x1_t_plus);

	%y(x1,x2)
y(1,1,nx,ny) = x2(1,ny)*x2_plus(1,nx);
y(2,1,nx,ny) = x2(1,ny)*x2_plus(2,nx);
y(3,1,nx,ny) = x2(1,ny)*x2_plus(3,nx);

y(1,2,nx,ny) = (x2(2,ny)-x2_t_minus)*x2_plus(1,nx)/(x2_t_plus-x2_t_minus)+...
               (x2(2,ny)-x2_t_plus)*x2_minus(1,nx)/(x2_t_minus-x2_t_plus);
y(2,2,nx,ny) = (x2(2,ny)-x2_t_minus)*x2_plus(2,nx)/(x2_t_plus-x2_t_minus)+...
               (x2(2,ny)-x2_t_plus)*x2_minus(2,nx)/(x2_t_minus-x2_t_plus);
y(3,2,nx,ny) = (x2(2,ny)-x2_t_minus)*x2_plus(3,nx)/(x2_t_plus-x2_t_minus)+...
               (x2(2,ny)-x2_t_plus)*x2_minus(3,nx)/(x2_t_minus-x2_t_plus);
           
y(1,3,nx,ny) = (x2(3,ny)-P2)*x2_plus(1,nx)/(x2_t_plus-P2)+...
               (x2(3,ny)-x2_t_plus)*P2/(P2-x2_t_plus);
y(2,3,nx,ny) = (x2(3,ny)-P2)*x2_plus(2,nx)/(x2_t_plus-P2)+...
               (x2(3,ny)-x2_t_plus)*P2/(P2-x2_t_plus);
y(3,3,nx,ny) = (x2(3,ny)-P2)*x2_plus(3,nx)/(x2_t_plus-P2)+...
               (x2(3,ny)-x2_t_plus)*P2/(P2-x2_t_plus);
    %}
%%%%%%%%%%%%%%%%%%%%%% derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dx/dx1
dx_x1(1,1,nx,ny) = x1_minus(1,ny)/x1_t_minus;
dx_x1(1,2,nx,ny) = x1_minus(2,ny)/x1_t_minus;
dx_x1(1,3,nx,ny) = x1_minus(3,ny)/x1_t_minus;

dx_x1(2,1,nx,ny) = (x1_plus(1,ny)-x1_minus(1,ny))/(x1_t_plus-x1_t_minus);
dx_x1(2,2,nx,ny) = (x1_plus(2,ny)-x1_minus(2,ny))/(x1_t_plus-x1_t_minus);
dx_x1(2,3,nx,ny) = (x1_plus(3,ny)-x1_minus(3,ny))/(x1_t_plus-x1_t_minus);

dx_x1(3,1,nx,ny) = (x1_plus(1,ny)-P1)/(x1_t_plus-P1);
dx_x1(3,2,nx,ny) = (x1_plus(2,ny)-P1)/(x1_t_plus-P1);
dx_x1(3,3,nx,ny) = (x1_plus(3,ny)-P1)/(x1_t_plus-P1);

    %dy/dx2
dy_x2(1,1,nx,ny) = x2_minus(1,nx)/x2_t_minus;
dy_x2(2,1,nx,ny) = x2_minus(2,nx)/x2_t_minus;
dy_x2(3,1,nx,ny) = x2_minus(3,nx)/x2_t_minus;

dy_x2(1,2,nx,ny) = (x2_plus(1,nx)-x2_minus(1,nx))/(x2_t_plus-x2_t_minus);
dy_x2(2,2,nx,ny) = (x2_plus(2,nx)-x2_minus(2,nx))/(x2_t_plus-x2_t_minus);
dy_x2(3,2,nx,ny) = (x2_plus(3,nx)-x2_minus(3,nx))/(x2_t_plus-x2_t_minus);

dy_x2(1,3,nx,ny) = (x2_plus(1,nx)-P2)/(x2_t_plus-P2);
dy_x2(2,3,nx,ny) = (x2_plus(2,nx)-P2)/(x2_t_plus-P2);
dy_x2(3,3,nx,ny) = (x2_plus(3,nx)-P2)/(x2_t_plus-P2);

    %dx/dx2
dx_x2(1,1,nx,ny) = x1(1,nx)*dx1_minus(1,ny)/x1_t_minus;
dx_x2(1,2,nx,ny) = x1(1,nx)*dx1_minus(2,ny)/x1_t_minus;
dx_x2(1,3,nx,ny) = x1(1,nx)*dx1_minus(3,ny)/x1_t_minus;

dx_x2(2,1,nx,ny) = (x1(2,nx)-x1_t_minus)*dx1_plus(1,ny)/(x1_t_plus-x1_t_minus)+...
                   (x1(2,nx)-x1_t_plus)*dx1_minus(1,ny)/(x1_t_minus-x1_t_plus);
dx_x2(2,2,nx,ny) = (x1(2,nx)-x1_t_minus)*dx1_plus(2,ny)/(x1_t_plus-x1_t_minus)+...
                   (x1(2,nx)-x1_t_plus)*dx1_minus(2,ny)/(x1_t_minus-x1_t_plus);
dx_x2(2,3,nx,ny) = (x1(2,nx)-x1_t_minus)*dx1_plus(3,ny)/(x1_t_plus-x1_t_minus)+...
                   (x1(2,nx)-x1_t_plus)*dx1_minus(3,ny)/(x1_t_minus-x1_t_plus);  
               
dx_x2(3,1,nx,ny) = (x1(3,nx)-P1)*dx1_plus(1,ny)/(x1_t_plus-P1);
dx_x2(3,2,nx,ny) = (x1(3,nx)-P1)*dx1_plus(2,ny)/(x1_t_plus-P1);
dx_x2(3,3,nx,ny) = (x1(3,nx)-P1)*dx1_plus(3,ny)/(x1_t_plus-P1);             

    %dy/dx1
dy_x1(1,1,nx,ny) = x2(1,ny)*dx2_minus(1,nx)/x2_t_minus;
dy_x1(2,1,nx,ny) = x2(1,ny)*dx2_minus(2,nx)/x2_t_minus;
dy_x1(3,1,nx,ny) = x2(1,ny)*dx2_minus(3,nx)/x2_t_minus;  

dy_x1(1,2,nx,ny) = (x2(2,ny)-x2_t_minus)*dx2_plus(1,nx)/(x2_t_plus-x2_t_minus)+...
                   (x2(2,ny)-x2_t_plus)*dx2_minus(1,nx)/(x2_t_minus-x2_t_plus);
dy_x1(2,2,nx,ny) = (x2(2,ny)-x2_t_minus)*dx2_plus(2,nx)/(x2_t_plus-x2_t_minus)+...
                   (x2(2,ny)-x2_t_plus)*dx2_minus(2,nx)/(x2_t_minus-x2_t_plus);
dy_x1(3,2,nx,ny) = (x2(2,ny)-x2_t_minus)*dx2_plus(3,nx)/(x2_t_plus-x2_t_minus)+...
                   (x2(2,ny)-x2_t_plus)*dx2_minus(3,nx)/(x2_t_minus-x2_t_plus); 
               
dy_x1(1,3,nx,ny) = (x2(3,ny)-P2)*dx2_plus(1,nx)/(x2_t_plus-P2);
dy_x1(2,3,nx,ny) = (x2(3,ny)-P2)*dx2_plus(2,nx)/(x2_t_plus-P2);
dy_x1(3,3,nx,ny) = (x2(3,ny)-P2)*dx2_plus(3,nx)/(x2_t_plus-P2);  

end
end
