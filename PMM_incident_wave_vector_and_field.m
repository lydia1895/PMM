function [alpha0,beta0,gamma0,k0,Ex0,Ey0] =...
        PMM_incident_wave_vector_and_field(lambda,theta,phi,delta,n1)
    
    k0 = 2*pi/lambda;
    
    alpha0 = k0*n1*sin(theta)*cos(phi);
    beta0  = k0*n1*sin(theta)*sin(phi);
    gamma0 = k0*n1*cos(theta);
    
    
    %incident wave
    
    TETM = [0; cos(delta); sin(delta)];
    TETMmatrix = [sin(theta)*cos(phi) cos(theta)*cos(phi) -sin(phi); ...
        sin(theta)*sin(phi) cos(theta)*sin(phi) cos(phi);...
        cos(phi) -sin(theta) 0];
    E = TETMmatrix*TETM;
    Ex = E(1);
    Ey = E(2);
    
    k1 = k0*n1;
    kz1v = gamma0;
    A1_nul = ( k1^2 - alpha0^2)/(k0*kz1v);
    B1_nul = ( k1^2 - beta0^2)/(k0*kz1v);
    C1_nul = alpha0*beta0/(k0*kz1v);
    
    norm = A1_nul*abs(Ey)^2 + B1_nul*abs(Ex)^2 +...
        C1_nul*( Ex*conj(Ey)+Ey*conj(Ex) );
    
    Ex0=Ex/sqrt(norm);
    Ey0=Ey/sqrt(norm);
    
    