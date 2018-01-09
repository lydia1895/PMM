function [ integral_der ] = PMM_integral_deriv( ni, nui, alphai, betai )
%Integral for x=[-1,1] that we use for <Cn,d(Cm)/dx>
%Integrate[(1-x)^alpha*(1+x)^beta*Gegenbauer(n,nu), {x,-1,1}]

%for hypergeometric function 3F2(a,b,z) with coefficients a=[a1 a2 a3], b=[b1 b2]
a = [-ni, ni+2*nui, alphai+1];
b = [nui+0.5 alphai+betai+2];

integral_der = (2^(alphai+betai+1))*gamma(alphai+1)*gamma(betai+1)*gamma(ni+2*nui)/...
    (factorial(ni)*gamma(2*nui)*gamma(alphai+betai+2))*hypergeom(a,b,1);

end

