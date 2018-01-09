function [integral1] = int_exp_gegenbauer(a,n,La)

integral1 = pi*(2^(1-La))*((1j)^n)*gamma(2*La+n)*a^(-La)*besselj(La+n,a)/...
    (factorial(n)*gamma(La));