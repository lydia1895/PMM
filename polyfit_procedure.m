function [g_fit] = polyfit_procedure(ksi1,ksi2,g)

[XOut, YOut, ZOut] = prepareSurfaceData(ksi1, ksi2, g);

g_fit = fit([XOut, YOut],ZOut,'poly55');
%plot(g_sqrt_inv_12,[XOut,YOut],ZOut,'Style','Residuals')
