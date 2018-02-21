

clc
clear all

load('eig_nonpeturbed.mat','H_nonperturbed_sorted','gammasqr_nonperturbed_sorted');
load('eig_peturbed.mat','H_perturbed_sorted','gammasqr_perturbed_sorted');

theta_btw_perturbed_nonperturbed=zeros(450,1);
theta_btw_perturbed_neighbour=zeros(450,1);
theta_btw_nonperturbed_neighbour=zeros(450,1);

for nvector=1:450
    u(:) = H_perturbed_sorted(1,1,:,nvector);
    v(:) = H_nonperturbed_sorted(1,1,:,nvector);
    %cos_theta = dot(u,v)/(norm(u)*norm(v));
    cos_theta = (sum(u.*conj(v)))/( sqrt ( dot(u,u) )*sqrt( dot(v,v) ) );
    theta_btw_perturbed_nonperturbed(nvector) = acosd(abs(cos_theta));
end

for nvector=1:449
    u(:) = H_perturbed_sorted(1,1,:,nvector);
    v(:) = H_perturbed_sorted(1,1,:,nvector+1);
    %cos_theta = dot(u,v)/(norm(u)*norm(v));
    cos_theta = (sum(u.*conj(v)))/( sqrt ( dot(u,u) )*sqrt( dot(v,v) ) );
    theta_btw_perturbed_neighbour(nvector) = acosd(abs(cos_theta));
    
    u(:) = H_nonperturbed_sorted(1,1,:,nvector);
    v(:) = H_nonperturbed_sorted(1,1,:,nvector+1);
    %cos_theta = dot(u,v)/(norm(u)*norm(v));
    cos_theta =(sum(u.*conj(v)))/( sqrt ( dot(u,u)) )*sqrt( dot(v,v)  );
    theta_btw_nonperturbed_neighbour(nvector) = acosd(abs(cos_theta));
    cos_nonperturbed(nvector) = cos_theta
end
abs(cos_nonperturbed)

gammasqr_nonperturbed(:,1) = gammasqr_nonperturbed_sorted(1,1,:);
gammasqr_perturbed(:,1) = gammasqr_perturbed_sorted(1,1,:);

save('angles_perturbation.mat','theta_btw_perturbed_nonperturbed',...
'theta_btw_perturbed_neighbour','theta_btw_nonperturbed_neighbour',...
'gammasqr_perturbed_sorted','gammasqr_nonperturbed_sorted');

dlmwrite('theta_btw_perturbed_nonperturbed.txt',...
    theta_btw_perturbed_nonperturbed);%,'delimiter','\t','precision',3)
dlmwrite('theta_btw_perturbed_neighbour.txt',...
    theta_btw_perturbed_neighbour)%,'delimiter','\t','precision',3)
dlmwrite('theta_btw_nonperturbed_neighbour.txt',...
    theta_btw_nonperturbed_neighbour)%,'delimiter','\t','precision',3)
dlmwrite('gammasqr_nonperturbed_sorted.txt',...
    gammasqr_nonperturbed_sorted)%,'delimiter','\t','precision',3)

matrix = cat(2,gammasqr_nonperturbed,theta_btw_perturbed_nonperturbed,...
    theta_btw_perturbed_neighbour,theta_btw_nonperturbed_neighbour);


dlmwrite('angles.txt',matrix,'delimiter','\t','precision',4)



