function [sr,dsr_dudg] = source(p,udg,param,time)

[ng,nc] = size(udg);

mu = param{1};
x = p(:,1);
y = p(:,2);

sr = zeros(ng,nc);
sr(:,1) = -mu*(3*y - (pi^2*sin((pi*y)/2))/8);

dsr_dudg = zeros(ng,3,nc); 



   