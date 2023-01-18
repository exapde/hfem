function [sr,dsr_dudg] = source(p,udg,param,time)

[ng,nc] = size(udg);

mu     = param{1};
lambda = param{2};
rho    = param{4};
A      = param{5};

x = p(:,1);
y = p(:,2);

sr = zeros(ng,nc);

% Linear Elastic case :
sr(:,3) = A*(2*mu-rho)*pi^2*sin(pi*time)*sin(pi*x).*sin(pi*y);
    
dsr_dudg = zeros(ng,3,nc);



   