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
% sr(:,3) = (2*mu-rho)*pi^2*sin(pi*time)*sin(pi*x).*sin(pi*y);

% Saint Venant-Kirchhoff
sr(:,1) =-A^2*pi^3*sin(pi*time)^2*cos(pi*x).*sin(pi*x).*(lambda*cos(2*pi*y) - mu + 2*mu*cos(2*pi*y));
sr(:,2) =-A^2*pi^3*sin(pi*time)^2*cos(pi*y).*sin(pi*y).*(lambda*cos(2*pi*x) - mu + 2*mu*cos(2*pi*x));
sr(:,3) = A*pi^2*sin(pi*time)*sin(pi*x).*sin(pi*y).*(2*mu - rho ...
        - 2*A^2*lambda*pi^2*sin(pi*time)^2 - 4*A^2*mu*pi^2*sin(pi*time)^2 ...
        + 4*A^2*lambda*pi^2*sin(pi*time)^2*sin(pi*x).^2 ...
        + 4*A^2*lambda*pi^2*sin(pi*time)^2*sin(pi*y).^2 ...
        + 8*A^2*mu*pi^2*sin(pi*time)^2*sin(pi*x).^2 ...
        + 8*A^2*mu*pi^2*sin(pi*time)^2*sin(pi*y).^2 ...
        - 6*A^2*lambda*pi^2*sin(pi*time)^2*sin(pi*x).^2 .*sin(pi*y).^2 ...
        - 12*A^2*mu*pi^2*sin(pi*time)^2*sin(pi*x).^2 .*sin(pi*y).^2);
    
dsr_dudg = zeros(ng,3,nc);



   