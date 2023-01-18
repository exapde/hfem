function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

sr = ones(ng,nch);

% Vx = param{2}(1);
% Vy = param{2}(2);
% sr = (Vy.*x.*exp(Vx.*x).*exp(Vy.*y) + Vx.*y.*exp(Vx.*x).*exp(Vy.*y) + Vy.*x.*exp(Vx.*x).*exp(Vy) - Vx.*y.*exp(Vx.*x).*exp(Vy) - Vy.*x.*exp(Vy.*y).*exp(Vx) + Vx.*y.*exp(Vy.*y).*exp(Vx) - Vy.*x.*exp(Vx).*exp(Vy) - Vx.*y.*exp(Vx).*exp(Vy))./(exp(Vx) + exp(Vy) - exp(Vx).*exp(Vy) - 1);
    
sr_udg = zeros(ng,nch,nc); 

% if nc==3
%     xg = p(:,1);
%     yg = p(:,2);
%     %sr = 2*pi^2*sin(pi*xg).*sin(pi*yg);
% elseif nc==4
%     xg = p(:,1);
%     yg = p(:,2);
%     zg = p(:,3);
%     sr = 3*pi^2*sin(pi*xg).*sin(pi*yg).*sin(pi*zg);
% end
