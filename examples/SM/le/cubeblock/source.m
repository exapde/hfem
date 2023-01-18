function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

if nc==3
    xg = p(:,1);
    yg = p(:,2);
    sr = 2*pi^2*sin(pi*xg).*sin(pi*yg);
elseif nc==4
    xg = p(:,1);
    yg = p(:,2);
    zg = p(:,3);
    sr = 3*pi^2*sin(pi*xg).*sin(pi*yg).*sin(pi*zg);
    %sr = 0*ones(ng,nch);
end

sr_udg = zeros(ng,nch,nc); 



