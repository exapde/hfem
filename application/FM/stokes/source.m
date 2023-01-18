function [sr,dsr_dudg] = source(pg,udg,param,time)

[ng,nc] = size(udg);
nd  = size(pg,2);
nch = nd;

sr = zeros(ng,nch);
dsr_dudg = zeros(ng,nch,nc); 

if nd==3
    x = pg(:,1);
    y = pg(:,2);
    z = pg(:,3);
    sr(:,nch)=-9*pi^2*cos(pi*z).*sin(pi*x).*sin(pi*y);
end