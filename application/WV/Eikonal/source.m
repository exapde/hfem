function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

x  = p(:,1);
y  = p(:,2);
s  = 2*sqrt(y.^2.*(1-x.^2).^2+x.^2.*(1-y.^2).^2);
v  = sqrt(udg(:,2).^2+udg(:,3).^2);
sr = 1 - v;

sr_udg = zeros(ng,nch,nc);
sr_udg(:,1,2) = -udg(:,2)./v;
sr_udg(:,1,3) = -udg(:,3)./v;
ind = find(abs(v)<1e-10);
sr_udg(ind,1,2) = -1;
sr_udg(ind,1,3) = -1;





