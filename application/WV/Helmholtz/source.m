function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

k = param{1};
sr = (k.^2).*udg(:,1) + 100*exp(-100*sqrt((p(:,1)-0.5).^2+(p(:,2)-0.5).^2));

sr_udg = zeros(ng,nch,nc); 
sr_udg(:,1,1) = k.^2;
