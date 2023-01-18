function [sr,sr_udg] = source2(p,udg,param,time)

[ng,nc] = size(udg);
nch = 2;

eps0 = param{3};

sr = zeros(ng,nch);
sr(:,1) = p(:,5)/eps0;

sr_udg = zeros(ng,nch,nc); 



