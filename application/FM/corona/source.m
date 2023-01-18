function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 2;

eps0 = param{3};

sr = zeros(ng,nch);
sr(:,1) = udg(:,2)/eps0;

sr_udg = zeros(ng,nch,nc); 
sr_udg(:,1,2) = 1/eps0;



