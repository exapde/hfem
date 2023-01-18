function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

sr = 10*ones(ng,nch);

sr_udg = zeros(ng,nch,nc); 
