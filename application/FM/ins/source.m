function [sr,dsr_dudg] = source(pg,udg,param,time)

[ng,nc] = size(udg);
nd  = size(pg,2);
nch = nd;

sr = zeros(ng,nch);
dsr_dudg = zeros(ng,nch,nc); 

