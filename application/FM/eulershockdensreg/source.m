function [sr,sr_udg] = source(p,udg,param,time)

[ng,nch] = size(udg);
nc = nch;

sr = zeros(ng,nc);

sr_udg = zeros(ng,nc,nc); 

