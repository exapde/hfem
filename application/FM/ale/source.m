function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);

sr = zeros(ng,nc);

sr_udg = zeros(ng,nc,nc); 
