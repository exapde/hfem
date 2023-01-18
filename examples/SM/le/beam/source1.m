function [sr,dsr_dudg] = source(p,udg,param,time)

[ng,nc] = size(udg);

sr = zeros(ng,nc);

dsr_dudg = zeros(ng,3,nc); 



   