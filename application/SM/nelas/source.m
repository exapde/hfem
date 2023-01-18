function [sr,dsr_dudg] = source(pg,udg,param,time)

[ng,nc] = size(udg);

if nc == 7
    nch = 2;
else 
    nch = 3;
end

sr = zeros(ng,nch);
dsr_dudg = zeros(ng,nch,nc); 




   