function [sr,dsr_dudg] = source(p,udg,param,time)

[ng,nc] = size(udg);

if nc == 6
    nch = 2;
else 
    nch = 3;
end

sr = zeros(ng,nc);
dsr_dudg = zeros(ng,nch,nc); 




   