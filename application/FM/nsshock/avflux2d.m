function [f,f_udg] = avflux2d(pg,udg,param,time)


[ng,nc] = size(udg);
nch  = 4;

sa = sensor2d1(pg,udg,param,time);
ind1 = find(sa>30);
ind2 = setdiff((1:ng)',ind1);

f = zeros(ng,nch,2);
f_udg = zeros(ng,nch,2,nc);

[f(ind1,:,:),f_udg(ind1,:,:,:)] = avflux2d1(pg(ind1,:),udg(ind1,:),param,time);
[f(ind2,:,:),f_udg(ind2,:,:,:)] = avflux2d2(pg(ind2,:),udg(ind2,:),param,time);

