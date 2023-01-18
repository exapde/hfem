function [f,f_udg] = avflux2d(pg,udg,param,time)


s = sensor2d(pg,udg,param,time);
ind1 = find(s<=30);
ind2 = setdiff(1:size(s,1),ind1);

[ng,nc] = size(udg);
f = zeros(ng,1);
f_udg = zeros(ng,nc);

[f(ind1),f_udg(ind1,:)] = avflux2d1(pg(ind1,:),udg(ind1,:),param,time);
[f(ind2),f_udg(ind2,:)] = avflux2d2(pg(ind2,:),udg(ind2,:),param,time);

