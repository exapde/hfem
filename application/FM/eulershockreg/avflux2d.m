function [f,f_udg] = avflux2d(pg,udg,param,time)

% cpar = 20.0;
% [ng,nc] = size(udg);
% nch  = 4;
% 
% rmin = param{8};
% r = udg(:,1);
% 
% inda = find(cpar*r/rmin>30);
% sa = sensor2d1(pg(inda,:),udg(inda,:),param,time);
% ia = sa>30;
% ind1 = inda(ia);
% ind2 = setdiff(inda,ind1);
% 
% indb = setdiff((1:ng)',inda);
% sb = sensor2d2(pg(indb,:),udg(indb,:),param,time);
% ib = sb>30;
% ind3 = indb(ib);
% ind4 = setdiff(indb,ind3);
% 
% f = zeros(ng,nch,2);
% f_udg = zeros(ng,nch,2,nc);
% 
% [f(ind1,:,:),f_udg(ind1,:,:,:)] = avflux2d1(pg(ind1,:),udg(ind1,:),param,time);
% [f(ind2,:,:),f_udg(ind2,:,:,:)] = avflux2d2(pg(ind2,:),udg(ind2,:),param,time);
% [f(ind3,:,:),f_udg(ind3,:,:,:)] = avflux2d5(pg(ind3,:),udg(ind3,:),param,time);
% [f(ind4,:,:),f_udg(ind4,:,:,:)] = avflux2d6(pg(ind4,:),udg(ind4,:),param,time);
% 
% ind = [ind1; ind2; ind3; ind4];
% if length(ind)~=ng || max(abs(sort(ind)-(1:ng)'))>0 
%     error('here');
% end

cpar = 20.0;
[ng,nc] = size(udg);
nch  = 4;

gam  = param{1};
gam1 = gam-1;
rmin = param{8};
pmin = param{9};
r  = udg(:,1);
ru = udg(:,2);
rv = udg(:,3);
rE = udg(:,4);
r1 = 1./r;
u  = ru.*r1;
v  = rv.*r1;
q  = 0.5*(u.*u+v.*v);
p  = gam1*(rE-r.*q);

inda = find(r>rmin);
pa   = p(inda);
sa   = sensor2d1(pg(inda,:),udg(inda,:),param,time);
ia1  = (pa>pmin) & (sa>30);  % r > rmin and p > pmin
ia2  = (pa>pmin) & (sa<=30); % r > rmin and p > pmin
ia3  = (pa<=pmin) & (sa>30); % r > rmin and p <= pmin
ia4  = (pa<=pmin) & (sa<=30);% r > rmin and p <= pmin
ind1 = inda(ia1);
ind2 = inda(ia2);
ind3 = inda(ia3);
ind4 = inda(ia4);


indb = setdiff((1:ng)',inda);
rt = (rmin/cpar)*log(cpar+exp(cpar*(r/rmin)));
r1 = 1./rt;
u  = ru.*r1;
v  = rv.*r1;
q  = 0.5*(u.*u+v.*v);
p  = gam1*(rE-rt.*q);

pb   = p(indb);
sb   = sensor2d2(pg(indb,:),udg(indb,:),param,time);
ib1  = (pb>pmin) & (sb>30);  % r <= rmin and p > pmin
ib2  = (pb>pmin) & (sb<=30); % r <= rmin and p > pmin
ib3  = (pb<=pmin) & (sb>30); % r <= rmin and p <= pmin
ib4  = (pb<=pmin) & (sb<=30);% r <= rmin and p <= pmin
ind5 = indb(ib1);
ind6 = indb(ib2);
ind7 = indb(ib3);
ind8 = indb(ib4);

f = zeros(ng,nch,2);
f_udg = zeros(ng,nch,2,nc);

[f(ind1,:,:),f_udg(ind1,:,:,:)] = avflux2d1(pg(ind1,:),udg(ind1,:),param,time);


if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    any(isnan(f_udg(:)))
    any(isinf(f_udg(:)))
    error('here');
end

[f(ind2,:,:),f_udg(ind2,:,:,:)] = avflux2d2(pg(ind2,:),udg(ind2,:),param,time);


if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    any(isnan(f_udg(:)))
    any(isinf(f_udg(:)))
    error('here');
end

[f(ind3,:,:),f_udg(ind3,:,:,:)] = avflux2d3(pg(ind3,:),udg(ind3,:),param,time);


if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    any(isnan(f_udg(:)))
    any(isinf(f_udg(:)))
    error('here');
end

[f(ind4,:,:),f_udg(ind4,:,:,:)] = avflux2d4(pg(ind4,:),udg(ind4,:),param,time);


if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    any(isnan(f_udg(:)))
    any(isinf(f_udg(:)))
    error('here');
end

[f(ind5,:,:),f_udg(ind5,:,:,:)] = avflux2d5(pg(ind5,:),udg(ind5,:),param,time);


if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    any(isnan(f_udg(:)))
    any(isinf(f_udg(:)))
    error('here');
end

[f(ind6,:,:),f_udg(ind6,:,:,:)] = avflux2d6(pg(ind6,:),udg(ind6,:),param,time);


if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    any(isnan(f_udg(:)))
    any(isinf(f_udg(:)))
    error('here');
end

[f(ind7,:,:),f_udg(ind7,:,:,:)] = avflux2d7(pg(ind7,:),udg(ind7,:),param,time);


if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    any(isnan(f_udg(:)))
    any(isinf(f_udg(:)))
    error('here');
end

[f(ind8,:,:),f_udg(ind8,:,:,:)] = avflux2d8(pg(ind8,:),udg(ind8,:),param,time);


if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    any(isnan(f_udg(:)))
    any(isinf(f_udg(:)))
    error('here');
end

ind = [ind1; ind2; ind3; ind4; ind5; ind6; ind7; ind8];
if length(ind)~=ng || max(abs(sort(ind)-(1:ng)'))>0 
    error('here');
end

