function [f,f_udg] = flux2d(pg,udg,param,time)

cpar = 20.0;
[ng,nc] = size(udg);
nch  = 4;

gam  = param{1};
gam1 = gam-1;
vis  = param{7};
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
ia   = p(inda)>pmin;
ind1 = inda(ia);           % r > rmin and p > pmin
ind2 = setdiff(inda,ind1); % r > rmin and p <= pmin

indb = setdiff((1:length(r))',inda);
rt = (rmin/cpar)*log(cpar+exp(cpar*(r/rmin)));
r1 = 1./rt;
u  = ru.*r1;
v  = rv.*r1;
q  = 0.5*(u.*u+v.*v);
p  = gam1*(rE-rt.*q);
ib   = p(indb)>pmin;
ind3 = indb(ib);          % r <= rmin and p > pmin
ind4 = setdiff(indb,ind3);% r <= rmin and p <= pmin

f = zeros(ng,nch,2);
f_udg = zeros(ng,nch,2,nc);
[f(ind1,:,:),f_udg(ind1,:,:,:)] = flux2d1(pg(ind1,:),udg(ind1,:),param,time);
[f(ind2,:,:),f_udg(ind2,:,:,:)] = flux2d2(pg(ind2,:),udg(ind2,:),param,time);
[f(ind3,:,:),f_udg(ind3,:,:,:)] = flux2d3(pg(ind3,:),udg(ind3,:),param,time);
[f(ind4,:,:),f_udg(ind4,:,:,:)] = flux2d4(pg(ind4,:),udg(ind4,:),param,time);

if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    any(isnan(f_udg(:)))
    any(isinf(f_udg(:)))
    error('here');
end

[fav,fav_udg] = avflux2d(pg,udg,param,time);
f = f + vis*bsxfun(@times,pg(:,4),fav);
f_udg = f_udg + vis*bsxfun(@times,pg(:,4),fav_udg);

% check 
if any(isnan(f_udg(:))) || any(isinf(f_udg(:)))
    any(isnan(f_udg(:)))
    any(isinf(f_udg(:)))
    error('here');
end
ind = [ind1; ind2; ind3; ind4];
if length(ind)~=ng || max(abs(sort(ind)-(1:ng)'))>0 
    error('here');
end


