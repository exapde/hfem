function [sr,sr_udg] = source(p,udg,param,time)

[ng,nch] = size(udg);
nc = nch;

gam  = param{1};
gam1 = gam - 1.0;

r    = udg(:,1);
ru   = udg(:,2);
rv   = udg(:,3);
rE   = udg(:,4);

r1   = 1./r;
uv   = ru.*r1;
vv   = rv.*r1;
af   = 0.5*(uv.*uv+vv.*vv);
pr    = gam1*(rE-r.*af);

sr = zeros(ng,nc);
sr(:,3) = pr;

sr_udg = zeros(ng,nc,nc); 
sr_udg(:,3,1) = gam1*af;
sr_udg(:,3,2) = -gam1*uv;
sr_udg(:,3,3) = -gam1*vv;
sr_udg(:,3,4) = gam1;

