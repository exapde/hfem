function [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param,time)
%FHAT flux function
%   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Lambda

[ng,nc] = size(u);

gam  = param{1};
gam1 = gam - 1.0;
                                             
nx   = nl(:,1);              
ny   = nl(:,2);

[f,fm,fq] = flux(p,m,q,param);

r    = m(:,1);
ru   = m(:,2);
rv   = m(:,3);
rE   = m(:,4);

zer  = zeros(ng,1);
one  = ones(ng,1);

r1   = 1./r;
r1m1 = -1./(r.^2);
r1m2 = zer;
r1m3 = zer;
r1m4 = zer;

uv   = ru.*r1;
uvm1 =          ru.*r1m1;
uvm2 =     r1 + ru.*r1m2;
uvm3 =          ru.*r1m3;
uvm4 =          ru.*r1m4;

vv   = rv.*r1;
vvm1 =          rv.*r1m1;
vvm2 =          rv.*r1m2;
vvm3 =     r1 + rv.*r1m3;
vvm4 =          rv.*r1m4;

E    = rE.*r1;
Em1  =          rE.*r1m1;
Em2  =          rE.*r1m2;
Em3  =          rE.*r1m3;
Em4  =     r1 + rE.*r1m4;

af   = 0.5*(uv.*uv+vv.*vv);
afm1 = uv.*uvm1 + vv.*vvm1;
afm2 = uv.*uvm2 + vv.*vvm2;
afm3 = uv.*uvm3 + vv.*vvm3;
afm4 = uv.*uvm4 + vv.*vvm4;

p    = gam1*(rE -r.*af);
pm1  = gam1*(   -   af - r.*afm1);
pm2  = gam1*(          - r.*afm2);
pm3  = gam1*(          - r.*afm3);
pm4  = gam1*(one       - r.*afm4);

h    = E   + p.*r1;
hm1  = Em1 + pm1.*r1 + p.*r1m1;
hm2  = Em2 + pm2.*r1 + p.*r1m2;
hm3  = Em3 + pm3.*r1 + p.*r1m3;
hm4  = Em4 + pm4.*r1 + p.*r1m4;

c2   = gam* p.*r1;
c2m1 = gam*(pm1.*r1 + p.*r1m1);
c2m2 = gam*(pm2.*r1 + p.*r1m2);
c2m3 = gam*(pm3.*r1 + p.*r1m3);
c2m4 = gam*(pm4.*r1 + p.*r1m4);

c    = sqrt(c2);
cm1  = 0.5*c2m1./c;
cm2  = 0.5*c2m2./c;
cm3  = 0.5*c2m3./c;
cm4  = 0.5*c2m4./c;

un   = uv.*nx   + vv.*ny;
unm1 = uvm1.*nx + vvm1.*ny;
unm2 = uvm2.*nx + vvm2.*ny;
unm3 = uvm3.*nx + vvm3.*ny;
unm4 = uvm4.*nx + vvm4.*ny;
                                         
rlam   = 1.0*(abs(un)+c);
rlamm1 = 1.0*(sign(un).*unm1+cm1);
rlamm2 = 1.0*(sign(un).*unm2+cm2);
rlamm3 = 1.0*(sign(un).*unm3+cm3);
rlamm4 = 1.0*(sign(un).*unm4+cm4);

An      = zeros(ng,nc,nc);
Anm     = zeros(ng,nc,nc,nc);

An(:,:,1)  = [rlam, zer, zer, zer];
Anm(:,:,1,1) = [rlamm1, zer, zer, zer];
Anm(:,:,1,2) = [rlamm2, zer, zer, zer];
Anm(:,:,1,3) = [rlamm3, zer, zer, zer];
Anm(:,:,1,4) = [rlamm4, zer, zer, zer];

An(:,:,2)  = [zer, rlam, zer, zer];
Anm(:,:,2,1)  = [zer, rlamm1, zer, zer];
Anm(:,:,2,2)  = [zer, rlamm2, zer, zer];
Anm(:,:,2,3)  = [zer, rlamm3, zer, zer];
Anm(:,:,2,4)  = [zer, rlamm4, zer, zer];

An(:,:,3)  = [zer, zer, rlam, zer];
Anm(:,:,3,1)  = [zer, zer, rlamm1, zer];
Anm(:,:,3,2)  = [zer, zer, rlamm2, zer];
Anm(:,:,3,3)  = [zer, zer, rlamm3, zer];
Anm(:,:,3,4)  = [zer, zer, rlamm4, zer];

An(:,:,4)  = [zer, zer, zer, rlam];
Anm(:,:,4,1)  = [zer, zer, zer, rlamm1];
Anm(:,:,4,2)  = [zer, zer, zer, rlamm2];
Anm(:,:,4,3)  = [zer, zer, zer, rlamm3];
Anm(:,:,4,4)  = [zer, zer, zer, rlamm4];

fh = multiprod(f(:,:,1),nx,2)+multiprod(f(:,:,2),ny,2) + multiprod(An,(u-m),[2,3],2);
fhu = An;
fhq = zeros(ng,nc,nc,2);
fhm = multiprod(permute(fm(:,:,1,:),[1,2,4,3]),nx,[2,3],2)+multiprod(permute(fm(:,:,2,:),[1,2,4,3]),ny,[2,3],2) + multiprod(Anm,(u-m),[2,3],2) - An;


