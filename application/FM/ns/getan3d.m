function [An,Anm] = getan3d(nl,uh,param,absolute)
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

[ng,nc] = size(uh);

gam   = param{1};
gam1 = gam - 1.0;

nx   = nl(:,1);              
ny   = nl(:,2);
nz   = nl(:,3);

r    = uh(:,1);
ru   = uh(:,2);
rv   = uh(:,3);
rw   = uh(:,4);
rE   = uh(:,5);

zer  = zeros(ng,1);
one  = ones(ng,1);

r1   = 1./r;
r1m1 = -1./(r.^2);
r1m2 = zer;
r1m3 = zer;
r1m4 = zer;
r1m5 = zer;

uv   = ru.*r1;
uvm1 =          ru.*r1m1;
uvm2 =     r1 + ru.*r1m2;
uvm3 =          ru.*r1m3;
uvm4 =          ru.*r1m4;
uvm5 =          ru.*r1m5;   

vv   = rv.*r1;
vvm1 =          rv.*r1m1;
vvm2 =          rv.*r1m2;
vvm3 =     r1 + rv.*r1m3;
vvm4 =          rv.*r1m4;
vvm5 =          rv.*r1m5;

wv   = rw.*r1;
wvm1 =          rw.*r1m1;
wvm2 =          rw.*r1m2;
wvm3 =          rw.*r1m3;
wvm4 =    r1 +  rw.*r1m4;
wvm5 =          rw.*r1m5;

E    = rE.*r1;
Em1  =          rE.*r1m1;
Em2  =          rE.*r1m2;
Em3  =          rE.*r1m3;
Em4  =          rE.*r1m4;
Em5  =     r1 + rE.*r1m5;

af   = 0.5*(uv.*uv+vv.*vv+wv.*wv);
afm1 = uv.*uvm1 + vv.*vvm1 + wv.*wvm1;
afm2 = uv.*uvm2 + vv.*vvm2 + wv.*wvm2;
afm3 = uv.*uvm3 + vv.*vvm3 + wv.*wvm3;
afm4 = uv.*uvm4 + vv.*vvm4 + wv.*wvm4;
afm5 = uv.*uvm5 + vv.*vvm5 + wv.*wvm5;

p    = gam1*(rE -r.*af);
pm1  = gam1*(   -   af - r.*afm1);
pm2  = gam1*(          - r.*afm2);
pm3  = gam1*(          - r.*afm3);
pm4  = gam1*(          - r.*afm4);
pm5  = gam1*(one       - r.*afm5);

h    = E   + p.*r1;
hm1  = Em1 + pm1.*r1 + p.*r1m1;
hm2  = Em2 + pm2.*r1 + p.*r1m2;
hm3  = Em3 + pm3.*r1 + p.*r1m3;
hm4  = Em4 + pm4.*r1 + p.*r1m4;
hm5  = Em5 + pm5.*r1 + p.*r1m5;

c2   = gam* p.*r1;
c2m1 = gam*(pm1.*r1 + p.*r1m1);
c2m2 = gam*(pm2.*r1 + p.*r1m2);
c2m3 = gam*(pm3.*r1 + p.*r1m3);
c2m4 = gam*(pm4.*r1 + p.*r1m4);
c2m5 = gam*(pm5.*r1 + p.*r1m5);

c    = sqrt(c2);
cm1  = 0.5*c2m1./c;
cm2  = 0.5*c2m2./c;
cm3  = 0.5*c2m3./c;
cm4  = 0.5*c2m4./c;
cm5  = 0.5*c2m5./c;

un   = uv.*nx   + vv.*ny + wv.*nz;
unm1 = uvm1.*nx + vvm1.*ny + wvm1.*nz;
unm2 = uvm2.*nx + vvm2.*ny + wvm2.*nz;
unm3 = uvm3.*nx + vvm3.*ny + wvm3.*nz;
unm4 = uvm4.*nx + vvm4.*ny + wvm4.*nz;                             
unm5 = uvm5.*nx + vvm5.*ny + wvm5.*nz;                             

if absolute
    rlam1   = abs(un+c);
    rlam1m1 = sign(un+c).*(unm1+cm1);
    rlam1m2 = sign(un+c).*(unm2+cm2);
    rlam1m3 = sign(un+c).*(unm3+cm3);
    rlam1m4 = sign(un+c).*(unm4+cm4);
    rlam1m5 = sign(un+c).*(unm5+cm5);
    
    epslm= param{2};
    if epslm>0
        rlam = 0.5*(rlam1.*rlam1./(epslm*c)+epslm*c);
        rlamm1 = rlam1.*rlam1m1./(epslm*c) + 0.5*(1 - rlam1.*rlam1./(epslm*c).^2).*(epslm*cm1);
        rlamm2 = rlam1.*rlam1m2./(epslm*c) + 0.5*(1 - rlam1.*rlam1./(epslm*c).^2).*(epslm*cm2);
        rlamm3 = rlam1.*rlam1m3./(epslm*c) + 0.5*(1 - rlam1.*rlam1./(epslm*c).^2).*(epslm*cm3);
        rlamm4 = rlam1.*rlam1m4./(epslm*c) + 0.5*(1 - rlam1.*rlam1./(epslm*c).^2).*(epslm*cm4);
        rlamm5 = rlam1.*rlam1m5./(epslm*c) + 0.5*(1 - rlam1.*rlam1./(epslm*c).^2).*(epslm*cm5);
        ic = rlam1 < epslm*c;
        rlam1 = ic.*rlam + (1-ic).*rlam1;
        rlam1m1 = ic.*rlamm1 + (1-ic).*rlam1m1;
        rlam1m2 = ic.*rlamm2 + (1-ic).*rlam1m2;
        rlam1m3 = ic.*rlamm3 + (1-ic).*rlam1m3;
        rlam1m4 = ic.*rlamm4 + (1-ic).*rlam1m4;
        rlam1m5 = ic.*rlamm5 + (1-ic).*rlam1m5;
    end
    
    rlam2   = abs(un-c);
    rlam2m1 = sign(un-c).*(unm1-cm1);
    rlam2m2 = sign(un-c).*(unm2-cm2);
    rlam2m3 = sign(un-c).*(unm3-cm3);
    rlam2m4 = sign(un-c).*(unm4-cm4);
    rlam2m5 = sign(un-c).*(unm5-cm5);

    if epslm>0
        rlam = 0.5*(rlam2.*rlam2./(epslm*c)+epslm*c);
        rlamm1 = rlam2.*rlam2m1./(epslm*c) + 0.5*(1 - rlam2.*rlam2./(epslm*c).^2).*(epslm*cm1);
        rlamm2 = rlam2.*rlam2m2./(epslm*c) + 0.5*(1 - rlam2.*rlam2./(epslm*c).^2).*(epslm*cm2);
        rlamm3 = rlam2.*rlam2m3./(epslm*c) + 0.5*(1 - rlam2.*rlam2./(epslm*c).^2).*(epslm*cm3);
        rlamm4 = rlam2.*rlam2m4./(epslm*c) + 0.5*(1 - rlam2.*rlam2./(epslm*c).^2).*(epslm*cm4);
        rlamm5 = rlam2.*rlam2m5./(epslm*c) + 0.5*(1 - rlam2.*rlam2./(epslm*c).^2).*(epslm*cm5);
        ic = rlam2 < epslm*c;
        rlam2 = ic.*rlam + (1-ic).*rlam2;
        rlam2m1 = ic.*rlamm1 + (1-ic).*rlam2m1;
        rlam2m2 = ic.*rlamm2 + (1-ic).*rlam2m2;
        rlam2m3 = ic.*rlamm3 + (1-ic).*rlam2m3;
        rlam2m4 = ic.*rlamm4 + (1-ic).*rlam2m4;
        rlam2m5 = ic.*rlamm5 + (1-ic).*rlam2m5;
    end
    
    rlam3   = abs(un);
    rlam3m1 = sign(un).*unm1;
    rlam3m2 = sign(un).*unm2;
    rlam3m3 = sign(un).*unm3;
    rlam3m4 = sign(un).*unm4; 
    rlam3m5 = sign(un).*unm5; 
    
    if epslm>0
        rlam = 0.5*(rlam3.*rlam3./(epslm*c)+epslm*c);
        rlamm1 = rlam3.*rlam3m1./(epslm*c) + 0.5*(1 - rlam3.*rlam3./(epslm*c).^2).*(epslm*cm1);
        rlamm2 = rlam3.*rlam3m2./(epslm*c) + 0.5*(1 - rlam3.*rlam3./(epslm*c).^2).*(epslm*cm2);
        rlamm3 = rlam3.*rlam3m3./(epslm*c) + 0.5*(1 - rlam3.*rlam3./(epslm*c).^2).*(epslm*cm3);
        rlamm4 = rlam3.*rlam3m4./(epslm*c) + 0.5*(1 - rlam3.*rlam3./(epslm*c).^2).*(epslm*cm4);
        rlamm5 = rlam3.*rlam3m5./(epslm*c) + 0.5*(1 - rlam3.*rlam3./(epslm*c).^2).*(epslm*cm5);
        ic = rlam3 < epslm*c;
        rlam3 = ic.*rlam + (1-ic).*rlam3;
        rlam3m1 = ic.*rlamm1 + (1-ic).*rlam3m1;
        rlam3m2 = ic.*rlamm2 + (1-ic).*rlam3m2;
        rlam3m3 = ic.*rlamm3 + (1-ic).*rlam3m3;
        rlam3m4 = ic.*rlamm4 + (1-ic).*rlam3m4;
        rlam3m5 = ic.*rlamm5 + (1-ic).*rlam3m5;
    end
else
    rlam1   = un+c;
    rlam1m1 = unm1+cm1;
    rlam1m2 = unm2+cm2;
    rlam1m3 = unm3+cm3;
    rlam1m4 = unm4+cm4;
    rlam1m5 = unm5+cm5;
    
    rlam2   = un-c;
    rlam2m1 = unm1-cm1;
    rlam2m2 = unm2-cm2;
    rlam2m3 = unm3-cm3;
    rlam2m4 = unm4-cm4;    
    rlam2m5 = unm5-cm5;    
    
    rlam3   = un;
    rlam3m1 = unm1;
    rlam3m2 = unm2;
    rlam3m3 = unm3;
    rlam3m4 = unm4;
    rlam3m5 = unm5;
end
%[min((un(:))) max((un(:))) min((c(:))) max((c(:)))]

s1      = 0.5*(rlam1   + rlam2);
s1m1    = 0.5*(rlam1m1 + rlam2m1);
s1m2    = 0.5*(rlam1m2 + rlam2m2);
s1m3    = 0.5*(rlam1m3 + rlam2m3);
s1m4    = 0.5*(rlam1m4 + rlam2m4);
s1m5    = 0.5*(rlam1m5 + rlam2m5);

s2      = 0.5*(rlam1   - rlam2);
s2m1    = 0.5*(rlam1m1 - rlam2m1);
s2m2    = 0.5*(rlam1m2 - rlam2m2);
s2m3    = 0.5*(rlam1m3 - rlam2m3);
s2m4    = 0.5*(rlam1m4 - rlam2m4);
s2m5    = 0.5*(rlam1m5 - rlam2m5);

An      = zeros(ng,nc,nc);
Anm     = zeros(ng,nc,nc,nc);

cc1   = gam1*(s1-rlam3).*af./c2-(s2.*un./c);
cc1m1 = gam1*((s1m1-rlam3m1).*af./c2 + (s1-rlam3).*afm1./c2 - (s1-rlam3).*af.*c2m1./c2.^2) - s2m1.*un./c - s2.*unm1./c + s2.*un.*cm1./c.^2;
cc1m2 = gam1*((s1m2-rlam3m2).*af./c2 + (s1-rlam3).*afm2./c2 - (s1-rlam3).*af.*c2m2./c2.^2) - s2m2.*un./c - s2.*unm2./c + s2.*un.*cm2./c.^2;
cc1m3 = gam1*((s1m3-rlam3m3).*af./c2 + (s1-rlam3).*afm3./c2 - (s1-rlam3).*af.*c2m3./c2.^2) - s2m3.*un./c - s2.*unm3./c + s2.*un.*cm3./c.^2;
cc1m4 = gam1*((s1m4-rlam3m4).*af./c2 + (s1-rlam3).*afm4./c2 - (s1-rlam3).*af.*c2m4./c2.^2) - s2m4.*un./c - s2.*unm4./c + s2.*un.*cm4./c.^2;
cc1m5 = gam1*((s1m5-rlam3m5).*af./c2 + (s1-rlam3).*afm5./c2 - (s1-rlam3).*af.*c2m5./c2.^2) - s2m5.*un./c - s2.*unm5./c + s2.*un.*cm5./c.^2;

cc2   = gam1*s2.*af./c-(s1-rlam3).*un;
cc2m1 = gam1*(s2m1.*af./c + s2.*afm1./c - s2.*af.*cm1./c.^2) - (s1m1-rlam3m1).*un - (s1-rlam3).*unm1;
cc2m2 = gam1*(s2m2.*af./c + s2.*afm2./c - s2.*af.*cm2./c.^2) - (s1m2-rlam3m2).*un - (s1-rlam3).*unm2;
cc2m3 = gam1*(s2m3.*af./c + s2.*afm3./c - s2.*af.*cm3./c.^2) - (s1m3-rlam3m3).*un - (s1-rlam3).*unm3;
cc2m4 = gam1*(s2m4.*af./c + s2.*afm4./c - s2.*af.*cm4./c.^2) - (s1m4-rlam3m4).*un - (s1-rlam3).*unm4;
cc2m5 = gam1*(s2m5.*af./c + s2.*afm5./c - s2.*af.*cm5./c.^2) - (s1m5-rlam3m5).*un - (s1-rlam3).*unm5;

An(:,:,1)  = [rlam3+cc1, cc1.*uv+cc2.*nx, cc1.*vv+cc2.*ny, cc1.*wv+cc2.*nz, cc1.*h+cc2.*un];
Anm(:,:,1,1) = [rlam3m1+cc1m1, cc1m1.*uv+cc1.*uvm1+cc2m1.*nx, cc1m1.*vv+cc1.*vvm1+cc2m1.*ny, cc1m1.*wv+cc1.*wvm1+cc2m1.*nz, cc1m1.*h+cc1.*hm1+cc2m1.*un+cc2.*unm1];
Anm(:,:,1,2) = [rlam3m2+cc1m2, cc1m2.*uv+cc1.*uvm2+cc2m2.*nx, cc1m2.*vv+cc1.*vvm2+cc2m2.*ny, cc1m2.*wv+cc1.*wvm2+cc2m2.*nz, cc1m2.*h+cc1.*hm2+cc2m2.*un+cc2.*unm2];
Anm(:,:,1,3) = [rlam3m3+cc1m3, cc1m3.*uv+cc1.*uvm3+cc2m3.*nx, cc1m3.*vv+cc1.*vvm3+cc2m3.*ny, cc1m3.*wv+cc1.*wvm3+cc2m3.*nz, cc1m3.*h+cc1.*hm3+cc2m3.*un+cc2.*unm3];
Anm(:,:,1,4) = [rlam3m4+cc1m4, cc1m4.*uv+cc1.*uvm4+cc2m4.*nx, cc1m4.*vv+cc1.*vvm4+cc2m4.*ny, cc1m4.*wv+cc1.*wvm4+cc2m4.*nz, cc1m4.*h+cc1.*hm4+cc2m4.*un+cc2.*unm4];
Anm(:,:,1,5) = [rlam3m5+cc1m5, cc1m5.*uv+cc1.*uvm5+cc2m5.*nx, cc1m5.*vv+cc1.*vvm5+cc2m5.*ny, cc1m5.*wv+cc1.*wvm5+cc2m5.*nz, cc1m5.*h+cc1.*hm5+cc2m5.*un+cc2.*unm5];

cc1   = -gam1*(s1-rlam3).*uv./c2+(s2.*nx./c);
cc1m1 = -gam1*((s1m1-rlam3m1).*uv./c2 + (s1-rlam3).*uvm1./c2 - (s1-rlam3).*uv.*c2m1./c2.^2) + s2m1.*nx./c - s2.*nx.*cm1./c.^2;
cc1m2 = -gam1*((s1m2-rlam3m2).*uv./c2 + (s1-rlam3).*uvm2./c2 - (s1-rlam3).*uv.*c2m2./c2.^2) + s2m2.*nx./c - s2.*nx.*cm2./c.^2;
cc1m3 = -gam1*((s1m3-rlam3m3).*uv./c2 + (s1-rlam3).*uvm3./c2 - (s1-rlam3).*uv.*c2m3./c2.^2) + s2m3.*nx./c - s2.*nx.*cm3./c.^2;
cc1m4 = -gam1*((s1m4-rlam3m4).*uv./c2 + (s1-rlam3).*uvm4./c2 - (s1-rlam3).*uv.*c2m4./c2.^2) + s2m4.*nx./c - s2.*nx.*cm4./c.^2;
cc1m5 = -gam1*((s1m5-rlam3m5).*uv./c2 + (s1-rlam3).*uvm5./c2 - (s1-rlam3).*uv.*c2m5./c2.^2) + s2m5.*nx./c - s2.*nx.*cm5./c.^2;

cc2   = -gam1*s2.*uv./c + (s1-rlam3).*nx;
cc2m1 = -gam1*(s2m1.*uv./c + s2.*uvm1./c - s2.*uv.*cm1./c.^2) + (s1m1-rlam3m1).*nx;
cc2m2 = -gam1*(s2m2.*uv./c + s2.*uvm2./c - s2.*uv.*cm2./c.^2) + (s1m2-rlam3m2).*nx;
cc2m3 = -gam1*(s2m3.*uv./c + s2.*uvm3./c - s2.*uv.*cm3./c.^2) + (s1m3-rlam3m3).*nx;
cc2m4 = -gam1*(s2m4.*uv./c + s2.*uvm4./c - s2.*uv.*cm4./c.^2) + (s1m4-rlam3m4).*nx;
cc2m5 = -gam1*(s2m5.*uv./c + s2.*uvm5./c - s2.*uv.*cm5./c.^2) + (s1m5-rlam3m5).*nx;

An(:,:,2)  = [cc1, rlam3+cc1.*uv+cc2.*nx, cc1.*vv+cc2.*ny, cc1.*wv+cc2.*nz, cc1.*h+cc2.*un];
Anm(:,:,2,1)  = [cc1m1, rlam3m1+cc1m1.*uv+cc1.*uvm1+cc2m1.*nx, cc1m1.*vv+cc1.*vvm1+cc2m1.*ny, cc1m1.*wv+cc1.*wvm1+cc2m1.*nz, cc1m1.*h+cc1.*hm1+cc2m1.*un+cc2.*unm1];
Anm(:,:,2,2)  = [cc1m2, rlam3m2+cc1m2.*uv+cc1.*uvm2+cc2m2.*nx, cc1m2.*vv+cc1.*vvm2+cc2m2.*ny, cc1m2.*wv+cc1.*wvm2+cc2m2.*nz, cc1m2.*h+cc1.*hm2+cc2m2.*un+cc2.*unm2];
Anm(:,:,2,3)  = [cc1m3, rlam3m3+cc1m3.*uv+cc1.*uvm3+cc2m3.*nx, cc1m3.*vv+cc1.*vvm3+cc2m3.*ny, cc1m3.*wv+cc1.*wvm3+cc2m3.*nz, cc1m3.*h+cc1.*hm3+cc2m3.*un+cc2.*unm3];
Anm(:,:,2,4)  = [cc1m4, rlam3m4+cc1m4.*uv+cc1.*uvm4+cc2m4.*nx, cc1m4.*vv+cc1.*vvm4+cc2m4.*ny, cc1m4.*wv+cc1.*wvm4+cc2m4.*nz, cc1m4.*h+cc1.*hm4+cc2m4.*un+cc2.*unm4];
Anm(:,:,2,5)  = [cc1m5, rlam3m5+cc1m5.*uv+cc1.*uvm5+cc2m5.*nx, cc1m5.*vv+cc1.*vvm5+cc2m5.*ny, cc1m5.*wv+cc1.*wvm5+cc2m5.*nz, cc1m5.*h+cc1.*hm5+cc2m5.*un+cc2.*unm5];

cc1   = -gam1*(s1-rlam3).*vv./c2+(s2.*ny./c);
cc1m1 = -gam1*((s1m1-rlam3m1).*vv./c2 + (s1-rlam3).*vvm1./c2 - (s1-rlam3).*vv.*c2m1./c2.^2) + s2m1.*ny./c - s2.*ny.*cm1./c.^2;
cc1m2 = -gam1*((s1m2-rlam3m2).*vv./c2 + (s1-rlam3).*vvm2./c2 - (s1-rlam3).*vv.*c2m2./c2.^2) + s2m2.*ny./c - s2.*ny.*cm2./c.^2;
cc1m3 = -gam1*((s1m3-rlam3m3).*vv./c2 + (s1-rlam3).*vvm3./c2 - (s1-rlam3).*vv.*c2m3./c2.^2) + s2m3.*ny./c - s2.*ny.*cm3./c.^2;
cc1m4 = -gam1*((s1m4-rlam3m4).*vv./c2 + (s1-rlam3).*vvm4./c2 - (s1-rlam3).*vv.*c2m4./c2.^2) + s2m4.*ny./c - s2.*ny.*cm4./c.^2;
cc1m5 = -gam1*((s1m5-rlam3m5).*vv./c2 + (s1-rlam3).*vvm5./c2 - (s1-rlam3).*vv.*c2m5./c2.^2) + s2m5.*ny./c - s2.*ny.*cm5./c.^2;

cc2   = -gam1*s2.*vv./c+(s1-rlam3).*ny;
cc2m1 = -gam1*(s2m1.*vv./c + s2.*vvm1./c - s2.*vv.*cm1./c.^2) + (s1m1-rlam3m1).*ny;
cc2m2 = -gam1*(s2m2.*vv./c + s2.*vvm2./c - s2.*vv.*cm2./c.^2) + (s1m2-rlam3m2).*ny;
cc2m3 = -gam1*(s2m3.*vv./c + s2.*vvm3./c - s2.*vv.*cm3./c.^2) + (s1m3-rlam3m3).*ny;
cc2m4 = -gam1*(s2m4.*vv./c + s2.*vvm4./c - s2.*vv.*cm4./c.^2) + (s1m4-rlam3m4).*ny;
cc2m5 = -gam1*(s2m5.*vv./c + s2.*vvm5./c - s2.*vv.*cm5./c.^2) + (s1m5-rlam3m5).*ny;

An(:,:,3)  = [cc1, cc1.*uv+cc2.*nx, rlam3+cc1.*vv+cc2.*ny, cc1.*wv+cc2.*nz, cc1.*h+cc2.*un];
Anm(:,:,3,1)  = [cc1m1, cc1m1.*uv+cc1.*uvm1+cc2m1.*nx, rlam3m1+cc1m1.*vv+cc1.*vvm1+cc2m1.*ny, cc1m1.*wv+cc1.*wvm1+cc2m1.*nz, cc1m1.*h+cc1.*hm1+cc2m1.*un+cc2.*unm1];
Anm(:,:,3,2)  = [cc1m2, cc1m2.*uv+cc1.*uvm2+cc2m2.*nx, rlam3m2+cc1m2.*vv+cc1.*vvm2+cc2m2.*ny, cc1m2.*wv+cc1.*wvm2+cc2m2.*nz, cc1m2.*h+cc1.*hm2+cc2m2.*un+cc2.*unm2];
Anm(:,:,3,3)  = [cc1m3, cc1m3.*uv+cc1.*uvm3+cc2m3.*nx, rlam3m3+cc1m3.*vv+cc1.*vvm3+cc2m3.*ny, cc1m3.*wv+cc1.*wvm3+cc2m3.*nz, cc1m3.*h+cc1.*hm3+cc2m3.*un+cc2.*unm3];
Anm(:,:,3,4)  = [cc1m4, cc1m4.*uv+cc1.*uvm4+cc2m4.*nx, rlam3m4+cc1m4.*vv+cc1.*vvm4+cc2m4.*ny, cc1m4.*wv+cc1.*wvm4+cc2m4.*nz, cc1m4.*h+cc1.*hm4+cc2m4.*un+cc2.*unm4];
Anm(:,:,3,5)  = [cc1m5, cc1m5.*uv+cc1.*uvm5+cc2m5.*nx, rlam3m5+cc1m5.*vv+cc1.*vvm5+cc2m5.*ny, cc1m5.*wv+cc1.*wvm5+cc2m5.*nz, cc1m5.*h+cc1.*hm5+cc2m5.*un+cc2.*unm5];

cc1   = -gam1*(s1-rlam3).*wv./c2+(s2.*nz./c);
cc1m1 = -gam1*((s1m1-rlam3m1).*wv./c2 + (s1-rlam3).*wvm1./c2 - (s1-rlam3).*wv.*c2m1./c2.^2) + s2m1.*nz./c - s2.*nz.*cm1./c.^2;
cc1m2 = -gam1*((s1m2-rlam3m2).*wv./c2 + (s1-rlam3).*wvm2./c2 - (s1-rlam3).*wv.*c2m2./c2.^2) + s2m2.*nz./c - s2.*nz.*cm2./c.^2;
cc1m3 = -gam1*((s1m3-rlam3m3).*wv./c2 + (s1-rlam3).*wvm3./c2 - (s1-rlam3).*wv.*c2m3./c2.^2) + s2m3.*nz./c - s2.*nz.*cm3./c.^2;
cc1m4 = -gam1*((s1m4-rlam3m4).*wv./c2 + (s1-rlam3).*wvm4./c2 - (s1-rlam3).*wv.*c2m4./c2.^2) + s2m4.*nz./c - s2.*nz.*cm4./c.^2;
cc1m5 = -gam1*((s1m5-rlam3m5).*wv./c2 + (s1-rlam3).*wvm5./c2 - (s1-rlam3).*wv.*c2m5./c2.^2) + s2m5.*nz./c - s2.*nz.*cm5./c.^2;

cc2   = -gam1*s2.*wv./c+(s1-rlam3).*nz;
cc2m1 = -gam1*(s2m1.*wv./c + s2.*wvm1./c - s2.*wv.*cm1./c.^2) + (s1m1-rlam3m1).*nz;
cc2m2 = -gam1*(s2m2.*wv./c + s2.*wvm2./c - s2.*wv.*cm2./c.^2) + (s1m2-rlam3m2).*nz;
cc2m3 = -gam1*(s2m3.*wv./c + s2.*wvm3./c - s2.*wv.*cm3./c.^2) + (s1m3-rlam3m3).*nz;
cc2m4 = -gam1*(s2m4.*wv./c + s2.*wvm4./c - s2.*wv.*cm4./c.^2) + (s1m4-rlam3m4).*nz;
cc2m5 = -gam1*(s2m5.*wv./c + s2.*wvm5./c - s2.*wv.*cm5./c.^2) + (s1m5-rlam3m5).*nz;

An(:,:,4)  = [cc1, cc1.*uv+cc2.*nx, cc1.*vv+cc2.*ny, rlam3+cc1.*wv+cc2.*nz, cc1.*h+cc2.*un];
Anm(:,:,4,1)  = [cc1m1, cc1m1.*uv+cc1.*uvm1+cc2m1.*nx, cc1m1.*vv+cc1.*vvm1+cc2m1.*ny, rlam3m1+cc1m1.*wv+cc1.*wvm1+cc2m1.*nz, cc1m1.*h+cc1.*hm1+cc2m1.*un+cc2.*unm1];
Anm(:,:,4,2)  = [cc1m2, cc1m2.*uv+cc1.*uvm2+cc2m2.*nx, cc1m2.*vv+cc1.*vvm2+cc2m2.*ny, rlam3m2+cc1m2.*wv+cc1.*wvm2+cc2m2.*nz, cc1m2.*h+cc1.*hm2+cc2m2.*un+cc2.*unm2];
Anm(:,:,4,3)  = [cc1m3, cc1m3.*uv+cc1.*uvm3+cc2m3.*nx, cc1m3.*vv+cc1.*vvm3+cc2m3.*ny, rlam3m3+cc1m3.*wv+cc1.*wvm3+cc2m3.*nz, cc1m3.*h+cc1.*hm3+cc2m3.*un+cc2.*unm3];
Anm(:,:,4,4)  = [cc1m4, cc1m4.*uv+cc1.*uvm4+cc2m4.*nx, cc1m4.*vv+cc1.*vvm4+cc2m4.*ny, rlam3m4+cc1m4.*wv+cc1.*wvm4+cc2m4.*nz, cc1m4.*h+cc1.*hm4+cc2m4.*un+cc2.*unm4];
Anm(:,:,4,5)  = [cc1m5, cc1m5.*uv+cc1.*uvm5+cc2m5.*nx, cc1m5.*vv+cc1.*vvm5+cc2m5.*ny, rlam3m5+cc1m5.*wv+cc1.*wvm5+cc2m5.*nz, cc1m5.*h+cc1.*hm5+cc2m5.*un+cc2.*unm5];

cc1   = gam1*(s1-rlam3)./c2;
cc1m1 = gam1*((s1m1-rlam3m1)./c2 - (s1-rlam3).*c2m1./c2.^2);
cc1m2 = gam1*((s1m2-rlam3m2)./c2 - (s1-rlam3).*c2m2./c2.^2);
cc1m3 = gam1*((s1m3-rlam3m3)./c2 - (s1-rlam3).*c2m3./c2.^2);
cc1m4 = gam1*((s1m4-rlam3m4)./c2 - (s1-rlam3).*c2m4./c2.^2);
cc1m5 = gam1*((s1m5-rlam3m5)./c2 - (s1-rlam3).*c2m5./c2.^2);

cc2   = gam1*s2./c;
cc2m1 = gam1*(s2m1./c - s2.*cm1./c.^2);
cc2m2 = gam1*(s2m2./c - s2.*cm2./c.^2);
cc2m3 = gam1*(s2m3./c - s2.*cm3./c.^2);
cc2m4 = gam1*(s2m4./c - s2.*cm4./c.^2);
cc2m5 = gam1*(s2m5./c - s2.*cm5./c.^2);

An(:,:,5)  = [cc1, cc1.*uv+cc2.*nx, cc1.*vv+cc2.*ny, cc1.*wv+cc2.*nz, rlam3+cc1.*h+cc2.*un];
Anm(:,:,5,1)  = [cc1m1, cc1m1.*uv+cc1.*uvm1+cc2m1.*nx, cc1m1.*vv+cc1.*vvm1+cc2m1.*ny, cc1m1.*wv+cc1.*wvm1+cc2m1.*nz, rlam3m1+cc1m1.*h+cc1.*hm1+cc2m1.*un+cc2.*unm1];
Anm(:,:,5,2)  = [cc1m2, cc1m2.*uv+cc1.*uvm2+cc2m2.*nx, cc1m2.*vv+cc1.*vvm2+cc2m2.*ny, cc1m2.*wv+cc1.*wvm2+cc2m2.*nz, rlam3m2+cc1m2.*h+cc1.*hm2+cc2m2.*un+cc2.*unm2];
Anm(:,:,5,3)  = [cc1m3, cc1m3.*uv+cc1.*uvm3+cc2m3.*nx, cc1m3.*vv+cc1.*vvm3+cc2m3.*ny, cc1m3.*wv+cc1.*wvm3+cc2m3.*nz, rlam3m3+cc1m3.*h+cc1.*hm3+cc2m3.*un+cc2.*unm3];
Anm(:,:,5,4)  = [cc1m4, cc1m4.*uv+cc1.*uvm4+cc2m4.*nx, cc1m4.*vv+cc1.*vvm4+cc2m4.*ny, cc1m4.*wv+cc1.*wvm4+cc2m4.*nz, rlam3m4+cc1m4.*h+cc1.*hm4+cc2m4.*un+cc2.*unm4];
Anm(:,:,5,5)  = [cc1m5, cc1m5.*uv+cc1.*uvm5+cc2m5.*nx, cc1m5.*vv+cc1.*vvm5+cc2m5.*ny, cc1m5.*wv+cc1.*wvm5+cc2m5.*nz, rlam3m5+cc1m5.*h+cc1.*hm5+cc2m5.*un+cc2.*unm5];

% Bn = [                                                                                                            0,                                                                                            nx,                                                                                            ny,                                                                                            nz,                               0;...
%  -(3*nx*ru^2 + nx*rv^2 + nx*rw^2 - gam*nx*ru^2 - gam*nx*rv^2 - gam*nx*rw^2 + 2*ny*ru*rv + 2*nz*ru*rw)/(2*r^2),                                                       (3*nx*ru + ny*rv + nz*rw - gam*nx*ru)/r,                                                               (ny*ru)/r - (nx*rv*(gam - 1))/r,                                                               (nz*ru)/r - (nx*rw*(gam - 1))/r,                    nx*(gam - 1);...
%  -(ny*ru^2 + 3*ny*rv^2 + ny*rw^2 - gam*ny*ru^2 - gam*ny*rv^2 - gam*ny*rw^2 + 2*nx*ru*rv + 2*nz*rv*rw)/(2*r^2),                                                               (nx*rv)/r - (ny*ru*(gam - 1))/r,                                                       (nx*ru + 3*ny*rv + nz*rw - gam*ny*rv)/r,                                                               (nz*rv)/r - (ny*rw*(gam - 1))/r,                    ny*(gam - 1);...
%  -(nz*ru^2 + nz*rv^2 + 3*nz*rw^2 - gam*nz*ru^2 - gam*nz*rv^2 - gam*nz*rw^2 + 2*nx*ru*rw + 2*ny*rv*rw)/(2*r^2),                                                               (nx*rw)/r - (nz*ru*(gam - 1))/r,                                                               (ny*rw)/r - (nz*rv*(gam - 1))/r,                                                       (nx*ru + ny*rv + 3*nz*rw - gam*nz*rw)/r,                    nz*(gam - 1);...
%           ((gam - 1)*(nx*ru + ny*rv + nz*rw)*(ru^2 + rv^2 + rw^2))/r^3 - (gam*rE*(nx*ru + ny*rv + nz*rw))/r^2, (gam*nx*rE)/r - ((gam - 1)*(3*nx*ru^2 + 2*ny*ru*rv + 2*nz*ru*rw + nx*rv^2 + nx*rw^2))/(2*r^2), (gam*ny*rE)/r - ((gam - 1)*(ny*ru^2 + 2*nx*ru*rv + 3*ny*rv^2 + 2*nz*rv*rw + ny*rw^2))/(2*r^2), (gam*nz*rE)/r - ((gam - 1)*(nz*ru^2 + 2*nx*ru*rw + nz*rv^2 + 2*ny*rv*rw + 3*nz*rw^2))/(2*r^2), (gam*(nx*ru + ny*rv + nz*rw))/r];
%  


