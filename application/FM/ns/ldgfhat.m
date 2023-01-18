function fh = ldgfhat(nl,p,udg1,udg2,uh,param,time)
% fhg = fhat(nlg1, pg1, udg1, udg2, uhg, app.arg, app.time);
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
%      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q
% 
% f1 = flux(p,udg1,param);
% f2 = flux(p,udg2,param);
% f  = 0.5*(f1+f2);
% 
% 
% An = getan(nl,0.5*(udg1+udg2),param,1);
% fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + mapContractK(An,udg1-udg2,2,3,1,2,[],1),[2 1]);    

fh = euleri_roe(udg1,udg2,nl,p,param,time);

fv1 = viscousflux(p,udg1,param,time);
fv2 = viscousflux(p,udg2,param,time);
fv  = 0.5*(fv1+fv2);

tau = 0.1;
[~,nch] = size(uh);
fh = fh + permute(mapContractK(fv,nl,2,3,1,2,[],1),[2 1]) + tau*(udg1(:,1:nch)-udg2(:,1:nch));    


% function fh = ldgfhat(nl,p,udg1,udg2,uh,param,time)
% % fhg = fhat(nlg1, pg1, udg1, udg2, uhg, app.arg, app.time);
% %FHAT flux function
% %   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
% %
% %      NL(N,ND)              Normal N points
% %      P(N,ND)               Coordinates for N points
% %      U(N,NC)               Unknown vector for N points with NC components
% %      Q(N,NC,ND)            Flux vector for N points with NC components in the
% %                            coordinate directions
% %      M(N,NC)               Hybrid unkowns for N points with NC components
% %      PARAM                 Parameter list
% %      FH(N,NC):              Volume flux at N points
% %      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
% %      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
% %      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q
% 
% [ng,nc] = size(udg1);
% nch = 1;
% nq = nc-nch;
% nd = nq;
% 
% kappa = param{1};
% c     = param{2};
% tau   = param{end};
% 
% u1 = udg1(:,1);
% q1x = udg1(:,2);
% q1y = udg1(:,3);
% u2 = udg2(:,1);
% q2x = udg2(:,2);
% q2y = udg2(:,3);
% u  = 0.5*(u1+u2);
% qx = 0.5*(q1x+q2x);
% qy = 0.5*(q1y+q2y);
% 
% % ng x nch
% fh = kappa.*(qx.*nl(:,1)+qy.*nl(:,2)) + (nl*c(:)).*u + tau.*(u1(:,1)-u2(:,1));
% 


function An = getan(nl,m,param,absolute)

[ng,nc] = size(m);

gam   = param{1};
epslm = param{2};

gam1 = gam - 1.0;

nx   = nl(:,1);              
ny   = nl(:,2);

r    = m(:,1);
ru   = m(:,2);
rv   = m(:,3);
rE   = m(:,4);

r1   = 1./r;
uv   = ru.*r1;
vv   = rv.*r1;
E    = rE.*r1;
af   = 0.5*(uv.*uv+vv.*vv);
p    = gam1*(rE -r.*af);
h    = E   + p.*r1;
c2   = gam* p.*r1;
c    = sqrt(c2);
un   = uv.*nx   + vv.*ny;
if absolute
    rlam1   = abs(un+c);
    if epslm>0
        rlam = 0.5*(rlam1.*rlam1./(epslm*c)+epslm*c);
        ic = rlam1 < epslm*c;
        rlam1 = ic.*rlam + (1-ic).*rlam1;
    end
else
    rlam1   = un+c;
end

if absolute
    rlam2   = abs(un-c);
    if epslm>0
        rlam = 0.5*(rlam2.*rlam2./(epslm*c)+epslm*c);
        ic = rlam2 < epslm*c;
        rlam2 = ic.*rlam + (1-ic).*rlam2;
    end
else
    rlam2   = un-c;
end

if absolute
    rlam3   = abs(un);
    if epslm>0
        rlam = 0.5*(rlam3.*rlam3./(epslm*c)+epslm*c);
        ic = rlam3 < epslm*c;
        rlam3 = ic.*rlam + (1-ic).*rlam3;
    end
else
    rlam3   = un;
end

s1      = 0.5*(rlam1   + rlam2);
s2      = 0.5*(rlam1   - rlam2);

An      = zeros(ng,nc,nc);
cc1   = gam1*(s1-rlam3).*af./c2-(s2.*un./c);
cc2   = gam1*s2.*af./c-(s1-rlam3).*un;
An(:,:,1)  = [rlam3+cc1, cc1.*uv+cc2.*nx, cc1.*vv+cc2.*ny, cc1.*h+cc2.*un];
cc1   = -gam1*(s1-rlam3).*uv./c2+(s2.*nx./c);
cc2   = -gam1*s2.*uv./c + (s1-rlam3).*nx;
An(:,:,2)  = [cc1, rlam3+cc1.*uv+cc2.*nx, cc1.*vv+cc2.*ny, cc1.*h+cc2.*un];
cc1   = -gam1*(s1-rlam3).*vv./c2+(s2.*ny./c);
cc2   = -gam1*s2.*vv./c+(s1-rlam3).*ny;
An(:,:,3)  = [cc1, cc1.*uv+cc2.*nx, rlam3+cc1.*vv+cc2.*ny, cc1.*h+cc2.*un];
cc1   = gam1*(s1-rlam3)./c2;
cc2   = gam1*s2./c;
An(:,:,4)  = [cc1, cc1.*uv+cc2.*nx, cc1.*vv+cc2.*ny, rlam3+cc1.*h+cc2.*un];

function fn = euleri_roe(up,um,np,p,param,time)
%EULERI Calculate Interface Roe flux for the Euler equations.
%   FN=EULERI_ROE(UL,UR,N,P,PARAM,TIME)
%
%      UL(np,4):     np left (or plus) states
%      UR(np,4):     np right (or minus) states
%      NP(np,2):     np normal vectors (pointwing outwars the p element) 
%      P:            Not used
%      PARAM{1}:     Cell array containing the value of gamma
%      TIME:         Not used
%      FN(np,4):     np normal fluxes (f plus)    
%                          
% - Written by: J. Peraire
%

gam = param{1};
gam1  = gam - 1.0;
                                             
nx   = np(:,1);              
ny   = np(:,2);

rr   = um(:,1);            
rum  = um(:,2);
rvr  = um(:,3);
rEr  = um(:,4);

rl   = up(:,1);
rup  = up(:,2);
rvl  = up(:,3);
rEl  = up(:,4);

rr1  = 1./rr;
um   = rum.*rr1;
vr   = rvr.*rr1;
Er   = rEr.*rr1;
u2r  = um.*um+vr.*vr;
pr   = gam1*(rEr-0.5*rr.*u2r);
hr   = Er+pr.*rr1;
unr  = um.*nx+vr.*ny;

rl1  = 1./rl;
up   = rup.*rl1;
vl   = rvl.*rl1;
El   = rEl.*rl1;
u2l  = up.*up+vl.*vl;
pl   = gam1*(rEl-0.5*rl.*u2l);
hl   = El+pl.*rl1;
unl  = up.*nx+vl.*ny;
                                         
fn = 0.5*[(rr.*unr+rl.*unl), ...        
          (rum.*unr+rup.*unl)+nx.*(pr+pl), ...
          (rvr.*unr+rvl.*unl)+ny.*(pr+pl), ...
          (rr.*hr.*unr+rl.*hl.*unl)];

di   = sqrt(rr.*rl1);      
d1   = 1./(di+1);
ui   = (di.*um+up).*d1;
vi   = (di.*vr+vl).*d1;
hi   = (di.*hr+hl).*d1;
ci2  = gam1*(hi-0.5*(ui.*ui+vi.*vi));
ci   = sqrt(ci2);
af   = 0.5*(ui.*ui+vi.*vi);
uni  = ui.*nx+vi.*ny;

dr    = rr-rl;
dru   = rum-rup;
drv   = rvr-rvl;
drE   = rEr-rEl;

rlam1 = abs(uni+ci);
rlam2 = abs(uni-ci);
rlam3 = abs(uni);

s1    = 0.5*(rlam1+rlam2);
s2    = 0.5*(rlam1-rlam2);
al1x  = gam1*(af.*dr-ui.*dru-vi.*drv+drE);
al2x  = -uni.*dr+dru.*nx+drv.*ny;
cc1   = ((s1-rlam3).*al1x./ci2)+(s2.*al2x./ci);
cc2   = (s2.*al1x./ci)+(s1-rlam3).*al2x;
      
fn(:,1)  = fn(:,1) - 0.5*(rlam3.*dr+cc1);
fn(:,2)  = fn(:,2) - 0.5*(rlam3.*dru+cc1.*ui+cc2.*nx);
fn(:,3)  = fn(:,3) - 0.5*(rlam3.*drv+cc1.*vi+cc2.*ny);
fn(:,4)  = fn(:,4) - 0.5*(rlam3.*drE+cc1.*hi+cc2.*uni);

function fv = viscousflux(pg,udg,param,time)

[ng,nc] = size(udg);
nch = nc/3;
zero = zeros(ng,1);
%one  = ones(ng,1);

gam  = param{1};
gam1 = gam-1.0;
Re   = param{3};
Re1  = 1/Re;
Pr   = param{4};
Minf = param{5};
M2   = Minf^2;

c23  = 2/3;
fc = 1/(gam1*M2*Re*Pr);
                                             
r    = udg(:,1);
ru   = udg(:,2);
rv   = udg(:,3);
rE   = udg(:,4);

rx   = udg(:,5);
rux  = udg(:,6);
rvx  = udg(:,7);
rEx  = udg(:,8);

ry   = udg(:,9);
ruy  = udg(:,10);
rvy  = udg(:,11);
rEy  = udg(:,12);

r1   = 1./r;
u    = ru.*r1;
v    = rv.*r1;
%E    = rE.*r1;
q    = 0.5*(u.*u+v.*v);
p    = gam1*(rE-r.*q);
                        
ux  = (rux - rx.*u).*r1;
vx  = (rvx - rx.*v).*r1;
qx  = u.*ux + v.*vx;
px  = gam1*(rEx - rx.*q - r.*qx);
Tx  = gam*M2*(px.*r - p.*rx).*r1.^2;

uy  = (ruy - ry.*u).*r1;
vy  = (rvy - ry.*v).*r1;
qy  = u.*uy + v.*vy;
py  = gam1*(rEy - ry.*q - r.*qy);
Ty  = gam*M2*(py.*r - p.*ry).*r1.^2;

txx = Re1*c23*(2*ux - vy);
txy = Re1*(uy + vx);
tyy = Re1*c23*(2*vy - ux);

fv = zeros(ng,nch,2);
fv(:,:,1) = [zero, txx, txy, u.*txx + v.*txy + fc*Tx];
fv(:,:,2) = [zero, txy, tyy, u.*txy + v.*tyy + fc*Ty];

