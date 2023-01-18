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