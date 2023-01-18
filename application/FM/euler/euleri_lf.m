function fn = euleri_lf(up,um,np,p,param,time)
%EULERI Calculate Interface Lax-Friedrichs flux for the Euler equations.
%   FN=EULERI_LF(UL,UR,N,P,PARAM,TIME)
%
%      UL(np,4):     np left (or plus) states
%      UR(np,4):     np right (or minus) states
%      N(np,2):      np normal vectors 
%      P:            Not used
%      PARAM{1}:     Cell array containing the value of gamma
%      TIME:         Not used
%      FN(np,4):     np normal fluxes (f plus)    
%                          
% - Written by: J. Peraire
%

gam = param{1};

uvl = up(:,2)./up(:,1);
vvl = up(:,3)./up(:,1);
pl = (gam-1)*(up(:,4) - 0.5*(up(:,2).*uvl + up(:,3).*vvl));
cl = sqrt(gam*pl./up(:,1));
u2l = sqrt(uvl.*uvl + vvl.*vvl);

fxl = [ up(:,2), up(:,2).*uvl+pl,    up(:,3).*uvl, uvl.*(up(:,4)+pl)];
fyl = [ up(:,3),    up(:,2).*vvl, up(:,3).*vvl+pl, vvl.*(up(:,4)+pl)];

uvr = um(:,2)./um(:,1);
vvr = um(:,3)./um(:,1);
pr = (gam-1)*(um(:,4) - 0.5*(um(:,2).*uvr + um(:,3).*vvr));
cr = sqrt(gam*pr./um(:,1));
u2r = sqrt(uvr.*uvr + vvr.*vvr);

fxr = [ um(:,2), um(:,2).*uvr+pr,    um(:,3).*uvr, uvr.*(um(:,4)+pr)];
fyr = [ um(:,3),    um(:,2).*vvr, um(:,3).*vvr+pr, vvr.*(um(:,4)+pr)];

fav = 0.5*diag(np(:,1))*(fxl+fxr) + 0.5*diag(np(:,2))*(fyl+fyr);

fn = fav + 0.25*diag(u2l + cl + u2r + cr)*(up - um);