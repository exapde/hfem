function fh = ldgfbou(ib,ui,nl,p,udg,uh,param,time)
% fhg = fbou(ib, uinf, nlg1, pg1, udg1, uhg, app.arg, app.time);
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

% if ib == 1                 % Far field
%     um = repmat(ui,size(udg,1),1);
% elseif  ib == 2            % Inviscid wall
%     un = udg(:,2).*nl(:,1)+udg(:,3).*nl(:,2);
%     um = [udg(:,1),udg(:,2)-2*un.*nl(:,1),udg(:,3)-2*un.*nl(:,2),udg(:,4)];
% end 
% 
% %fn = ldgfhat(up,um,np,p,param,time);
% fh = ldgfhat(nl,p,udg,um,uh,param,time);

up = udg;
%np = nl;
if ib == 1                 % Far field
    um = repmat(ui,size(up,1),1);
    fh = ldgfhatbou(nl,p,up,um,uh,param,time);
% elseif  ib == 2            % Reflect
%     un = up(:,2).*np(:,1)+up(:,3).*np(:,2);
%     um = [up(:,1),up(:,2)-2*un.*np(:,1),up(:,3)-2*un.*np(:,2),up(:,4)];
%     fh = ldgfhat(nl,p,up,um,uh,param,time);
%     uv = up(:,2,:)./up(:,1,:);
%     vv = up(:,3,:)./up(:,1,:);
%     gam = param{1};
%     p = (gam-1)*(up(:,4,:) - 0.5*(up(:,2,:).*uv + up(:,3,:).*vv));
%     fh = 0*up;
%     fh(:,2) = p.*nl(:,1);
%     fh(:,3) = p.*nl(:,2);
%     [f1 fh]
%     pause
elseif  ib == 2 % adiabatic wall    
    um = [up(:,1),0*up(:,2),0*up(:,3),up(:,4)];
    fh = ldgfhatbou(nl,p,up,um,uh,param,time);
end 
%fh = ldgfhat(nl,p,up,um,uh,param,time);

%fh = euleri_roe(up,um,np,p,param,time);

% %[ng,nc] = size(udg);
% nch = 1;
% 
% tau   = param{end};
% u     = udg(:,nch);
% switch ib
%     case 1  % Dirichlet                
%         c  = param{2};
%         qn = udg(:,2).*nl(:,1) + udg(:,3).*nl(:,2);        
%         fh = qn + (nl*c(:)).*u + tau*(u-ui);                
%     case 2  
%         fh = ui;        
%     case 3  % Prescribed flux
% %         x = p(:,1);
% %         y = p(:,2);
% %         ui = sin(x)*sin(y);            
%     otherwise
%         error('unknown boundary type');
% end
% 
% 
% % [ng,nc] = size(udg);
% % nch = 1;
% % nq = nc-nch;
% % nd = nq;
% % 
% % kappa = param{1};
% % tau   = param{end};
% % 
% % u1 = udg1(:,1);
% % q1x = udg1(:,2);
% % q1y = udg1(:,3);
% % u2 = udg2(:,1);
% % q2x = udg2(:,2);
% % q2y = udg2(:,3);
% % qx = 0.5*(q1x+q2x);
% % qy = 0.5*(q1y+q2y);
% % 
% % % ng x nch
% % fh = kappa.*(qx.*nl(:,1)+qy.*nl(:,2)) + tau.*(u1(:,1)-u2(:,1));
% % 
% % 
% % 

function fh = ldgfhatbou(nl,p,udg1,udg2,uh,param,time)
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

tau = 0.1;
[~,nch] = size(uh);
fh = fh + tau*(udg1(:,1:nch)-udg2(:,1:nch));    

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
