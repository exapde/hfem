function [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time)
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

[ng,nch] = size(uh);

u = udg(:,1:nch);

caseFlux = 1;
if caseFlux == 1
    [f,f_uhdg] = flux(p,uh,param);
elseif caseFlux==2
    [f,f_udg] = flux(p,u,param);
elseif caseFlux==3
    [f1,f1_uhdg] = flux(p,uh,param);
    [f2,f2_udg] = flux(p,u,param);
    f = 0.5*(f1+f2);
    f1_uhdg = 0.5*f1_uhdg;
    f2_udg = 0.5*f2_udg;
end

caseStab = 4;
if caseStab==1 % Roe
    [An,An_uh] = getan(nl,uh,param,1);
elseif caseStab==2 % Upwind 
    [Am,Amm] = getan(nl,uh,param,0);
    [Ap,Apm] = getan(nl,uh,param,1);
    An = 0.5*(Ap+Am);
    An_uh = 0.5*(Apm+Amm);
elseif caseStab==3 % Local Lax-Friedrich
    [An,An_uh] = getan(nl,uh,param,2);    
elseif caseStab==4 % global Lax-Friedrich
    An  = zeros(ng,nch,nch);
    for i=1:nch
        An(:,i,i)=param{end};
    end
    An_uh = zeros(ng,nch,nch,nch);
elseif caseStab==5
    [An,An_uh] = getan2(nl,uh,param);
end

fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);    

if caseFlux == 1
    fh_udg = An;
    fh_uh = permute(mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1)+mapContractK(An_uh,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;
elseif caseFlux==2
    fh_udg = permute(mapContractK(f_udg,nl,[2 4],3,1,2,[],1),[3 1 2])+An;
    fh_uh  = permute(mapContractK(An_uh,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;
elseif caseFlux==3
    fh_udg = permute(mapContractK(f2_udg,nl,[2 4],3,1,2,[],1),[3 1 2])+An;
    fh_uh  = permute(mapContractK(f1_uhdg,nl,[2 4],3,1,2,[],1)+mapContractK(An_uh,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;
end


function [An,Anm] = getan2(nl,m,param)
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

[ng,nc] = size(m);

gam   = param{1};
%epslm = param{2};

gam1 = gam - 1.0;


nx   = nl(:,1);              
ny   = nl(:,2);   

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

% b1 = 5;
% b2 = 0;
% bet   = b1*(abs(un) - b2*un);
% betm1 = b1*(sign(un).*unm1-b2*unm1);
% betm2 = b1*(sign(un).*unm2-b2*unm2);
% betm3 = b1*(sign(un).*unm3-b2*unm3);
% betm4 = b1*(sign(un).*unm4-b2*unm4);    

% b1 = 1;
% b2 = 0;
% bet   = b1*(abs(un+c) - b2*(un+c));
% betm1 = b1*(sign(un+c).*(unm1+cm1)-b2*(unm1+cm1));
% betm2 = b1*(sign(un+c).*(unm2+cm2)-b2*(unm2+cm2));
% betm3 = b1*(sign(un+c).*(unm3+cm3)-b2*(unm3+cm3));
% betm4 = b1*(sign(un+c).*(unm4+cm4)-b2*(unm4+cm4));    

% [y,dydx] = smoothramp1(un-c);
% b1 = 2;
% bet = b1*y;
% betm1 = b1*(dydx.*(unm1-cm1));
% betm2 = b1*(dydx.*(unm2-cm2));
% betm3 = b1*(dydx.*(unm3-cm3));
% betm4 = b1*(dydx.*(unm4-cm4));

% 0.5*(abs(un+c) + (un+c))
[y,dydx] = smoothramp1(un+c);
b1 = 2;
bet = b1*y;
betm1 = b1*(dydx.*(unm1+cm1));
betm2 = b1*(dydx.*(unm2+cm2));
betm3 = b1*(dydx.*(unm3+cm3));
betm4 = b1*(dydx.*(unm4+cm4));

%[y 0.5*(abs(un)+un)]
% [y,dydx] = smoothramp1(un);
% b1 = 4;
% bet = b1*y;
% betm1 = b1*(dydx.*(unm1));
% betm2 = b1*(dydx.*(unm2));
% betm3 = b1*(dydx.*(unm3));
% betm4 = b1*(dydx.*(unm4));

lam   = 1.0*(abs(un)+c);
lamm1 = 1.0*(sign(un).*unm1+cm1);
lamm2 = 1.0*(sign(un).*unm2+cm2);
lamm3 = 1.0*(sign(un).*unm3+cm3);
lamm4 = 1.0*(sign(un).*unm4+cm4);    

% if epslm>0
%     lam = 0.5*(lam.*lam./(epslm*c)+epslm*c);
%     lamm1 = lam.*lamm1./(epslm*c) + 0.5*(1 - lam.*lam./(epslm*c).^2).*(epslm*cm1);
%     lamm2 = lam.*lamm2./(epslm*c) + 0.5*(1 - lam.*lam./(epslm*c).^2).*(epslm*cm2);
%     lamm3 = lam.*lamm3./(epslm*c) + 0.5*(1 - lam.*lam./(epslm*c).^2).*(epslm*cm3);
%     lamm4 = lam.*lamm4./(epslm*c) + 0.5*(1 - lam.*lam./(epslm*c).^2).*(epslm*cm4);
% end    

An      = zeros(ng,nc,nc);
Anm     = zeros(ng,nc,nc,nc);

An(:,:,1)  = [bet, zer, zer, zer];
Anm(:,:,1,1) = [betm1, zer, zer, zer];
Anm(:,:,1,2) = [betm2, zer, zer, zer];
Anm(:,:,1,3) = [betm3, zer, zer, zer];
Anm(:,:,1,4) = [betm4, zer, zer, zer];

An(:,:,2)  = [zer, bet, zer, zer];
Anm(:,:,2,1)  = [zer, betm1, zer, zer];
Anm(:,:,2,2)  = [zer, betm2, zer, zer];
Anm(:,:,2,3)  = [zer, betm3, zer, zer];
Anm(:,:,2,4)  = [zer, betm4, zer, zer];

An(:,:,3)  = [zer, zer, bet, zer];
Anm(:,:,3,1)  = [zer, zer, betm1, zer];
Anm(:,:,3,2)  = [zer, zer, betm2, zer];
Anm(:,:,3,3)  = [zer, zer, betm3, zer];
Anm(:,:,3,4)  = [zer, zer, betm4, zer];

An(:,:,4)  = [zer, zer, zer, bet];
Anm(:,:,4,1)  = [zer, zer, zer, betm1];
Anm(:,:,4,2)  = [zer, zer, zer, betm2];
Anm(:,:,4,3)  = [zer, zer, zer, betm3];
Anm(:,:,4,4)  = [zer, zer, zer, betm4];

% An(:,:,1)  = [lam, zer, zer, zer];
% Anm(:,:,1,1) = [lamm1, zer, zer, zer];
% Anm(:,:,1,2) = [lamm2, zer, zer, zer];
% Anm(:,:,1,3) = [lamm3, zer, zer, zer];
% Anm(:,:,1,4) = [lamm4, zer, zer, zer];
% 
% An(:,:,2)  = [zer, lam, zer, zer];
% Anm(:,:,2,1)  = [zer, lamm1, zer, zer];
% Anm(:,:,2,2)  = [zer, lamm2, zer, zer];
% Anm(:,:,2,3)  = [zer, lamm3, zer, zer];
% Anm(:,:,2,4)  = [zer, lamm4, zer, zer];
% 
% An(:,:,3)  = [zer, zer, lam, zer];
% Anm(:,:,3,1)  = [zer, zer, lamm1, zer];
% Anm(:,:,3,2)  = [zer, zer, lamm2, zer];
% Anm(:,:,3,3)  = [zer, zer, lamm3, zer];
% Anm(:,:,3,4)  = [zer, zer, lamm4, zer];

% An(:,:,4)  = [zer, zer, zer, lam];
% Anm(:,:,4,1)  = [zer, zer, zer, lamm1];
% Anm(:,:,4,2)  = [zer, zer, zer, lamm2];
% Anm(:,:,4,3)  = [zer, zer, zer, lamm3];
% Anm(:,:,4,4)  = [zer, zer, zer, lamm4];

function [y,dydx] = smoothramp1(x) 

c = 10;
y = (x.*(atan(c*x)/pi + 1/2) + 1/2 - atan(c)/pi);
dydx = atan(c*x)/pi + (c*x)./(pi*(c^2*x.^2 + 1)) + 1/2;

function [y,dydx] = smoothramp2(x) 

c = 10;
y = (-x.*(atan(-c*x)/pi + 1/2) + 1/2 - atan(c)/pi);
dydx = atan(c*x)/pi + (c*x)./(pi*(c^2*x.^2 + 1)) - 1/2;






