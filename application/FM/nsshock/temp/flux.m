function [f,f_udg] = flux(pg,udg,param,time)
% FLUX Volume flux function
%   [f,fu,fq] = flux(p,u,q,param)
% 
%      P(N,ND)              Coordinates for N points
%      U(N,NC)              Unknown vector for N points with NC components
%      Q(N,NC,ND)           Flux vector for N points with NC components in the
%                           coordinate directions
%      PARAM                Parameter list
%      F(N,NC,ND):          Volume flux at N points
%      FU(N,NC,ND,NC):      Jacobian of the flux flux vector w.r.t. U
%      FQ(N,NC,ND,NC,ND):   Jacobian of the flux flux vector w.r.t. Q

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
E    = rE.*r1;
q    = 0.5*(u.*u+v.*v);
p    = gam1*(rE-r.*q);
h    = E+p.*r1;
                                        
[eps,eps_udg] = avflux(pg,udg,param,time);

wc=2;
if wc==1
    f = zeros(ng,nch,2);
    f(:,:,1) = [ru+eps.*rx, ru.*u+p+eps.*rux, rv.*u+eps.*rvx,   ru.*h+eps.*(rEx)];
    f(:,:,2) = [rv+eps.*ry, ru.*v+eps.*ruy,   rv.*v+p+eps.*rvy, rv.*h+eps.*(rEy)];

    f_u = zeros(ng,nch,2,nch);
    f_u(:,:,1,1) = -[-eps_udg(:,1).*rx, 0.5*((3-gam)*u.*u-gam1*v.*v)-eps_udg(:,1).*rux, u.*v-eps_udg(:,1).*rvx, gam*E.*u-2*gam1*u.*q-eps_udg(:,1).*rEx];
    f_u(:,:,1,2) = -[-ones(ng,1)-eps_udg(:,2).*rx, (gam-3)*u-eps_udg(:,2).*rux, -v-eps_udg(:,2).*rvx, -gam*E+0.5*gam1*(3*u.*u+v.*v)-eps_udg(:,2).*rEx];
    f_u(:,:,1,3) = -[-eps_udg(:,3).*rx, gam1*v-eps_udg(:,3).*rux, -u-eps_udg(:,3).*rvx, gam1*u.*v-eps_udg(:,3).*rEx];
    f_u(:,:,1,4) = -[-eps_udg(:,4).*rx, -gam1*ones(ng,1)-eps_udg(:,4).*rux, zeros(ng,1)-eps_udg(:,4).*rvx, -gam*u-eps_udg(:,4).*rEx];

    f_u(:,:,2,1) = -[-eps_udg(:,1).*ry, u.*v-eps_udg(:,1).*ruy, 0.5*((3-gam)*v.*v-gam1*u.*u)-eps_udg(:,1).*rvy, gam*E.*v-2*gam1*v.*q-eps_udg(:,1).*rEy];
    f_u(:,:,2,2) = -[-eps_udg(:,2).*ry, -v-eps_udg(:,2).*ruy, gam1*u-eps_udg(:,2).*rvy, gam1*u.*v-eps_udg(:,2).*rEy];
    f_u(:,:,2,3) = -[-ones(ng,1)-eps_udg(:,3).*ry, -u-eps_udg(:,3).*ruy, (gam-3)*v-eps_udg(:,3).*rvy,  -gam*E+0.5*gam1*(3*v.*v+u.*u)-eps_udg(:,3).*rEy];
    f_u(:,:,2,4) = -[-eps_udg(:,4).*ry, zeros(ng,1)-eps_udg(:,4).*ruy, -gam1*ones(ng,1)-eps_udg(:,4).*rvy, -gam*v-eps_udg(:,4).*rEy];

    f_q = zeros(ng,nch,2,nch,2);
    f_q(:,:,1,1,1) = [eps+eps_udg(:,5).*rx, eps_udg(:,5).*rux, eps_udg(:,5).*rvx, eps_udg(:,5).*rEx];
    f_q(:,:,1,2,1) = [eps_udg(:,6).*rx, eps+eps_udg(:,6).*rux, eps_udg(:,6).*rvx, eps_udg(:,6).*rEx];
    f_q(:,:,1,3,1) = [eps_udg(:,7).*rx, eps_udg(:,7).*rux, eps+eps_udg(:,7).*rvx, eps_udg(:,7).*rEx];
    f_q(:,:,1,4,1) = [eps_udg(:,8).*rx, eps_udg(:,8).*rux, eps_udg(:,8).*rvx, eps+eps_udg(:,8).*rEx];
    f_q(:,:,1,1,2) = [eps_udg(:,9).*rx, eps_udg(:,9).*rux , eps_udg(:,9).*rvx , eps_udg(:,9).*rEx];
    f_q(:,:,1,2,2) = [eps_udg(:,10).*rx, eps_udg(:,10).*rux , eps_udg(:,10).*rvx , eps_udg(:,10).*rEx];
    f_q(:,:,1,3,2) = [eps_udg(:,11).*rx, eps_udg(:,11).*rux , eps_udg(:,11).*rvx , eps_udg(:,11).*rEx];
    f_q(:,:,1,4,2) = [eps_udg(:,12).*rx, eps_udg(:,12).*rux , eps_udg(:,12).*rvx , eps_udg(:,12).*rEx];

    f_q(:,:,2,1,1) = [eps_udg(:,5).*ry, eps_udg(:,5).*ruy , eps_udg(:,5).*rvy , eps_udg(:,5).*rEy];
    f_q(:,:,2,2,1) = [eps_udg(:,6).*ry, eps_udg(:,6).*ruy , eps_udg(:,6).*rvy , eps_udg(:,6).*rEy];
    f_q(:,:,2,3,1) = [eps_udg(:,7).*ry, eps_udg(:,7).*ruy , eps_udg(:,7).*rvy , eps_udg(:,7).*rEy];
    f_q(:,:,2,4,1) = [eps_udg(:,8).*ry, eps_udg(:,8).*ruy , eps_udg(:,8).*rvy , eps_udg(:,8).*rEy];
    f_q(:,:,2,1,2) = [eps+eps_udg(:,9).*ry, eps_udg(:,9).*ruy , eps_udg(:,9).*rvy , eps_udg(:,9).*rEy];
    f_q(:,:,2,2,2) = [eps_udg(:,10).*ry, eps+eps_udg(:,10).*ruy , eps_udg(:,10).*rvy , eps_udg(:,10).*rEy];
    f_q(:,:,2,3,2) = [eps_udg(:,11).*ry, eps_udg(:,11).*ruy , eps+eps_udg(:,11).*rvy , eps_udg(:,11).*rEy];
    f_q(:,:,2,4,2) = [eps_udg(:,12).*ry, eps_udg(:,12).*ruy , eps_udg(:,12).*rvy , eps+eps_udg(:,12).*rEy];

elseif wc==2    
    ux  = (rux - rx.*u).*r1;
    vx  = (rvx - rx.*v).*r1;
    qx  = u.*ux + v.*vx;
    px  = gam1*(rEx - rx.*q - r.*qx);

    rHx     = rEx+px;
    rHx_r   = -gam1.*(r.*((ru.^2.*rx)./r.^4 + (rv.^2.*rx)./r.^4 - (2.*ru.*(rux - (ru.*rx)./r))./r.^3 - (2.*rv.*(rvx - (rv.*rx)./r))./r.^3) - rx.*(ru.^2./r.^3 + rv.^2./r.^3) + (ru.*(rux - (ru.*rx)./r))./r.^2 + (rv.*(rvx - (rv.*rx)./r))./r.^2);
    rHx_ru  = -gam1.*(r.*((rux - (ru.*rx)./r)./r.^2 - (ru.*rx)./r.^3) + (ru.*rx)./r.^2);
    rHx_rv  = -gam1.*(r.*((rvx - (rv.*rx)./r)./r.^2 - (rv.*rx)./r.^3) + (rv.*rx)./r.^2);
    rHx_rE  = 0;
    rHx_rx  = -gam1.*(ru.^2./(2.*r.^2) + rv.^2./(2.*r.^2) - r.*(ru.^2./r.^3 + rv.^2./r.^3));
    rHx_rux = -(gam1.*ru)./r;
    rHx_rvx = -(gam1.*rv)./r;
    rHx_rEx = gam1 + 1;
    rHx_ry  = 0;
    rHx_ruy = 0;
    rHx_rvy = 0;
    rHx_rEy = 0;

    uy  = (ruy - ry.*u).*r1;
    vy  = (rvy - ry.*v).*r1;
    qy  = u.*uy + v.*vy;
    py  = gam1*(rEy - ry.*q - r.*qy);

    rHy     = rEy+py;
    rHy_r   = -gam1.*(r.*((ru.^2.*ry)./r.^4 + (rv.^2.*ry)./r.^4 - (2.*ru.*(ruy - (ru.*ry)./r))./r.^3 - (2.*rv.*(rvy - (rv.*ry)./r))./r.^3) - ry.*(ru.^2./r.^3 + rv.^2./r.^3) + (ru.*(ruy - (ru.*ry)./r))./r.^2 + (rv.*(rvy - (rv.*ry)./r))./r.^2);
    rHy_ru  = -gam1.*(r.*((ruy - (ru.*ry)./r)./r.^2 - (ru.*ry)./r.^3) + (ru.*ry)./r.^2);
    rHy_rv  = -gam1.*(r.*((rvy - (rv.*ry)./r)./r.^2 - (rv.*ry)./r.^3) + (rv.*ry)./r.^2);
    rHy_rE  = 0;
    rHy_rx  = 0;
    rHy_rux = 0;
    rHy_rvx = 0;
    rHy_rEx = 0;
    rHy_ry  = -gam1.*(ru.^2./(2.*r.^2) + rv.^2./(2.*r.^2) - r.*(ru.^2./r.^3 + rv.^2./r.^3));
    rHy_ruy = -(gam1.*ru)./r;
    rHy_rvy = -(gam1.*rv)./r;
    rHy_rEy = gam1 + 1;

    f = zeros(ng,nch,2);
    f(:,:,1) = [ru+eps.*rx, ru.*u+p+eps.*rux, rv.*u+eps.*rvx,   ru.*h+eps.*(rHx)];
    f(:,:,2) = [rv+eps.*ry, ru.*v+eps.*ruy,   rv.*v+p+eps.*rvy, rv.*h+eps.*(rHy)];

    f_u = zeros(ng,nch,2,nch);
    f_u(:,:,1,1) = -[-eps_udg(:,1).*rx, 0.5*((3-gam)*u.*u-gam1*v.*v)-eps_udg(:,1).*rux, u.*v-eps_udg(:,1).*rvx, gam*E.*u-2*gam1*u.*q-eps_udg(:,1).*rHx-eps.*rHx_r];
    f_u(:,:,1,2) = -[-ones(ng,1)-eps_udg(:,2).*rx, (gam-3)*u-eps_udg(:,2).*rux, -v-eps_udg(:,2).*rvx, -gam*E+0.5*gam1*(3*u.*u+v.*v)-eps_udg(:,2).*rHx-eps.*rHx_ru];
    f_u(:,:,1,3) = -[-eps_udg(:,3).*rx, gam1*v-eps_udg(:,3).*rux, -u-eps_udg(:,3).*rvx, gam1*u.*v-eps_udg(:,3).*rHx-eps.*rHx_rv];
    f_u(:,:,1,4) = -[-eps_udg(:,4).*rx, -gam1*ones(ng,1)-eps_udg(:,4).*rux, zeros(ng,1)-eps_udg(:,4).*rvx, -gam*u-eps_udg(:,4).*rHx-eps.*rHx_rE];

    f_u(:,:,2,1) = -[-eps_udg(:,1).*ry, u.*v-eps_udg(:,1).*ruy, 0.5*((3-gam)*v.*v-gam1*u.*u)-eps_udg(:,1).*rvy, gam*E.*v-2*gam1*v.*q-eps_udg(:,1).*rHy-eps.*rHy_r];
    f_u(:,:,2,2) = -[-eps_udg(:,2).*ry, -v-eps_udg(:,2).*ruy, gam1*u-eps_udg(:,2).*rvy, gam1*u.*v-eps_udg(:,2).*rHy-eps.*rHy_ru];
    f_u(:,:,2,3) = -[-ones(ng,1)-eps_udg(:,3).*ry, -u-eps_udg(:,3).*ruy, (gam-3)*v-eps_udg(:,3).*rvy,  -gam*E+0.5*gam1*(3*v.*v+u.*u)-eps_udg(:,3).*rHy-eps.*rHy_rv];
    f_u(:,:,2,4) = -[-eps_udg(:,4).*ry, zeros(ng,1)-eps_udg(:,4).*ruy, -gam1*ones(ng,1)-eps_udg(:,4).*rvy, -gam*v-eps_udg(:,4).*rHy-eps.*rHy_rE];

    f_q = zeros(ng,nch,2,nch,2);
    f_q(:,:,1,1,1) = [eps+eps_udg(:,5).*rx, eps_udg(:,5).*rux, eps_udg(:,5).*rvx, eps_udg(:,5).*rHx+eps.*rHx_rx];
    f_q(:,:,1,2,1) = [eps_udg(:,6).*rx, eps+eps_udg(:,6).*rux, eps_udg(:,6).*rvx, eps_udg(:,6).*rHx+eps.*rHx_rux];
    f_q(:,:,1,3,1) = [eps_udg(:,7).*rx, eps_udg(:,7).*rux, eps+eps_udg(:,7).*rvx, eps_udg(:,7).*rHx+eps.*rHx_rvx];
    f_q(:,:,1,4,1) = [eps_udg(:,8).*rx, eps_udg(:,8).*rux, eps_udg(:,8).*rvx, eps_udg(:,8).*rHx+eps.*rHx_rEx];
    f_q(:,:,1,1,2) = [eps_udg(:,9).*rx, eps_udg(:,9).*rux, eps_udg(:,9).*rvx , eps_udg(:,9).*rHx+eps.*rHx_ry];
    f_q(:,:,1,2,2) = [eps_udg(:,10).*rx, eps_udg(:,10).*rux , eps_udg(:,10).*rvx , eps_udg(:,10).*rHx+eps.*rHx_ruy];
    f_q(:,:,1,3,2) = [eps_udg(:,11).*rx, eps_udg(:,11).*rux , eps_udg(:,11).*rvx , eps_udg(:,11).*rHx+eps.*rHx_rvy];
    f_q(:,:,1,4,2) = [eps_udg(:,12).*rx, eps_udg(:,12).*rux , eps_udg(:,12).*rvx , eps_udg(:,12).*rHx+eps.*rHx_rEy];

    f_q(:,:,2,1,1) = [eps_udg(:,5).*ry, eps_udg(:,5).*ruy , eps_udg(:,5).*rvy , eps_udg(:,5).*rHy+eps.*rHy_rx];
    f_q(:,:,2,2,1) = [eps_udg(:,6).*ry, eps_udg(:,6).*ruy , eps_udg(:,6).*rvy , eps_udg(:,6).*rHy+eps.*rHy_rux];
    f_q(:,:,2,3,1) = [eps_udg(:,7).*ry, eps_udg(:,7).*ruy , eps_udg(:,7).*rvy , eps_udg(:,7).*rHy+eps.*rHy_rvx];
    f_q(:,:,2,4,1) = [eps_udg(:,8).*ry, eps_udg(:,8).*ruy , eps_udg(:,8).*rvy , eps_udg(:,8).*rHy+eps.*rHy_rEx];
    f_q(:,:,2,1,2) = [eps+eps_udg(:,9).*ry, eps_udg(:,9).*ruy , eps_udg(:,9).*rvy , eps_udg(:,9).*rHy+eps.*rHy_ry];
    f_q(:,:,2,2,2) = [eps_udg(:,10).*ry, eps+eps_udg(:,10).*ruy , eps_udg(:,10).*rvy , eps_udg(:,10).*rHy+eps.*rHy_ruy];
    f_q(:,:,2,3,2) = [eps_udg(:,11).*ry, eps_udg(:,11).*ruy , eps+eps_udg(:,11).*rvy , eps_udg(:,11).*rHy+eps.*rHy_rvy];
    f_q(:,:,2,4,2) = [eps_udg(:,12).*ry, eps_udg(:,12).*ruy , eps_udg(:,12).*rvy , eps_udg(:,12).*rHy+eps.*rHy_rEy];    
end


u_r  = -u.*r1;
u_ru =  r1;
v_r  = -v.*r1;
v_rv =  r1;

ux  = (rux - rx.*u).*r1;
vx  = (rvx - rx.*v).*r1;
%Ex  = (rEx - rx.*E).*r1;
qx  = u.*ux + v.*vx;
px  = gam1*(rEx - rx.*q - r.*qx);
Tx  = gam*M2*(px.*r - p.*rx).*r1.^2;

uy  = (ruy - ry.*u).*r1;
vy  = (rvy - ry.*v).*r1;
%Ey  = (rEy - ry.*E).*r1;
qy  = u.*uy + v.*vy;
py  = gam1*(rEy - ry.*q - r.*qy);
Ty  = gam*M2*(py.*r - p.*ry).*r1.^2;

txx = Re1*c23*(2*ux - vy);
txy = Re1*(uy + vx);
tyy = Re1*c23*(2*vy - ux);

% Tx  = gam*M2*(px.*r - p.*rx).*r1.^2;
% Ty  = gam*M2*(py.*r - p.*ry).*r1.^2;

fv = zeros(ng,nch,2);
fv(:,:,1) = [zero, txx, txy, u.*txx + v.*txy + fc*Tx];
fv(:,:,2) = [zero, txy, tyy, u.*txy + v.*tyy + fc*Ty];

txx_r  =  Re1*c23*((4*ru.*rx-2*rv.*ry)-r.*(2*rux-rvy)).*r1.^3;
txx_ru = -Re1*c23*2*rx.*r1.^2;
txx_rv =  Re1*c23*ry.*r1.^2;
txx_rE =  zero;

txx_rx  = -Re1*c23*2*ru.*r1.^2;
txx_rux =  Re1*c23*2*r1;
txx_rvx =  zero;
txx_rEx =  zero;

txx_ry  =  Re1*c23*rv.*r1.^2;
txx_ruy =  zero;
txx_rvy = -Re1*c23.*r1;
txx_rEy =  zero;

txy_r  =  Re1*(2*(ru.*ry+rv.*rx)-r.*(ruy+rvx)).*r1.^3;
txy_ru = -Re1*ry.*r1.^2;
txy_rv = -Re1*rx.*r1.^2;
txy_rE =  zero;

txy_rx  = -Re1*rv.*r1.^2;
txy_rux =  zero;
txy_rvx =  Re1*r1;
txy_rEx =  zero;

txy_ry  = -Re1*ru.*r1.^2;
txy_ruy =  Re1*r1;
txy_rvy =  zero;
txy_rEy =  zero;

tyy_r  =  Re1*c23*((4*rv.*ry-2*ru.*rx)-r.*(2*rvy-rux)).*r1.^3;
tyy_ru =  Re1*c23*rx.*r1.^2;
tyy_rv = -Re1*c23*2*ry.*r1.^2;
tyy_rE =  zero;

tyy_rx  =  Re1*c23*ru.*r1.^2;
tyy_rux = -Re1*c23*r1;
tyy_rvx =  zero;
tyy_rEx =  zero;

tyy_ry  = -Re1*c23*2*rv.*r1.^2;
tyy_ruy =  zero;
tyy_rvy =  Re1*c23*2*r1;
tyy_rEy =  zero;

Tx_r  = -M2*gam*gam1*(rEx.*r.^2-2*rux.*r.*ru-2*rvx.*r.*rv-2*rE.*rx.*r+3*rx.*(ru.^2+rv.^2)).*r1.^4;
Tx_ru = -M2*gam*gam1*(r.*rux-2*ru.*rx).*r1.^3;
Tx_rv = -M2*gam*gam1*(r.*rvx-2*rv.*rx).*r1.^3;
Tx_rE = -M2*gam*gam1*rx.*r1.^2;

Tx_rx  =  M2*gam*gam1*(ru.^2+rv.^2-r.*rE).*r1.^3;
Tx_rux = -M2*gam*gam1*ru.*r1.^2;
Tx_rvx = -M2*gam*gam1*rv.*r1.^2;
Tx_rEx =  M2*gam*gam1*r1;

Ty_r  = -M2*gam*gam1*(rEy.*r.^2-2*ruy.*r.*ru-2*rvy.*r.*rv-2*rE.*ry.*r+3*ry.*(ru.^2+rv.^2)).*r1.^4;
Ty_ru = -M2*gam*gam1*(r.*ruy-2*ru.*ry).*r1.^3;
Ty_rv = -M2*gam*gam1*(r.*rvy-2*rv.*ry).*r1.^3;
Ty_rE = -M2*gam*gam1*ry.*r1.^2;

Ty_ry  = Tx_rx;
Ty_ruy = Tx_rux;
Ty_rvy = Tx_rvx;
Ty_rEy = Tx_rEx;

% fv(:,:,1) = [zero, txx, txy, u.*txx + v.*txy + fc*Tx];
% fv(:,:,2) = [zero, txy, tyy, u.*txy + v.*tyy + fc*Ty];

fv_u = zeros(ng,nch,2,nch);
fv_u(:,:,1,1) = [zero, txx_r , txy_r , u_r.*txx  + u.*txx_r  + v_r.*txy  + v.*txy_r  + fc*Tx_r ];
fv_u(:,:,1,2) = [zero, txx_ru, txy_ru, u_ru.*txx + u.*txx_ru             + v.*txy_ru + fc*Tx_ru];
fv_u(:,:,1,3) = [zero, txx_rv, txy_rv,             u.*txx_rv + v_rv.*txy + v.*txy_rv + fc*Tx_rv];
fv_u(:,:,1,4) = [zero, txx_rE, txy_rE,             u.*txx_rE             + v.*txy_rE + fc*Tx_rE];

fv_u(:,:,2,1) = [zero, txy_r , tyy_r , u_r.*txy  + u.*txy_r  + v_r.*tyy  + v.*tyy_r  + fc*Ty_r];
fv_u(:,:,2,2) = [zero, txy_ru, tyy_ru, u_ru.*txy + u.*txy_ru             + v.*tyy_ru + fc*Ty_ru];
fv_u(:,:,2,3) = [zero, txy_rv, tyy_rv,             u.*txy_rv + v_rv.*tyy + v.*tyy_rv + fc*Ty_rv];
fv_u(:,:,2,4) = [zero, txy_rE, tyy_rE,             u.*txy_rE             + v.*tyy_rE + fc*Ty_rE];

fv_q = zeros(ng,nch,2,nch,2);
fv_q(:,:,1,1,1) = [zero, txx_rx , txy_rx , u.*txx_rx  + v.*txy_rx  + fc*Tx_rx ];
fv_q(:,:,1,2,1) = [zero, txx_rux, txy_rux, u.*txx_rux + v.*txy_rux + fc*Tx_rux];
fv_q(:,:,1,3,1) = [zero, txx_rvx, txy_rvx, u.*txx_rvx + v.*txy_rvx + fc*Tx_rvx];
fv_q(:,:,1,4,1) = [zero, txx_rEx, txy_rEx, u.*txx_rEx + v.*txy_rEx + fc*Tx_rEx];
fv_q(:,:,1,1,2) = [zero, txx_ry , txy_ry , u.*txx_ry  + v.*txy_ry             ];
fv_q(:,:,1,2,2) = [zero, txx_ruy, txy_ruy, u.*txx_ruy + v.*txy_ruy            ];
fv_q(:,:,1,3,2) = [zero, txx_rvy, txy_rvy, u.*txx_rvy + v.*txy_rvy            ];
fv_q(:,:,1,4,2) = [zero, txx_rEy, txy_rEy, u.*txx_rEy + v.*txy_rEy            ];

fv_q(:,:,2,1,1) = [zero, txy_rx , tyy_rx , u.*txy_rx  + v.*tyy_rx             ];
fv_q(:,:,2,2,1) = [zero, txy_rux, tyy_rux, u.*txy_rux + v.*tyy_rux            ];
fv_q(:,:,2,3,1) = [zero, txy_rvx, tyy_rvx, u.*txy_rvx + v.*tyy_rvx            ];
fv_q(:,:,2,4,1) = [zero, txy_rEx, tyy_rEx, u.*txy_rEx + v.*tyy_rEx            ];
fv_q(:,:,2,1,2) = [zero, txy_ry , tyy_ry , u.*txy_ry  + v.*tyy_ry  + fc*Ty_ry ];
fv_q(:,:,2,2,2) = [zero, txy_ruy, tyy_ruy, u.*txy_ruy + v.*tyy_ruy + fc*Ty_ruy];
fv_q(:,:,2,3,2) = [zero, txy_rvy, tyy_rvy, u.*txy_rvy + v.*tyy_rvy + fc*Ty_rvy];
fv_q(:,:,2,4,2) = [zero, txy_rEy, tyy_rEy, u.*txy_rEy + v.*tyy_rEy + fc*Ty_rEy];

f = f+fv;
f_udg = cat(4,f_u+fv_u,reshape(f_q+fv_q,ng,nch,2,2*nch));