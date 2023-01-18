function [f,f_udg] = flux(pg,udg,param,time)
%FLUX Volume flux function
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
nch  = 4;
zero = zeros(ng,1);

gam  = param{1};
gam1 = gam - 1.0;

r    = udg(:,1);
ru   = udg(:,2);
rv   = udg(:,3);
rE   = udg(:,4);

if min(r(:))<=1e-3
    error('negative density');
end

rx   = udg(:,5);
rux  = udg(:,6);
rvx  = udg(:,7);
rEx  = udg(:,8);

ry   = udg(:,9);
ruy  = udg(:,10);
rvy  = udg(:,11);
rEy  = udg(:,12);

r1   = 1./r;
uv   = ru.*r1;
vv   = rv.*r1;
E    = rE.*r1;
af   = 0.5*(uv.*uv+vv.*vv);
p    = gam1*(rE-r.*af);
h    = E+p.*r1;

[eps,eps_udg] = avflux(pg,udg,param,time);

wc=2;
if wc==1
    f = zeros(ng,nch,2);
    f(:,:,1) = [ru+eps.*rx, ru.*uv+p+eps.*rux, rv.*uv+eps.*rvx,   ru.*h+eps.*(rEx)];
    f(:,:,2) = [rv+eps.*ry, ru.*vv+eps.*ruy,   rv.*vv+p+eps.*rvy, rv.*h+eps.*(rEy)];

    f_u = zeros(ng,nch,2,nch);
    f_u(:,:,1,1) = -[-eps_udg(:,1).*rx, 0.5*((3-gam)*uv.*uv-gam1*vv.*vv)-eps_udg(:,1).*rux, uv.*vv-eps_udg(:,1).*rvx, gam*E.*uv-2*gam1*uv.*af-eps_udg(:,1).*rEx];
    f_u(:,:,1,2) = -[-ones(ng,1)-eps_udg(:,2).*rx, (gam-3)*uv-eps_udg(:,2).*rux, -vv-eps_udg(:,2).*rvx, -gam*E+0.5*gam1*(3*uv.*uv+vv.*vv)-eps_udg(:,2).*rEx];
    f_u(:,:,1,3) = -[-eps_udg(:,3).*rx, gam1*vv-eps_udg(:,3).*rux, -uv-eps_udg(:,3).*rvx, gam1*uv.*vv-eps_udg(:,3).*rEx];
    f_u(:,:,1,4) = -[-eps_udg(:,4).*rx, -gam1*ones(ng,1)-eps_udg(:,4).*rux, zeros(ng,1)-eps_udg(:,4).*rvx, -gam*uv-eps_udg(:,4).*rEx];

    f_u(:,:,2,1) = -[-eps_udg(:,1).*ry, uv.*vv-eps_udg(:,1).*ruy, 0.5*((3-gam)*vv.*vv-gam1*uv.*uv)-eps_udg(:,1).*rvy, gam*E.*vv-2*gam1*vv.*af-eps_udg(:,1).*rEy];
    f_u(:,:,2,2) = -[-eps_udg(:,2).*ry, -vv-eps_udg(:,2).*ruy, gam1*uv-eps_udg(:,2).*rvy, gam1*uv.*vv-eps_udg(:,2).*rEy];
    f_u(:,:,2,3) = -[-ones(ng,1)-eps_udg(:,3).*ry, -uv-eps_udg(:,3).*ruy, (gam-3)*vv-eps_udg(:,3).*rvy,  -gam*E+0.5*gam1*(3*vv.*vv+uv.*uv)-eps_udg(:,3).*rEy];
    f_u(:,:,2,4) = -[-eps_udg(:,4).*ry, zeros(ng,1)-eps_udg(:,4).*ruy, -gam1*ones(ng,1)-eps_udg(:,4).*rvy, -gam*vv-eps_udg(:,4).*rEy];

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
    q   = 0.5*(uv.*uv+vv.*vv);

    ux  = (rux - rx.*uv).*r1;
    vx  = (rvx - rx.*vv).*r1;
    qx  = uv.*ux + vv.*vx;
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

    uy  = (ruy - ry.*uv).*r1;
    vy  = (rvy - ry.*vv).*r1;
    qy  = uv.*uy + vv.*vy;
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
    f(:,:,1) = [ru+eps.*rx, ru.*uv+p+eps.*rux, rv.*uv+eps.*rvx,   ru.*h+eps.*(rHx)];
    f(:,:,2) = [rv+eps.*ry, ru.*vv+eps.*ruy,   rv.*vv+p+eps.*rvy, rv.*h+eps.*(rHy)];

    f_u = zeros(ng,nch,2,nch);
    f_u(:,:,1,1) = -[-eps_udg(:,1).*rx, 0.5*((3-gam)*uv.*uv-gam1*vv.*vv)-eps_udg(:,1).*rux, uv.*vv-eps_udg(:,1).*rvx, gam*E.*uv-2*gam1*uv.*af-eps_udg(:,1).*rHx-eps.*rHx_r];
    f_u(:,:,1,2) = -[-ones(ng,1)-eps_udg(:,2).*rx, (gam-3)*uv-eps_udg(:,2).*rux, -vv-eps_udg(:,2).*rvx, -gam*E+0.5*gam1*(3*uv.*uv+vv.*vv)-eps_udg(:,2).*rHx-eps.*rHx_ru];
    f_u(:,:,1,3) = -[-eps_udg(:,3).*rx, gam1*vv-eps_udg(:,3).*rux, -uv-eps_udg(:,3).*rvx, gam1*uv.*vv-eps_udg(:,3).*rHx-eps.*rHx_rv];
    f_u(:,:,1,4) = -[-eps_udg(:,4).*rx, -gam1*ones(ng,1)-eps_udg(:,4).*rux, zeros(ng,1)-eps_udg(:,4).*rvx, -gam*uv-eps_udg(:,4).*rHx-eps.*rHx_rE];

    f_u(:,:,2,1) = -[-eps_udg(:,1).*ry, uv.*vv-eps_udg(:,1).*ruy, 0.5*((3-gam)*vv.*vv-gam1*uv.*uv)-eps_udg(:,1).*rvy, gam*E.*vv-2*gam1*vv.*af-eps_udg(:,1).*rHy-eps.*rHy_r];
    f_u(:,:,2,2) = -[-eps_udg(:,2).*ry, -vv-eps_udg(:,2).*ruy, gam1*uv-eps_udg(:,2).*rvy, gam1*uv.*vv-eps_udg(:,2).*rHy-eps.*rHy_ru];
    f_u(:,:,2,3) = -[-ones(ng,1)-eps_udg(:,3).*ry, -uv-eps_udg(:,3).*ruy, (gam-3)*vv-eps_udg(:,3).*rvy,  -gam*E+0.5*gam1*(3*vv.*vv+uv.*uv)-eps_udg(:,3).*rHy-eps.*rHy_rv];
    f_u(:,:,2,4) = -[-eps_udg(:,4).*ry, zeros(ng,1)-eps_udg(:,4).*ruy, -gam1*ones(ng,1)-eps_udg(:,4).*rvy, -gam*vv-eps_udg(:,4).*rHy-eps.*rHy_rE];

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
f_udg = cat(4,f_u,reshape(f_q,ng,nch,2,2*nch));


