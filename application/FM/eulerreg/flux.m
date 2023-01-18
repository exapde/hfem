function [f,f_udg] = flux(p,udg,param,time)
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

% [ng,nch] = size(udg);
% 
% gam  = param{1};
% gam1 = gam - 1.0;
%                                              
% r    = udg(:,1);
% ru   = udg(:,2);
% rv   = udg(:,3);
% rE   = udg(:,4);
% 
% r1   = 1./r;
% uv   = ru.*r1;
% vv   = rv.*r1;
% E    = rE.*r1;
% af   = 0.5*(uv.*uv+vv.*vv);
% p    = gam1*(rE-r.*af);
% h    = E+p.*r1;
%                                         
% f = zeros(ng,nch,2);
% f(:,:,1) = [ru, ru.*uv+p, rv.*uv,   ru.*h];
% f(:,:,2) = [rv, ru.*vv,   rv.*vv+p, rv.*h];
% 
% f_udg = zeros(ng,nch,2,nch);
% f_udg(:,:,1,1) = -[zeros(ng,1), 0.5*((3-gam)*uv.*uv-gam1*vv.*vv), uv.*vv, gam*E.*uv-2*gam1*uv.*af];
% f_udg(:,:,1,2) = -[-ones(ng,1), (gam-3)*uv, -vv, -gam*E+0.5*gam1*(3*uv.*uv+vv.*vv)];
% f_udg(:,:,1,3) = -[zeros(ng,1), gam1*vv, -uv, gam1*uv.*vv];
% f_udg(:,:,1,4) = -[zeros(ng,1), -gam1*ones(ng,1), zeros(ng,1), -gam*uv];
% f_udg(:,:,2,1) = -[zeros(ng,1), uv.*vv, 0.5*((3-gam)*vv.*vv-gam1*uv.*uv), gam*E.*vv-2*gam1*vv.*af];
% f_udg(:,:,2,2) = -[zeros(ng,1), -vv, gam1*uv, gam1*uv.*vv];
% f_udg(:,:,2,3) = -[-ones(ng,1), -uv, (gam-3)*vv,  -gam*E+0.5*gam1*(3*vv.*vv+uv.*uv) ];
% f_udg(:,:,2,4) = -[zeros(ng,1), zeros(ng,1), -gam1*ones(ng,1), -gam*vv];

[ng,nc] = size(udg);
nch = 4;
nd = 2;
one = ones(ng,1);
param1 = param{1};
param4 = param{4};
u1 = udg(:,1);
u2 = udg(:,2);
u3 = udg(:,3);
u4 = udg(:,4);
zero = zeros(ng,1);
t2 = 1.0./param4;
t3 = t2.*u1.*2.0e1;
t4 = exp(t3);
t5 = t4+1.0;
t6 = log(t5);
t7 = 1.0./param4.^2;
t8 = 1.0./t6.^2;
t9 = u2.^2;
t10 = 1.0./t6;
t11 = t7.*t8.*t9.*2.0e2;
t12 = u3.^2;
t13 = t7.*t8.*t12.*2.0e2;
t14 = t11+t13;
t18 = param4.*t6.*t14.*(1.0./2.0e1);
t15 = -t18+u4;
t16 = param1-1.0;
t17 = t2.*t10.*u2.*u3.*2.0e1;
t19 = t15.*t16;
t20 = t2.*t10.*u4.*2.0e1;
t21 = t2.*t10.*t15.*t16.*2.0e1;
t22 = t20+t21;
f = [u2;t19+t2.*t9.*t10.*2.0e1;t17;t22.*u2;u3;t17;t19+t2.*t10.*t12.*2.0e1;t22.*u3];
if nargout > 1
    t23 = 1.0./t5;
    t24 = 1.0./param4.^3;
    t25 = 1.0./t6.^3;
    t26 = t4.*t14.*t23;
    t27 = t4.*t9.*t23.*t24.*t25.*8.0e3;
    t28 = t4.*t12.*t23.*t24.*t25.*8.0e3;
    t29 = t27+t28;
    t31 = param4.*t6.*t29.*(1.0./2.0e1);
    t30 = t26-t31;
    t32 = t2.*t10.*t16.*t30.*2.0e1;
    t33 = t4.*t7.*t8.*t23.*u4.*4.0e2;
    t34 = t4.*t7.*t8.*t15.*t16.*t23.*4.0e2;
    t35 = t32+t33+t34;
    t36 = t2.*t10.*u3.*2.0e1;
    t37 = t2.*t10.*u2.*2.0e1;
    t38 = one.*t16;
    t39 = t2.*t10.*2.0e1;
    t40 = t2.*t10.*t16.*2.0e1;
    t41 = t39+t40;
    f_udg = [zero;-t16.*t30-t4.*t7.*t8.*t9.*t23.*4.0e2;t4.*t7.*t8.*t23.*u2.*u3.*-4.0e2;-t35.*u2;zero;t4.*t7.*t8.*t23.*u2.*u3.*-4.0e2;-t16.*t30-t4.*t7.*t8.*t12.*t23.*4.0e2;-t35.*u3;one;t2.*t10.*u2.*4.0e1-t2.*t10.*t16.*u2.*2.0e1;t36;t20+t21-t7.*t8.*t9.*t16.*4.0e2;zero;t36;t2.*t10.*t16.*u2.*-2.0e1;t7.*t8.*t16.*u2.*u3.*-4.0e2;zero;t2.*t10.*t16.*u3.*-2.0e1;t37;t7.*t8.*t16.*u2.*u3.*-4.0e2;one;t37;t2.*t10.*u3.*4.0e1-t2.*t10.*t16.*u3.*2.0e1;t20+t21-t7.*t8.*t12.*t16.*4.0e2;zero;t38;zero;t41.*u2;zero;zero;t38;t41.*u3];
end
f = reshape(f,ng,nch,nd);
f_udg = reshape(f_udg,ng,nch,nd,nc);
