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

[ng,nc] = size(udg);
nch = 1;
nq = nc-nch;
nd = nq;

kappa = param{1};

% x = p(:,1);
% y = p(:,2);
% kappa = sin(x).*sin(y)

if nd==2
    u = udg(:,1);
    qx = udg(:,2);
    qy = udg(:,3);

    f = zeros(ng,nch,nd);
    f(:,1,1) = kappa*qx;
    f(:,1,2) = kappa*qy;

    f_udg = zeros(ng,nch,nd,nc);
    f_udg(:,1,1,2) = kappa;
    f_udg(:,1,2,3) = kappa;
else
    u = udg(:,1);
    qx = udg(:,2);
    qy = udg(:,3);
    qz = udg(:,4);

    f = zeros(ng,nch,nd);
    f(:,1,1) = kappa*qx;
    f(:,1,2) = kappa*qy;
    f(:,1,3) = kappa*qz;

    f_udg = zeros(ng,nch,nd,nc);
    f_udg(:,1,1,2) = kappa;
    f_udg(:,1,2,3) = kappa;    
    f_udg(:,1,3,4) = kappa;    
end
