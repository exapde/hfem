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
nch = 2;
nd = (size(udg,2)-nch)/nch;

K = param{1};
D  = param{2};
Wind = param{9};

f = zeros(ng,nch,nd);
f_udg = zeros(ng,nch,nd,nc);

%phi = udg(:,1);
rho = udg(:,2);
phix = udg(:,3);
rhox = udg(:,4);
phiy = udg(:,5);
rhoy = udg(:,6);

if size(p,2)>nd
    vx = Wind*p(:,nd+1);
    vy = Wind*p(:,nd+2);
else
    vx = 0;
    vy = 0;
end

f(:,1,1) = phix;
f(:,1,2) = phiy;
f(:,2,1) = K*rho.*(phix+vx) + D*rhox;
f(:,2,2) = K*rho.*(phiy+vy) + D*rhoy;

f_udg(:,1,1,1) = 0;
f_udg(:,1,1,2) = 0;
f_udg(:,1,1,3) = 1;
f_udg(:,1,1,4) = 0;
f_udg(:,1,1,5) = 0;
f_udg(:,1,1,6) = 0;

f_udg(:,1,2,1) = 0;
f_udg(:,1,2,2) = 0;
f_udg(:,1,2,3) = 0;
f_udg(:,1,2,4) = 0;
f_udg(:,1,2,5) = 1;
f_udg(:,1,2,6) = 0;

f_udg(:,2,1,1) = 0;
f_udg(:,2,1,2) = K*(phix+vx);
f_udg(:,2,1,3) = K*rho;
f_udg(:,2,1,4) = D;
f_udg(:,2,1,5) = 0;
f_udg(:,2,1,6) = 0;

f_udg(:,2,2,1) = 0;
f_udg(:,2,2,2) = K*(phiy+vy);
f_udg(:,2,2,3) = 0;
f_udg(:,2,2,4) = 0;
f_udg(:,2,2,5) = K*rho;
f_udg(:,2,2,6) = D;

