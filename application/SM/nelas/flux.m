function [f,f_udg] = flux(pg,udg,param,time)
%FLUX Volume flux function
%   [f,fu,fq,fp] = flux(x,u,q,p,param)
%
%      X(N,ND)              Coordinates for N points
%      U(N,NC)              Unknown vector for N points with NC components
%      Q(N,NC,ND)           Flux vector for N points with NC components in the
%                           coordinate directions
%      P(N,1)               Pressure vector for N points with 1 component
%      PARAM                Parameter list
%      F(N,NC,ND):          Volume flux at N points
%      FU(N,NC,ND,NC):      Jacobian of the flux flux vector w.r.t. U
%      FQ(N,NC,ND,NC,ND):   Jacobian of the flux flux vector w.r.t. Q
%      FP(N,NC,ND):         Jacobian of the flux flux vector w.r.t. P

% 2D plane strain neo-Hookean materials

nc = size(udg,2);

if nc==7
    [f,f_udg] = flux2d(pg,udg,param,time);
else
    [f,f_udg] = flux3d(pg,udg,param,time);
end

return;


[ng,nc] = size(udg);
nch = 2;
 
if size(pg,2)<=2
    mu = param{1};
else
    mu = pg(:,3);
end

u1  = udg(:,1);
u2  = udg(:,2);
pres= udg(:,3);
F11 = udg(:,4);
F21 = udg(:,5);
F12 = udg(:,6);
F22 = udg(:,7);

wcase=1;
[P,dPdF,dPdp,g,dgdF] = Deriv(F11,F12,F21,F22,pres,mu,wcase);
dgdF(:,[2 3]) = dgdF(:,[3 2]);

f = zeros(ng,nch,2);
% f(:,:,1) = [P(:,1) P(:,2)];
% f(:,:,2) = [P(:,3) P(:,4)];
f(:,:,1) = [P(:,1) P(:,3)];
f(:,:,2) = [P(:,2) P(:,4)];

f_p = zeros(ng,nch,2);
% f_p(:,:,1) = [dPdp(:,1) dPdp(:,2)];
% f_p(:,:,2) = [dPdp(:,3) dPdp(:,4)];
f_p(:,:,1) = [dPdp(:,1) dPdp(:,3)];
f_p(:,:,2) = [dPdp(:,2) dPdp(:,4)];

f_u = zeros(ng,nch,2,nch);

f_q = zeros(ng,nch,2,2*nch);
% f_q(:,:,1,1) = [dPdF(:,1) dPdF(:,2)];
% f_q(:,:,1,2) = [dPdF(:,3) dPdF(:,4)];
% f_q(:,:,1,3) = [dPdF(:,5) dPdF(:,6)];
% f_q(:,:,1,4) = [dPdF(:,7) dPdF(:,8)];
% f_q(:,:,2,1) = [dPdF(:,9) dPdF(:,10)];
% f_q(:,:,2,2) = [dPdF(:,11) dPdF(:,12)];
% f_q(:,:,2,3) = [dPdF(:,13) dPdF(:,14)];
% f_q(:,:,2,4) = [dPdF(:,15) dPdF(:,16)];
f_q(:,:,1,1) = [dPdF(:,1) dPdF(:,9)];
f_q(:,:,1,2) = [dPdF(:,5) dPdF(:,13)];
f_q(:,:,1,3) = [dPdF(:,3) dPdF(:,11)];
f_q(:,:,1,4) = [dPdF(:,7) dPdF(:,15)];
f_q(:,:,2,1) = [dPdF(:,2) dPdF(:,10)];
f_q(:,:,2,2) = [dPdF(:,6) dPdF(:,14)];
f_q(:,:,2,3) = [dPdF(:,4) dPdF(:,12)];
f_q(:,:,2,4) = [dPdF(:,8) dPdF(:,16)];

 
f_udg = cat(4,cat(4,f_u,f_p),f_q);

%g = detF-1;
g_udg = [zeros(ng,3) dgdF];

% u1  = udg(:,1);
% u2  = udg(:,2);
% % F11 = udg(:,3);
% % F12 = udg(:,5);
% % F21 = udg(:,4);
% % F22 = udg(:,6);
% ux = udg(:,3);
% vx = udg(:,4);
% uy = udg(:,5);
% vy = udg(:,6);
% 
% 
% wcase=1;
% [P,dPdF] = Deriv(ux,uy,vx,vy,mu,kappa,wcase);
% 
% f = zeros(ng,nch,2);
% f(:,:,1) = [P(:,1) P(:,3)];
% f(:,:,2) = [P(:,2) P(:,4)];
% 
% f_u = zeros(ng,nch,2,nch);
% 
% f_q = zeros(ng,nch,2,2*nch);
% f_q(:,:,1,1) = [dPdF(:,1) dPdF(:,9)];
% f_q(:,:,1,2) = [dPdF(:,5) dPdF(:,13)];
% f_q(:,:,1,3) = [dPdF(:,3) dPdF(:,11)];
% f_q(:,:,1,4) = [dPdF(:,7) dPdF(:,15)];
% f_q(:,:,2,1) = [dPdF(:,2) dPdF(:,10)];
% f_q(:,:,2,2) = [dPdF(:,6) dPdF(:,14)];
% f_q(:,:,2,3) = [dPdF(:,4) dPdF(:,12)];
% f_q(:,:,2,4) = [dPdF(:,8) dPdF(:,16)];
%  
% f_udg = cat(4,f_u,f_q);
% 
% 
% 
