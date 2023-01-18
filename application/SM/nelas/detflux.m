function [J,J_udg] = detflux(p,udg,param,time)
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

[ng,nc] = size(udg);
nch = 2;

if size(p,2)<=2
    mu = param{1};
else
    mu = p(:,3);
end

if size(p,2)<=3
    lambda = 1;
else
    lambda = p(:,4);
end


u1  = udg(:,1);
u2  = udg(:,2);
pres= udg(:,3);
F11 = udg(:,4);
F21 = udg(:,5);
F12 = udg(:,6);
F22 = udg(:,7);

J = lambda.*((F11.*F22 - F12.*F21) - 1);
J_udg = [0*F11, 0*F11, 0*F11, lambda.*F22, -lambda.*F12, -lambda.*F21,  lambda.*F11];  


