function [fx,fy] = eulerv(u,p,param,time)
%EULERV Calculate the volume flux for the Euler equations.
%   [FX,FY]=EULERV(U,P,PARAM,TIME)
%
%      U(np,4):      np left (or plus) states
%      P:            Not used
%      PARAM{1}:     Cell array containing the value of gamma
%      TIME:         Not used
%      FX(np,4):     np fluxes in the x direction (f plus)  
%      FY(np,4):     np fluxes in the y direction (f plus)  
%                          
% - Written by: J. Peraire
%
gam = param{1};

uv = u(:,2)./u(:,1);
vv = u(:,3)./u(:,1);

p = (gam-1)*(u(:,4) - 0.5*(u(:,2).*uv + u(:,3).*vv));

fx = [ u(:,2), u(:,2).*uv+p,   u(:,3).*uv, uv.*(u(:,4)+p)];
fy = [ u(:,3),   u(:,2).*vv, u(:,3).*vv+p, vv.*(u(:,4)+p)];


