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

nc = size(udg,2);

if nc==12
    vis  = param{7};
    [f,f_udg] = flux2d(pg,udg,param,time);
    [fav,fav_udg] = avflux2d(pg,udg,param,time);
    f = f + vis*bsxfun(@times,pg(:,4),fav);
    f_udg = f_udg + vis*bsxfun(@times,pg(:,4),fav_udg);
end
