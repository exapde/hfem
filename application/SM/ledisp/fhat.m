function [fh,fh_udg,fh_uh] = fhat(nl,pg,udg,uh,param,time)
%FHAT flux function
%   [fh,fhu,fhq,fhp,fhm] = fhat(nl,x,u,q,p,m,param)
%
%      NL(N,ND)              Normal N points
%      X(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      P(N,1)                Pressure vector for N points with 1 component
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):             Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHP(N,NC,1):          Jacobian of the flux flux vector w.r.t. P
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. M

[~,nc] = size(udg);

if nc==6
    [fh,fh_udg,fh_uh] = fhat2d(nl,pg,udg,uh,param,time);    
else
    [fh,fh_udg,fh_uh] = fhat3d(nl,pg,udg,uh,param,time);    
end

