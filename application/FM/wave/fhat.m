function [fh,fh_udg,fh_uh] = fhat(nl,pg,udg,uh,param,time,signe)
%FHAT flux function
%   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q
% 

[ng,nc] = size(udg);
nd = size(pg,2);
nch = nc/(nd+1);

c2 = param{1};

u  = udg(:,nch);
q  = reshape(udg(:,nch+1:nc),[ng nch nd]);

% tau  = c2;
% tau = 10;
tau = param{end};
fh = c2*reshape(mapContractK(nl,q,[],2,1,3,2,1),[ng nch]) + tau*(u-uh);

fh_uh = -tau*ones(ng,nch);

fh_u = tau*ones(ng,nch);            
fh_q = c2*reshape(nl,[ng,nch,nd]);
fh_udg = cat(3,fh_u,fh_q);