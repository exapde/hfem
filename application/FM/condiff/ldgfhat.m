function fh = ldgfhat(nl,p,udg1,udg2,uh,param,time)
% fhg = fhat(nlg1, pg1, udg1, udg2, uhg, app.arg, app.time);
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

[ng,nc] = size(udg1);
nch = 1;
nq = nc-nch;
nd = nq;

kappa = param{1};
c     = param{2};
tau   = param{end};

u1 = udg1(:,1);
q1x = udg1(:,2);
q1y = udg1(:,3);
u2 = udg2(:,1);
q2x = udg2(:,2);
q2y = udg2(:,3);
u  = 0.5*(u1+u2);
qx = 0.5*(q1x+q2x);
qy = 0.5*(q1y+q2y);

% ng x nch
fh = kappa.*(qx.*nl(:,1)+qy.*nl(:,2)) + (nl*c(:)).*u + tau.*(u1(:,1)-u2(:,1));



