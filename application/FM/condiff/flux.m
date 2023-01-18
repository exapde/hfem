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
c     = param{2};

u     = udg(:,nch);
q     = reshape(udg(:,nch+1:nc),[ng nch nd]);

f = bsxfun(@times,u,reshape(c,[1,1,nd]))+ kappa*reshape(q,[ng nch nd]);

f_u = bsxfun(@times,ones(ng,1,1),reshape(c,[1,nch,nd]));
f_q = bsxfun(@times,ones(ng,1,1,1),kappa*reshape(eye(nd),[1,nch,nd,nd]));

f_udg = cat(4,f_u,f_q);


