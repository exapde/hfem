function [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time,signe)
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

[ng,nc] = size(udg);
nch = 1;
nq = nc-nch;
nd = nq;
onesng = ones(ng,1);

kappa = param{1};
c     = param{2};
tau   = param{3};
fhform = param{4};
Aform = param{5};

u     = udg(:,nch);
q     = reshape(udg(:,nch+1:nc),[ng nch nd]);
An    = nl(:,1)*c(:,1)+nl(:,2).*c(:,2);

switch (fhform)
    case 1
        fh1 = An.*uh;
        fh1_uh = An;
        fh1_u = 0*onesng;
    case 2
        fh1 = An.*u;
        fh1_uh = 0*onesng;
        fh1_u = An;        
    case 3
        fh1 = 0.5*An.*(uh+u);
        fh1_uh = 0.5*An;
        fh1_u = 0.5*An;
    otherwise
        error('Not yet implemented');
end

epsilon = 1e-10;
switch (Aform)
    case 1
        fh2 = tau.*(u-uh);
        fh2_uh = -tau*onesng;
        fh2_u = tau*onesng;
    case 2
        fh2 = (0.5*abs(An)+kappa+epsilon).*(u-uh);
        fh2_uh = -(0.5*abs(An)+kappa+epsilon);
        fh2_u = (0.5*abs(An)+kappa+epsilon);
    case 3
        fh2 = (0.5*(abs(An)+An)+kappa+epsilon).*(u-uh);
        fh2_uh = -(0.5*(abs(An)+An)+kappa+epsilon);
        fh2_u = (0.5*(abs(An)+An)+kappa+epsilon);
    case 4
        fh2 = (0.5*(abs(An)-An)+kappa+epsilon).*(u-uh);    
        fh2_uh = -(0.5*(abs(An)-An)+kappa+epsilon);
        fh2_u = (0.5*(abs(An)-An)+kappa+epsilon);
    otherwise
        error('Not yet implemented');
end

fh = fh1 + fh2 +  kappa*reshape(mapContractK(nl,q,[],2,1,3,2,1),[ng nch]);

fh_uh = fh1_uh + fh2_uh;

fh_u = fh1_u + fh2_u;
fh_q = kappa*reshape(nl,[ng,nch,nd]);
fh_udg = cat(3,fh_u,fh_q);


% %tau = 0.5*abs(An);
% % tau  = tauc + kappa;
% %tau = param{end}*ones(ng,1);
% tau = 0.5*(abs(An))+kappa;
% %tau = 1/2*ones(ng,1);
% fh = An.*uh + kappa*reshape(mapContractK(nl,q,[],2,1,3,2,1),[ng nch]) + tau.*(u-uh);
% 
% fh_uh = An-tau;
% 
% fh_u = tau;            
% fh_q = kappa*reshape(nl,[ng,nch,nd]);
% fh_udg = cat(3,fh_u,fh_q);
