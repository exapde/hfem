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

[ng,nc] = size(udg);


if nc==6
    %[fh,fh_udg,fh_uh] = fhat2d(nl,pg,udg,uh,param,time);
    nch = 2;
else
    [fh,fh_udg,fh_uh] = fhat3d(nl,pg,udg,uh,param,time);
    nch = 3;
end

% u = udg(:,1:nch);
% q = reshape(udg(:,nch+1:end),ng,nch,[]);
% 
% [f,f_uhdg] = flux(pg,[uh,udg(:,nch+1:end)],param,time);
% 
% %[An,Anm] = getan(nl,uh,param,2);
% An  = zeros(ng,nch,nch);
% for i=1:nch
%     An(:,i,i)=param{end};
% end
% Anm = zeros(ng,nch,nch,nch);
% 
% fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
% %fh = permute(mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
% 
% fn_udgh = mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1);
% 
% fh_u = An;
% fh_q = permute(fn_udgh(:,nch+1:nc,:),[3 1 2]);
% fh_udg = cat(3,fh_u,fh_q);
% 
% fh_uh = permute(fn_udgh(:,1:nch,:)+mapContractK(Anm,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;    
% 
