function [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time)
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
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Lambda

[ng,nch] = size(uh);
nc = size(udg,2);

[f,f_uhdg] = flux(p,[uh,udg(:,nch+1:end)],param,time);

tau = param{6};    
fh = permute(mapContractK(f,nl,2,3,1,2,[],1),[2 1]) + tau*(udg(:,1:nch)-uh);

fh_udg = zeros(ng,nch,nc);
temp = permute(mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1),[3 1 2]);
fh_udg(:,:,nch+1:end) = temp(:,:,nch+1:end);    
fh_uh = temp(:,:,1:nch);
for i=1:nch
    fh_udg(:,i,i) = tau;
    fh_uh(:,i,i) = fh_uh(:,i,i) - tau;
end    

% [ng,nch] = size(uh);
% 
% u = udg(:,1:nch);
% q = reshape(udg(:,nch+1:3*nch),ng,nch,2);
% 
% [f,f_uhdg] = flux(p,[uh,udg(:,nch+1:3*nch)],param,time);
% 
% [An,Anm] = getan(nl,uh,param,2);
% 
% fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
% %fh = permute(mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
% 
% fn_udgh = mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1);
% 
% fh_u = An;
% fh_q = permute(fn_udgh(:,nch+1:3*nch,:),[3 1 2]);
% fh_udg = cat(3,fh_u,fh_q);
% 
% fh_uh = permute(fn_udgh(:,1:nch,:)+mapContractK(Anm,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;    