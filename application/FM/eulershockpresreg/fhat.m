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
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Lambda

[ng1,nch] = size(uh);
nc = size(udg,2);

u = udg(:,1:nch);
% if min(u(:,1))<=1e-3
%     error('negative density');
% end

if nch == 4    
    [f,f_uhdg] = flux(p,[uh,udg(:,nch+1:end)],param,time);

%     [An,Anm] = getan(nl,uh,param,2);
% 
%     %fh = multiprod(f(:,:,1),nl(:,1),2)+multiprod(f(:,:,2),nl(:,2),2) + multiprod(An,(u-uh),[2,3],2);
%     fh = zeros(ng1,nch);
%     for i = 1:nch
%        fh(:,i) = f(:,i,1).*nl(:,1) + f(:,i,2).*nl(:,2) + sum(squeeze(An(:,i,:)).*(u - uh),2);    
%     end
%     fh_udg = An;
% 
%     fh_uh = multiprod(permute(f_uhdg(:,:,1,1:nch),[1,2,4,3]),nl(:,1),[2,3],2) +   ...
%             multiprod(permute(f_uhdg(:,:,2,1:nch),[1,2,4,3]),nl(:,2),[2,3],2) + multiprod(Anm,(u-uh),[2,3],2) - An;
    tau = param{6};    
    fh = permute(mapContractK(f,nl,2,3,1,2,[],1),[2 1]) + tau*(u-uh);
    fh_udg = zeros(ng1,nch,nc);
    temp = permute(mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1),[3 1 2]);
    fh_udg(:,:,nch+1:end) = temp(:,:,nch+1:end);    
    fh_uh = temp(:,:,1:nch);
    for i=1:nch
        fh_udg(:,i,i) = tau;
        fh_uh(:,i,i) = fh_uh(:,i,i) - tau;
    end    
    
else        
    umuh = (u-uh)';    

    [f,f_uhdg] = flux(p,[uh,udg(:,nch+1:end)],param,time);

    [An,Anm] = getan(nl,uh,param,2);

    n1d = spdiags(nl(:,1),0,size(nl,1),size(nl,1));
    n2d = spdiags(nl(:,2),0,size(nl,1),size(nl,1));

    fh = n1d*f(:,:,1) + n2d*f(:,:,2);
    temp = permute(An,[2,3,1]);
    for ig = 1:ng1
        fh(ig,:) = fh(ig,:) + (temp(:,:,ig)*(umuh(:,ig)))';
    end

    temp = reshape(n1d*reshape(f_uhdg(:,:,1,:),ng1,[]) + n2d*reshape(f_uhdg(:,:,2,:),ng1,[]),ng1,nch,[]);

    fh_u = An;
    fh_q = temp(:,:,nch+1:end);    
    fh_uh = -An + temp(:,:,1:nch);

    temp = reshape(permute(Anm,[2,4,3,1]),nch*nch,nch,ng1);
    for ig = 1:ng1
        fh_uh(ig,:,:) = fh_uh(ig,:,:) + reshape(temp(:,:,ig)*umuh(:,ig),1,nch,nch);
    end

    fh_udg = zeros(ng1,nch,nc);
    fh_udg(:,:,1:nch) = fh_u;
    fh_udg(:,:,nch+1:end) = fh_q;        
end
