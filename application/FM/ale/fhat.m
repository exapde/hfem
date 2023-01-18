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

u = udg(:,1:nch);
q = reshape(udg(:,nch+1:3*nch),ng1,nch,2);
umuh = (u-uh)';

[f,f_uhdg] = flux(p,[uh,udg(:,nch+1:3*nch)],param,time);

[An,Anm] = getan(p,nl,uh,param,2,time);

n1d = spdiags(nl(:,1),0,size(nl,1),size(nl,1));
n2d = spdiags(nl(:,2),0,size(nl,1),size(nl,1));

%fh = multiprod(f(:,:,1),nl(:,1),2)+multiprod(f(:,:,2),nl(:,2),2) + multiprod(An,(u-uh),[2,3],2);
%fh = zeros(ng1,nch);

fh = n1d*f(:,:,1) + n2d*f(:,:,2);
temp = permute(An,[2,3,1]);
for ig = 1:ng1
    fh(ig,:) = fh(ig,:) + (temp(:,:,ig)*umuh(:,ig))';
end

% fh_uh = multiprod(permute(f_uhdg(:,:,1,1:nch),[1,2,4,3]),nl(:,1),[2,3],2) +   ...
%         multiprod(permute(f_uhdg(:,:,2,1:nch),[1,2,4,3]),nl(:,2),[2,3],2) + multiprod(Anm,(u-uh),[2,3],2) - An;

%n1d=full(n1d);n2d=full(n2d);
temp = reshape(n1d*reshape(f_uhdg(:,:,1,:),ng1,[]) + n2d*reshape(f_uhdg(:,:,2,:),ng1,[]),ng1,nch,[]);

fh_u = An;
fh_q = temp(:,:,nch+1:3*nch);
fh_uh = -An + temp(:,:,1:nch);

temp = reshape(permute(Anm,[2,4,3,1]),nch*nch,nch,ng1);
for ig = 1:ng1
    fh_uh(ig,:,:) = fh_uh(ig,:,:) + reshape(temp(:,:,ig)*umuh(:,ig),1,nch,nch);
end
    
fh_udg = cat(3,fh_u,fh_q);


