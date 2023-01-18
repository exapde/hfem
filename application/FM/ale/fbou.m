function [fh,fh_udg,fh_uh] = fbou(ib,uinf,nl,p,udg,uh,param,time)
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
%      FH(N,NC):             Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Lambda

[ng1,nch] = size(uh);

gam = param{1};
gam1 = gam-1.0;
Minf = param{5};
M2   = Minf^2;

one = ones(ng1,1);
zero = zeros(ng1,1);

u = udg(:,1:nch);
%q = reshape(udg(:,nch+1:3*nch),ng1,nch,2);

p = mapping(p,time);
switch ib
    case 1  % Far field
        %uinf = repmat(ui,ng1,1).*repmat(p(:,12),[1 nch]);
        [an,anm] = getan(p,nl,uh,param,0,time);
        [An,Anm] = getan(p,nl,uh,param,1,time);
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_u = an+An;
        fh_q = zeros(ng1,nch,nch,2);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-2*An;
    case 2  % Adiabatic Wall
        [fh,fh_udg,fh_uh] = fhatwall(nl,p,udg,uh,param,time,[]);
        fh_u = fh_udg(:,:,1:nch);
        fh_q = fh_udg(:,:,nch+1:3*nch);
                
        fh(:,2:3) = uh(:,2:3) - p(:,[5 6]).*repmat(p(:,12),[1 2]);
        fh_u(:,2:3,:) = 0;
        fh_q(:,2:3,:) = 0;
        fh_uh(:,2:3,[1,4]) = 0;
        %fh_uh(:,2:3,2:3) = multiprod(one,eye(2,2),2,[1,2]);
        fh_uh(:,2,2) = 1;
        fh_uh(:,3,3) = 1;
    case 3 % Isothermal wall
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time,[]);
        fh_u = fh_udg(:,:,1:nch);
        fh_q = fh_udg(:,:,nch+1:3*nch);
        
        fh(:,2:3) = uh(:,2:3) - p(:,[5 6]).*repmat(p(:,12),[1 2]);
        fh(:,4) = gam*gam1*M2*uh(:,4)./uh(:,1) - ui(4);
        fh_u(:,2:4,:) = 0;
        fh_q(:,2:4,:) = 0;
        fh_uh(:,2:3,[1,4]) = 0;
        fh_uh(:,2,2) = 1;
        fh_uh(:,3,3) = 1;
        fh_uh(:,4,:) = gam*gam1*M2*[-uh(:,4)./uh(:,1).^2, zero, zero, 1./uh(:,1)];
    otherwise
        error('unknown boundary type');
end

fh_udg = cat(3,fh_u,reshape(fh_q,ng1,nch,2*nch));



