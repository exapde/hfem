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

[ng,nch] = size(uh);

gam = param{1};

u = udg(:,1:nch);

switch ib
    case 1  % Far field
        [an,anm] = getan(nl,uh,param,2);
        [An,Anm] = getan(nl,uh,param,3);
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_udg = an+An;
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-2*An;
    case 2  % Wall
        un = u(:,2).*nl(:,1) + u(:,3).*nl(:,2);        
        uinf = u;
        uinf(:,2) = uinf(:,2) - nl(:,1).*un;
        uinf(:,3) = uinf(:,3) - nl(:,2).*un;
        
        fh = uinf - uh;
                
        fh_udg = zeros(ng,nch,nch);
        fh_udg(:,2,2:3) = [ones(ng,1)-nl(:,1).*nl(:,1), -nl(:,1).*nl(:,2)]; 
        fh_udg(:,3,2:3) = [ -nl(:,1).*nl(:,2), ones(ng,1)-nl(:,2).*nl(:,2)];
        fh_udg(:,4,4) = ones(ng,1); 
        fh_udg(:,1,1) = ones(ng,1); 
       
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,4,4) = -1; 

    case 3 % Prescribed pressure
        p  = uinf(:,1);
        uinf = [ u(:,1), u(:,2), u(:,3), p./(gam-1) + 0.5*u(:,2).*u(:,2)./u(:,1) + 0.5*u(:,3).*u(:,3)./u(:,1)];
        uinfu = bsxfun(@times,ones(ng,1,1),reshape(eye(nch),[1 nch nch]));
        uinfu(:,4,:) = [-0.5*(u(:,2).*u(:,2)+u(:,3).*u(:,3))./(u(:,1).*u(:,1)), u(:,2)./u(:,1), u(:,3)./u(:,1), zeros(ng,1)];
        
        [an,anm] = getan(nl,uh,param,2);
        [An,Anm] = getan(nl,uh,param,3);
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_udg = an+An-permute(mapContractK(an-An,uinfu,2,3,1,2,3,1),[3 1 2]);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2]) - 2*An;
        
    case 4 % Prescribed pressure
        p  = (uinf(:,4)-0.5)*(gam-1);
        uinf = [ u(:,1), u(:,2), u(:,3), p./(gam-1) + 0.5*u(:,2).*u(:,2)./u(:,1) + 0.5*u(:,3).*u(:,3)./u(:,1)];
        uinfu = bsxfun(@times,ones(ng,1,1),reshape(eye(nch),[1 nch nch]));
        uinfu(:,4,:) = [-0.5*(u(:,2).*u(:,2)+u(:,3).*u(:,3))./(u(:,1).*u(:,1)), u(:,2)./u(:,1), u(:,3)./u(:,1), zeros(ng,1)];
        
        [an,anm] = getan(nl,uh,param,2);
        [An,Anm] = getan(nl,uh,param,3);
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_udg = an+An-permute(mapContractK(an-An,uinfu,2,3,1,2,3,1),[3 1 2]);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2]) - 2*An;        
    case 5                    
        uinf(:,1) = u(:,1);
        uinf(:,4) = u(:,4);
        
        fh = uinf - uh;
                
        fh_udg = zeros(ng,nch,nch);
        fh_udg(:,4,4) = ones(ng,1); 
        fh_udg(:,1,1) = ones(ng,1); 
       
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,4,4) = -1;         
    otherwise
        error('unknown boundary type');
end



