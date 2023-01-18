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
nc = size(udg,2);

% gam = param{1};
% one = ones(ng,1);
% zero = zeros(ng,1);

u = udg(:,1:nch);

switch ib
    case 1  % Far field: Roe
        [an,anm] = getan(nl,uh,param,0);
        [An,Anm] = getan(nl,uh,param,1);
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_u = an+An;
        fh_q = zeros(ng,nch,nch,2);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-2*An;
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
    case 2  % Slip Adiabatic Wall
        un = u(:,2).*nl(:,1) + u(:,3).*nl(:,2);        
        uinf = u;
        uinf(:,2) = uinf(:,2) - nl(:,1).*un;
        uinf(:,3) = uinf(:,3) - nl(:,2).*un;
        
        fh = uinf - uh;
                
        fh_udg = zeros(ng,nch,nc);
        fh_udg(:,2,2:3) = [ones(ng,1)-nl(:,1).*nl(:,1), -nl(:,1).*nl(:,2)]; 
        fh_udg(:,3,2:3) = [ -nl(:,1).*nl(:,2), ones(ng,1)-nl(:,2).*nl(:,2)];
        fh_udg(:,4,4) = ones(ng,1); 
        fh_udg(:,1,1) = ones(ng,1); 
       
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,4,4) = -1; 
        
%         un = u(:,2).*nl(:,1) + u(:,3).*nl(:,2);        
%         
%         uinf = u;
%         uinfu = zeros(ng,nch,nch);
%         uinfu(:,1,:) = [ones(ng,1) zeros(ng,1) zeros(ng,1) zeros(ng,1)]; 
%         
%         uinf(:,2) = uinf(:,2) - nl(:,1).*un;
%         uinfu(:,2,:) = [zeros(ng,1) ones(ng,1)-nl(:,1).*nl(:,1) -nl(:,1).*nl(:,2) zeros(ng,1)]; 
%         uinf(:,3) = uinf(:,3) - nl(:,2).*un;
%         uinfu(:,3,:) = [zeros(ng,1) -nl(:,1).*nl(:,2) ones(ng,1)-nl(:,2).*nl(:,2) zeros(ng,1)];
%         
%         uinfu(:,4,:) = [zeros(ng,1) zeros(ng,1) zeros(ng,1) ones(ng,1)]; 
%                                                 
%         [An,Anm] = getan(nl,uh,param,2);   
% 
%         fh = permute(mapContractK(An,uinf-uh,2,3,1,2,[],1),[2 1]);
%         fh_u = permute(mapContractK(An,uinfu,2,3,1,2,3,1),[3 1 2]); 
%         fh_q = zeros(ng,nch,nch,2);
%         fh_uh = permute(mapContractK(Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-An;    
%         fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
    case 3 % supersonic inflow                                              
        fh = uinf - uh;                                
        
        fh_udg = zeros(ng,nch,nc);
        
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,4,4) = -1; 
    case 4 % supersonic outflow                                              
        fh = u - uh;                                
        
        fh_udg = zeros(ng,nch,nc);
        fh_udg(:,1,1) = 1; 
        fh_udg(:,2,2) = 1; 
        fh_udg(:,3,3) = 1; 
        fh_udg(:,4,4) = 1;        
        
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,4,4) = -1; 
    otherwise
        error('unknown boundary type');
end



