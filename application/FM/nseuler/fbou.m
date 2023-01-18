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
nd = nch-2;

gam = param{1};
gam1 = gam-1.0;
Minf = param{5};
M2   = Minf^2;


u = udg(:,1:nch);
q = reshape(udg(:,nch+1:3*nch),ng,nch,2);

switch ib
    case 1  % Far field
        [an,anm] = getan(nl,uh,param,0);
        [An,Anm] = getan(nl,uh,param,1);
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_u = an+An;
        fh_q = zeros(ng,nch,nch,2);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-2*An;
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
    case 2  % Adiabatic Wall     
        uinf = u;
        uinf(:,2:3) = 0;
        
        [fh4,fh_udg4,fh_uh4] = fhat(nl,p,udg,uh,param,time);
        
        fh = uinf - uh;
        fh(:,4) = fh4(:,4);
        
        fh_u = zeros(ng,nch,nch);
        fh_u(:,1,1) = ones(ng,1);        
        fh_u(:,4,:) = fh_udg4(:,4,1:nch);
   
        fh_q = zeros(ng,nch,nch,2);
        fh_q(:,4,:) = fh_udg4(:,4,nch+1:3*nch);
        
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;                
        fh_uh(:,4,:) = fh_uh4(:,4,:);
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));                
    case 3 % Isothermal wall
        %[fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        %fh_u = fh_udg(:,:,1:nch);
        %fh_q = fh_udg(:,:,nch+1:3*nch);
        
        %one = ones(ng,1);        
        zero = zeros(ng,1);
        fh(:,1)   = u(:,1)-uh(:,1);        
        fh(:,2:3) = -uh(:,2:3);
        fh(:,4) = gam*gam1*M2*uh(:,4)./uh(:,1) - uinf(:,end);
        
        fh_u = zeros(ng,nch,nch);
        fh_u(:,1,1) = ones(ng,1);  
        
        fh_q = zeros(ng,nch,nch,2);        
        
        fh_uh = zeros(ng,nch,nch);
        fh_uh(:,1,1) = -1;
        fh_uh(:,2,2) = -1;
        fh_uh(:,3,3) = -1;
        fh_uh(:,4,:) = gam*gam1*M2*[-uh(:,4)./uh(:,1).^2, zero, zero, 1./uh(:,1)];
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
    case 4 % Prescribed pressure
        %p  = uinf(:,1);        
        pinf = (uinf(:,4)-0.5)*(gam-1);
        
        uinf = [ u(:,1), u(:,2), u(:,3), pinf./(gam-1) + 0.5*u(:,2).*u(:,2)./u(:,1) + 0.5*u(:,3).*u(:,3)./u(:,1)];
        uinfu = bsxfun(@times,ones(ng,1,1),reshape(eye(nch),[1 nch nch]));
        uinfu(:,4,:) = [-0.5*(u(:,2).*u(:,2)+u(:,3).*u(:,3))./(u(:,1).*u(:,1)), u(:,2)./u(:,1), u(:,3)./u(:,1), zeros(ng,1)];                 
        
        [an,anm] = getan(nl,uh,param,0);
        [An,Anm] = getan(nl,uh,param,1);
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_u = an+An-permute(mapContractK(an-An,uinfu,2,3,1,2,3,1),[3 1 2]);
        fh_q = zeros(ng,nch,nch,2);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2]) - 2*An;
        
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));    
    case 5 % periodic
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time); 
    case 6 % Dirichlet bc   
        %max(abs(uinf(:)))
        fh = uh-uinf;
        fh_udg = zeros(ng,nch,3*nch);
        fh_uh = bsxfun(@times,ones(ng,1,1),reshape(eye(nch),[1 nch nch]));      
    case 7
        pinf = (uinf(:,4)-0.5)*(gam-1);        
        uinf = [ u(:,1), u(:,2), u(:,3), pinf];
        uinfu = bsxfun(@times,ones(ng,1,1),reshape(eye(nch),[1 nch nch]));
        uinfu(:,4,:) = 0; 

        uhinf = [uh(:,1), uh(:,2), uh(:,3), (gam-1)*(uh(:,4) - 0.5*uh(:,2).*uh(:,2)./uh(:,1) - 0.5*uh(:,3).*uh(:,3)./uh(:,1))];
        uhinfu = bsxfun(@times,ones(ng,1,1),reshape(eye(nch),[1 nch nch]));
        uhinfu(:,4,:) = (gam-1)*[0.5*(uh(:,2).*uh(:,2)+uh(:,3).*uh(:,3))./(uh(:,1).*uh(:,1)), -uh(:,2)./uh(:,1), -uh(:,3)./uh(:,1), ones(ng,1)];                 
        
        [An,Anm] = getan(nl,uh,param,1);
        fh = permute(mapContractK(An,uinf-uhinf,2,3,1,2,[],1),[2 1]);
        fh_u = permute(mapContractK(An,uinfu,2,3,1,2,3,1),[3 1 2]);
        fh_q = zeros(ng,nch,nch,2);
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
        fh_uh = permute(mapContractK(Anm,uinf-uhinf,[2 4],3,1,2,[],1),[3 1 2])-permute(mapContractK(An,uhinfu,2,3,1,2,3,1),[3 1 2]);                
    case 10 % Prescribing everything except pressure (taken from FM/rans)
        uinf(:,4) = u(:,4);
        uinfu = zeros(ng,nch,nch);      
        uinfu(:,4,4) = 1;
        
        [An,Anm] = getan(nl,uh,param,1);
        fh = permute(mapContractK(An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_u = permute(mapContractK(An,uinfu,2,3,1,2,3,1),[3 1 2]);
        fh_q = zeros(ng,nch,nch,2);
        fh_udg = cat(3,fh_u,reshape(fh_q,ng,nch,2*nch));
        fh_uh = permute(mapContractK(Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-An;                
    otherwise
        error('unknown boundary type');
end

if ib==2 || ib == 3
    alpha = param{6};
    [fheu,fheu_udg,fheu_uh] = fboueuler(nl,udg,uh);
    fh     = alpha*fh + (1-alpha)*fheu;
    fh_udg = alpha*fh_udg + (1-alpha)*fheu_udg;
    fh_uh  = alpha*fh_uh + (1-alpha)*fheu_uh;
end

function [fh,fh_udg,fh_uh] = fboueuler(nl,udg,uh)

[ng,nch] = size(uh);
nc = size(udg,2);

u = udg(:,1:nch);
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
