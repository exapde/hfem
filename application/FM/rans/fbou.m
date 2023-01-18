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
q = reshape(udg(:,nch+1:3*nch),ng1,nch,2);

switch ib
    case 1  % Far field        
        [an,anm] = getan(nl,uh,param,3);
        [An,Anm] = getan(nl,uh,param,2);
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_u = an+An;
        fh_q = zeros(ng1,nch,nch,2);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-2*An;
    case 2  % Adiabatic Wall
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time,[]);
        fh_u = fh_udg(:,:,1:nch);
        fh_q = fh_udg(:,:,nch+1:3*nch);
                
        fh(:,[2 3 5]) = uh(:,[2 3 5]);
        fh_u(:,[2 3 5],:) = 0;
        fh_q(:,[2 3 5],:) = 0;
        fh_uh(:,[2 3 5],:) = 0;
        %fh_uh(:,[2 3 5],[2 3 5]) = multiprod(one,eye(3,3),2,[1,2]);
        fh_uh(:,2,2) = 1;
        fh_uh(:,3,3) = 1;
        fh_uh(:,5,5) = 1;                                
    case 3 % Isothermal wall
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time,[]);
        fh_u = fh_udg(:,:,1:nch);
        fh_q = fh_udg(:,:,nch+1:3*nch);
        
        fh(:,[2 3 5]) = uh(:,[2 3 5]);
        fh(:,4) = gam*gam1*M2*uh(:,4)./uh(:,1) - ui(4);
        fh_u(:,2:5,:) = 0;
        fh_q(:,2:5,:) = 0;
        fh_uh(:,[2 3 5],[1,4]) = 0;
        fh_uh(:,2,2) = 1;
        fh_uh(:,3,3) = 1;
        fh_uh(:,5,5) = 1;
        fh_uh(:,4,:) = gam*gam1*M2*[-uh(:,4)./uh(:,1).^2, zero, zero, 1./uh(:,1), zero];
    case 4 % Prescribed pressure
        p  = repmat(uinf(1,nch-1),ng1,1);        
        uinf = [ u(:,1), u(:,2), u(:,3), p./(gam-1) + 0.5*u(:,2).*u(:,2)./u(:,1) + 0.5*u(:,3).*u(:,3)./u(:,1), u(:,5)];
        uinfu = zeros(ng1,nch,nch);
        uinfu(:,1,1) = 1;
        uinfu(:,2,2) = 1;
        uinfu(:,3,3) = 1;
        uinfu(:,4,:) = [-0.5*(u(:,2).*u(:,2)+u(:,3).*u(:,3))./(u(:,1).*u(:,1)), u(:,2)./u(:,1), u(:,3)./u(:,1), zero, zero];
        uinfu(:,5,5) = 1;
        
        [an,anm] = getan(nl,uh,param,0);
        [An,Anm] = getan(nl,uh,param,1);
        fh = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fh_u = an+An - permute(mapContractK(an-An,uinfu,2,3,1,2,3,1),[3 1 2]);
        fh_q = zeros(ng1,nch,nch,2);
        fh_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-2*An;
    case 5  % Slip Adiabatic Wall
        un = u(:,2).*nl(:,1) + u(:,3).*nl(:,2);        
        
        uinf(:,1:nch-1) = u(:,1:nch-1);        
        uinfu = zeros(ng1,nch,nch);
        uinfu(:,1,1) = 1;         
        uinf(:,2) = uinf(:,2) - nl(:,1).*un;
        uinfu(:,2,:) = [zeros(ng1,1) ones(ng1,1)-nl(:,1).*nl(:,1) -nl(:,1).*nl(:,2) zeros(ng1,1) zero]; 
        uinf(:,3) = uinf(:,3) - nl(:,2).*un;
        uinfu(:,3,:) = [zeros(ng1,1) -nl(:,1).*nl(:,2) ones(ng1,1)-nl(:,2).*nl(:,2) zeros(ng1,1) zero];                
        uinfu(:,4,4) = 1;
        %uinfu(:,5,5) = 1;
                        
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time,[]);        
        fh_u = fh_udg(:,:,1:nch);
        fh_q = fh_udg(:,:,nch+1:3*nch);
                
        ind = [1 2 3 5];        
        fh(:,ind) = uinf(:,ind) - uh(:,ind);        
        
        fh_u(:,ind,:) = uinfu(:,ind,:);
        fh_q(:,ind,:) = 0;        
        fh_uh(:,ind,:) = 0;
        for j=1:length(ind)
            fh_uh(:,ind(j),ind(j)) = -1;
        end

    case 6 %%% symmetric BC
        fh = zeros(ng1,nch);
        fh(:,1) = q(:,1,1).*nl(:,1) + q(:,1,2).*nl(:,2); %%% drho/dn = 0
        t1620 = 1.0./u(:,1);
        t1621 = gam1;
        t1622 = 1.0./u(:,1).^2;
        t1623 = u(:,2).^2;
        t1624 = u(:,3).^2;
        t1625 = t1623+t1624;
        t1626 = 1.0./uh(:,1).^3;
        fh(:,2) = nl(:,1).*t1621.*(q(:,4,1)-uh(:,2).*q(:,2,1).*t1620-uh(:,3).*q(:,3,1).*t1620+q(:,1,1).*t1622.*t1625.*(1.0./2.0)) + ... %%% dP/dn = 0
                  nl(:,2).*t1621.*(q(:,4,2)-uh(:,2).*q(:,2,2).*t1620-uh(:,3).*q(:,3,2).*t1620+q(:,1,2).*t1622.*t1625.*(1.0./2.0));
        fh(:,3) = q(:,5,1).*nl(:,1) + q(:,5,2).*nl(:,2); %%% d(nu_tilde)/dn = 0
        fh(:,4) = uh(:,2).*nl(:,1) + uh(:,3).*nl(:,2);   %%% non-penetration
        fh(:,5) = (q(:,2,1).*nl(:,2)-q(:,3,1).*nl(:,1)).*nl(:,1) + (q(:,2,2).*nl(:,2)-q(:,3,2).*nl(:,1)).*nl(:,2); %%% d(rho u_tang1ential)/dn = 0
        
        fh_u = zeros(ng1,nch,nch);
        fh_uh = zeros(ng1,nch,nch);
        
        fh_uh(:,2,1) = nl(:,1).*t1621.*(uh(:,2).*q(:,2,1).*t1622+uh(:,3).*q(:,3,1).*t1622-q(:,1,1).*t1625.*t1626)+nl(:,2).*t1621.*(uh(:,2).*q(:,2,2).*t1622+uh(:,3).*q(:,3,2).*t1622-q(:,1,2).*t1625.*t1626);
        fh_uh(:,2,2) = -nl(:,1).*t1621.*(q(:,2,1).*t1620-uh(:,2).*q(:,1,1).*t1622)-nl(:,2).*t1621.*(q(:,2,2).*t1620-uh(:,2).*q(:,1,2).*t1622);
        fh_uh(:,2,3) = -nl(:,1).*t1621.*(q(:,3,1).*t1620-uh(:,3).*q(:,1,1).*t1622)-nl(:,2).*t1621.*(q(:,3,2).*t1620-uh(:,3).*q(:,1,2).*t1622);
        fh_uh(:,2,4) = 0;
        fh_uh(:,2,5) = 0;
        
        fh_uh(:,4,2) = nl(:,1);
        fh_uh(:,4,3) = nl(:,2);
       
        
        fh_q = zeros(ng1,nch,2*nch);
        fh_q(:,1,1) = nl(:,1);
        fh_q(:,1,6) = nl(:,2);
        
        fh_q(:,2,1) = nl(:,1).*t1621.*t1622.*t1625.*(1.0./2.0);
        fh_q(:,2,2) = -nl(:,1).*uh(:,2).*t1620.*t1621;
        fh_q(:,2,3) = -nl(:,1).*uh(:,3).*t1620.*t1621;
        fh_q(:,2,4) = nl(:,1).*t1621;
        fh_q(:,2,5) = 0;
        fh_q(:,2,6) = nl(:,2).*t1621.*t1622.*t1625.*(1.0./2.0);
        fh_q(:,2,7) = -nl(:,2).*uh(:,2).*t1620.*t1621;
        fh_q(:,2,8) = -nl(:,2).*uh(:,3).*t1620.*t1621;
        fh_q(:,2,9) = nl(:,2).*t1621;
        fh_q(:,2,10)= 0;
        
        fh_q(:,3,5) = nl(:,1);
        fh_q(:,3,10) = nl(:,2);
        
        fh_q(:,5,2) = nl(:,2).*nl(:,1);
        fh_q(:,5,3) = -nl(:,1).^2;
        fh_q(:,5,7) = nl(:,2).^2;
        fh_q(:,5,8) = -nl(:,1).*nl(:,2);        
    otherwise
        error('unknown boundary type');
end

fh_udg = cat(3,fh_u,reshape(fh_q,ng1,nch,2*nch));



