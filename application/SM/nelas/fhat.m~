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

nc = size(udg,2);

if nc==7
    [fh,fh_udg,fh_uh] = fhat2d(nl,pg,udg,uh,param,time);
else
    [fh,fh_udg,fh_uh] = fhat3d(nl,pg,udg,uh,param,time);
end

return;

[ngf,nch] = size(uh);
nc = size(udg,2);
nd = size(pg,2);
tau = param{end};

[f,f_udg] = flux(pg,udg,param,time);

fh = tau*(udg(:,1:nch) - uh);
for ic=1:nch
    fh(:,ic) = fh(:,ic) + sum(reshape(f(:,ic,:),[ngf nd]).*nl,2);
end

fh_udg = zeros(ngf,nch,nc);
fh_uh = zeros(ngf,nch,nch);    
for d=1:nd
    fh_udg(:,d,d) = tau;
    fh_uh(:,d,d) = -tau;
    for ic=nch+1:nc
        fh_udg(:,d,ic) = sum(reshape(f_udg(:,d,:,ic),[ngf nd]).*nl,2);    
    end
end

return;

[ng1,nc] = size(udg);
nch = 2;

if size(p,2)<=2
    mu = param{1};
else
    mu = p(:,3);
end

tau    = param{end};

u1  = udg(:,1);
u2  = udg(:,2);
pres= udg(:,3);
F11 = udg(:,4);
F21 = udg(:,5);
F12 = udg(:,6);
F22 = udg(:,7);

if tau==0
    beta = 16; gamma=16;
    tau  = beta*mu.*(F11.^2 + F12.^2 + F21.^2 + F22.^2 + 1) + gamma*abs(pres);
else    
    beta=0; gamma=0;
end

wcase=1;
[P,dPdF,dPdp] = Deriv(F11,F12,F21,F22,pres,mu,wcase);
%dPdp(:,[2 3]) = dPdp(:,[3 2]);
%dgdF(:,[2 3]) = dgdF(:,[3 2]);

fh(:,1) = P(:,1).*nl(:,1) + P(:,2).*nl(:,2) + tau*(udg(:,1)-uh(:,1));
fh(:,2) = P(:,3).*nl(:,1) + P(:,4).*nl(:,2) + tau*(udg(:,2)-uh(:,2));

if nargout > 1
    fh_u  = zeros(ng1,nch,nch);
    fh_p  = zeros(ng1,nch,1);
    fh_q  = zeros(ng1,nch,2*nch);
    fh_uh = zeros(ng1,nch,nch);
    
    fh_u(:,1,1) = tau;
    fh_u(:,2,2) = tau;
        
    fh_p(:,1,1) = dPdp(:,1).*nl(:,1) + dPdp(:,2).*nl(:,2);
    fh_p(:,2,1) = dPdp(:,3).*nl(:,1) + dPdp(:,4).*nl(:,2);
    
    fh_q(:,1,1) = dPdF(:,1).*nl(:,1)+dPdF(:,2).*nl(:,2);
    fh_q(:,1,2) = dPdF(:,5).*nl(:,1)+dPdF(:,6).*nl(:,2);
    fh_q(:,1,3) = dPdF(:,3).*nl(:,1)+dPdF(:,4).*nl(:,2);
    fh_q(:,1,4) = dPdF(:,7).*nl(:,1)+dPdF(:,8).*nl(:,2);
    
    fh_q(:,2,1) = dPdF(:,9).*nl(:,1) +dPdF(:,10).*nl(:,2);
    fh_q(:,2,2) = dPdF(:,13).*nl(:,1)+dPdF(:,14).*nl(:,2);
    fh_q(:,2,3) = dPdF(:,11).*nl(:,1)+dPdF(:,12).*nl(:,2);
    fh_q(:,2,4) = dPdF(:,15).*nl(:,1)+dPdF(:,16).*nl(:,2);
        
    fh_udg = cat(3,cat(3,fh_u,fh_p),fh_q);
        
    fh_uh(:,1,1) = -tau;
    fh_uh(:,2,2) = -tau;
end



