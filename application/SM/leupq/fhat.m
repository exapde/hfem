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

%f = kappa*reshape(q,[ng nch nd]);

% kappa = param{1};
% tau   = param{end};
% 
% nx   = nl(:,1);
% ny   = nl(:,2);
% 
% [f,f_udg] = flux(p,udg,param,time);
% 
% fh(:,1) = f(:,1,1).*nx + f(:,1,2).*ny + tau*(udg(:,1)-uh(:,1));
% fh(:,2) = f(:,2,1).*nx + f(:,2,2).*ny + tau*(udg(:,2)-uh(:,2));
% 
% if nargout > 1            
%     for i = 1:nc
%         fh_udg(:,1,i) = f_udg(:,1,1,i).*nx + f_udg(:,1,2,i).*ny;
%         fh_udg(:,2,i) = f_udg(:,2,1,i).*nx + f_udg(:,2,2,i).*ny;
%     end    
%     
% %     fh_udg(:,1,:) = squeeze(f_udg(:,1,1,:)).*repmat(nx,[1 nc]) + squeeze(f_udg(:,1,2,:)).*repmat(ny,[1 nc]);
% %     fh_udg(:,2,:) = squeeze(f_udg(:,2,1,:)).*repmat(nx,[1 nc]) + squeeze(f_udg(:,2,2,:)).*repmat(ny,[1 nc]);
%        
%     fh_udg(:,1,1) = tau;
%     fh_udg(:,2,2) = tau;           
%                 
%     fh_uh = zeros(ng1,nch,nch);    
%     fh_uh(:,1,1) = -tau;
%     fh_uh(:,2,2) = -tau;        
% end
