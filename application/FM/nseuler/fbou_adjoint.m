function [fh,fh_udg,fh_uh] = fbou_adjoint(ib,ui,nl,p,udg,uh,param,time)
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
%      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q

[ng,nc] = size(udg);
nch = size(uh,2);

alfa = 0.0; % angle of attack
switch ib
    case 1  % Drag        
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);   
        beta = [cos(alfa),sin(alfa)];
        fh = fh(:,2)*beta(1)+fh(:,3)*beta(2);
        fh_udg = fh_udg(:,2,:)*beta(1)+fh_udg(:,3,:)*beta(2);
        fh_uh = fh_uh(:,2,:)*beta(1)+fh_uh(:,3,:)*beta(2);
        fh_udg = reshape(fh_udg,[ng nc]);
        fh_uh = reshape(fh_uh,[ng nch]);
%         fh = uh(:,1);
%         fh_udg = zeros(ng,nc);
%         fh_uh = zeros(ng,nch);
% %         fh_uh(:,1) = 1;
%         r    = uh(:,1);
%         ru   = uh(:,2);
%         rv   = uh(:,3);
%         rE   = uh(:,4);
%         r1   = 1./r;
%         u    = ru.*r1;
%         v    = rv.*r1;        
%         q    = 0.5*(u.*u+v.*v);
%         gam1 = param{1} - 1;
%         p    = gam1*(rE-r.*q);
%         fh   = p(:,1);        
%         fh_uh(:,1) = -gam1*(ru.^2./(2*r.^2) + rv.^2./(2*r.^2) - r.*(ru.^2./r.^3 + rv.^2./r.^3));
%         fh_uh(:,2) = -(gam1*ru)./r;
%         fh_uh(:,3) = -(gam1*rv)./r;
%         fh_uh(:,4) = gam1;
% 
    case 2  % Lift                
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);    
        beta = [-sin(alfa),cos(alfa)];
%         fh = fh(:,2).*nl(:,1) + fh(:,3).*nl(:,2);
%         fh_udg = bsxfun(@times,fh_udg(:,2,:),nl(:,1))+bsxfun(@times,fh_udg(:,3,:),nl(:,2));
%         fh_uh = bsxfun(@times,fh_uh(:,2,:),nl(:,1))+bsxfun(@times,fh_uh(:,3,:),nl(:,2));                
        fh = fh(:,2)*beta(1)+fh(:,3)*beta(2);
        fh_udg = fh_udg(:,2,:)*beta(1)+fh_udg(:,3,:)*beta(2);
        fh_uh = fh_uh(:,2,:)*beta(1)+fh_uh(:,3,:)*beta(2);
        fh_udg = reshape(fh_udg,[ng nc]);
        fh_uh = reshape(fh_uh,[ng nch]);
    case 3  % zero                
        fh = zeros(ng,1);        
        fh_udg = zeros(ng, nc);
        fh_uh = zeros(ng, nch);        
    otherwise
        error('unknown boundary type');
end

