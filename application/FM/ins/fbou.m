function [fh,fh_udg,fh_uh] = fbou(ib,ui,nl,pg,udg,uh,param,time)
%FBOU boundary flux function

%      IC                    Boundary number
%      IB                    Boundary type
%      FUNC                  Function to compute boundary data
%      NL(N,ND)              Normal N points
%      X(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      P(N,1)                Pressure vector for N points with 1 component
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHP(N,NC,1):          Jacobian of the flux flux vector w.r.t. P
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. M

[ng,nch] = size(uh);
nd = size(pg,2);

%mu = param{1};
tau = param{end};
uinf = ui(:,1:nch);

switch (ib)
    case 1 % Dirichlet bc
        fh = tau*(uinf-uh);
        fh_udg = zeros(ng,nch,(nd+1)*nch+1);
        fh_uh = -bsxfun(@times,ones(ng,1,1),reshape(tau*eye(nch),[1 nch nch]));
    case 2 % homogeneous total stress bc        
        [fh,fh_udg,fh_uh] = fhatbou1(nl,pg,udg,uh,param,time);           
    case 3 % homogeneous total gradient bc        
        [fh,fh_udg,fh_uh] = fhat(nl,pg,udg,uh,param,time);                   
    otherwise
        error('unknown boundary type');
end

% total stress bc
function [fh,fh_udg,fh_uh] = fhatbou1(nl,pg,udg,uh,param,time)

[ngf,nch] = size(uh);
nc = size(udg,2);
nd = size(pg,2);
mu  = param{1};
tau = param{end};

p = udg(:,nch+1);
q = udg(:,nch+2:end);
q = reshape(q,[ngf nch nd]);

% stress tensor
f = mu*(q + permute(q,[1 3 2]));
for d=1:nd
    f(:,d,d) = f(:,d,d) + p;    
end

f_udg = zeros(ngf,nch,nd,nc);
h_udg = zeros(ngf,nch,nd,nc);
for d=1:nd
    f_udg(:,d,d,nch+1) = 1;
    for j=1:nch
        k = (nch+1) + (d-1)*nch+j;
        f_udg(:,j,d,k) = mu;        
        h_udg(:,d,j,k) = mu;        
    end
end
f_udg = f_udg + h_udg;

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


% % % total stress bc
% function [fh,fh_udg,fh_uh] = fhatbou2(nl,p,udg,uh,param,time)
% 
% [ng,nch] = size(uh);
% [ng,nc] = size(udg);
% nd = nch;
% nq = nc-nch-1;
% 
% mu  = param{1};
% tau = param{end};
% 
% u = udg(:,1:nch);
% pr = udg(:,nch+1);
% q = reshape(udg(:,nch+2:nc),[ng nch nd]);
% 
% str = mu*(q+permute(q,[1 3 2]));
% 
% fh = permute(mapContractK(str,nl,2,3,1,2,[],1),[2 1]) + bsxfun(@times,pr,nl) + tau*(u-uh);
% 
% if nargout > 1
%     fh_u = bsxfun(@times,ones(ng,1,1),reshape(tau*eye(nch),[1 nch nch]));
%     fh_p = nl;
%     
%     fh_q  = zeros(ng,nch,nq); 
%     fh_q(:,1,1) = 2*mu*nl(:,1);
%     fh_q(:,1,2) = mu*nl(:,2);
%     fh_q(:,1,3) = mu*nl(:,2);
%     fh_q(:,2,2) = mu*nl(:,1);
%     fh_q(:,2,3) = mu*nl(:,1);
%     fh_q(:,2,4) = 2*mu*nl(:,2);
%         
%     fh_udg = cat(3,cat(3,fh_u,fh_p),fh_q);
%         
%     fh_uh = -fh_u;
% end
% 
% 
