function [fh,fh_udg,fh_uh] = fbou(ib,ui,nl,pg,udg,uh,param,time)
%FBOU boundary flux function

%      IC                    Boundary number
%      IB                    Boundary type
%      FUNC                  Function to compute boundary data
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
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

[ngf,nch] = size(uh);
nc = size(udg,2);

tau   = param{end};
ginf = repmat(ui(1,1:nch),ngf,1);

switch (ib)
    case 1 % Dirichlet bc        
        fh = tau*(ginf-uh);
        fh_udg = zeros(ngf,nch,nc);
        fh_uh = -bsxfun(@times,ones(ngf,1,1),reshape(tau*eye(nch),[1 nch nch]));
    case 2 % Neumman bc        
        [fh,fh_udg,fh_uh] = fhat(nl,pg,udg,uh,param,time);
        fh = fh + ginf;
    case 3 % mixed bc
        uinf = repmat(ui(1,1),ngf,1);
        finf = repmat(ui(1,2),ngf,1);
        
        fh   = zeros(ngf,nch);        
        fh_uh = zeros(ngf,nch,nch);
        fh_udg = zeros(ngf,nch,nc);
        
        fh(:,1) = tau*(uinf - uh(:,1));        
        fh_uh(:,1,1) = -tau;
        
        [gh,gh_udg,gh_uh] = fhat(nl,pg,udg,uh,param,time);
        fh(:,2) = gh(:,2) + finf;
        fh_uh(:,2,:)  = gh_uh(:,2,:);
        fh_udg(:,2,:) = gh_udg(:,2,:);        
    case 4 % mixed bc
        uinf = repmat(ui(1,1),ngf,1);
        finf = repmat(ui(1,2),ngf,1);
        
        fh   = zeros(ngf,nch);        
        fh_uh = zeros(ngf,nch,nch);
        fh_udg = zeros(ngf,nch,nc);
        
        fh(:,2) = tau*(uinf - uh(:,2));        
        fh_uh(:,2,2) = -tau;
        
        [gh,gh_udg,gh_uh] = fhat(nl,pg,udg,uh,param,time);
        fh(:,1) = gh(:,1) + finf;
        fh_uh(:,1,:)  = gh_uh(:,1,:);
        fh_udg(:,1,:) = gh_udg(:,1,:);            
    otherwise
        error('unknown boundary type');
end

