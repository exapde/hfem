function [fh,fh_udg,fh_uh] = fbou(ib,ui,nl,p,udg,uh,param,time)
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

% 2D plane strain neo-Hookean materials

[ng1,nc] = size(udg);
nch = 2;

mu     = param{1};
tau    = 1e2;

switch (ib)
    case 1 % Dirichlet bc                        
        fh = multiprod(tau,p-uh,[2,3],2);
        fh_udg = zeros(ng1,nch,nc);
        fh_uh = zeros(ng1,nch,nch);
        fh_uh(:,1,1) = -tau;
        fh_uh(:,2,2) = -tau;  
    case 2 % Neumman bc 
        finf = repmat(ui(1,1:nch),ng1,1);
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);    
        fh = fh + finf;        
    case 3 % mixed Dirichlet and stress bcs
        uinf = p(:,1); %repmat(ui(1,1),ng1,1);
        finf = repmat(ui(1,1),ng1,1);
        
        fh   = zeros(ng1,2);        
        fh_uh = zeros(ng1,nch,nch);
        fh_udg = zeros(ng1,nch,nc);
        
        fh(:,1) = tau*(uinf - uh(:,1));        
        fh_uh(:,1,1) = -tau;
        
        [gh,gh_udg,gh_uh] = fhat(nl,p,udg,uh,param,time);
        fh(:,2) = gh(:,2) + finf;
        fh_uh(:,2,:)  = gh_uh(:,2,:);
        fh_udg(:,2,:) = gh_udg(:,2,:);        
    case 4 % mixed Dirichlet and stress bcs
        uinf = p(:,2); %repmat(ui(1,1),ng1,1);
        finf = repmat(ui(1,1),ng1,1);
        
        fh   = zeros(ng1,2);        
        fh_uh = zeros(ng1,nch,nch);
        fh_udg = zeros(ng1,nch,nc);
        
        fh(:,2) = tau*(uinf - uh(:,2));        
        fh_uh(:,2,2) = -tau;
        
        [gh,gh_udg,gh_uh] = fhat(nl,p,udg,uh,param,time);
        fh(:,1) = gh(:,1) + finf;
        fh_uh(:,1,:)  = gh_uh(:,1,:);
        fh_udg(:,1,:) = gh_udg(:,1,:);     
    otherwise
        error('unknown boundary type');
end

