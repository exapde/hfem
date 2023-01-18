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

% 3D plane strain neo-Hookean materials

[ng1,nc] = size(udg);
nch = 3;

mu     = param{1};
kappa  = param{2};
beta   = 1;
tau    = param{end};

if tau==0
    tau = mu;
end


switch (ib)
    case 1 % Dirichlet bc                        
        fh = zeros(ng1,nch);
        fh_udg = zeros(ng1,nch,nc);
        fh_uh = zeros(ng1,nch,nch);
        fh(:,1) = tau*(p(:,1)-uh(:,1));
        fh(:,2) = tau*(p(:,2)-uh(:,2));
        fh(:,3) = tau*(p(:,3)-uh(:,3));
        fh_uh(:,1,1) = -tau;
        fh_uh(:,2,2) = -tau;
        fh_uh(:,3,3) = -tau;
    case 2 % Neumman bc (Tractions equivalent to a moment)
        xm = mean(p(:,1));
        zm = mean(p(:,3));
        nx = mean(nl(:,1));
        nz = mean(nl(:,3));
        % Distance of the Gauss points from the plate's mid-surface
        td = (p(:,3)-zm)*nx - (p(:,1)-xm)*nz;
        % Select points for different polynomial orders
        switch size(udg,1)
            case 4
                np1 = 2; np2 = 3;
            case 9
                np1 = 3; np2 = 7;
            case 16
                np1 = 4; np2 = 13;
        end
        % Computing normal in the current configuration
        v1 = udg(np1,1:3) - udg(1,1:3);
        v2 = udg(np2,1:3) - udg(1,1:3);
        n  = cross(v1,v2);
        n = 1./norm(n) * n;
        % Applying linear traction equivalent to moment
        sigmamax = -ui(1,2);
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        fh(:,1) = fh(:,1) + sigmamax*n(1)*td ;
        fh(:,3) = fh(:,3) + sigmamax*n(3)*td ;
    case 3 % Neumman bc
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);            
        fh = fh+ui;
    otherwise
        error('unknown boundary type');
end
