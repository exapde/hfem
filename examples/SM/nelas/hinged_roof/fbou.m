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

% Coordinates center of the axis rotation
R = 254;
xc = R*sin(1);
yc = R*cos(1);

switch (ib)
    case 1 % Dirichlet                      
        fh = zeros(ng1,nch);
        fh_udg = zeros(ng1,nch,nc);
        fh_uh = zeros(ng1,nch,nch);
        fh(:,1) = tau*(p(:,1)-uh(:,1));
        fh(:,2) = tau*(p(:,2)-uh(:,2));
        fh(:,3) = tau*(p(:,3)-uh(:,3));
        fh_uh(:,1,1) = -tau;
        fh_uh(:,2,2) = -tau;
        fh_uh(:,3,3) = -tau;
    case 2 % Symmetry : Block displacement along X axis
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        fh_udg(:,1,:) = zeros(ng1,1,nc);
        fh(:,1) = tau*(p(:,1)-uh(:,1));
        fh_uh(:,1,1) = -tau;
    case 3 % Symmetry : Block displacement along Y axis
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        fh_udg(:,2,:) = zeros(ng1,1,nc);
        fh(:,2) = tau*(p(:,2)-uh(:,2));
        fh_uh(:,2,2) = -tau;
    case 4 % Symmetry : Block displacement along Z axis
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        fh_udg(:,3,:) = zeros(ng1,1,nc);
        fh(:,3) = tau*(p(:,3)-uh(:,3));
        fh_uh(:,3,3) = -tau;
    case 5 % Neumman bc
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);            
        fh = fh+ui;
    case 6 % Block displacements in 0xy plane
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);            
        fh_udg(:,1,:) = zeros(ng1,1,nc);
        fh_udg(:,2,:) = zeros(ng1,1,nc);
        fh(:,1) = tau*(p(:,1)-uh(:,1));
        fh(:,2) = tau*(p(:,2)-uh(:,2));
        fh_uh(:,1,1) = -tau;
        fh_uh(:,2,2) = -tau;
    case 7 % Hinged Boundary condition
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);  
        %r0 = sqrt((xc - p(:,1)).^2 +(yc - p(:,2)).^2);
        %r  = sqrt((xc -uh(:,1)).^2 +(yc -uh(:,2)).^2);
        %fh(:,1) = tau*(r0 - r);
        %fh_udg(:,1,:) = zeros(ng1,1,nc);
        %fh(:,2) = tau*(r0 - r);
        %fh_uh(:,1,1) = -tau * (uh(:,1)-xc) ./ r(:);
        %fh_uh(:,1,2) = -tau * (uh(:,2)-yc) ./ r(:);
        %fh_uh(:,2,1) = -tau * (uh(:,1)-xc) ./ r(:);
        %fh_uh(:,2,2) = -tau * (uh(:,2)-yc) ./ r(:);
        % Block displacements in the Z direction
        fh(:,3) = tau*(p(:,3)-uh(:,3));
        fh_udg(:,3,:) = zeros(ng1,1,nc);
        fh_uh(:,3,3) = -tau;
    otherwise
        error('unknown boundary type');
end
