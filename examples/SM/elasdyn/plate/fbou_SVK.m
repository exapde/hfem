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
nch = 3;

mu     = param{1};
lambda = param{2};
tau    = param{3};
rho    = param{4};
A      = param{5};

% Gaussian excitation parameters
% t0 = 20.; sigma = 5.;
t0 = 5.; sigma = 1.;
dgauss = -0.01*(time-t0)/(sigma^2) * exp(-0.5*((time-t0)/sigma)^2);

if tau==0
    tau = mu;
end


switch (ib)
    case 1 % Dirichlet bc                        
        fh = zeros(ng1,nch);
        fh_udg = zeros(ng1,nch,nc);
        fh_uh = zeros(ng1,nch,nch);
        fh(:,1) = tau*(0*p(:,1)-uh(:,1));
        fh(:,2) = tau*(0*p(:,2)-uh(:,2));
        fh(:,3) = tau*(0*p(:,3)-uh(:,3));
        fh_uh(:,1,1) = -tau;
        fh_uh(:,2,2) = -tau;  
        fh_uh(:,3,3) = -tau;  
    case 2 % Neumman bc         
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);            
        fh = fh+ui;          
    case 3 % Neumman bc         
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);            
        %fh = fh+sin(time)*ui;
        %fh = fh+dgauss*ui;
        x = p(:,1);
        y = p(:,2);        
        s = 0*fh;
        % Saint Venant-Kirchhoff
        s(:,1) = A*mu*pi*sin(pi*time)*cos(pi*x).*sin(pi*y).*nl(:,3);
        s(:,2) = A*mu*pi*sin(pi*time)*cos(pi*y).*sin(pi*x).*nl(:,3);
        s(:,3) = 0.5* (A^2*pi^2*sin(pi*time)^2*(lambda + 2*mu)...
          * (sin(pi*x).^2 - 2*sin(pi*x).^2 .*sin(pi*y).^2 + sin(pi*y).^2)).*nl(:,3);
        fh = fh + s;
    case 4 % Imposed time varying Dirichlet BC
        fh = zeros(ng1,nch);
        fh_udg = zeros(ng1,nch,nc);
        fh_uh = zeros(ng1,nch,nch);
        fh(:,1) = tau*(dgauss-uh(:,1));
        fh(:,2) = tau*(0*p(:,2)-uh(:,2));
        fh(:,3) = tau*(0*p(:,3)-uh(:,3));
        %fh(:,3) = tau*(dgauss-uh(:,3));
        fh_uh(:,1,1) = -tau;
        fh_uh(:,2,2) = -tau;
        fh_uh(:,3,3) = -tau;
    otherwise
        error('unknown boundary type');
end
