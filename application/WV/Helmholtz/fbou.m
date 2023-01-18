function [fh,fh_udg,fh_uh] = fbou(ib,ui,nl,p,udg,uh,param,time)
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
nch = 1;
nq = nc-nch;
nd = nq;

k     = param{1};
%kappa = param{2};
tau   = param{3};

u     = udg(:,nch);
%q     = reshape(udg(:,nch+1:nc),[ng nch nd]);

switch ib
    case 1  % Dirichlet
        fh = tau.*(ui-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = -tau;
    case 2  % "Extrapolate  m = u 
        fh = tau.*(u-uh);
        fh_u = tau;
        fh_q = zeros(ng,nch,nq);
        fh_udg = cat(3,fh_u,fh_q);
        fh_uh = -tau;
    case 3  % Prescribed flux
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        fh = fh + ui;
    case 4  % sound-soft
        kx = 1;
        ky = 0;
        ui = exp(sqrt(-1)*k*(kx*p(:,1)+ky*p(:,2)));        
        fh = tau.*(ui+uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = tau;    
    case 5  % sound-hard
        t = pi/4;
        kx = cos(t);
        ky = sin(t);
        ui = sqrt(-1)*k*(kx*nl(:,1)+ky*nl(:,2)).*exp(sqrt(-1)*k*(kx*p(:,1)+ky*p(:,2)));        
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);        
        fh = fh - ui;    
    case 6 % absorbing            
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        fh = fh + sqrt(-1)*k*uh;
        fh_uh = fh_uh + sqrt(-1)*k;
    case 7        
        kx = 1;
        ky = 0;
        ui = exp(sqrt(-1)*k*(kx*p(:,1)+ky*p(:,2)));        
        fh = tau.*(ui-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh = -tau;            
    otherwise
        error('unknown boundary type');
end

