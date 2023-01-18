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
nch = 1;
nq = nc-nch;

tau   = param{end};
u     = udg(:,nch);
switch ib
    case 3  % Dirichlet
        fh = uh;
        fh_udg = zeros(ng,nc);
        fh_uh = ones(ng,nch);
    case 2  % "Extrapolate  m = u 
        fh = tau.*(u-uh);
        fh_u = tau;
        fh_q = zeros(ng,nch,nq);
        fh_udg = cat(3,fh_u,fh_q);
        fh_uh = -tau;
    case 1  % Prescribed flux
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);    
        fh_udg = reshape(fh_udg,[ng nc]);
        fh_uh = reshape(fh_uh,[ng nch]);
    otherwise
        error('unknown boundary type');
end

