function [fh,fh_udg,fh_uh] = irkfhat(nl,p,udg,uh,param,time)
%FHAT flux function
%   [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time,signe)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      UDG(N,NC)             Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NCH):            Volume flux at N points
%      FH_UDG(N,NCH,NC):     Jacobian of the flux vector w.r.t. [U Q]
%      FH_UH(N,NCH):         Jacobian of the flux vector w.r.t. UH
%--------------------------------------------------------------------------

% iflux = irk.iflux;
% ifhat = irk.ifhat;
% ifbou = irk.ifbou;
% isource = irk.isource;

irk = param{2};
Ns = irk.Ns;          % number of snapshots 
nch0 = irk.nch0;      % number of primary field vars in single snapshot
dt = irk.dt;          % time step size
t = irk.t;            % time coefficients

[ng,nc] = size(udg);
nch = Ns*nch0;
nd1 = nc/nch;
fh = zeros(ng,nch);
fh_udg = zeros(ng,nch,nc);
fh_uh = zeros(ng,nch);
for i=1:Ns
    % subset of U variables describing this snapshot
    uvars = (i-1)*nch0+1:i*nch0; 
    
    % subset of UDG variables describing this snapshot
    udgvars = zeros(nch0,nd1);
    for j=1:nch0; udgvars(j,:) = (i-1)*nch0+j:nch:nc; end
    udgvars= sort(udgvars(:))';
    
    % time coordinate of this snapshot    
    timei = time + t(i)*dt;    
    
    [fh(:,uvars),fh_udg(:,uvars,udgvars),fh_uh(:,uvars,uvars)] ...
        = fhat(nl,p,udg(:,udgvars),uh(:,uvars),param{1},timei);
end

