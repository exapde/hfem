function [sr,sr_udg] = irksource(p,udg,param,time)
%SOURCE Source function for time spectral method
%   [sr,sr_udg] = source(p,udg,param,time)
%
%      P(N,ND)               Coordinates for N points
%      UDG(N,NC)             Unknown vector for N points with NC components
%      PARAM                 Parameter list
%      SR(N,NCH):            Source term at N points
%      SR_UDG(N,NCH,NC):     Jacobian of the source term w.r.t. [U Q]
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
D = irk.D;            % IRK matrix

[ng,nc] = size(udg);
nch = Ns*nch0;
nd1 = nc/nch;
sr = zeros(ng,nch);
sr_udg = zeros(ng,nch,nc); 
for i=1:Ns
    % subset of U variables describing this snapshot
    uvars = (i-1)*nch0+1:i*nch0; 
    
    % subset of UDG variables describing this snapshot
    udgvars = zeros(nch0,nd1);
    for j=1:nch0; udgvars(j,:) = (i-1)*nch0+j:nch:nc; end
    udgvars= sort(udgvars(:))';
    
    % time coordinate of this snapshot    
    timei = time + t(i)*dt;    
    
    % source term S
    [sr(:,uvars),sr_udg(:,uvars,udgvars)] = source(p,udg(:,udgvars),param{1},timei);
    
    % source term -D*U 
    for j=1:nch0
        sr(:,uvars(j)) = sr(:,uvars(j)) - udg(:,j:nch0:nch)*(D(i,:)');
        sr_udg(:,uvars(j),j:nch0:nch) = sr_udg(:,uvars(j),j:nch0:nch) ...
          - permute(repmat(D(i,:),ng,1).*ones(size(udg(:,j:nch0:nch))),[1 3 2]);           
    end        
end

