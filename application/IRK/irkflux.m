function [f,f_udg] = irkflux(p,udg,param,time)
%FLUX Volume flux function for time-spectral method
%   [f,fu,fq] = flux(p,u,q,param)
%
%      P(N,ND)              Coordinates for N points
%      U(N,NC)              Unknown vector for N poiNss with NC components
%      Q(N,NC,ND)           Flux vector for N points with NC components in the
%                           coordinate directions
%      PARAM                Parameter list
%      F(N,NCH,ND):         Volume flux at N points
%      F_U(N,NCH,ND):       Jacobian of the flux vector w.r.t. U
%      F_Q(N,NCH,ND,ND):    Jacobian of the flux vector w.r.t. Q
%      F_UDG(N,NCH,ND,NC):  Jacobian of the flux vector w.r.t. [U Q]

% iflux = irk.iflux;
% ifhat = irk.ifhat;
% ifbou = irk.ifbou;
% isource = irk.isource;

irk = param{2};
Ns = irk.Ns;          % number of snapshots 
nch0 = irk.nch0;      % number of primary field vars in single snapshot
dt = irk.dt;          % time step size
t = irk.t;            % time coefficients
nd = irk.nd;

[ng,nc] = size(udg);
nch = Ns*nch0;
nd1 = nc/nch;
f = zeros(ng,nch,nd);
f_udg = zeros(ng,nch,nd,nc);
for i=1:Ns
    % subset of U variables describing this snapshot
    uvars = (i-1)*nch0+1:i*nch0; 
    
    % subset of UDG variables describing this snapshot
    if nd1==1
        udgvars = uvars;
    else
        udgvars = zeros(nch0,nd1);
        for j=1:nch0; udgvars(j,:) = (i-1)*nch0+j:nch:nc; end
        udgvars= sort(udgvars(:))';
    end
    
    % time coordinate of this snapshot    
    timei = time + t(i)*dt;
    
    % Assemble fluxes
    [f(:,uvars,:),f_udg(:,uvars,:,udgvars)] = flux(p,udg(:,udgvars),param{1},timei);
end

