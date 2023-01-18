function [f,f_udg] = flux(pg,udg,param,time)
%FLUX Volume flux function
%   [f,fu,fq,fp] = flux(x,u,q,p,param)
%
%      X(N,ND)              Coordinates for N points
%      U(N,NC)              Unknown vector for N points with NC components
%      Q(N,NC,ND)           Flux vector for N points with NC components in the
%                           coordinate directions
%      P(N,1)               Pressure vector for N points with 1 component
%      PARAM                Parameter list
%      F(N,NC,ND):          Volume flux at N points
%      FU(N,NC,ND,NC):      Jacobian of the flux flux vector w.r.t. U
%      FQ(N,NC,ND,NC,ND):   Jacobian of the flux flux vector w.r.t. Q
%      FP(N,NC,ND):         Jacobian of the flux flux vector w.r.t. P
% 

[ng,nc] = size(udg);
nd = size(pg,2);
nch = (nc-1)/(nd+1);
 
kappa = param{1};

%u = udg(:,1:nch);
p = udg(:,nch+1);
q = udg(:,nch+2:end);

f = kappa*reshape(q,[ng nch nd]);
for d=1:nd
    f(:,d,d) = f(:,d,d) + p;    
end

f_udg = zeros(ng,nch,nd,nc);
for d=1:nd        
    f_udg(:,d,d,nch+1) = 1;
    for j=1:nch
        k = (nch+1) + (d-1)*nch+j;
        f_udg(:,j,d,k) = kappa;        
    end    
end
