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

[ng,nc] = size(udg);
nch = size(pg,2);
ncs = nc - nch - 1;

mu     = param{1};
lambda = 1;

f   = zeros(ng,nch,nch);
f_u = zeros(ng,nch,nch,nch);
f_p = zeros(ng,nch,nch,1);
f_s = zeros(ng,nch,nch,ncs);
if nch==2
    f(:,1,1) = mu*udg(:,nch+2) + lambda*udg(:,nch+1);
    f(:,1,2) = mu*udg(:,nch+3);
    f(:,2,1) = mu*udg(:,nch+3);
    f(:,2,2) = mu*udg(:,nch+4) + lambda*udg(:,nch+1);

    f_p(:,1,1,1) = lambda;
    f_p(:,2,2,1) = lambda;    

    f_s(:,1,1,1) = mu;
    f_s(:,1,2,2) = mu;    
    f_s(:,2,1,2) = mu;    
    f_s(:,2,2,3) = mu;    
elseif nch==3
    f(:,1,1) = mu*udg(:,nch+2) + lambda*udg(:,nch+1);
    f(:,1,2) = mu*udg(:,nch+3);
    f(:,1,3) = mu*udg(:,nch+4);
    f(:,2,1) = mu*udg(:,nch+3);
    f(:,2,2) = mu*udg(:,nch+5) + lambda*udg(:,nch+1);
    f(:,2,3) = mu*udg(:,nch+6);
    f(:,3,1) = mu*udg(:,nch+4);
    f(:,3,2) = mu*udg(:,nch+6);
    f(:,3,3) = mu*udg(:,nch+7) + lambda*udg(:,nch+1);

    f_p(:,1,1,1) = lambda;
    f_p(:,2,2,1) = lambda;   
    f_p(:,3,3,1) = lambda;   

    f_s(:,1,1,1) = mu;
    f_s(:,1,2,2) = mu;    
    f_s(:,1,3,3) = mu;    
    
    f_s(:,2,1,2) = mu;    
    f_s(:,2,2,4) = mu;        
    f_s(:,2,3,5) = mu;        
    
    f_s(:,3,1,3) = mu;    
    f_s(:,3,2,5) = mu;        
    f_s(:,3,3,6) = mu;        
else
    error('Only can handle 2D and 3D');
end

f_udg = cat(4,cat(4,f_u,f_p),f_s);


