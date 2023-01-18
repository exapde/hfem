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

u = udg(:,1:nch);
p = udg(:,nch+1);
q = udg(:,nch+2:end);

f = kappa*reshape(q,[ng nch nd]);
for d=1:nd
    f(:,d,d) = f(:,d,d) + p;    
    for j=1:nd
        f(:,j,d) = f(:,j,d) + u(:,d).*u(:,j);    
    end
end

f_udg = zeros(ng,nch,nd,nc);
for d=1:nd        
    f_udg(:,d,d,nch+1) = 1;
    for j=1:nch
        k = (nch+1) + (d-1)*nch+j;
        f_udg(:,j,d,k) = kappa;        
    end
    for j=1:nd
        f_udg(:,j,d,j) = u(:,d);           
        f_udg(:,j,d,d) = f_udg(:,j,d,d) + u(:,j);    
    end
end

% f_u = zeros(ng,nch,2,nch);
% f_u(:,:,1,1) = [2*u1,        u2];
% f_u(:,:,1,2) = [zeros(ng,1), u1];
% f_u(:,:,2,1) = [u2, zeros(ng,1)];
% f_u(:,:,2,2) = [u1,        2*u2];

% u1   = udg(:,1);
% u2   = udg(:,2);
% pr   = udg(:,3);
% q11  = udg(:,4);
% q12  = udg(:,5);
% q21  = udg(:,6);
% q22  = udg(:,7);
% 
% f = zeros(ng,nch,2);
% f(:,:,1) = [u1.*u1+kappa*q11+pr, u1.*u2+kappa*q12];
% f(:,:,2) = [u1.*u2+kappa*q21, u2.*u2+kappa*q22+pr];
% 
% f_p = zeros(ng,nch,2);
% f_p(:,:,1) = [ones(ng,1),  zeros(ng,1)];
% f_p(:,:,2) = [zeros(ng,1),  ones(ng,1)];
% 
% f_u = zeros(ng,nch,2,nch);
% f_u(:,:,1,1) = [2*u1,        u2];
% f_u(:,:,1,2) = [zeros(ng,1), u1];
% f_u(:,:,2,1) = [u2, zeros(ng,1)];
% f_u(:,:,2,2) = [u1,        2*u2];
% 
% f_q = zeros(ng,nch,2,nch,2);
% f_q(:,:,1,1,1) = [kappa*ones(ng,1)  zeros(ng,1)];
% f_q(:,:,1,2,1) = [zeros(ng,1)  kappa*ones(ng,1)];
% f_q(:,:,1,1,2) = [zeros(ng,1)       zeros(ng,1)];
% f_q(:,:,1,2,2) = [zeros(ng,1)       zeros(ng,1)];
% f_q(:,:,2,1,1) = [zeros(ng,1)       zeros(ng,1)];
% f_q(:,:,2,2,1) = [zeros(ng,1)       zeros(ng,1)];
% f_q(:,:,2,1,2) = [kappa*ones(ng,1)  zeros(ng,1)];
% f_q(:,:,2,2,2) = [zeros(ng,1)  kappa*ones(ng,1)];
% 
% f_q =reshape(f_q,ng,nch,2,2*nch);
%  
% f_udg = cat(4,cat(4,f_u,f_p),f_q);
% 
