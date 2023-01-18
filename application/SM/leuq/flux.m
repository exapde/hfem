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

if nc==6
    % 2D plane strain neo-Hookean materials
    [f,f_udg] = flux2d(pg,udg,param,time);
else
    %[f1,f1_udg] = flux3d(pg,udg,param,time);
    
    mu = param{1};
    kappa = param{2};
    ncu = 3;
    nch = 3;
    nd = 3;
    f = zeros(ng,ncu*nd);
    
    q11 = udg(:,4);
    q21 = udg(:,5);
    q31 = udg(:,6);
    q12 = udg(:,7);
    q22 = udg(:,8);
    q32 = udg(:,9);
    q13 = udg(:,10);
    q23 = udg(:,11);
    q33 = udg(:,12);
    
    tq = (q11 + q22 + q33);
    f(:,1) = 2*mu*q11 + kappa*tq;
    f(:,2) = mu*(q12 + q21);
    f(:,3) = mu*(q13 + q31);
    f(:,4) = mu*(q12 + q21);
    f(:,5) = 2*mu*q22 + kappa*tq;
    f(:,6) = mu*(q23 + q32);
    f(:,7) = mu*(q13 + q31);
    f(:,8) = mu*(q23 + q32);
    f(:,9) = 2*mu*q33 + kappa*tq;
    
    f_udg = zeros(ng,ncu*nd,nc);
    f_udg(:,1,4) = kappa+2*mu;
    f_udg(:,1,8) = kappa;
    f_udg(:,1,12) = kappa;    
    f_udg(:,2,5) = mu;
    f_udg(:,2,7) = mu;        
    f_udg(:,3,6) = mu;        
    f_udg(:,3,10) = mu;        
    f_udg(:,4,5) = mu;
    f_udg(:,4,7) = mu;    
    f_udg(:,5,4) = kappa;
    f_udg(:,5,8) = kappa+2*mu;
    f_udg(:,5,12) = kappa;    
    f_udg(:,6,9) = mu;
    f_udg(:,6,11) = mu;        
    f_udg(:,7,6) = mu;        
    f_udg(:,7,10) = mu;        
    f_udg(:,8,9) = mu;
    f_udg(:,8,11) = mu;    
    f_udg(:,9,4) = kappa;
    f_udg(:,9,8) = kappa;
    f_udg(:,9,12) = kappa+2*mu;
    
    f = reshape(f,ng,nch,nd);
    f_udg = reshape(f_udg,ng,nch,nd,nc);
    
%     max(abs(f(:)-f1(:)))
%     max(abs(f_udg(:)-f1_udg(:)))
end


return;

