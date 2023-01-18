function [fh,fh_udg,fh_uh] = fboucorona(ib,ui,nl,p,udg,uh,param,time)
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
nch = size(uh,2);
nq = nc-nch;

T1star = param{6};
tau   = param{end};
u     = udg(:,1:nch);

switch ib
    case 1  % Dirichlet     
        if isempty(time)
            time = T1star;
        end        
        ui(:,1) = min((time/T1star),1.0)*ui(:,1);
        fh = tau.*(ui(:,1:nch)-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end
    case 2  % "Extrapolate  m = u 
        fh = tau*(u-uh);
        fh_u = zeros(ng,nch,nch);        
        fh_q = zeros(ng,nch,nq);
        fh_uh = zeros(ng,nch,nch);
        for i = 1:nch
            fh_u(:,i,i) = tau;
            fh_uh(:,i,i) = -tau;
        end
        fh_udg = cat(3,fh_u,fh_q);        
    case 3  % dirichlet and extrapolation        
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = tau*(ui(:,1)-uh(:,1));
        fh(:,2) = tau*(u(:,2)-uh(:,2));   
        fh_udg(:,2,2) = tau;        
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end
    case 4 % fluxes
        [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time);
        fh = fh + ui;    
    case 5 % symmetry
        nx = nl(:,1);
        ny = nl(:,2);
        phi = udg(:,1);
        rho = udg(:,2);        
        phix = udg(:,3);
        rhox = udg(:,4);
        phiy = udg(:,5);
        rhoy = udg(:,6);
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = phix.*nx + phiy.*ny + tau*(phi - uh(:,1));
        fh(:,2) = rhox.*nx + rhoy.*ny + tau*(rho - uh(:,2));        
        fh_udg(:,1,1) = tau;
        fh_udg(:,1,3) = nx;
        fh_udg(:,1,5) = ny;
        fh_udg(:,2,2) = tau;
        fh_udg(:,2,4) = nx;
        fh_udg(:,2,6) = ny;        
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end
    case 6 % plate
        nx = nl(:,1);
        ny = nl(:,2);        
        rho = udg(:,2);                
        rhox = udg(:,4);        
        rhoy = udg(:,6);
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = tau*(ui(:,1) - uh(:,1));
        fh(:,2) = rhox.*nx + rhoy.*ny + tau*(rho - uh(:,2));        
        fh_udg(:,2,2) = tau;
        fh_udg(:,2,4) = nx;
        fh_udg(:,2,6) = ny;        
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end        
    case 7    
        beta1 = param{8};
        beta2 = param{9};
        nx = nl(:,1);
        ny = nl(:,2);
        phi = udg(:,1);
        rho = udg(:,2);        
        phix = udg(:,3);
        rhox = udg(:,4);
        phiy = udg(:,5);
        rhoy = udg(:,6);
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = phix.*nx + phiy.*ny + beta1*uh(:,1) + tau*(phi - uh(:,1));
        fh(:,2) = rhox.*nx + rhoy.*ny + beta2*uh(:,2) + tau*(rho - uh(:,2));        
        fh_udg(:,1,1) = tau;
        fh_udg(:,1,3) = nx;
        fh_udg(:,1,5) = ny;
        fh_udg(:,2,2) = tau;
        fh_udg(:,2,4) = nx;
        fh_udg(:,2,6) = ny;        
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end        
        fh_uh(:,1,1) = fh_uh(:,1,1) + beta1;
        fh_uh(:,2,2) = fh_uh(:,2,2) + beta2;
    case 8    
        beta1 = param{8};        
        nx = nl(:,1);
        ny = nl(:,2);
        phi = udg(:,1);               
        phix = udg(:,3);        
        phiy = udg(:,5);        
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = phix.*nx + phiy.*ny + beta1*uh(:,1) + tau*(phi - uh(:,1));        
        fh(:,2) = tau*(u(:,2)-uh(:,2));   
        fh_udg(:,1,1) = tau;
        fh_udg(:,1,3) = nx;
        fh_udg(:,1,5) = ny;    
        fh_udg(:,2,2) = tau;  
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end        
        fh_uh(:,1,1) = fh_uh(:,1,1) + beta1;      
    case 9    
        beta1 = param{8};        
        nx = nl(:,1);
        ny = nl(:,2);
        phi = udg(:,1);               
        phix = udg(:,3);        
        phiy = udg(:,5);        
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = phix.*nx + phiy.*ny + beta1*uh(:,1) + tau*(phi - uh(:,1));        
        fh(:,2) = tau*(ui(:,2)-uh(:,2));   
        fh_udg(:,1,1) = tau;
        fh_udg(:,1,3) = nx;
        fh_udg(:,1,5) = ny;            
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end        
        fh_uh(:,1,1) = fh_uh(:,1,1) + beta1;                
    otherwise
        error('unknown boundary type');
end

