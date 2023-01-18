function [fh,fh_udg,fh_uh] = fboucorona2(ib,ui,nl,p,udg,uh,param,time)
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

foc   = param{end-1};
tau   = param{end};
u     = udg(:,1:nch);

switch ib
    case 0
%         ex = p(:,5);
%         ey = p(:,6);
%         en = abs(ex.*nl(:,1) + ey.*nl(:,2));
%         e1 = en-foc(2);
%         e1(e1<0) = 0;
%         ui(:,2) = foc(1)*e1;
%         fh = tau.*(ui(:,1:nch)-uh);
%         fh_udg = zeros(ng,nch,nc);
%         fh_uh =  zeros(ng,nch,nch);
%         for i = 1:nch           
%             fh_uh(:,i,i) = -tau;
%         end                        
        ex = udg(:,3);
        ey = udg(:,5)+param{5};
        en = -(ex.*nl(:,1) + ey.*nl(:,2));
        e1 = en-foc(2);
        in1 = e1<0; 
        in2 = find(e1>=0);
        e1(in1) = 0;
        ui(:,2) = foc(1)*e1;        
       
%         [uh(:,2) ui(:,2) e1 en]
%         foc(2)
%         pause
        
        fh = tau.*(ui(:,1:nch)-uh);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end                     
        fh_udg(in2,2,3) = -tau*foc(1)*nl(in2,1);
        fh_udg(in2,2,5) = -tau*foc(1)*nl(in2,2);
    case 1  % Dirichlet     
        % compute fourier modes
%         [theta,~]=cart2pol(p(:,1,:),p(:,2,:));        
%         n = 10;
%         fom = zeros(ng,2*n+1);
%         fom(:,1) = 1;
%         for j=1:n            
%             fom(:,2*j) = cos(j*theta);
%             fom(:,2*j+1) = sin(j*theta);
%         end            
%         
%         % boundary condition for ui
%         nfm = length(foc);
%         ui(:,2) = foc(1)*fom(:,1); 
%         for j=2:nfm
%             ui(:,2) = ui(:,2) + foc(j)*fom(:,j);
%         end        
%         ui(:,2) = 0.665*(0.6*(1-p(:,1)/0.01).^2).*ui(:,2);
                
        %9.5*foc(1)
        %ui(:,2) =  49*exp(-abs(1+p(:,1)/0.01).^1.0/(0.25^2));
        x = p(:,1);
        y = p(:,2);
%         i1 = find(y<0);
%         g = zeros(size(x));
%         g(i1) = 49*exp(-abs(1+x(i1)/0.01).^1.0/(0.23^2));
%         i2 = find(y>=0);
%         g(i2) = 59*exp(-abs(1+x(i2)/0.01).^1.0/(0.23^2));
%         ui(:,2) =  g;
        
        t = cart2pol(x,y);
        i1 = find(t<0);
        t(i1) = 2*pi+t(i1);
        %ui(:,2) = 50.12*exp(-(t-0.962*pi).^2/(0.35^2));
%         ui(:,2) = 49.825*exp(-(t-0.969*pi).^2/(0.35^2));
%         ui(:,2) = 51.175*exp(-(t-0.969*pi).^2/(0.35^2));
%         ui(:,2) = 446*exp(-(t-0.973*pi).^2/(0.2495^2));
        ui(:,2) = foc(1)*exp(-(t-foc(2)).^2/(foc(3)^2));
%         a = foc(2); b = foc(3);
%         ui(:,2) = foc(1)*(exp(-(t-a).^2/(b^2))+ exp(-(t-2*pi).^2/(b^2)));
%         i1 = find(t>=1.5*pi);    
%         t(i1) = t(i1)-2*pi;
%         ui(:,2) = foc(1)*exp(-(t-foc(2)).^2/(foc(3)^2));
%         i1 = find(t<pi/2);    
%         t(i1) = t(i1)+2*pi;
%         ui(:,2) = foc(1)*exp(-(t-foc(2)).^2/(foc(3)^2));
        
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
%         vx = p(:,2+1);
%         vy = p(:,2+2)+1;
%         nx = vx./sqrt(vx.^2+vy.^2);
%         ny = vy./sqrt(vx.^2+vy.^2);        
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
    case 10 % plate
        fh =  zeros(ng,nch);
        fh_udg = zeros(ng,nch,nc);
        fh_uh =  zeros(ng,nch,nch);
        fh(:,1) = tau*(ui(:,1) - uh(:,1));
        fh(:,2) = tau*(ui(:,2) - uh(:,2));
        for i = 1:nch           
            fh_uh(:,i,i) = -tau;
        end                
    otherwise
        error('unknown boundary type');
end

