function [fb,fb_udg,fb_uh] = fbou3d(ib,uinf,nl,pg,udg,uh,param,time)
%FBOU3D.M
%    [FB,FB_UDG,FB_UH] = FBOU3D(IB,UINF,NL,PG,UDG,UH,PARAM,TIME)
[ng,nc] = size(udg);
nch = 5;
nd = 3;
gam = param{1};
gam1 = gam-1.0;
Minf = param{5};
Pr   = param{4};
M2   = Minf^2;
tau  = param{end};

u = udg(:,1:nch);
switch (ib)
    case 1          
%         fb = (uh-uinf);
%         fb_udg = zeros(ng,nch,nc);
%         fb_uh = zeros(ng,nch,nch);
%         for i=1:nch
%             fb_uh(:,i,i) = 1;
%         end
%     case 2
%         [fb,fb_udg,fb_uh] = fbat(nl,pg,udg,uh,param,time);
%         fb = fb - uinf;
%     case 1  % Far field
% % An      = zeros(ng,nc,nc);
% % Anm     = zeros(ng,nc,nc,nc);        
        [an,anm] = getan(nl,uh,param,0);
        [An,Anm] = getan(nl,uh,param,1);
        fb = permute(mapContractK(an+An,u-uh,2,3,1,2,[],1)-mapContractK(an-An,uinf-uh,2,3,1,2,[],1),[2 1]);
        fb_u = an+An;
        fb_q = zeros(ng,nch,nch,nd);
        fb_uh = permute(mapContractK(anm+Anm,u-uh,[2 4],3,1,2,[],1)-mapContractK(anm-Anm,uinf-uh,[2 4],3,1,2,[],1),[3 1 2])-2*An;        
        fb_udg = cat(3,fb_u,reshape(fb_q,ng,nch,nd*nch));
    case 2  % Adiabatic Wall     
        uinf = udg(:,1:nch);
        uinf(:,2:nd+1) = 0;
        
        [fbn,fbn_udg,fbn_uh] = fhat(nl,pg,udg,uh,param,time);
        
        fb = uinf - uh;
        fb(:,nch) = fbn(:,nch);
        
        fb_u = zeros(ng,nch,nch);
        fb_u(:,1,1) = ones(ng,1);        
        fb_u(:,nch,:) = fbn_udg(:,nch,1:nch);
   
        fb_q = zeros(ng,nch,nch*nd);
        fb_q(:,nch,:) = fbn_udg(:,nch,nch+1:(nd+1)*nch);
        
        fb_uh = zeros(ng,nch,nch);
        fb_uh(:,1,1) = -1;
        fb_uh(:,2,2) = -1;
        fb_uh(:,3,3) = -1;           
        if nd==3
            fb_uh(:,4,4) = -1;           
        end
        fb_uh(:,nch,:) = fbn_uh(:,nch,:);
        
        fb_udg = cat(3,fb_u,reshape(fb_q,ng,nch,nd*nch));        
    case 3 % Isothermal wall
        nd = 3;
        zero = zeros(ng,1);
        fb(:,1)   = udg(:,1)-uh(:,1);        
        fb(:,2:4) = -uh(:,2:4);
        fb(:,5) = gam*gam1*M2*uh(:,5)./uh(:,1) - uinf(:,end);
        
        fb_u = zeros(ng,nch,nch);
        fb_u(:,1,1) = ones(ng,1);  
        
        fb_q = zeros(ng,nch,nch,nd);        
        
        fb_uh = zeros(ng,nch,nch);
        fb_uh(:,1,1) = -1;
        fb_uh(:,2,2) = -1;
        fb_uh(:,3,3) = -1;
        fb_uh(:,4,4) = -1;
        fb_uh(:,5,:) = gam*gam1*M2*[-uh(:,5)./uh(:,1).^2, zero, zero, zero, 1./uh(:,1)];
        
        fb_udg = cat(3,fb_u,reshape(fb_q,ng,nch,nd*nch));    
    case 5 % periodic
        [fb,fb_udg,fb_uh] = fhat(nl,pg,udg,uh,param,time);         
    case 8        
        y = pg(:,2);                
        T0 = 0.75;
        T1 = 0.85;
        Ti = T0 + y*(T1-T0) + ((gam1*gam*M2*Pr)/(2*gam))*y.*(1-y);        
        fb(:,1) = 1./Ti-uh(:,1);        
        fb(:,2) = y./Ti-uh(:,2);
        fb(:,3) = -uh(:,3);
        fb(:,4) = -uh(:,4);
        
        r = uh(:,1);
        ru = uh(:,2);
        rv = uh(:,3);
        rw = uh(:,4);
        rE = uh(:,5);
        t = gam*(gam-1)*M2*(rE-0.5*(ru.*ru./r+rv.*rv./r+rw.*rw./r))./r;
        t_r = (M2*gam*(gam - 1)*(ru.^2 + rv.^2 + rw.^2 - r.*rE))./r.^3;
        t_ru = -(M2*gam*ru*(gam - 1))./r.^2;
        t_rv = -(M2*gam*rv*(gam - 1))./r.^2;
        t_rw = -(M2*gam*rw*(gam - 1))./r.^2;
        t_rE = (M2*gam*(gam - 1))./r;        
        fb(:,5) = Ti-t;                
        
        fb_u = zeros(ng,nch,nch);
        %fb_u(:,1,1) = ones(ng,1);  
        
        fb_q = zeros(ng,nch,nch,3);        
        
        fb_uh = zeros(ng,nch,nch);
        fb_uh(:,1,1) = -1;
        fb_uh(:,2,2) = -1;
        fb_uh(:,3,3) = -1;        
        fb_uh(:,4,4) = -1;        
        fb_uh(:,5,:) = -[t_r, t_ru, t_rv, t_rw, t_rE];                
        
        fb_udg = cat(3,fb_u,reshape(fb_q,ng,nch,3*nch));        
    otherwise
         error('unknown boundary type');
end
