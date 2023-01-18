function [fh,fh_udg,fh_uh] = fbou(ib,ui,nl,pg,udg,uh,param,time)
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

[ngf,nc] = size(udg);
nd = size(pg,2);
nch = nc/(nd+1);

if size(ui,1) == ngf
    uinf = ui(:,1:nch);
else
    uinf = repmat(ui(1,1:nch),ng,1);
end

u = udg(:,1:nch);
%q = reshape(udg(:,nch+1:3*nch),ng,nch,2);

c2 = param{1};
k  = param{2};
freq = param{3};
tau = param{end};

kx    = k(1);
ky    = k(2);
kk    = sqrt(kx^2+ky^2);

tau  = tau*bsxfun(@times,ones(ngf,1,1),reshape(eye(nch),[1 nch nch]));
switch ib
    case 1  % Dirichlet
        fh = tau.*(uinf-uh);
        fh_udg = zeros(ngf,nch,nc);        
        fh_uh = -bsxfun(@times,ones(ngf,1,1),reshape(c2*eye(nch),[1 nch nch]));
    case 2  % "Extrapolate  m = u         
%         fh = tau.*(u-uh);
%         fh_u = tau;        
%         fh_udg = cat(3,fh_u,zeros(ngf,nch,nd*nch));        
%         fh_uh = -bsxfun(@times,ones(ngf,1,1),reshape(c2*eye(nch),[1 nch nch]));
        
        [fh,fh_udg,fh_uh] = fhat(nl,pg,udg,uh,param,time);
        fh = fh - uh;
        fh_uh = fh_uh - 1;
        
        x = pg(:,1);
        ind = find(x>=-2.5 & x<=-2);        
        [fht,fht_udg,fht_uh] = fhat(nl(ind,:),pg(ind,:),udg(ind,:),uh(ind,:),param,time);                
        fht = fht + 8*sin(freq*2*pi*time).*exp(-pi*freq*time^2);
        fh(ind,:) = fht;
        fh_udg(ind,:,:) = fht_udg;
        fh_uh(ind,:,:) = fht_uh;        
    case 3  % Prescribed flux
        [fh,fh_udg,fh_uh] = fhat(nl,pg,udg,uh,param,time);
        u1 = (kx)*cos(kx*pg(:,1)+ky*pg(:,2) - kk*time); 
        u2 = (ky)*cos(kx*pg(:,1)+ky*pg(:,2) - kk*time); 
        finf = u1.*nl(:,1) + u2.*nl(:,2);
        fh = fh - finf;
    case 4 % 1st order abosrbing         
        [fh,fh_udg,fh_uh] = fhat(nl,pg,udg,uh,param,time);
        fh = fh - uh;
        fh_uh = fh_uh - 1;
    case 5        
        [fh,fh_udg,fh_uh] = fhat(nl,pg,udg,uh,param,time);
    otherwise
        error('unknown boundary type');
end

