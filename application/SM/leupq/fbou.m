function [fh,fh_udg,fh_uh] = fbou(ib,ui,nl,pg,udg,uh,param,time)
%FBOU boundary flux function

%      IC                    Boundary number
%      IB                    Boundary type
%      FUNC                  Function to compute boundary data
%      NL(N,ND)              Normal N points
%      X(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      P(N,1)                Pressure vector for N points with 1 component
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHP(N,NC,1):          Jacobian of the flux flux vector w.r.t. P
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. M

[ngf,nch] = size(uh);
nc = size(udg,2);

%kappa = param{1};
tau   = param{end};
ginf  = repmat(ui(1,1:nch),ngf,1);
switch (ib)
    case 1 % Dirichlet bc        
        fh = tau*(ginf-uh);
        fh_udg = zeros(ngf,nch,nc);
        fh_uh = -bsxfun(@times,ones(ngf,1,1),reshape(tau*eye(nch),[1 nch nch]));
    case 2 % stress bc        
        [fh,fh_udg,fh_uh] = fhatbou1(nl,pg,udg,uh,param,time);            
        fh = fh+ginf;        
    case 3 % gradient bc                
        [fh,fh_udg,fh_uh] = fhat(nl,pg,udg,uh,param,time);        
        fh = fh+ginf;
    case 4 % mixed Dirichlet and stress bcs
        uinf = repmat(ui(1,1),ngf,1);
        finf = repmat(ui(1,2),ngf,1);
        
        fh   = zeros(ngf,nch);        
        fh_uh = zeros(ngf,nch,nch);
        fh_udg = zeros(ngf,nch,nc);
        
        fh(:,1) = tau*(uinf - uh(:,1));        
        fh_uh(:,1,1) = -tau;
        
        [gh,gh_udg,gh_uh] = fhatbou1(nl,pg,udg,uh,param,time);
        fh(:,2) = gh(:,2) + finf;
        fh_uh(:,2,:)  = gh_uh(:,2,:);
        fh_udg(:,2,:) = gh_udg(:,2,:);        
    case 5 % mixed Dirichlet and stress bcs
        uinf = repmat(ui(1,1),ngf,1);
        finf = repmat(ui(1,2),ngf,1);
        
        fh   = zeros(ngf,nch);        
        fh_uh = zeros(ngf,nch,nch);
        fh_udg = zeros(ngf,nch,nc);
        
        fh(:,2) = tau*(uinf - uh(:,2));        
        fh_uh(:,2,2) = -tau;
        
        [gh,gh_udg,gh_uh] = fhabout1(nl,pg,udg,uh,param,time);
        fh(:,1) = gh(:,1) + finf;
        fh_uh(:,1,:)  = gh_uh(:,1,:);
        fh_udg(:,1,:) = gh_udg(:,1,:);     
    otherwise
        error('unknown boundary type');
end


% total stress bc
function [fh,fh_udg,fh_uh] = fhatbou1(nl,pg,udg,uh,param,time)

[ngf,nch] = size(uh);
nc = size(udg,2);
nd = size(pg,2);
mu  = param{1};
tau = param{end};

lambda  = param{2};
cp = lambda/(mu+lambda);
cp = 1;
%[mu lambda tau]

p = udg(:,nch+1);
q = udg(:,nch+2:end);
q = reshape(q,[ngf nch nd]);

% stress tensor
f = mu*(q + permute(q,[1 3 2]));
for d=1:nd
    f(:,d,d) = f(:,d,d) + cp*p;    
end

f_udg = zeros(ngf,nch,nd,nc);
h_udg = zeros(ngf,nch,nd,nc);
for d=1:nd
    f_udg(:,d,d,nch+1) = cp;
    for j=1:nch
        k = (nch+1) + (d-1)*nch+j;
        f_udg(:,j,d,k) = mu;        
        h_udg(:,d,j,k) = mu;        
    end
end
f_udg = f_udg + h_udg;

fh = tau*(udg(:,1:nch) - uh);
for ic=1:nch
    fh(:,ic) = fh(:,ic) + sum(reshape(f(:,ic,:),[ngf nd]).*nl,2);
end

fh_udg = zeros(ngf,nch,nc);
fh_uh = zeros(ngf,nch,nch);    
for d=1:nd
    fh_udg(:,d,d) = tau;
    fh_uh(:,d,d) = -tau;
    for ic=nch+1:nc
        fh_udg(:,d,ic) = sum(reshape(f_udg(:,d,:,ic),[ngf nd]).*nl,2);    
    end
end

function [fh,fh_udg,fh_uh] = fhatbou2(nl,p,udg,uh,param,time)


[ng1,nch] = size(uh);
ncs = 4;

mu  = param{1};
tau = mu;
lambda  = param{2};
cp = lambda/(mu+lambda);

fh(:,1) = mu*(2*udg(:,nch+2).*nl(:,1) + (udg(:,nch+3)+udg(:,nch+4)).*nl(:,2)) + cp*udg(:,nch+1).*nl(:,1) + tau*(udg(:,1)-uh(:,1));
fh(:,2) = mu*((udg(:,nch+3)+udg(:,nch+4)).*nl(:,1) + 2*udg(:,nch+5).*nl(:,2)) + cp*udg(:,nch+1).*nl(:,2) + tau*(udg(:,2)-uh(:,2));

if nargout > 1
    fh_u  = zeros(ng1,nch,nch);
    fh_p  = zeros(ng1,nch,1);
    fh_q  = zeros(ng1,nch,ncs);
    fh_uh = zeros(ng1,nch,nch);
    
    fh_u(:,1,1) = tau;
    fh_u(:,2,2) = tau;
    
    fh_p(:,1,1) = cp*nl(:,1);
    fh_p(:,2,1) = cp*nl(:,2);
    
    fh_q(:,1,1) = 2*mu*nl(:,1);
    fh_q(:,1,2) = mu*nl(:,2);
    fh_q(:,1,3) = mu*nl(:,2);
    fh_q(:,2,2) = mu*nl(:,1);
    fh_q(:,2,3) = mu*nl(:,1);
    fh_q(:,2,4) = 2*mu*nl(:,2);    
    
    fh_udg = cat(3,cat(3,fh_u,fh_p),fh_q);
        
    fh_uh(:,1,1) = -tau;
    fh_uh(:,2,2) = -tau;
end

function [fh,fh_udg,fh_uh] = fhatbou3(nl,p,udg,uh,param,time)


[ng1,nch] = size(uh);
ncs = 4;

mu  = param{1};
lambda = param{2};
tau = param{end};
eta = lambda/(mu+lambda);

fh(:,1) = mu*(2*udg(:,nch+2).*nl(:,1) + (udg(:,nch+3)+udg(:,nch+4)).*nl(:,2)) + eta*udg(:,nch+1).*nl(:,1) + tau*(udg(:,1)-uh(:,1));
fh(:,2) = mu*((udg(:,nch+3)+udg(:,nch+4)).*nl(:,1) + 2*udg(:,nch+5).*nl(:,2)) + eta*udg(:,nch+1).*nl(:,2) + tau*(udg(:,2)-uh(:,2));

if nargout > 1
    fh_u  = zeros(ng1,nch,nch);
    fh_p  = zeros(ng1,nch,1);
    fh_q  = zeros(ng1,nch,ncs);
    fh_uh = zeros(ng1,nch,nch);
    
    fh_u(:,1,1) = tau;
    fh_u(:,2,2) = tau;
    
    fh_p(:,1,1) = eta*nl(:,1);
    fh_p(:,2,1) = eta*nl(:,2);
    
    fh_q(:,1,1) = 2*mu*nl(:,1);
    fh_q(:,1,2) = mu*nl(:,2);
    fh_q(:,1,3) = mu*nl(:,2);
    fh_q(:,2,2) = mu*nl(:,1);
    fh_q(:,2,3) = mu*nl(:,1);
    fh_q(:,2,4) = 2*mu*nl(:,2);
        
    fh_udg = cat(3,cat(3,fh_u,fh_p),fh_q);
        
    fh_uh(:,1,1) = -tau;
    fh_uh(:,2,2) = -tau;
end


