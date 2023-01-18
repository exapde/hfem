function [RU,RH,RQ] = hdg_residual(master,app,mesh,UDG,UH,SH)

npv = size(UDG,1);                
ne  = size(mesh.dgnodes,3);
npf = size(master.perm,1);
nfe = size(master.perm,2);
ncu = app.ncu;
nc  = app.nc;
ns  = 500;

nb = ceil(ne/ns);          
nk = 1:ns:ne;
nm = [nk(1:end); [nk(2:end)-1,ne]];    

UDG = permute(UDG,[1 3 2]);
SH  = permute(SH,[1 3 2]);
mesh.dgnodes = permute(mesh.dgnodes,[1 3 2]);

UH = reshape(UH(:,mesh.elcon),[ncu npf*nfe ne]);
UH = permute(UH,[2 3 1]);

RU = zeros(npv,ncu,ne);
RQ = zeros(npv,nc-ncu,ne);
RH = zeros(ncu,nfe*npf,ne);
for j=1:nb
    id = nm(1,j):nm(2,j);    
    
    % compute volume integrals    
    [Ru, Rq] = volint(master,app,mesh.dgnodes(:,id,:),UDG(:,id,:),SH(:,id,:));          
    
    %RU(:,:,id) = permute(Ru,[1 3 2]);
    
    % compute face integrals        
    [Ru, Rq, Rh] = faceint(master,app,mesh.bf(:,id),mesh.dgnodes(:,id,:),UDG(:,id,:),UH(:,id,:),Ru,Rq);                
        
    RU(:,:,id) = Ru;
    RQ(:,:,id) = Rq;
    RH(:,:,id) = Rh;        
end    

nb = npf*ncu;
% assemble the residual
il = zeros(nfe*npf*ncu,nfe*npf*ncu,ne);
for i=1:ne
    con = repmat((mesh.elcon(:,i)'-1)*ncu,ncu,1)+repmat((1:ncu)',1,nfe*npf);
    con = reshape(con,nfe*npf*ncu,1);    
    il(:,:,i) = repmat(con ,1,nfe*nb);       
end

RH = sparse(reshape(il(:,1,:),nfe*nb*ne,1),ones(nfe*nb*ne,1),reshape(RH,nfe*nb*ne,1));                                                                      


function [Ru, Rq] = volint(master,app,dgnodes,UDG,SH)
% VOLINTND compute volume integrals 

% Obtain dgnodes, Jacobian matrix and determinant at Gauss points
[pg, Xx, jac] = volgeom(master.shapmv,dgnodes);

ne   = size(UDG,2);
nd   = master.nd;
npv  = master.npv;
ngv  = master.ngv;

nc   = app.nc;
ncu  = app.ncu;
arg  = app.arg;
time = app.time;
tdep = app.tdep;
fc_q = app.fc_q;
fc_u = app.fc_u;
localsolve = app.localsolve;
source = str2func(app.source);
flux   = str2func(app.flux);

% Shap functions and derivatives
shapvt = master.shapvt;
shapvg = reshape(master.shapvg,[npv ngv*(nd+1)]);
shapvgdotshapvl = reshape(master.shapvgdotshapvl,[npv*npv ngv (nd+1)]);

% DG solution at Gauss points
udgg = reshape(UDG,[npv ne*nc]);
udgg = shapvt(:,:,1)*udgg;
udgg = reshape(udgg,[ngv*ne nc]);

if localsolve == 0    
    Rq = [];
else
    % Mass matrix
    M = fc_q*reshape(shapvgdotshapvl(:,:,1)*reshape(jac,[ngv ne]),[npv npv ne]);
    
    % mass inverse times convection matrices
    C = zeros(npv,npv,nd,ne);
    for i=1:nd
        tmp = reshape(shapvgdotshapvl(:,:,2)*Xx(:,:,i,1),[npv npv ne]);
        for j=2:nd
            tmp = tmp + reshape(shapvgdotshapvl(:,:,j+1)*Xx(:,:,i,j),[npv npv ne]);
        end    
        for k=1:ne
            C(:,:,i,k) = tmp(:,:,k);        
        end
    end        
    
    Rq = reshape(mapContractK(M,SH(:,:,ncu+1:end)/fc_q-UDG(:,:,ncu+1:end),1,2,3,1,3,2),[npv ncu nd ne]);     
    for i=1:nd
        Rq(:,:,i,:) = Rq(:,:,i,:)+reshape(mapContractK(C(:,:,i,:),UDG(:,:,1:ncu),1,[2 3],4,1,3,2),[npv ncu 1 ne]);
    end           
end

% Fluxes and source at Gauss points
f = flux( pg, udgg, arg, time);
s = source( pg, udgg, arg, time); 
f     = reshape(f,[ngv ne ncu nd]);
s     = reshape(s(:,1:ncu),[ngv*ne ncu]);

% Update source term for time-dependent problems
if tdep    
    Stn = reshape(SH(:,:,1:ncu),[npv ne*ncu]);
    Stg = shapvt(:,:,1)*Stn;
    Stg = reshape(Stg,[ngv*ne ncu]);

    s = s + Stg - udgg(:,1:ncu)*fc_u;    
end

% compute wrk and wrl to time with shape functions
wrk = zeros(ngv*(nd+1),ne*ncu);
wrk(1:ngv,:) =  reshape(bsxfun(@times,s,jac),[ngv ne*ncu]);
for i=1:nd
    fk = bsxfun(@times,f(:,:,:,1),Xx(:,:,1,i));    
    for j=2:nd
        fk = fk + bsxfun(@times,f(:,:,:,j),Xx(:,:,j,i));        
    end
    wrk(i*ngv+1:(i+1)*ngv,:) = reshape(fk,[ngv ne*ncu]);   
end

% Volume residual
% [Phi Phi_xi Phi_eta] x [S.*jac; Fx.*Xx(:,:,1,1)+Fy.*Xx(:,:,2,1); Fx.*Xx(:,:,1,2)+Fy.*Xx(:,:,2,2)]
Ru = shapvg*wrk; % [npv ngv*(nd+1)] x [ngv*(nd+1) ne*ncu] 
Ru = reshape(Ru,[npv ne ncu]); 


function [Ru, Rq, Rh] = faceint(master,app,bf,dgnodes,UDG,UH,Ru,Rq)
% FACEINTND compute face integrals 

% Obtain dgnodes, Jacobian matrix and determinant at Gauss points
[pg, nlg, jac] = facegeom(master.shapmf,dgnodes,master.permgeom);

ne   = size(UDG,2);
nd   = master.nd;
npv  = master.npv;
npf  = master.npf;
ngf  = master.ngf;
nfe  = size(master.perm,2);
nc   = app.nc;
ncu  = app.ncu;
nch  = app.nch;
arg  = app.arg;
time = app.time;
bcm  = app.bcm;        
bcs  = app.bcs;
localsolve = app.localsolve;
fbou   = str2func(app.fbou);
fhat   = str2func(app.fhat);

% Shap functions 
perm            = master.perm(:,:,1);
shapft          = master.shapft(:,:,1);
shapfg          = master.shapfg(:,:,1);
shapfgdotshapfc = master.shapfgdotshapfc(:,:,1);

if localsolve==0    
else    
    tmp = zeros(npv,npf*nfe,ne);
    E = zeros(npv,npf*nfe,nd,ne);
    for i=1:nd    
        njc = reshape(nlg(:,i).*(jac),[ngf,nfe*ne]);
        wrk = reshape(shapfgdotshapfc*njc,[npf npf*nfe ne]);   
        for j=1:nfe
            tmp(perm(:,j),(j-1)*npf+1:j*npf,:) = wrk(1:npf,(j-1)*npf+1:j*npf,:);      
        end    
        for k=1:ne
            E(:,:,i,k) = tmp(:,:,k);        
        end
    end
    for i = 1:nd
        Rq(:,:,i,:) = Rq(:,:,i,:) - reshape(mapContractK(E(:,:,i,:),UH,1,[2 3],4,1,3,2),[npv ncu 1 ne]);
    end
end
    
% DG solution at Gauss points
udgn = reshape(UDG(perm,:,:),[npf nfe*ne*nc]);
udgg = shapft*udgn;
udgg = reshape(udgg,[ngf*nfe*ne nc]);

uh  = reshape(UH,[npf nfe*ne*nch]);
uhg = shapft*uh;
uhg = reshape(uhg,[ngf*nfe*ne nch]);

% Compute numerical fluxes 
FH = fhat(nlg, pg, udgg, uhg, arg, time);     

tm     = find(bf<0);
nbf    = length(tm(:));
BH     = zeros(ngf*nbf,nch);
an     = zeros(nbf*npf,1);
bn     = zeros(nbf*ngf,1);
j = 1;
for i = 1:length(bcm)    
    [I,J] = find(bf==-i);
    nfb = length(I);
    pgb  = zeros(ngf*nfb, size(pg,2));
    nlgb = zeros(ngf*nfb, nd);
    udgb = zeros(ngf*nfb, nc);
    uhgb = zeros(ngf*nfb, ncu);
    im = zeros(ngf*nfb,1);
    for k = 1:nfb
        i0 = (j-1)*ngf+1:j*ngf;
        i1 = (k-1)*ngf+1:k*ngf; 
        i2 = (J(k)-1)*nfe*ngf+(I(k)-1)*ngf+1:(J(k)-1)*nfe*ngf+I(k)*ngf; 
        bn(i0) = i2;
        im(i1)  = i0;
        pgb(i1,:)  = pg(i2,:);
        nlgb(i1,:) = nlg(i2,:);
        udgb(i1,:) = udgg(i2,:);
        uhgb(i1,:) = uhg(i2,:);
        i0 = (j-1)*npf+1:j*npf;
        i2 = (J(k)-1)*nfe*npf+(I(k)-1)*npf+1:(J(k)-1)*nfe*npf+I(k)*npf; 
        an(i0) = i2;
        j = j + 1;
    end  
    
    ib = bcm(i);
    bv = repmat(bcs(i,:),[ngf*nfb 1]);
    fh = fbou(ib,bv,nlgb,pgb,udgb,uhgb,arg,time);                                                     
    BH(im,:) = fh;
end

wrk = reshape(bsxfun(@times,FH,jac),[ngf nfe*ne*nch]);
Run = reshape(shapfg*wrk,[npf nfe ne nch]);   
Rut = zeros(npv,1,ne,ncu);
for is=1:nfe  % Adding face contributions - vector dependencies avoided
    IP = perm(:,is);
    Rut(IP,1,:,:) = Rut(IP,1,:,:) + Run(:,is,:,:);
end
Ru = permute(Ru-reshape(Rut,[npv ne ncu]),[1 3 2]);

wrk = reshape(bsxfun(@times,BH,jac(bn)),[ngf nbf*nch]);
Run = reshape(Run,[npf*nfe*ne nch]);
Run(an,:) = reshape(shapfg*wrk,[npf*nbf nch]);

Rh  = permute(reshape(Run,[npf*nfe ne nch]), [3 1 2]);

