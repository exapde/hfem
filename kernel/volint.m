function [Ru, Rq, BD, M, C, L, Q, Ju, Jq, wrl] = volint(master,app,dgnodes,UDG,SH)
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
flg_p = app.flg_p;
localsolve = app.localsolve;
adjoint = app.adjoint;
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
    M  = [];
    C = [];
    L = [];
    Q = [];
    Rq = [];
    L = shapvg(1:npv,1:ngv)*reshape(jac,[ngv ne]); % npv x ne
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
    
    if flg_p
        ncq = nc - ncu - 1;        
        [w, w_q] = press(udgg(:,ncu+2:end)); 
        L = shapvg(1:npv,1:ngv)*reshape(w.*jac,[ngv ne]); % npv x ne
        Q = shapvgdotshapvl(:,:,1)*reshape(bsxfun(@times,w_q,reshape(jac,[ngv*ne 1])),[ngv ne*ncq]);
        Q = permute(reshape(Q,[npv npv ne ncq]),[1 2 4 3]);
    else
        L = [];
        Q = [];
    end    
        
%     Rq(:,:,1,1)
%     Rq(:,:,2,1)
%     pause
end
% udgg(1:ngv,:)
% pause

% Fluxes and source at Gauss points
if (app.uqpk == 1)
    [f, f_udg] = flux(pg, udgg, arg, time);      
    f  = reshape(f,[ngv ne ncu*nd]);
    f_udg  = reshape(f_udg,[ngv ne ncu*nd*nc]);
    jac = reshape(jac, [ngv ne]);
    for ie=1:ne         
        % mass matrix
        Me = (1/fc_q)*M(:,:,ie);        
        
        % L2 projection of the flux                
        ft = Me\(shapvg(:,1:npv)*reshape(bsxfun(@times,reshape(f(:,ie,:),[ngv ncu*nd]),jac(:,ie)),[ngv ncu*nd])); % f at the nodal points
        % values of the flux at the Gauss points
        f(:,ie,:) = reshape(shapvt(:,:,1)*ft,[ngv 1 ncu*nd]); % f at the gauss points                
        
        % L2 projection of the flux derivatives
        ft = Me\(shapvg(:,1:npv)*reshape(bsxfun(@times,reshape(f_udg(:,ie,:),[ngv ncu*nd*nc]),jac(:,ie)),[ngv ncu*nd*nc])); % f at the nodal points
        % values of the flux derivarives at the Gauss points
        f_udg(:,ie,:) = reshape(shapvt(:,:,1)*ft,[ngv 1 ncu*nd*nc]); % f at the gauss points                
    end
    f  = reshape(f,[ngv*ne ncu*nd]);         
    f_udg  = reshape(f_udg,[ngv*ne ncu*nd*nc]);    
    jac = reshape(jac, [ngv*ne 1]);
else
    [f, f_udg] = flux( pg, udgg, arg, time);
end

[s, s_udg] = source( pg, udgg, arg, time); 
f     = reshape(f,[ngv ne ncu nd]);
f_udg = permute(reshape(f_udg,[ngv ne ncu nd nc]),[1 2 3 5 4]); 
s     = reshape(s(:,1:ncu),[ngv*ne ncu]);
s_udg = reshape(s_udg(:,1:ncu,1:nc),[ngv*ne ncu nc]);


if adjoint==1
    [~,JUDG] = source_adjoint(pg, udgg, arg, time);     
    JUDG = bsxfun(@times,reshape(JUDG(:,:,:,1),[ngv ne nc]),reshape(jac,[ngv ne 1]));        
    Juq = reshape(master.shapvg(:,:,1)*reshape(JUDG,[ngv ne*nc]),[npv ne nc]);                 
    Ju = Juq(:,:,1:ncu);
    Jq = Juq(:,:,ncu+1:end);
else
    Ju = [];
    Jq = [];
end

% Update source term for time-dependent problems
if tdep    
    Stn = reshape(SH(:,:,1:ncu),[npv ne*ncu]);
    Stg = shapvt(:,:,1)*Stn;
    Stg = reshape(Stg,[ngv*ne ncu]);

    % axis symmetry
    axisflag = isfield(app,'axisymmetry');
    if axisflag
        xr = app.axisymmetry.*abs(pg(:,2));
        xr = reshape(xr,[ngv*ne 1]);
    else
        xr = ones(ngv*ne,1);
    end
    
    % time-derivative coefficients
    dtcoef = isfield(app,'dtcoef');
    if dtcoef
        dtcoef = app.dtcoef;
    else
        dtcoef = ones(ncu,1);
    end
    
%    s = s + bsxfun(@times,xr,Stg-udgg(:,1:ncu)*fc_u);    
%     s = s + Stg - udgg(:,1:ncu)*fc_u;    
    for i=1:ncu
        s(:,i) = s(:,i) + dtcoef(i)*xr.*(Stg(:,i)-udgg(:,i)*fc_u);
        s_udg(:,i,i) = s_udg(:,i,i) - dtcoef(i)*fc_u*xr;
%         s_udg(:,i,i) = s_udg(:,i,i) - fc_u;
    end    
end

% compute wrk and wrl to time with shape functions
wrk = zeros(ngv*(nd+1),ne*ncu);
wrl = zeros(ngv*(nd+1),ne*ncu*nc);
wrk(1:ngv,:) =  reshape(bsxfun(@times,s,jac),[ngv ne*ncu]);
wrl(1:ngv,:) = -reshape(bsxfun(@times,s_udg,reshape(jac,[ngv*ne 1 1])),[ngv ne*ncu*nc]);
for i=1:nd
    fk = bsxfun(@times,f(:,:,:,1),Xx(:,:,1,i));
    fl = bsxfun(@times,f_udg(:,:,:,:,1),Xx(:,:,1,i));
    for j=2:nd
        fk = fk + bsxfun(@times,f(:,:,:,j),Xx(:,:,j,i));
        fl = fl + bsxfun(@times,f_udg(:,:,:,:,j),Xx(:,:,j,i));
    end
    wrk(i*ngv+1:(i+1)*ngv,:) = reshape(fk,[ngv ne*ncu]);
    wrl(i*ngv+1:(i+1)*ngv,:) = -reshape(fl,[ngv ne*ncu*nc]);
end

% Volume residual
% [Phi Phi_xi Phi_eta] x [S.*jac; Fx.*Xx(:,:,1,1)+Fy.*Xx(:,:,2,1); Fx.*Xx(:,:,1,2)+Fy.*Xx(:,:,2,2)]
Ru = shapvg*wrk; % [npv ngv*(nd+1)] x [ngv*(nd+1) ne*ncu] 
Ru = reshape(Ru,[npv ne ncu]); 

% volume matrices
BD = reshape(shapvgdotshapvl,[npv*npv ngv*(nd+1)])*wrl;
BD = reshape(BD,[npv npv ne ncu nc]);




