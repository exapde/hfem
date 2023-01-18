function [RU,UDG,UH] = ldg_residual(master,app,mesh,U,S)

npv = size(U,1);                
ne  = size(mesh.dgnodes,3);
npf = size(master.perm,1);
% nfe = size(master.perm,2);
ncu = app.ncu;
%nc  = app.nc;
ns  = 500;

nb = ceil(ne/ns);          
nk = 1:ns:ne;
nm = [nk(1:end); [nk(2:end)-1,ne]];    

if app.flg_q
    % compute uhat and q
    UH = ldg_uhat(mesh,master,app,U(:,1:ncu,:));
    %UH = permute(reshape(UH,[ncu npf mesh.nf]),[2 1 3]);
    QDG = getq(master, mesh, U, reshape(permute(UH,[2 1 3]),[ncu npf*mesh.nf]));
    UDG = cat(2,U,QDG);
else    
    UH   = zeros(npf,ncu,mesh.nf);
    UDG = U;
end

% compute face integral to form RU
RU = zeros(npv,ncu,ne);
RU = faceint(master,mesh,app,UDG,UH,RU);

% compute volume integral to form RU
UDG = permute(UDG,[1 3 2]);
S  = permute(S,[1 3 2]);
mesh.dgnodes = permute(mesh.dgnodes,[1 3 2]);
for j=1:nb
    id = nm(1,j):nm(2,j);    
    
    % compute volume integrals    
    ru = volint(master,app,mesh.dgnodes(:,id,:),UDG(:,id,:),S(:,id,:));                       
    RU(:,:,id) = RU(:,:,id) + ru;
end

UDG = permute(UDG,[1 3 2]);


function Ru = volint(master,app,dgnodes,UDG,SH)
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
fc_u = app.fc_u;
source = str2func(app.source);
flux   = str2func(app.flux);

% Shap functions and derivatives
shapvt = master.shapvt;
shapvg = reshape(master.shapvg,[npv ngv*(nd+1)]);

% DG solution at Gauss points
udgg = reshape(UDG,[npv ne*nc]);
udgg = shapvt(:,:,1)*udgg;
udgg = reshape(udgg,[ngv*ne nc]);

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
Ru = permute(reshape(Ru,[npv ne ncu]),[1 3 2]); 

function RU = faceint(master,mesh,app,UDG,UH,RU)
% FACEINTND compute face integrals 

npf = size(mesh.perm,1);
fbou   = str2func(app.fbou);
fhat   = str2func(app.fhat);

% Shap functions 
perm            = master.perm(:,:,1);
shapft          = master.shapft(:,:,1);
shapfg          = master.shapfg(:,:,1);
elcon = reshape(mesh.elcon,[npf mesh.nfe mesh.ne]);
for i = 1:mesh.nf
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i      
    if fi(2)>0           % face i is an interior face                
        kf = mesh.t2f(fi,:);         % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element
        i2 = kf(2,:)==i;  % obtain the index of face i in the 2nd element                                            
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        j2 = elcon(:,i2,fi(2)) - (i-1)*npf;                
        udg1 = shapft*UDG(perm(j1,i1),:,fi(1));                
        udg2 = shapft*UDG(perm(j2,i2),:,fi(2));                
        uhg  = shapft*UH(:,:,i);
        xdg1 = mesh.dgnodes(perm(j1,i1),:,fi(1));
        [pg1, nlg1, jac1] = facegeom(master.shapmf,xdg1);
        
        %FH = fhat(nlg, pg, udgg, uhg, arg, time);   
        fhg = fhat(nlg1, pg1, udg1, udg2, uhg, app.arg, app.time);
        cnt = shapfg*diag(jac1)*fhg;                
        RU(perm(j1,i1),:,fi(1)) = RU(perm(j1,i1),:,fi(1)) - cnt;
        RU(perm(j2,i2),:,fi(2)) = RU(perm(j2,i2),:,fi(2)) + cnt;          
    else % face i is a boundary face
        kf = mesh.t2f(fi(1),:); % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element                 
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        udg1 = shapft*UDG(perm(j1,i1),:,fi(1));
        uhg  = shapft*UH(:,:,i);
        xdg1 = mesh.dgnodes(perm(j1,i1),:,fi(1));
        [pg1, nlg1, jac1] = facegeom(master.shapmf,xdg1);                        
        b = -fi(2);
        ib = app.bcm(b);
        uinf = app.bcs(b,:);
        fhg = fbou(ib, uinf, nlg1, pg1, udg1, uhg, app.arg, app.time);
        cnt = shapfg*diag(jac1)*fhg;               
        RU(perm(j1,i1),:,fi(1)) = RU(perm(j1,i1),:,fi(1)) - cnt;        
    end        
end 

function [pg, nlg, jac] = facegeom(shapgeomft,pn)
% FACEGEOM computes dg nodes, Jacobian determinant and normal vectors at Gauss points  

%   [pg, nlg, jac] = facegeom(shapgeomft,dgnodes,perm)
%
%    SHAPGEOMFT :  Shape functions and derivatives at Gauss points
%    DGNODES    :  Geometry DG nodes 
%    PERM       :  Indices of the boundary nodes
%    PG         :  Physical nodes at Gauss points 
%    nlg        :  Normal vector at Gauss points
%    jac        :  Determinant of the Jacobian mapping 

nq    = size(pn,2);
ngf   = size(shapgeomft,1);
npf   = size(shapgeomft,2);
nd    = size(shapgeomft,3);

if nd>1
    dshapft  = reshape(permute(shapgeomft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
    dpg = dshapft*pn(:,1:nd);
    dpg = permute(reshape(dpg,[ngf nd-1 nd]), [1 3 4 5 2]);    
    dpg = reshape(dpg,[ngf,nd,nd-1]);    
end

shapgeomft   = shapgeomft(:,:,1);
pg = shapgeomft*pn;
pg = reshape(pg,[ngf nq]);

switch nd
    case 1
        jac = ones(1,1);
        nlg = [ones(1,1); -ones(1,1)];        
        nlg = nlg(:);
    case 2
        jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nlg   = [dpg(:,2),-dpg(:,1)];
        nlg   = bsxfun(@rdivide, nlg, jac);
    case 3
        nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
        nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
        nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
        jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
        nlg   = bsxfun(@rdivide, nlg, jac);
    otherwise
        error('Dimension is not implemented');
end



% function r = rinvexpl(master,mesh,app,u,time)
% %RINVEXPL Calculates the residual vector for explicit time stepping.
% %   R=RINVEXPL(MASTER,MESH,APP,U,TIME)
% %
% %      MASTER:       Master structure
% %      MESH:         Mesh structure
% %      APP:          Application structure
% %      U(NPL,NC,NT): Vector of unknowns
% %                    NPL = size(mesh.plocal,1)
% %                    NC = app.nc (number of equations in system)
% %                    NT = size(mesh.t,1)
% %      TIME:         Time
% %      R(NPL,NC,NT): Residual vector (=dU/dt) (already divided by mass
% %                    matrix)                             
% %
% % - Written by: J. Peraire
% %
% nt  = size(mesh.t,1);
% nf  = size(mesh.f,1);
% nc  = app.nc;
% npl = size(mesh.dgnodes,1);
% ng  = size(master.gwgh,1);
% ng1d = size(master.gw1d,1);
% 
% r = zeros(size(u));
% 
% % Interfaces
% sh1d = squeeze(master.sh1d(:,1,:));
% perm = master.perm;
% ni = find(mesh.f(:,4) < 0,1)-1;
% 
% % Interior First
% for i=1:ni
%     ipt  = sum(mesh.f(i,1:2));
%     el  = mesh.f(i,3);
%     er  = mesh.f(i,4);
%     
%     ipl = sum(mesh.t(el,:))-ipt;
%     isl = find(mesh.t(el,:) == ipl);
%     iol = 1; if mesh.t2f(el,isl) < 0, iol = 2; end
%     
%     ipr = sum(mesh.t(er,:))-ipt;
%     isr = find(mesh.t(er,:) == ipr);
%     ior = 1; if mesh.t2f(er,isr) < 0, ior = 2; end
%     
%     if app.pg
%         pl = mesh.dgnodes(perm(:,isl,iol),:,el);
%         plg = sh1d'*pl;
%     else
%         plg = [];
%     end
%     
%     if mesh.fcurved(i)
%        xxi = squeeze(master.sh1d(:,2,:))'*squeeze(mesh.dgnodes(perm(:,isl,iol),1,el)); 
%        yxi = squeeze(master.sh1d(:,2,:))'*squeeze(mesh.dgnodes(perm(:,isl,iol),2,el));  
%        dsdxi = sqrt(xxi.^2+yxi.^2);
%        nl = [yxi./dsdxi,-xxi./dsdxi];
%        dws = master.gw1d.*dsdxi;
%     else
%        dx = mesh.p(mesh.f(i,2),:)-mesh.p(mesh.f(i,1),:);
%        dsdxi = sqrt(sum(dx.^2));
%        nl = [dx(2),-dx(1)]/dsdxi;
%        nl = repmat(nl,ng1d,1);
%        dws = master.gw1d*dsdxi;
%     end
%     
%     ul = u(perm(:,isl,iol),:,el);
%     ulg = sh1d'*ul;
%         
%     ur = u(perm(:,isr,ior),:,er);
%     urg = sh1d'*ur;
%  
%     fng = app.finvi( ulg, urg, nl, plg, app.arg, time);
%     cnt = sh1d*diag(dws)*fng;
%    
%     r(perm(:,isl,iol),:,el) = r(perm(:,isl,iol),:,el) - cnt;
%     r(perm(:,isr,ior),:,er) = r(perm(:,isr,ior),:,er) + cnt;  
% end
%     
% % Now Boundary
% for i=ni+1:nf
%     ipt  = sum(mesh.f(i,1:2));
%     el  = mesh.f(i,3);
%     ib  = -mesh.f(i,4);
%     
%     ipl = sum(mesh.t(el,:))-ipt;
%     isl = find(mesh.t(el,:) == ipl);
%     iol = 1; if mesh.t2f(el,isl) < 0, iol = 2; end
% 
%     if app.pg
%         pl = mesh.dgnodes(perm(:,isl,iol),:,el);
%         plg = sh1d'*pl;
%     else
%         plg = [];
%     end
%     
%     if mesh.fcurved(i)
%        xxi = squeeze(master.sh1d(:,2,:))'*squeeze(mesh.dgnodes(perm(:,isl,iol),1,el)); 
%        yxi = squeeze(master.sh1d(:,2,:))'*squeeze(mesh.dgnodes(perm(:,isl,iol),2,el));  
%        dsdxi = sqrt(xxi.^2+yxi.^2);
%        nl = [yxi./dsdxi,-xxi./dsdxi];
%        dws = master.gw1d.*dsdxi;
%     else
%        dx = mesh.p(mesh.f(i,2),:)-mesh.p(mesh.f(i,1),:);
%        dsdxi = sqrt(sum(dx.^2));
%        nl = [dx(2),-dx(1)]/dsdxi;
%        nl = repmat(nl,ng1d,1);
%        dws = master.gw1d*dsdxi;
%     end
%     
%     ul = u(perm(:,isl,iol),:,el);
%     ulg = sh1d'*ul;
%     
%     fng = app.finvb( ulg, nl, app.bcm(ib), app.bcs(app.bcm(ib),:), plg, app.arg, time);
%         
%     r(perm(:,isl,iol),:,el) = r(perm(:,isl,iol),:,el) - sh1d*diag(dws)*fng;
% end
% 
% 
% % Volume integral
% shap = squeeze(master.shap(:,1,:));
% shapxi = squeeze(master.shap(:,2,:))*diag(master.gwgh);
% shapet = squeeze(master.shap(:,3,:))*diag(master.gwgh); 
% 
% for i=1:nt
%    if app.pg
%        pg = shap'*mesh.dgnodes(:,:,i);
%    else
%        pg = [];
%    end
%    
%    if mesh.tcurved(i)
%        xxi = squeeze(master.shap(:,2,:))'*squeeze(mesh.dgnodes(:,1,i));
%        xet = squeeze(master.shap(:,3,:))'*squeeze(mesh.dgnodes(:,1,i));
%        yxi = squeeze(master.shap(:,2,:))'*squeeze(mesh.dgnodes(:,2,i));
%        yet = squeeze(master.shap(:,3,:))'*squeeze(mesh.dgnodes(:,2,i));
%        jac = xxi.*yet - xet.*yxi;
%        shapx =   shapxi*diag(yet) - shapet*diag(yxi);
%        shapy = - shapxi*diag(xet) + shapet*diag(xxi);
%        M = shap*diag(master.gwgh.*jac)*shap';
%    else
%        xxi = mesh.p(mesh.t(i,2),1) - mesh.p(mesh.t(i,1),1);
%        xet = mesh.p(mesh.t(i,3),1) - mesh.p(mesh.t(i,1),1);
%        yxi = mesh.p(mesh.t(i,2),2) - mesh.p(mesh.t(i,1),2);
%        yet = mesh.p(mesh.t(i,3),2) - mesh.p(mesh.t(i,1),2);
%        jac = xxi*yet - xet*yxi;
%        shapx =   shapxi*yet - shapet*yxi;
%        shapy = - shapxi*xet + shapet*xxi; 
%        M = master.mass*jac;
%    end
% 
%    ug = shap'*squeeze(u(:,:,i));
%    
%    if ~isempty(app.src)
%        sr = app.src( ug, [], pg, app.arg, time);
%        r(:,:,i) =  r(:,:,i) + shap*diag(master.gwgh)*sr;
%    end   
%    
%    [fgx,fgy] = app.finvv( ug, pg, app.arg, time);
%    r(:,:,i) =  r(:,:,i) + shapx*fgx+shapy*fgy;
%    
%    r(:,:,i) = M\r(:,:,i);
% end
%    
%    
