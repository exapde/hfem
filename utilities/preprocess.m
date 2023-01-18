function [master,mesh,app] = preprocess(master,mesh,app,nHdgFaces)

if nargin < 3 
    hybrid = 'hdg';   
elseif isstruct(app) == 1
    hybrid = app.hybrid;    
elseif isstruct(app) == 0    
    hybrid = app;
end
if nargin<4; nHdgFaces = 0.1; end

% get dimensions
dim    = mesh.nd;
ngv    = master.ngv;
ngf    = master.ngf;
npv    = master.npv;
npf    = master.npf;
ne     = size(mesh.t,1);

%--------------- update MASTER structure ----------------%

% volume shape functions and their derivatives 
master.shapvgdotshapvl  = zeros(npv*npv,ngv,dim+1);      
for d=1:dim+1
    master.shapvt(:,:,d) = master.shapvl(:,:,d)';
    master.shapvg(:,:,d) = master.shapvl(:,:,d)*diag(master.gwvl);    
    for ii=1:npv
        for jj = 1:npv
            master.shapvgdotshapvl((ii-1)*npv+jj,:,d) = master.shapvg(jj,:,d).*master.shapvl(ii,:,1);                    
        end
    end            
end

% master.shapvndotshapvl  = zeros(npv*npv,npv,dim+1);      
% for d=1:dim+1
%     master.shapvnt(:,:,d) = master.shapnvl(:,:,d)';
%     master.shapvn(:,:,d) = master.shapnvl(:,:,d)*diag(master.gwnvl);    
%     for ii=1:npv
%         for jj = 1:npv
%             master.shapvndotshapvl((ii-1)*npv+jj,:,d) = master.shapvn(jj,:,d).*master.shapnvl(ii,:,1);                    
%         end
%     end
% end

% face shape functions and their derivatives 
master.shapfgdotshapfc  = zeros(npf*npf,ngf,dim);      
for d=1:dim
    master.shapft(:,:,d) = master.shapfc(:,:,d)';
    master.shapfg(:,:,d) = master.shapfc(:,:,d)*diag(master.gwfc);
    for ii=1:npf
        for jj = 1:npf
            master.shapfgdotshapfc((ii-1)*npf+jj,:,d) = master.shapfg(jj,:,d).*master.shapfc(ii,:,1);                    
        end
    end            
end

% master.shapfndotshapfc  = zeros(npf*npf,npf,dim);      
% for d=1:dim
%     master.shapfnt(:,:,d) = master.shapnfc(:,:,d)';
%     master.shapfn(:,:,d) = master.shapnfc(:,:,d)*diag(master.gwnfc);
%     for ii=1:npf
%         for jj = 1:npf
%             master.shapfndotshapfc((ii-1)*npf+jj,:,d) = master.shapfn(jj,:,d).*master.shapnfc(ii,:,1);                    
%         end
%     end            
% end

mesh.permgeom = mesh.perm(:,:,1);
master.shapmv = master.shapvt;
master.shapmf = master.shapft;
master.shapmh = master.shapmf;
master.permgeom = mesh.perm(:,:,1);

% --------------- update MESH structure ----------------%
ngrsiz = 800;
ngr    = ceil(ne/ngrsiz);
ngrne  = round(ne/ngr);
nk = 1:ngrne:ne;
nb = [nk(1:end); [nk(2:end)-1,ne]];
mesh.nb    = nb;

% Reorder faces
%[mesh.f,mesh.t2f,mesh.f2f,mesh.flev] = faceordering(mesh.f, mesh.t2f); 

mesh.bf = reshape(mesh.f(abs(mesh.t2f'),end),[size(master.perm,2) size(mesh.t,1)]);
mesh.f2f = mkf2f(mesh.f, mesh.t2f);

% if strcmp(hybrid,'edg')
%     hybridn = 1;
%     mesh.f(:,end+1) = 1;
%     %fprintf('\n --- EDG Algorithm: nsiz = %d \n\n', mesh.nsiz);
% elseif strcmp(hybrid,'hdg')
%     mesh.f(:,end+1) = 0;
%     hybridn = 0;
%     %fprintf('\n --- HDG Algorithm: nsiz = %d \n\n', mesh.nsiz);
% elseif strcmp(hybrid,'iedg')
%     hybridn = 2;
%     a = mesh.f(:,end);
%     i = (a>0);
%     mesh.f(:,end+1) = 0;
%     mesh.f(i,end)   = 1;
%     %fprintf('\n --- H-EDG Algorithm: nsiz = %d \n\n', mesh.nsiz);
% elseif strcmp(hybrid,'hedg')
%     hybridn = 2;
%     fhdg = hdgface(mesh,boundaryFaces);
%     mesh.f(:,end+1) = 1;
%     mesh.f(fhdg,end)= 0;
% elseif ~isempty(hybrid) 
%     mesh.f = feval(hybrid,mesh);
% else 
%     error('Hybrid flag is not valid');
% end

mesh.hybrid = hybrid;
[mesh.elcon,mesh.nsiz,mesh.edgnodes,hybridn] = elconnectivities(mesh,nHdgFaces);
mesh = block_crs(mesh,hybrid);

if mesh.porder==0
    mesh.dgnodes = mesh.edgnodes;
%     mesh.permgeom = master.perm(:,:,1);
%     master.permgeom = master.perm(:,:,1);
%    mastergeom = mkmaster(mesh,2);
    [plocal,~,plocfc,~,permnode,permedge,permface] = masternodes(1,dim,mesh.elemtype,mesh.nodetype);   
    shapvl = mkshape(1,plocal,master.gpvl,mesh.elemtype);
    shapfc = mkshape(1,plocfc,master.gpfc,mesh.elemtype);
    if dim==1
        permgeom = permnode;
    elseif dim==2
        permgeom = permedge;
    elseif dim==3
        permgeom = permface;
    end
    mesh.permgeom = permgeom;
    master.permgeom = permgeom;    
    for d=1:dim+1
        shapvt(:,:,d) = shapvl(:,:,d)';
    end    
    for d=1:dim
        shapft(:,:,d) = shapfc(:,:,d)';
    end    
    master.shapmv = shapvt;
    master.shapmf = shapft;
    master.shapmh = master.shapmf;
end

mesh.ncd = size(mesh.dgnodes,2);
mesh.nfe = size(mesh.t2f,2);
mesh.nve = size(mesh.t,2);
mesh.nvf = size(mesh.tlocfc,2);
mesh.nv  = mesh.np;
mesh.ndh = mesh.nsiz;
master.npe = master.npv;
master.nge = master.ngv;
master.nme = size(master.shapmv,2);
master.nmf = size(master.shapmf,2);

if nargin>=3 && nargout>=3
    if strcmp(hybrid,'edg') || strcmp(hybrid,'iedg') || strcmp(hybrid,'hedg')
        mesh.blkSize = app.ncu;
    elseif strcmp(hybrid,'hdg')    
        mesh.blkSize = app.ncu*master.npf;        
    end
    
    % --------------- update APP structure ----------------%    
    if isfield(app,'refine') == 0
        app.refine = 0;    
    end
    if isfield(app,'reflevel') == 0
        app.reflevel = 0;
    end 
    if app.refine==0
        app.reflevel = 0;
    end
    if isfield(app,'iterative') == 0
        app.iterative = 0;   
    end
    if isfield(app,'flux') == 0
        app.flux = 'flux';
    end
    if isfield(app,'fhat') == 0
        app.fhat = 'fhat';
    end
    if isfield(app,'fbou') == 0
        app.fbou = 'fbou';
    end
    if isfield(app,'source') == 0
        app.source = 'source';
    end
    if app.flg_q == 0
        app.localsolve = 0; % u
    elseif app.flg_q == 1 && app.flg_p == 0
        app.localsolve = 1; % uq
    elseif app.flg_q == 1 && app.flg_p == 1
        app.localsolve = 2; % upq     
    end    
    if strcmp(app.appname,'euler')
        appname = 0;    
    elseif strcmp(app.appname,'ns')
        appname = 1;        
    elseif strcmp(app.appname,'poisson')
        appname = 2;   
    elseif strcmp(app.appname,'ransSA')
        appname = 3;  
    else
        error('app.appname not implemented');
    end
    
    app.porder = master.porder;
    app.morder = mesh.porder;
    app.nb     = length(app.bcm);
    app.ndt    = length(app.dt);
    app.param  = reshape(cell2mat(app.arg),[],1); 
    app.flag   = [app.tdep app.wave app.alag app.adjoint app.linearproblem app.flg_q app.flg_p app.flg_g];
    app.factor = [app.fc_u app.fc_q app.fc_p app.dtfc app.alpha];
    app.problem = [hybridn appname app.linearSolver app.jacobianStep app.orderingStep];
    app.nparam = length(app.param);
    app.nflag = length(app.flag);
    app.nfactor = length(app.factor);
    app.nproblem = length(app.problem);    
end
