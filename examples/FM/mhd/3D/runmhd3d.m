clear

setapplicationpath('FM/mhd')

hybrid      = 'hdg';
elemtype    = 0;
nodetype    = 0;
porder      = 1;
nstage      = 3;
torder      = 3;
tau         = 1;
T           = 7;
dtps        = 1e-2;
ntime       = T/dtps;
gam         = 2;
app.arg     = {gam,tau};

app.source  = 'source_mhd3d';
app.flux    = 'flux_mhd3d';
app.fbou    = 'fbou_mhd3d';
app.fhat    = 'fhat_mhd3d';

app.iterative = 0;
app.hybrid  = hybrid;
app.localsolve = 0;
app.bcm     = [];
app.bcs     = [];
app.bcd     = [];
app.bcv     = [];

app.wave    = false;
app.tdep    = true;
app.alag    = false;
app.uqpk    = 0;

app.flg_q   = 0;
app.flg_p   = 0;
app.flg_g   = 0;
app.fc_q    = 0;
app.fc_u    = 1;
app.fc_p    = 0;

app.nd      = 3;
app.ncu     = 6 + app.nd;      % Number of components of U
app.nch     = app.ncu;         % Number of componets of UH
app.nq      = 0;               % Number of componets of Q
app.nc      = app.ncu;         % Number of componeents of UDG
app.nco     = 0;
app.ncq     = 0; 
app.ncp     = 0; 

app.time    = 0;
app.dt      = dtps*ones(1,ntime);
app.torder  = torder;
app.nstage  = nstage;
           
app.adjoint         = 0;               % 1 if adjoint problem. 0 otherwise
app.linearproblem   = 0;               % 0 if problem is linear. 1 if nonlinear
app.appname         = 'mhd';
app.linearSolver    = 1;
app.jacobianStep    = 0;
app.orderingStep    = 0;
app.morder          = [porder porder porder];
app.porder          = [porder porder porder];
app.nodetype        = nodetype;
app.pgauss          = 3*[porder porder porder];
app.pgaussR         = 3*[porder porder porder];
app.overlappinglevel= 1;
app.preconditioner  = 0;
app.quadtype        = [0 0 0];
check               = 0;
bndexpr = {'all(p(:,1)<1e-3)','all(p(:,1)>2*pi-1e-3)', ...
           'all(p(:,2)<1e-3)','all(p(:,2)>2*pi-1e-3)', ...
           'all(p(:,3)<1e-3)','all(p(:,3)>2*pi-1e-3)'}; 

% mesh
ngrid  = 24;
mesh   = mkmesh_cube(ngrid,ngrid,ngrid,porder,2*pi,2*pi,2*pi,elemtype,nodetype);
master = mkmaster(mesh,3*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% Periodic boundary conditions
periodicexpr = {1,'p(:,[2 3])',2,'p(:,[2 3])'; 3,'p(:,[1 3])',4,'p(:,[1 3])'; ... 
                      5,'p(:,[1 2])',6,'p(:,[1 2])'};

% Initialize solution
% 1 : scalar case
% 2 : Orszag-Tang vortex
% 3 : Rotor problem
UDG0    = initsol3d(mesh.dgnodes,1,gam);

app.ncd     = size(mesh.dgnodes,2);
elementtype = elemtype*ones(mesh.ne,1);

% GMRES Tolerance
app.restart     = 200;
app.gmrestol    = 1e-7;

% Number of processors :
nproc       = 12;
app.nfile   = nproc;

if nproc>1
    % parallel preprocessing
    apppar = digasopre(app,'mhdpar',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,[],[],elementtype,bndexpr,periodicexpr,nproc,0,check,mesh.perm); %% !!!!!!! Need to define master.perm inside digasopre
    apppar.fileout = 'mhdparout';
    return;
    % [UDGpar,UHpar] = getsolfrombinaryfile(apppar.fileout,nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
else
    % serial preprocessing
    appser = digasopre(app,'mhdser',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,[],[],elementtype,bndexpr,periodicexpr,nproc,0,check,mesh.perm);
    appser.fileout = 'mhdserout';
    return;
    % [UDGser,UHser] = getsolfrombinaryfile(appser.fileout,nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
end

