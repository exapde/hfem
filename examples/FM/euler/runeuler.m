clear

setapplicationpath('FM/euler')

% ------ Parameters --------- %
hybrid      = 'hdg';
elemtype    = 0;
nodetype    = 1;
porder      = 4;        % HDG order
nstage      = 3;        % n-stage DIRK
torder      = 3;        % t-order DIRK
T           = 5;        % Final time
ngrid       = 129;       % mesh discretization
dtps        = 1/(8*(ngrid-1));     % Time step
ntime       = T/dtps;   % number time iteration
ntime       = ceil(ntime);
gam         = 1.4;      % adiabatic constant depending on the physical properties of the fluid
icase       = 6;        % 1:scalar case, 2:Orszag-Tang vortex, 3:Rotor problem, 4:Smooth Alfven wave, 5:Smooth vortex problem


% app.source  = 'source_ns2d';
% app.flux    = 'flux_AV2_ns2d';
% app.fbou    = 'fbou_nsNEW';
% app.fhat    = 'fhat_ns2d';

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
app.flg_q   = 1;
app.flg_p   = 0;
app.flg_g   = 0;
app.fc_q    = 1;
app.fc_u    = 1;
app.fc_p    = 0;
app.ns = 1;

app.nd      = 2;            % Dimension
app.ncu     = 2 + app.nd;   % Number of components of U
app.nch     = app.ncu;      % Number of componets of UH
app.ncq     = app.nd*app.ncu; 
app.nc      = app.ncu+app.ncq;      % Number of componeents of UDG
app.nco     = 3;
app.ncp     = 0; 
app.AV = 10;
app.convStabMethod = 1;

app.time    = 0;
app.dt      = dtps*ones(1,ntime);
app.torder  = torder;
app.nstage  = nstage;
           
app.adjoint         = 0;               % 1 if adjoint problem. 0 otherwise
app.linearproblem   = 0;               % 0 if problem is linear. 1 if nonlinear
app.appname         = 'ns';
app.linearSolver    = 1;
app.jacobianStep    = 0;
app.orderingStep    = 0;
app.morder          = [porder porder];
app.porder          = [porder porder];
app.nodetype        = nodetype;
app.pgauss          = 2*[porder porder];
app.pgaussR         = 2*[porder porder];
app.overlappinglevel= 1;
app.preconditioner  = 0;
app.quadtype        = [0 0];
check               = 0;

% ----- Mesh depending of the test case ------- %
mesh = mkmesh_rect(ngrid,ngrid,porder,0,[-0.5,0.5,-0.5,0.5],elemtype);

master  = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% ------ Periodic boundary conditions ------- %
bndexpr = mesh.bndexpr;
periodicexpr = {1,'p(:,1)',3,'p(:,1)';4,'p(:,2)',2,'p(:,2)'};
mesh = get_hField(mesh);
mesh.dgnodes(:,3,:) = mesh.hField;
mesh.dgnodes(:,4,:) = getDistance2Boundary(mesh, 1);

% ------ Initialize solution ------------- %
%UDG0    = initsol(mesh.dgnodes,icase,gam);
UDG0 = initL2(master,mesh,icase,gam);
%UH0  = inituhat(master,mesh.elcon,UDG0,app.ncu);
UH0  = initL2f(master,mesh,icase,gam);
[Q,~,~] = getq(master,mesh,UDG0,UH0);
UDG0 = cat(2,UDG0,Q);
%  ----- Compute Tau and alpha1 ------ %
app.arg{1}   = gam;
app.arg      = {gam,0.0,1.0e30,0.72,0.4,2.0};

app.ncd     = size(mesh.dgnodes,2);
elementtype = elemtype*ones(mesh.ne,1);

% ------- GMRES Parameters --------- %
app.restart   = 200;
app.gmrestol  = 1e-12;
app.gmresiter = 2000;

% ------- Newton Parameters ---------- %
app.newtoniter = 100;  % def 10
app.newtontol  = 1e-8; % def 1e-7

% ------- Number of processors ------- %
nproc       = 4;
app.nfile   = nproc;

app.flag = 16;
if nproc>1 % parallel preprocessing
    apppar = digasopre(app,'euler',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,[],[],elementtype,bndexpr,periodicexpr,nproc,0,check,mesh.perm);
    apppar.fileout = 'eulerout';
    return;
else  % serial preprocessing
    appser = digasopre(app,'euler',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,[],[],elementtype,bndexpr,periodicexpr,nproc,0,check,mesh.perm);
    appser.fileout = 'eulerout';
    return;
end