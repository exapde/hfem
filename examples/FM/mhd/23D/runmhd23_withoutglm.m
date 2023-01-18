clear

setapplicationpath('FM/mhd')

hybrid      = 'hdg';
elemtype    = 0;
nodetype    = 0;
porder      = 1;
nstage      = 1;
torder      = 2;
T           = 1;
ngrid       = 12;
dtps        = 1/(4*ngrid);
ntime       = T/dtps;
gam         = 5/3;
icase       = 4; % 1:scalar case, 2:Orszag-Tang vortex, 3:Rotor problem, 4:Smooth Alfven wave, 5:Smooth vortex problem

app.source  = 'source_mhd';
app.flux    = 'flux_mhd';
app.fbou    = 'fbou_mhd';
app.fhat    = 'fhat_mhd';

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

app.nd      = 2;
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
app.morder          = [porder porder];
app.porder          = [porder porder];
app.nodetype        = nodetype;
app.pgauss          = 2*[porder porder];
app.pgaussR         = 2*[porder porder];
app.overlappinglevel= 1;
app.preconditioner  = 0;
app.quadtype        = [0 0 0];
check               = 0;

% mesh
if icase == 1 || icase == 2
   mesh = mkmesh_square(ngrid,ngrid,porder,0,2*pi,2*pi);
elseif icase == 3
   mesh = mkmesh_square(ngrid,ngrid,porder,0,1,1);
elseif icase == 4 
   mesh = mkmesh_square(ngrid,ngrid,porder,0,1/cos(pi/4),1/sin(pi/4));
elseif icase == 5
   mesh = mkmesh_rect(ngrid,ngrid,porder,0,[-5,5,-5,5]);
end
master  = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% Periodic boundary conditions
bndexpr = mesh.bndexpr;
periodicexpr = {1,'p(:,1)',3,'p(:,1)';4,'p(:,2)',2,'p(:,2)'};

% Initialize solution
UDG0    = initsol(mesh.dgnodes,icase,gam);
UDG0 = UDG0(:,1:8,:);
UH0  = inituhat(master,mesh.elcon,UDG0,app.ncu);

% Compute Tau and alpha1 values 
app.arg{1}   = gam;
[alpha1,tau] = glm(mesh,app,UDG0,UH0,ngrid);
app.arg      = {gam,tau};

app.ncd     = size(mesh.dgnodes,2);
%check = 0;
elementtype = elemtype*ones(mesh.ne,1);

% GMRES Parameters
app.restart   = 200;
app.gmrestol  = 1e-10;
%app.gmresiter = 2000;

% Newton Parameters
%app.newtoniter = 20;  % def 10
app.newtontol  = 1e-10; % def 1e-7

% Number of processors :
nproc       = 1;
app.nfile   = nproc;

if nproc>1
    % parallel preprocessing
    apppar = digasopre(app,'mhdpar',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,[],[],elementtype,bndexpr,periodicexpr,nproc,0,check,mesh.perm); 
    apppar.fileout = 'mhdparout';
    return;
%    [UDGpar,UHpar] = getsolfrombinaryfile(apppar.fileout,nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
else
    % serial preprocessing
    appser = digasopre(app,'mhdser',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,[],[],elementtype,bndexpr,periodicexpr,nproc,0,check,mesh.perm);
    appser.fileout = 'mhdserout';
    return;
    % [UDGser,UHser] = getsolfrombinaryfile(appser.fileout,nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
end

