%%% -------------------------------------------------------- %%%
%%% --------------- Orszag-Tang vortex problem ------------- %%%
%%% -------------------------------------------------------- %%%

setapplicationpath('FM/mhd')

porder      = 2;
nstage      = 3;
torder      = 3;
hybrid      = 'hdg';
elemtype    = 0;
nodetype    = 0;
T           = 3;
dtps        = 1e-3;
ntime       = T/dtps;

gam   = 5/3;
epslm = 0.0;
Minf  = 0.3;                  % Infinity conditions
pinf  = 1/(gam*Minf^2);
Re    = 1.85e6;
Pr    = 0.72;
%tau   = 100; 
%alpha1 = 2;
%alpha2 = 2;

% paramters files
app.source  = 'source23d';
app.flux    = 'flux23d';
app.fbou    = 'fbou23d';
app.fhat    = 'fhat23d';
app.arg     = {gam,100};


% Boundary conditions 
%app.arg     = {gam,epslm,Re,Pr,Minf,tau,alpha1,alpha2};
app.bcm     = [];           % 2: Wall, 1: Far-field
app.bcs     = [];
app.bcd     = [];           % 2: Slip wall, 1: Far-field
app.bcv     = [];

% various parameters
app.uqpk            = 0;
app.localsolve      = 0;             
app.hybrid          = hybrid;
app.iterative       = 0;              % 0 for direct solver, 1 for iterative solver.
app.wave            = false;
app.tdep            = true;
app.alag            = false;
app.adjoint         = 0;               % 1 if adjoint problem. 0 otherwise
app.linearproblem   = 0;               % 0 if problem is linear. 1 if nonlinear
app.appname         = 'mhd';
app.linearSolver    = 1;
app.jacobianStep    = 0;
app.orderingStep    = 0;

% Flags and factors depending of system
app.flg_q   = 0;              % Flag for the q equation
app.flg_p   = 0;              % Flag for the p equation
app.flg_g   = 0;              % Flag for the GCL equation
app.fc_q    = 0;              % Factor the q equation
app.fc_u    = 1;              % Factor the u equation
app.fc_p    = 0;              % Factor the p equation

% time parameters
app.torder  = torder;
app.nstage  = nstage;
%app.dt = [1e-4,1e-3,1e-2,0.1,1,10,100];
app.dt      = dtps*ones(1,ntime);
app.time    = 0;

% dim
app.nd   = 2;
app.ncu  = 7 + app.nd;               % Number of components of U
app.nch  = app.ncu;                  % Number of components of UH
app.nq   = 0;                        % Number of components of Q
app.nc   = app.ncu;                  % Number of components of UDG

% mesh
ngrid = 20;
mesh   = mkmesh_square(ngrid,ngrid,porder,0,2*pi,2*pi);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% Periodic boundary conditions
mesh = periodic(mesh,{1,'p(:,1)',3,'p(:,1)';4,'p(:,2)',2,'p(:,2)'},2);

% initialize solutions
% u = [rho, rho*ux, rho*uy, rho*e, Bx, By, phi] 
UDG = initsol(mesh.dgnodes,2,gam);
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

[alpha1,tau] = glm(mesh,app,UDG,UH,ngrid);
app.arg     = {gam,tau,alpha1};

% time step looping 
time = 0;

figure
for itime = 1:length(app.dt)
    dt = app.dt(itime);
    fprintf('Timestep :  %d,   Time :  %g\n', itime, time + dt);
    
        
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt,nstage,torder);    
    time = time + dt;  
    
    clf; scaplot(mesh,UDG(:,1,:),[],2,0); 
    axis on; axis square; axis tight; colorbar; colormap jet;
    drawnow
end
