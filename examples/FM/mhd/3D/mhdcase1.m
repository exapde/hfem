%%% -------------------------------------------------------- %%%
%%% --------------- Convergence case   --------------------- %%%
%%% -------------------------------------------------------- %%%

clear all
setapplicationpath('FM/mhd')

%tngrid = [6 12 24];
tngrid = 6;

porder      = 1;
nstage      = 3;
torder      = 3;
hybrid      = 'hdg';
elemtype    = 0;
nodetype    = 0;
T           = 7;
dtps        = 1e-2;
ntime       = T/dtps;
ngrid       = tngrid;
gam         = 2;
tau         = 1;
app.arg     = {gam,tau};

% paramters files
app.source  = 'source3d';
app.flux    = 'flux3d';
app.fbou    = 'fbou3d';
app.fhat    = 'fhat3d';


% Boundary conditions 
app.bcm     = [];           % 2: Wall, 1: Far-field
app.bcs     = [];
app.bcd     = [];           % 2: Slip wall, 1: Far-field
app.bcv     = [];

% various parameters
app.uqpk          = 0;
app.localsolve    = 0;          % U solve = 0, U,Q solve = 1          
app.hybrid        = hybrid;
app.iterative     = 0;          % 0 for direct solver, 1 for iterative solver.
app.wave          = false;
app.tdep          = true;
app.alag          = false;
app.adjoint       = 0;          % 1 if adjoint problem. 0 otherwise
app.linearproblem = 0;          % 0 if problem is linear. 1 if nonlinear
app.appname       = 'mhd';
app.linearSolver  = 1;
app.jacobianStep  = 0;
app.orderingStep  = 0;

% Flags and factors depending of system
app.flg_q = 0;              % Flag for the q equation
app.flg_p = 0;              % Flag for the p equation
app.flg_g = 0;              % Flag for the GCL equation
app.fc_q  = 0;              % Factor the q equation
app.fc_u  = 1;              % Factor the u equation
app.fc_p  = 0;              % Factor the p equation

% time parameters
app.torder = torder;
app.nstage = nstage;
app.dt     = dtps*ones(1,ntime);
app.time   = 0;

% dim
app.nd   = 3;
app.ncu  = 6 + app.nd;               % Number of components of U
app.nch  = app.ncu;                  % Number of components of UH
app.nq   = 0;                        % Number of components of Q
app.nc   = app.ncu;                  % Number of components of UDG

% mesh
mesh   = mkmesh_cube(ngrid,ngrid,ngrid,porder,2*pi,2*pi,2*pi,elemtype,nodetype);
master = mkmaster(mesh,3*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% Periodic boundary conditions
% mesh = periodic(mesh,{1,'p(:,1)',2,'p(:,1)'; 3,'p(:,2)',4,'p(:,2)'; ... 
%                       5,'p(:,3)',6,'p(:,3)'},2);                  
mesh = periodic(mesh,{1,'p(:,[2 3])',2,'p(:,[2 3])'; 3,'p(:,[1 3])',4,'p(:,[1 3])'; ... 
                      5,'p(:,[1 2])',6,'p(:,[1 2])'},2);
                  
% initialize solutions
% u = [rho, rho*ux, rho*uy, rho*uz, rho*e, Bx, By, Bz, phi] 
UDG = initsol(mesh.dgnodes,1,gam);
UH  = inituhat(master,mesh.elcon,UDG,app.ncu);


% % time step looping 
time = 0;
%TableError = zeros(length(ntime), 10);
for itime = 1:length(app.dt)
    dt = app.dt(itime);
    fprintf('Timestep :  %d,   Time :  %g\n', itime, time + dt);
            
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt,nstage,torder);    
    time = time + dt;  

%     % --- Evaluate L2 errors in time ---
%     UDGnum = zeros(size(UDG));
%     UDGnum(:,1,:) = 2 + sin(mesh.dgnodes(:,1,:) + mesh.dgnodes(:,2,:) - ...
%          (UDG(:,2,:)+UDG(:,3,:))*time./UDG(:,1,:));
%    % UDGnum(:,1,:) = UDG(:,1,:);
%     UDGnum(:,2,:) = UDG(:,2,:)./UDG(:,1,:);  
%     UDGnum(:,3,:) = UDG(:,3,:)./UDG(:,1,:); 
%    % UDGnum(:,4,:) = UDG(:,4,:) - UDG(:,1,:);
%     vv = UDGnum(:,2,:).*UDGnum(:,2,:) + UDGnum(:,3,:).*UDGnum(:,3,:);
%     bb = UDG(:,5,:).*UDG(:,5,:) + UDG(:,6,:).*UDG(:,6,:);
%     UDGnum(:,4,:) = (gam-1)*(UDG(:,4,:) - UDGnum(:,1,:).*vv/2 -bb/2);
%     err = calerror(UDGnum,mesh,master,@exactsol,time);
%     TableError(itime,:) = [time err'];

end

% maxL2rho = max(TableError(:,2))
