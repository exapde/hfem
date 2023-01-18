%%% -------------------------------------------------------- %%%
%%% --------------- Convergence case   --------------------- %%%
%%% -------------------------------------------------------- %%%

clear all
setapplicationpath('FM/mhd')

%tngrid = [6 12 24];
tngrid = 12;

porder      = 1;
nstage      = 1;
torder      = 2;
hybrid      = 'hdg';
elemtype    = 0;
nodetype    = 0;
T           = 1;
dtps        = 1/(4*tngrid);
ntime       = T/dtps;
ngrid       = tngrid;
gam         = 5/3;
tau         = 1;
app.arg     = {gam,tau};

% paramters files
% app.source  = 'source23d';
% app.flux    = 'flux23d';
% app.fbou    = 'fbou23d';
% app.fhat    = 'fhat23d';
app.source  = 'source2did';
app.flux    = 'flux2did';
app.fbou    = 'fbou2d';
app.fhat    = 'fhat2did';

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
app.nd   = 2;
%app.ncu  = 7 + app.nd;               % Number of components of U
app.ncu  = 6 + app.nd;
app.nch  = app.ncu;                  % Number of components of UH
app.nq   = 0;                        % Number of components of Q
app.nc   = app.ncu;                  % Number of components of UDG

% mesh
mesh   = mkmesh_square(ngrid,ngrid,porder,0,1/cos(pi/4),1/sin(pi/4));
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% Periodic boundary conditions
mesh = periodic(mesh,{1,'p(:,1)',3,'p(:,1)';4,'p(:,2)',2,'p(:,2)'},2);
                  
% initialize solutions
% u = [rho, rho*ux, rho*uy, rho*uz, rho*e, Bx, By, Bz, phi] 
UDG = initsol(mesh.dgnodes,4,gam);
UDG = UDG(:,1:8,:);
UH  = inituhat(master,mesh.elcon,UDG,app.ncu);

%[alpha1,tau] = glm(mesh,app,UDG,UH,ngrid);
%app.arg     = {gam,tau,alpha1};

% % time step looping 
time = 0;

TableError = zeros(length(ntime), 10);
for itime = 1:length(app.dt)
    dt = app.dt(itime);
    fprintf('Timestep :  %d,   Time :  %g\n', itime, time + dt);
    
%     [alpha1,alpha2,tau] = glm(mesh,app,UDG,UH,ngrid);
%     app.arg{2} =  tau;
            
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt,nstage,torder);    
    time = time + dt;  

%     % --- Evaluate L2 errors in time ---
%     err = calerror(UDG,mesh,master,@exactsolAlfven);
%     TableError(itime,:) = [time err'];
    
         figure(1)
         clf; scaplot(mesh,UDG(:,1,:),[],2,0); 
         axis on; axis square; axis tight; colorbar; colormap jet;
         title('density');
         drawnow
         
         figure(2)
         clf; scaplot(mesh,UDG(:,2,:),[],2,0); 
         axis on; axis square; axis tight; colorbar; colormap jet;
         title('rho*ux');
         drawnow       
        
         figure(3)
         clf; scaplot(mesh,UDG(:,3,:),[],2,0); 
         axis on; axis square; axis tight; colorbar; colormap jet;
         title('rho*uy');
         drawnow
         
         figure(4)
         clf; scaplot(mesh,UDG(:,4,:),[],2,0); 
         axis on; axis square; axis tight; colorbar; colormap jet;
         title('rho*uz');
         drawnow
         
         figure(5)
         clf; scaplot(mesh,UDG(:,5,:),[],2,0); 
         axis on; axis square; axis tight; colorbar; colormap jet;
         title('e');
         drawnow 
         
         figure(6)
         clf; scaplot(mesh,UDG(:,6,:),[],2,0); 
         axis on; axis square; axis tight; colorbar; colormap jet;
         title('Bx');
         drawnow
         
         figure(7)
         clf; scaplot(mesh,UDG(:,7,:),[],2,0); 
         axis on; axis square; axis tight; colorbar; colormap jet;
         title('By');
         drawnow  
        
         figure(8)
         clf; scaplot(mesh,UDG(:,8,:),[],2,0); 
         axis on; axis square; axis tight; colorbar; colormap jet;
         title('Bz');
         drawnow
         
%          figure(9)
%          clf; scaplot(mesh,UDG(:,9,:),[],2,0); 
%          axis on; axis square; axis tight; colorbar; colormap jet;
%          title('psi');
%          drawnow           
end

% maxL2rho  = max(TableError(:,2))
% maxL2rhou = max(TableError(:,3))
% maxL2rhov = max(TableError(:,4))
% maxL2rhow = max(TableError(:,5))
% maxL2rhoe = max(TableError(:,6))
% maxL2bx = max(TableError(:,7))
% maxL2by = max(TableError(:,8))
% maxL2bz = max(TableError(:,9))
% maxL2phi = max(TableError(:,10))
