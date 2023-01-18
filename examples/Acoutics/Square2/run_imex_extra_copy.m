clear; close all
setapplicationpath('FM/wave');

% Global IMEX data
global IMEXdata;

porder = 5;
nstage = 3; % Number of implicit stages
torder = 3;
hybrid = 'hdg';
nref = 1; % nref = 1 is no refinements
% dt_choose = [0.13, 0.03, 0.0125, 0.007]; % Will choose the appropriate regions for each porder
dt_choose = 0.004; %0.0011 for fully explicit
dt = dt_choose;

% Time Discretization parameters
% Final time to be one Period
Tend = 1.1*sqrt(2);
% Number and value of time steps
% ntime = 10000;
ntime = Tend/dt;
% dt = Tend/ntime;

% Acoustic Wave Speed c
c  = 1;
c2 = c;
tau = 1;
Kmax = 1;

% Name of application files
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

% Boundary Conditions (6 is for IMEX BC)
app.bcm = [1;1;1;1;6];
app.bcs = [0;0;0;0;0]; 
app.bcd = app.bcm;
app.bcv = app.bcs;

% Application Flags and Parameters
app.linear = 1;
app.localsolve=1;
app.hybrid=hybrid;
app.iterative=0;
app.tdep = true;
app.wave = true;
app.alag = false;
app.uqpk = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

% Dimensions of Variables
app.nd   = 2;                      % Number of dimensions
app.ncu  = 1;                      % Number of components of U
app.nch  = app.ncu;                % Number of components of UH
app.ncq  = app.nch*app.nd;         % Number of components of Q
app.nc   = app.ncu+app.ncq;        % Number of components of UDG
app.ncp = 0;                       % Number of components of pressure

% hdg_wave parameters
app.arg = {c2,tau};

% Time integration parameters
app.nstage = nstage;
app.torder = torder;
app.dt = dt;

%%%%%%%%%%%%%% IMEX MESH %%%%%%%%%%%%%%%
% Creation of the Mesh and Preprocessing
mesh=mkmesh_square_imex(porder);
mesh2 = mesh;
% Refine the Mesh
mesh = refinemesh(mesh,porder,nref);
% Determine if an element is explicit or implicit
% mesh = create_mesh_imex(mesh,dt_choose(porder)/nref);
mesh = create_mesh_imex(mesh,dt_choose,c,Kmax);
% Find implicit/explicit boundary
mesh = find_imex_bound(mesh);
% Find implicit faces
mesh = create_mesh_imface(mesh);
% Find IMEX Boundary Number
mesh.imexB = length(mesh.bndexpr) + 1;
% Split Implicit-Explicit Mesh
[meshIM,meshEX] = split_IMEX_mesh(mesh);

meshIM.imexB = mesh.imexB;
meshEX.imexB = mesh.imexB;

% Create master Element
master = mkmaster(mesh,2*porder);
master2 = mkmaster(mesh2,2*porder);
% mesh.fimex = []; mesh2.fimex = []; meshEX = mesh2; meshIM = mesh; % TO BE SUPPRESSED !!!!
meshIM.imexB = mesh.imexB;
meshEX.imexB = mesh.imexB;

% Preprocess all the meshes
[  ~   ,meshIM] = preprocess(master,meshIM,hybrid);
[master,meshEX] = preprocess(master,meshEX,hybrid);
[master,mesh]   = preprocess(master,mesh,hybrid);
mesh.fimex = [];

%%%%%%%%%%%% MESH FOR POSTPROCESSING %%%%%%%%%%%%%%%%%
% Creation of the Mesh and Preprocessing
mesh1 = mkmesh_square_imex(porder+1);
% Refine the Mesh
mesh1 = refinemesh(mesh1,porder+1,nref);
% Find IMEX Boundary Number
mesh1.imexB = length(mesh.bndexpr) + 1;
% Split Implicit-Explicit Mesh
mesh1.imex = mesh.imex;
[meshIM1,meshEX1] = split_IMEX_mesh(mesh1);

meshIM1.imexB = mesh.imexB;
meshEX1.imexB = mesh.imexB;

% Create master Element
master1 = mkmaster(mesh1,2*(porder+1));

% Preprocess all the meshes
[  ~    ,meshIM1] = preprocess(master1,meshIM1,hybrid);
[master1,meshEX1] = preprocess(master1,meshEX1,hybrid);
[master1,mesh1]   = preprocess(master1,mesh1,hybrid);
mesh1.fimex = [];
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);


%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%
% L2-projection of the initial sol on implicit mesh
SOLimp = l2eprojection(meshIM,master,@exactsol,[],0,4);
UDGimp(:,:,:) = SOLimp(:,1:3,:);
PDGimp(:,:,:) = SOLimp(:,4,:);
UHimp = inituhat(master,meshIM.elcon,UDGimp,app.ncu);
% UHimp = zeros(1,1);

% L2-projection of the initial sol on explicit mesh
SOLexp = l2eprojection(meshEX,master,@exactsol,[],0,4);
UDGexp(:,:,:) = SOLexp(:,1:3,:);
PDGexp(:,:,:) = SOLexp(:,4,:);
UHexp = inituhat(master,meshEX.elcon,UDGexp,app.ncu);
% UHexp = zeros(1,1);

% Inverse of Mass Matrix for Explicit mesh
[Minv]=massinv(master, meshEX);

% Initialize time
time=0;
app.time = time;

for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
    
    % IMEX SOLVER
    [UDGexp,UDGimp,UHexp,UHimp,PDGexp,PDGimp] = hdg_solve_IMEX(master,meshEX,meshIM,app,UDGexp,UDGimp,UHexp,UHimp,PDGexp,PDGimp,Minv,time,dt,nstage);
    
    % Purely implicit, old method
    %[UDGimp,UHimp,PDGimp] = hdg_solve_dirk(master,mesh,app,UDGimp,UHimp,PDGimp,time,dt,nstage,torder);
    
    % Purely explicit, old method
    % UDGexp(:,4,:) = SOLexp(:,4,:);
    % UDGexp = ehdgrk(master,meshEX,app,Minv,UDGexp); 

    time = time + dt;
    app.time = time; 
    
end

% Calculate UHexp
IMEXdata = getVhatIMEX(meshIM,UHimp);
[UHexp,QHexp] = mkfhat(meshEX,app,UDGexp);

% UHexp = reshape(UHexp, 1, size(UHexp,1)*size(UHexp,2));

% The gradients on the implicit and explicit sides have opposite signs
UDGimp(:,2:3,:) = -UDGimp(:,2:3,:);

% ustar = postprocess(master,mesh,master1,mesh1,UDG(:,4,:),UDG(:,2:3,:));

% Get Errors for implicit and explicit meshes
% erru = calerror(UDG,mesh,master,@exactsol,itime*dt);
UDGexpAll = cat(2, UDGexp, PDGexp);
UDGimpAll = cat(2, UDGimp, PDGimp);

erruexp = calerror(UDGexpAll,meshEX,master,@exactsol,itime*dt);
erruimp = calerror(UDGimpAll,meshIM,master,@exactsol,itime*dt);

% Get Exact solutions on each meshes
uex   = exactsol(mesh.dgnodes,itime*dt);
uexIM = exactsol(meshIM.dgnodes,itime*dt);
uexEX = exactsol(meshEX.dgnodes,itime*dt);

% Plot fields and errors fields
figure(1); clf; scaplot(mesh,uex(:,1,:),[],2,1); axis off;
figure(2); clf; scaplot(meshEX,UDGexp(:,1,:),[],2,1); axis off;
figure(3); clf; scaplot(meshEX,abs(uexEX(:,1,:)-UDGexp(:,1,:)),[],2,1); axis off;
figure(4); clf; scaplot(meshIM,UDGimp(:,1,:),[],2,1); axis off;
figure(5); clf; scaplot(meshIM,abs(UDGimp(:,1,:)-uexIM(:,1,:)),[],2,1); axis off;


% Post Process U and V

% Ustar for Implicit and Explicit Domains
ustarIM = postprocessnd(master,meshIM,master1,meshIM1,cat(2,PDGimp,-UDGimp(:,2:3,:)));
% ustarIM2 = postprocess(master,meshIM,master1,meshIM1,PDGimp,-UDGimp(:,2:3,:));
ustarEX = postprocessnd(master,meshEX,master1,meshEX1,cat(2,PDGexp,-UDGexp(:,2:3,:)));

% Vstar for Implicit and Explicit Domains
vstarIM = postprocess930(master,meshIM,master1,meshIM1,UDGimp(:,1,:),UHimp);
vstarEX = postprocess930exp(master,meshEX,master1,meshEX1,UDGexp(:,1,:),UHexp);

% Post Processed Errors

% Vstar Errors
erruimp_vstar = calerror(vstarIM,meshIM1,master1,@exactsol,itime*dt);
erruexp_vstar = calerror(vstarEX,meshEX1,master1,@exactsol,itime*dt);

% Ustar Errors
erruimp_ustar = calerror(cat(2,ustarIM,ustarIM,ustarIM,ustarIM) ,meshIM1,master1,@exactsol,itime*dt);
erruexp_ustar = calerror(cat(2,ustarEX,ustarEX,ustarEX,ustarEX) ,meshEX1,master1,@exactsol,itime*dt);

% Set Ustar Errors to be last column
erruimp_ustar = erruimp_ustar(end);
erruexp_ustar = erruexp_ustar(end);


figure(6); clf; scaplot(meshIM1,ustarIM,[],2,1); axis off;