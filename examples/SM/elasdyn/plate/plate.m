%setapplicationpath('SM/LE_uq')
%setapplicationpath('SM/NeoHookean_uq')
%setapplicationpath('SM/SVK_uq')

is_SVK = true;

if is_SVK  % Saint Venant-Kirchhoff
    setapplicationpath('SM/SVK_uq')
    app.fbou   = 'fbou_SVK';
    app.source = 'source_SVK';
else      % Linear elastic model
    setapplicationpath('SM/LE_uq')
    app.fbou   = 'fbou';
    app.source = 'source';
end

porder = 2;
nstage = 3;
torder = 3;
Cst_ht = 0.05;
elemtype = 1;
nodetype = 1;
hybrid = 'hdg';

mu = 1.;
tau = 2.;
lambda = 1.5;
rho = 1;
h = 0.0001;
ampl = 0.4;


Tend = 1.;
%ntime = 10;
%dt = Tend/ntime;

app.flux = 'flux';
app.fhat = 'fhat';
app.localsolve=1;
ui = [0, 0, 0];
app.bcm  = [1,1,1,1,3,3];  
app.bcs  = [ui; ui; ui; ui; ui; ui; ui];
app.bcd  = [1];  
app.bcv  = [ui];
app.arg = {mu,lambda,tau,rho,ampl};

app.iterative = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.wave = true;
app.tdep = true;
app.alag = false;
app.uqpk = true;

app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_q = 1;
app.fc_u = 1;
app.fc_p = 0;

app.nd   = 3;
app.ncu  = 3;               % Number of components of U
app.nch  = app.ncu;                % Number of componets of UH
app.ncq  = app.nch*app.nd;         % Number of componets of Q
app.nc   = app.ncu+app.ncq;         % Number of componeents of UDG
app.ncp = 0;

app.time = [];
app.dtfc = [];
app.alpha = [];


% Parameters for convergence study :
n_ref   = 5;        % Number of refinements
ngrid0  = 2;        % Initial nb of element in each directions
refcoef = sqrt(2.); % Refinement coefficient

% Initialization of errors results
L2errors.errors = zeros(n_ref,app.nc);
L2errors.orders = zeros(n_ref,app.nc);
L2errors.Uerror = zeros(n_ref,2);
L2errors.Qerror = zeros(n_ref,2);
L2errors.h      = zeros(1,n_ref);
L2errors.dt     = zeros(1,n_ref);

for i_ref = 1:n_ref
    
    fprintf('########################################### \n');
    fprintf('     REFINEMENT STEP : %d\n',i_ref);
    fprintf('########################################### \n');
    ngrid = round(ngrid0*refcoef^(i_ref-1))+1;
    fprintf('Number of elements in each direction :  %d\n', ngrid-1);
    dx = 1./(ngrid-1);
    
    % Time step computation
    dt = (Cst_ht*dx^(porder+1))^(1./torder);
    ntime = round(Tend/dt);
    dt = Tend/ntime;
    fprintf('Computed time step : %d  --- Number of time steps : %d\n', [dt ntime]);
    
    mesh = mkmesh_cube(ngrid,ngrid,2,porder,1,1,h,elemtype,nodetype);
    master = mkmaster(mesh,2*porder);
    [master,mesh] = preprocess(master,mesh,hybrid);
    
    x = mesh.dgnodes(:,1,:);
    y = mesh.dgnodes(:,2,:);
    z = mesh.dgnodes(:,3,:);
    UDG = initu(mesh,{0;0;0;0;0;0;0;0;0;0;0;0});
    UDG(:,1,:) = 0;
    UDG(:,2,:) = 0;
    UDG(:,3,:) = ampl*pi*sin(pi*x).*sin(pi*y);
    UDG(:,4,:) = -1;
    UDG(:,5,:) = 0;
    UDG(:,6,:) = 0;
    UDG(:,7,:) = 0;
    UDG(:,8,:) = -1;
    UDG(:,9,:) = 0;
    UDG(:,10,:) = 0;
    UDG(:,11,:) = 0;
    UDG(:,12,:) = -1;
    UH = inituhat(master,mesh.elcon,UDG,app.ncu);
    PDG = initu(mesh,{0,0,0});
    PDG(:,1,:) = x;
    PDG(:,2,:) = y;
    PDG(:,3,:) = z;
    
    mesh1=mesh;
    time=0;
    for itime = 1:ntime
        fprintf('Timestep : %d,  Time : %d\n', itime, time+dt);
        
        [UDG,UH,PDG] = hdg_solve_dirk(master,mesh,app,UDG,UH,PDG,time,dt,nstage,torder);
        time = time + dt;
        
        mesh1.dgnodes = PDG(:,1:3,:);
        %figure(1); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;
        %pause(0.05);
    end

    % Computation of the errors
    err = calerror(UDG,mesh,master,@exactsol,time);
    L2errors.h (i_ref)  = 1./(ngrid-1);   % Mesh size
    L2errors.dt(i_ref)  = dt;             % Time step
    L2errors.errors(i_ref,:) = err;       % Errors
    % Convergence Orders
    if i_ref== 1 
        L2errors.orders(i_ref,:) = 0.;
    else
        L2errors.orders(i_ref,:) = log(L2errors.errors(i_ref,:)./L2errors.errors(i_ref-1,:))...
                                 /(log(L2errors.h(i_ref)/L2errors.h(i_ref-1)));
    end
    
    % POSTPROCESSING
    % mesh1   = mkmesh_cube(ngrid,ngrid,ngrid,porder+1,1,1,1,elemtype,nodetype);
    % master1 = mkmaster(mesh1,2*(porder+1));
    % [master1,mesh1] = preprocess(master1,mesh1,hybrid);
    % UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);
    
end


% Computing error and convergence orders for vector U and tensor Q
L2errors.Uerror(:,1) = sqrt(L2errors.errors(:,1).^2 + L2errors.errors(:,2).^2 + L2errors.errors(:,3).^2);
L2errors.Uerror(2:end,2) = log(L2errors.Uerror(2:end,1)./L2errors.Uerror(1:end-1,1))...
                        ./(log(L2errors.h(2:end)./L2errors.h(1:end-1))');
L2errors.Qerror(:,1) = sqrt(L2errors.errors(:,4).^2 + L2errors.errors(:,5).^2 + L2errors.errors(:,6).^2 ...
                          + L2errors.errors(:,7).^2 + L2errors.errors(:,8).^2 + L2errors.errors(:,9).^2 ...
                          + L2errors.errors(:,10).^2+ L2errors.errors(:,11).^2+ L2errors.errors(:,12).^2);
L2errors.Qerror(2:end,2) = log(L2errors.Qerror(2:end,1)./L2errors.Qerror(1:end-1,1))...
                        ./(log(L2errors.h(2:end)./L2errors.h(1:end-1))');

% Convergence plot - Velocities
figure(2); clf; loglog(1./L2errors.h,L2errors.Uerror(:,1),'-s'); grid on;
title('Convergence for error on displacements u'); xlabel('1/h'); ylabel('||u-u_h||')
% Convergence plot - Deformation Gradient
figure(3); clf; loglog(1./L2errors.h,L2errors.Qerror(:,1),'-s'); grid on;
title('Convergence for error on the gradient of deformation'); xlabel('1/h'); ylabel('||F-F_h||')

