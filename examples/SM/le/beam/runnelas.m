setapplicationpath('SM/LE_uq')

hybrid = 'hdg';
elemtype = 1;
nodetype = 1;
ngrid = 9;
porder= 2;

mu = 1;
lambda = 1;
tau = 2;
ui = [0, 0, 0];

app.source = 'source1';
app.flux = 'flux';
app.fbou = 'fbou1';
app.fhat = 'fhat';

app.iterative = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {mu,lambda,tau};
app.bcm  = [1,2,2,2,2,2];  
app.bcs  = [ui; [0 0 -0.001]; ui; ui; ui; ui];
app.bcd  = [1,1,1,1,1,1];  
app.bcv  = [ui];

app.wave = false;
app.tdep = false;
app.alag = false;
app.uqpk = false;

app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_q = 1;
app.fc_u = 1;
app.fc_p = 0;

app.nd   = 3;
app.ncu  = 3;               % Number of components of U
app.nch  = 3;                % Number of componets of UH
app.ncq  = 9;         % Number of componets of Q
app.nc   = 12;         % Number of componeents of UDG
app.ncp = 0;

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh   = mkmesh_cube(ngrid,ngrid,ngrid,porder,5,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*(porder));
[master,mesh] = preprocess(master,mesh,hybrid);

UDG0 = initu(mesh,{0;0;0;0;0;0;0;0;0;0;0;0});
% x = mesh.dgnodes(:,1,:);
% y = mesh.dgnodes(:,2,:);
% z = mesh.dgnodes(:,3,:);
% UDG0(:,1,:) = x;
% UDG0(:,2,:) = y;
% UDG0(:,3,:) = z;
% UDG0(:,4,:) = -1;
% UDG0(:,5,:) = 0;
% UDG0(:,6,:) = 0;
% UDG0(:,7,:) = 0;
% UDG0(:,8,:) = -1;
% UDG0(:,9,:) = 0;
% UDG0(:,10,:) = 0;
% UDG0(:,11,:) = 0;
% UDG0(:,12,:) = -1;
UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);

fprintf('\nSteady State Solve\n');
app.linear = 1;
[UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,[]);

% original configuration
figure(1); clf; meshplot(mesh,1); axis equal; axis tight; axis on;

% deformed configuration
mesh1=mesh;
mesh1.dgnodes = UDG(:,1:3,:);
figure(2); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;


% preprocessing for c++ code
app.linearproblem = 1;
app.appname='leuq';
app.overlappinglevel=1;
app.preconditioner=0;
app.morder = [porder porder porder];
app.porder = [porder porder porder];
app.nodetype = nodetype;
app.pgauss = 2*[porder porder porder];
app.pgaussR = 2*[porder porder porder];
app.quadtype = [0 0 0];
app.dt = 0;
app.nco = 0;
app.ncd = size(mesh.dgnodes,2);
check = 0;
app.nfile=1;
elementtype = elemtype*ones(mesh.ne,1);
bndexpr=mesh.bndexpr;
% GMRES Tolerance
app.restart=200;
app.gmrestol=1e-7;
app.gmresiter = 2000;

% Number of processors :
nproc = 1;
app.nfile = nproc;

if nproc>1
    % parallel preprocessing
    apppar = digasopre(app,'lepar',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,bndexpr,[],nproc,0,check);
    apppar.fileout = 'leparout';
    % run this before using getsolfrombinaryfile
    % mpirun -np 4 digasopar lepar leparout
    [UDGpar,UHpar] = getsolfrombinaryfile(apppar.fileout,nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
    mesh1=mesh;
    mesh1.dgnodes = UDGpar(:,1:3,:);
else
    % serial preprocessing
    appser = digasopre(app,'leuqser',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,bndexpr,[],nproc,0,check);
    appser.fileout = 'leuqseroutsol';

    % run this before using getsolfrombinaryfile
    % digasoser leuqser leuqserout
    [UDGser,UHser] = getsolfrombinaryfile(appser.fileout,nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
    mesh1=mesh;
    mesh1.dgnodes = UDGser(:,1:3,:);
end

% deformed configuration
figure(3); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;
