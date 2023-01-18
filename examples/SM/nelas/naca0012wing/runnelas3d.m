%setapplicationpath('SM/nelasuq')
setapplicationpath('SM/NeoHookean_uq')
%setapplicationpath('SM/SVK_uq')

hybrid = 'hdg';
porder = 2;
foilthickness = 0.01;
n = 20;
nelem2d = n;
zz = linspace(0,4,n); 

mu = 1;
lambda = 1; % mu/lambda
tau = 1e6;
bodyforce = [0 0 0];
ui = [0, 0, 0];

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.iterative = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {mu,lambda,tau,[0,0,0]};
app.bcm  = [1,2,2,2];  
app.bcs  = [ui; [0 4e-4 0]; ui; ui];
app.bcd  = [1,2,2,2]; 
app.bcv  = [ui;ui;ui;ui];

app.wave = false;
app.tdep = false;
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

mesh = mkmesh_naca0012wing(porder,nelem2d,foilthickness,zz);
master = mkmaster(mesh,2*(porder));
[master,mesh] = preprocess(master,mesh,hybrid);

fprintf('\nSteady State Solve\n');
[UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,[]);
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

% original configuration
%figure(1); clf; meshplot(mesh,1);  axis equal; axis tight; axis on;

% deformed configuration
mesh1=mesh;
mesh1.dgnodes(:,1,:)=UDG(:,1,:);
mesh1.dgnodes(:,2,:)=UDG(:,2,:);
mesh1.dgnodes(:,3,:)=UDG(:,3,:);
figure(3); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;
% figure(3); scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% figure(4); scaplot(mesh,UDG(:,2,:),[],2,1); axis equal; axis tight;

