setapplicationpath('SM/LE_uq')

hybrid = 'hdg';
elemtype = 0;
nodetype = 0;
porder = 2;

mu = 1;
lambda = 1;
tau = 1e3;
bodyforce = [0 0 0];
ui = [0, 0, 0];
% Parameters for the imposed wing movement
ttime  = 1/4.;
Fpitch = 1.; Froll = 1.;
Apitch = pi/12.; % Angular Amplitude for the pitch rotation  around y axis
Aroll  = 0.2;   % Displac Amplitude for the roll translation along y axis

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.iterative = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {mu,lambda,tau,ttime,Fpitch,Froll};
app.bcm  = [4,9,1];
app.bcs  = [ui; [0. Apitch Aroll]; ui];
app.bcd  = [1,1,1];
app.bcv  = [ui;ui;ui];

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

load naca3c.mat; 
bndexpr={'p(:,2)<1e-6','sqrt(p(:,1).^2+p(:,3).^2)<3 & p(:,2)<10','true'};
mesh = mkmesh(msh.p',double(msh.t')+1,3,bndexpr,0,0);
ind = [4     7     9    10    13    15    16    18    19    20     3     6     8    12    14    17     2     5    11     1];
mesh.dgnodes = msh.p1(ind,:,:);
% generate new mesh
mesh = higherordermesh(mesh,porder);
master = mkmaster(mesh,2*(porder));
[master,mesh] = preprocess(master,mesh,hybrid);

% plot dimensions of the wing
faceswing = find(mesh.f(:,5)==-2);
nodeswing = mesh.f(faceswing,1:3);
nodeswing = unique(nodeswing(:));
xminwing = min(mesh.p(nodeswing,1)); xmaxwing = max(mesh.p(nodeswing,1));
yminwing = min(mesh.p(nodeswing,2)); ymaxwing = max(mesh.p(nodeswing,2));
zminwing = min(mesh.p(nodeswing,3)); zmaxwing = max(mesh.p(nodeswing,3));
fprintf('Wing dimensions : %d \n',xminwing,xmaxwing,yminwing,ymaxwing,zminwing,zmaxwing);
clear('faceswing','nodeswing');

UDG0 = initu(mesh,{0;0;0;0;0;0;0;0;0;0;0;0});
x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
z = mesh.dgnodes(:,3,:);
UDG0(:,1,:) = x;
UDG0(:,2,:) = y;
UDG0(:,3,:) = z;
UDG0(:,4,:) = -1;
UDG0(:,5,:) = 0;
UDG0(:,6,:) = 0;
UDG0(:,7,:) = 0;
UDG0(:,8,:) = -1;
UDG0(:,9,:) = 0;
UDG0(:,10,:) = 0;
UDG0(:,11,:) = 0;
UDG0(:,12,:) = -1;
UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);

% original configuration
figure(1); clf; meshplot(mesh,1);  
axis equal; axis tight; axis on;
view([0,0]);%view([-75,30]);
set(gca,'FontSize',17);
set(gca, 'LooseInset', get(gca, 'TightInset'));

%%%% MATLAB RESOLUTION %%%%
% fprintf('\nSteady State Solve\n');
% [UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,[]);

% preprocessing for c++ code
app.linearproblem = 1; % Should be equal to one
app.appname='LE_uq';
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


% Number of processors :
nproc = 12;
app.nfile = nproc;

if nproc>1
    % parallel preprocessing
    apppar = digasopre(app,'leuqpar',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,bndexpr,[],nproc,0,check);
    apppar.fileout = 'leuqparout';
    return;
    % run this before using plotsol 
    % mpirun -np 4 digasopar leuqpar leuqparout
    [UDGpar,UHpar] = getsolfrombinaryfile(apppar.fileout,nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
    mesh1=mesh;
    mesh1.dgnodes = UDGpar(:,1:3,:);
else
    % serial preprocessing
    appser = digasopre(app,'leuqser',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,bndexpr,[],nproc,0,check);
    appser.fileout = 'leuqserout';
    % run this before using getsolfrombinaryfile
    % digasoser leuqser leuqserout
    [UDGser,UHser] = getsolfrombinaryfile(appser.fileout,nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
    mesh1=mesh;
    mesh1.dgnodes = UDGser(:,1:3,:);
end

% deformed configuration
figure(3); clf; meshplot(mesh1,1); 
axis equal; axis tight; axis on; box on;
view([0,0]); %view([-75,30]);
set(gca,'FontSize',17);
set(gca, 'LooseInset', get(gca, 'TightInset'));

% Just modifying application with digasopre (does not work so far...)
apppar = digasopre(app,'leuqpar');


% Write outputs for Paraview visu
UDG0(:,1:3,:)   = UDGpar(:,1:3,:)-UDG0(:,1:3,:);
UDG0(:,4:end,:) = UDGpar(:,4:end,:);
writeINPfile('output_12prc', mesh1, UDG0);


