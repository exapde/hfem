setapplicationpath('SM/LE_uq')

% Elements type
hybrid = 'hdg';
elemtype = 1;
porder = 2;

% Mesh parameters
thickness = 0.03;
Ri = 6.;
R0 = 10.;
n  = 51;
nr = 9;
nt = 2;
nodetype = 0;
Sur = thickness*(R0-Ri);

% Mechanical parameters
Eyoung = 21.e6; nu = 0.;
lambda = Eyoung*nu/((1+nu)*(1-2*nu));
mu0 = Eyoung/(2*(1+nu));
lambda = lambda/mu0;
mu = 1;
tau = 1;

% Applied Total Load (traction)
nloads = 10;
applied_f = 3.2;
applied_trac = applied_f/Sur;
applied_trac = applied_trac/mu0;


app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.iterative = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {mu,lambda,tau};

% Boundary conditions
bodyforce = [0 0 0];
ui = [0, 0, 0];
app.bcm  = [3,3,3,1,3,3];  
app.bcs  = [ui; ui; [0, 0, applied_trac]; ui; ui; ui];
app.bcd  = [1];  
app.bcv  = [ui];

% Application flags
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
app.nch  = app.ncu;                % Number of componets of UH
app.ncq  = app.nch*app.nd;         % Number of componets of Q
app.nc   = app.ncu+app.ncq;         % Number of componeents of UDG
app.ncp = 0;

app.time = [];
app.dtfc = [];
app.alpha = [];

% Mesh generation
mesh = mkmesh_annular(n,nr,nt,R0,Ri,thickness,porder,elemtype);

master = mkmaster(mesh,2*(porder));
[master,mesh] = preprocess(master,mesh,hybrid);

% Initialization of results sets
res.WB1 = zeros(nloads,3); res.WB2 = zeros(nloads,3);
res.WA1 = zeros(nloads,3); res.WA2 = zeros(nloads,3);
res.Nit = zeros(nloads,1);

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

fprintf('\nSteady State Solve\n');

percent = 0.05;
UDG = UDG0;
UH  = UH0;
for iload = 1:1
    fprintf('########################################### \n'); 
    fprintf('Load Increment :  %d\n', iload); 
    fprintf('########################################### \n'); 
    app.bcs  = [ui; ui; [0, 0, applied_trac*percent]; ui; ui; ui];
    [UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
    % Height of end points of the slit
    %res.WA1(iload,:) = UDG(1,1:3,1);
    %res.WA2(iload,:) = UDG(5,1:3,1);
    %res.WB1(iload,:) = UDG(2,1:3,nr-1);
    %res.WB2(iload,:) = UDG(6,1:3,nr-1);
    res.WA1(iload,:) = UDG(1,1:3,1);
    res.WA2(iload,:) = UDG(19,1:3,1);
    res.WB1(iload,:) = UDG(3,1:3,nr-1);
    res.WB2(iload,:) = UDG(21,1:3,nr-1);
end

% original configuration
figure(1); clf; meshplot(mesh,1);  axis equal; axis tight; axis on;

% deformed configuration
mesh1=mesh;
mesh1.dgnodes(:,1:3,:)=UDG(:,1:3,:);
figure(3); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;
% figure(3); scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% figure(4); scaplot(mesh,UDG(:,2,:),[],2,1); axis equal; axis tight;

