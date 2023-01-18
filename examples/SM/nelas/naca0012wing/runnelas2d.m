setapplicationpath('SM/NeoHookean_uq')

hybrid = 'hdg';
elemtype = 1;
nodetype = 0;
ngrid = 10;
porder=5;

mu = 1;
lambda = 1;
tau = 1e5;
ui = [0, 0];

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.iterative = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {mu,lambda,tau,[0,0e-3]};
app.bcm  = [1,2];  
app.bcs  = [ui; [0 5e-4]];
app.bcd  = [1,2]; 
app.bcv  = [ui;ui];

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

app.nd   = 2;
app.ncu  = 2;               % Number of components of U
app.nch  = app.ncu;                % Number of componets of UH
app.ncq  = app.nch*app.nd;         % Number of componets of Q
app.nc   = 6;         % Number of componeents of UDG
app.ncp = 0;

app.time = [];
app.dtfc = [];
app.alpha = [];

thick = 12;
th = (pi:-pi/200:pi/2)';
xt = (cos(th)+1)*1.0089304129;  
xt = xt(end:-1:1);
yt=naca(xt,thick);  
xb = flipud(xt);   
yb=-naca(xb,thick);
yt(1) = 0; yb(end) = 0;
xf =[xt; xb(2:end)];
yf =[yt; yb(2:end)];
xf(end) = xf(1);
yf(end) = yf(1);

chord = max(xf);
xf = xf/chord;
yf = yf/chord;
mesh = mkmesh_naca0012foil( xf, yf, porder, 50, 0.01);
master = mkmaster(mesh,2*(porder));
[master,mesh] = preprocess(master,mesh,hybrid);

UDG0 = initu(mesh,{0;0;0;0;0;0});
x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
UDG0(:,1,:) = x;
UDG0(:,2,:) = y;
UDG0(:,3,:) = -ones(size(x));
UDG0(:,4,:) = 0*x;
UDG0(:,5,:) = 0*x;
UDG0(:,6,:) = -ones(size(x));
UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);

fprintf('\nSteady State Solve\n');
[UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,[]);

% original configuration
figure(1); clf; meshplot(mesh,1); axis equal; axis tight; axis on;

% deformed configuration
mesh1=mesh;
mesh1.dgnodes(:,1,:)=UDG(:,1,:);
mesh1.dgnodes(:,2,:)=UDG(:,2,:);
figure(2); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;
% figure(3); scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% figure(4); scaplot(mesh,UDG(:,2,:),[],2,1); axis equal; axis tight;


