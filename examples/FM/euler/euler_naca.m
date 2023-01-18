setapplicationpath('FM/euler')

hybrid = 'hdg';

m      = 20;
n      = 40;
gridNum = 1;
porder = 3;
torder = 1;
nstage = 1;

gam = 1.4;
epslm = 0.00;
Minf = 0.4;                  % Infinity conditions
pinf = 1/(gam*Minf^2);
alpha = 2*(pi/180);
Re = 1000;  % Irrelevant value
Pr = 0.73; % Irrelevant value
tau = 2;
ui = [1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.iterative=0;
app.hybrid = hybrid;
app.localsolve=0;
app.arg = {gam,epslm,Re,Pr,Minf,tau};

app.uqpk = 0;
app.bcm  = [2,1];  % 2: Slip wall, 1: Far-field
app.bcs  = [ui; ui];

app.bcd  = [1,3];  % 2: Slip wall, 1: Far-field
app.bcv  = [ui; ui];

app.wave = false;
app.tdep = false;
app.alag = false;
app.flg_q = 0;
app.flg_p = 0;
app.flg_g = 0;

app.np = 2;
app.nd = 2;
app.nch = 2+app.nd;                % Number of componets of UH
app.ncq = 0;
app.ncu = app.nch;                        % Number od components with time derivative
app.nc  = app.ncu;                 % Number of componeents of UDG
app.ncp = 0;

app.dtfc = 0;
app.alpha = 0;

app.adjoint = 0;
app.linearproblem = 0;
app.appname = 'euler';
app.linearSolver = 1;
app.jacobianStep = 0;
app.orderingStep = 0;
app.dt = [1e-4,1e-3,1e-2,0.1,1,10,100];
app.torder = torder;
app.nstage = nstage;
app.time = 0;
app.fc_q = 0;
app.fc_u = 0;
app.fc_p = 0;
app.ns   = 1;

mesh = mkmesh_naca0012(porder,1,gridNum);
%mesh = mkmesh_trefftz(m,n,porder,[0.05,0.05,1.98]);
master = mkmaster(mesh,2*porder);
[master,mesh,app] = preprocess(master,mesh,app);

UDG0 = initu(mesh,{ui(1),ui(2),ui(3),ui(4)});
UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);

fprintf('\nSteady State Solve\n');
[UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,[]);

figure(1); clf;
scaplot(mesh,eulereval(UDG(:,1:4,:),'M',1.4),[],3,0); 
hold on;
axis off; axis equal; axis tight;    
colormap('jet');
colorbar('FontSize',14);
axis([-0.2 1.2 -0.6 0.6]);

return;

