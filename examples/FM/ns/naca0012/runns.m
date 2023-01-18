
setapplicationpath('FM/ns')


hybrid = 'hdg';
elemtype = 0;
nodetype = 0;
gridNum=2;
porder=2;

gam = 1.4;
epslm = 0.05;
Minf = 0.5;                  % Infinity conditions
pinf = 1/(gam*Minf^2);
alpha = 2.0*(pi/180);
Re = 500;
Pr = 0.72;
tau = 2;
ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];

dt = [1e-3 1e-2 1e-1 1e0 1e1];

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.iterative = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {gam,epslm,Re,Pr,Minf,tau};
app.bcm  = [2,4,1];  % 2: Wall, 1: Far-field
app.bcs  = [ui; ui; ui];
% app.bcm  = [2,1];
% app.bcs  = [ui; ui];

app.bcd  = [2,4,1];  % 2: Wall, 1: Far-field
app.bcv  = [ui;ui;ui];
% app.bcd  = [1,3];  % 2: Wall, 1: Far-field
% app.bcv  = [ui;ui];

app.wave = false;
app.tdep = false;
app.alag = false;

app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_q = 1;
app.fc_u = 1;
app.fc_p = 0;

app.nd   = 2;
app.ncu  = 4;               % Number of components of U
app.nch  = app.ncu;                % Number of componets of UH
app.ncq  = app.nch*app.nd;         % Number of componets of Q
app.nc   = 12;         % Number of componeents of UDG
app.ncp = 0;

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh = mkmesh_naca0012(porder,elemtype,nodetype,gridNum);
master = mkmaster(mesh,2*(porder));
[master,mesh] = preprocess(master,mesh,hybrid);

UDG0 = initu(mesh,{ui(1),ui(2),ui(3),ui(4); 0,0,0,0; 0,0,0,0});
UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);

fprintf('\nSteady State Solve\n');
app.uqpk=0;
app.time = [];
app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
[UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,[]);

%app.tdep = true;
%app.dt = 0.01*ones(100,1);

close all;scaplot(mesh,UDG(:,2,:),[],1,0); colormap(jet);
%return;

% preprocessing for c++ code
app.appname='ns';
app.overlappinglevel=1;
app.preconditioner=0;
app.morder = [porder porder];
app.porder = [porder porder];
app.nodetype = nodetype;
app.pgauss = 2*[porder porder];
app.pgaussR = 2*[porder porder];
app.quadtype = [0 0 0];
app.dt = 0;
app.nco = 0;
app.ncd = size(mesh.dgnodes,2);
check = 0;
app.nfile=1;
elementtype = elemtype*ones(mesh.ne,1);
bndexpr={'sqrt((p(:,1)-.5).^2+p(:,2).^2)<2','p(:,1)>8','true'};

% % parallel preprocessing
nproc = 4;
apppar = digasopre(app,'ns',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,bndexpr,[],nproc,0,check);
apppar.fileout = 'nsout';

% run this before using plotsol 
% mpirun -np 4 ./../../main/digasopar ns nsout
 
% serial preprocessing
% nproc = 1;
% appser = digasopre(app,'nsser',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,bndexpr,[],nproc,0,check);
% appser.fileout = 'nsserout';

% run this before using plotsol 
% digasoser nsser nsserout

