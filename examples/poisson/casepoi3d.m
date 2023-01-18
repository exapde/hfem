setapplicationpath('FM/poi');

% nproc = 1;
fileName = 'poi3d';

porder = 2;
ngrid  = 17;
elemtype = 0;
nodetype = 1;
nstage = 0;
torder = 0;
hybrid = 'hdg';

kappa = 1;
tau = 1;

app.iterative = 0;
app.getdqdg = 1;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {kappa,tau};
app.bcm = [1;1;1;1;1;1];
app.bcs = [0;0;0;0;0;0];
app.bcd  = [1;1;1;1;1;1];  % [nothing,nothing,inlet,outlet,nothing,nothing,airfoil]
app.bcv  = [0;0;0;0;0;0];
app.wave = false;
app.tdep = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_u = 1;
app.fc_q = 1;
app.fc_p = 0;
app.nd   = 3;
app.ncu  = 1;               % Number of components of U
app.nch  = app.ncu;                % Number of componets of UH
app.ncq  = app.ncu*app.nd;         % Number of componets of Q
app.nc   = app.ncu+app.ncq;         % Number of componeents of UDG
app.ncp = 0; 
app.time = 0.0;
app.dtfc = 0;
app.alpha = 0;
app.ns   = 1;  
app.adjoint = 0;
app.linearproblem = 1;
app.appname = 'poisson';
app.linearSolver = 1;
app.jacobianStep = 0;
app.orderingStep = 0;
app.dt = 0;
app.torder = torder;
app.nstage = nstage;

mesh   = mkmesh_cube(ngrid,ngrid,ngrid,porder,1,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*porder);
[master,mesh,app] = preprocess(master,mesh,app);

UDG = initu(mesh,{0;0;0;0});
UH=inituhat(master,mesh.elcon,UDG,1);

return;

% writeBinaryFile(fileName,mesh,master,app,UDG,UH,[],1);
% 

filename = 'solpoi3d';
fileID = fopen([filename,'.bin'],'r');
data = fread(fileID,'double');
fclose(fileID);

N = master.npv*(mesh.nd+1)*mesh.ne;
UDG = reshape(data(1:N),[master.npv (mesh.nd+1) mesh.ne]);

x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));
z = (mesh.dgnodes(:,3,:));
u = sin(pi*x).*sin(pi*y).*sin(pi*z);
v = UDG(:,1,:);
max(abs(u(:)-v(:)))

u = pi*cos(pi*x).*sin(pi*y).*sin(pi*z);
v = -UDG(:,2,:);
max(abs(u(:)-v(:)))

u = pi*sin(pi*x).*cos(pi*y).*sin(pi*z);
v = -UDG(:,3,:);
max(abs(u(:)-v(:)))

u = pi*sin(pi*x).*sin(pi*y).*cos(pi*z);
v = -UDG(:,4,:);
max(abs(u(:)-v(:)))

