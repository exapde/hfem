setapplicationpath('FM/poi');

% nproc = 1;
fileName = 'poi2d';

porder = 4;
ngrid  = 9;
elemtype = 1;
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
app.bcm = [1;1;1;1];
app.bcs = [0;0;0;0];
app.bcd  = [1;1;1;1];  % [nothing,nothing,inlet,outlet,nothing,nothing,airfoil]
app.bcv  = [0;0;0;0];
app.wave = false;
app.tdep = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_u = 1;
app.fc_q = 1;
app.fc_p = 0;
app.nd   = 2;
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
app.linear = 1;
app.appname = 'poisson';
app.linearSolver = 1;
app.jacobianStep = 0;
app.orderingStep = 0;
app.dt = 0;
% app.torder = torder;
% app.nstage = nstage;

mesh   = mkmesh_square(ngrid,ngrid,porder,0,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*porder);
[master,mesh,app] = preprocess(master,mesh,app);

UDG0 = initu(mesh,{0;0;0});
UH0 = inituhat(master,mesh.elcon,UDG0,1);

[UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,0*UDG0);

x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));
u = sin(pi*x).*sin(pi*y);
v = UDG(:,1,:);
max(abs(u(:)-v(:)))

% preprocessing for c++ code
app.linearproblem = 1;
app.appname='poisson';
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
bndexpr=mesh.bndexpr;

% % parallel preprocessing
% nproc = 4;
% apppar = digasopre(app,'ns',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,bndexpr,[],nproc,0,check);
% apppar.fileout = 'nsout';

% run this before using plotsol 
% mpirun -np 4 digasopar ns nsout
 
% serial preprocessing
nproc = 1;
appser = digasopre(app,'poiser',mesh.p,mesh.t'-1,mesh.dgnodes,UDG0,UH0,[],elementtype,bndexpr,[],nproc,0,check);
appser.fileout = 'poiseroutsol';

% run this before using getsolfrombinaryfile 
!./../../../main/digasoser poiser poiserout
[UDGser,UHser] = getsolfrombinaryfile(appser.fileout,appser.nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);
UDGser = reshape(UDGser,size(UDG));

v = UDGser(:,1,:);
max(abs(u(:)-v(:)))

% deformed configuration
% mesh1=mesh;
% mesh1.dgnodes = UDGser(:,1:3,:);
% figure(3); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;



%writeBinaryFile(fileName,mesh,master,app,UDG,UH,[],1);
% app.flag(5) = 2;
% app.factor(5) = 1e-2;
% writeBinaryFile(fileName,mesh,master,app,UDG,UH,[],1);
% 
% 
% [RU,RH,RQ] = hdg_residual(master,app,mesh.dgnodes,mesh.bf,UDG,UH,0*UDG);
% 
% 
% filename = 'solpoi2d';
% fileID = fopen([filename,'.bin'],'r');
% data = fread(fileID,'double');
% fclose(fileID);
% 
% N = master.npv*(mesh.nd+1)*mesh.ne;
% UDG = reshape(data(1:N),[master.npv (mesh.nd+1) mesh.ne]);
% 
% x = (mesh.dgnodes(:,1,:));
% y = (mesh.dgnodes(:,2,:));
% u = sin(pi*x).*sin(pi*y);
% v = UDG(:,1,:);
% max(abs(u(:)-v(:)))
% 
% figure(1); scaplot(mesh,UDG(:,3,:),[],0,1); axis equal; axis tight;
% 
