setapplicationpath('FM/condiff');

porder = 3;
ngrid  = 5;
elemtype = 1;
nodetype = 1;
hybrid = 'hdg';

kappa = 1;
c = [0,0,0]; 
tau = 1;

app.localsolve=1;
app.arg = {kappa,c,tau};
app.bcm = [1;1;1;1;1;1];
app.bcs = [0;0;0;0;0;0]; %[1,0,0;1,0,0];

app.tdep = false;
app.wave = false;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_q = 1;
app.fc_u = 0;
app.fc_p = 0;

app.nd  = 3;
app.nch = 1;                       % Number of componets of UH
app.nc  = 4;    % Number of componeents of UDG
app.ncu = 1;

app.time = [];
app.dtfc = [];
app.alpha = [];

mesh   = mkmesh_cube(ngrid,ngrid,ngrid,porder,1,1,1,elemtype,nodetype);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{0;0;0;0});
%UH = inituhat(master,mesh,app,UDG);
UH=inituhat(master,mesh.elcon,UDG,1);

% HDG solver
tic
%[UDG,UH] = hdg_solve(master,mesh,UDG,UH,[],app);
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
%[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,0*UDG,[]);
toc

x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));
z = (mesh.dgnodes(:,3,:));
u = sin(pi*x).*sin(pi*y).*sin(pi*z);
v = UDG(:,1,:);
max(abs(u(:)-v(:)))

% mesh1   = mkmesh_cube(ngrid,ngrid,ngrid,porder+1,1,1,1,elemtype,nodetype);
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);
% 
% x = (mesh1.dgnodes(:,1,:));
% y = (mesh1.dgnodes(:,2,:));
% z = (mesh1.dgnodes(:,3,:));
% u = sin(pi*x).*sin(pi*y).*sin(pi*z);
% v = UDGstar(:,1,:);
% max(abs(u(:)-v(:)))
% 
