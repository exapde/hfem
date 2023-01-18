setapplicationpath('FM/wave');

porder = 2;
ngrid  = 9;
nstage =porder+1;
torder =porder+2;
hybrid = 'hdg';

dt = (1/(ngrid-1))/(2*(2*porder+1));
ntime = round(1/dt);

app.tau = 1;
app.kappa = 1;
app.localsolve=1;
app.arg = {1,1};
app.bcm = [1;1;1;1];
app.bcs = [0;0;0;0]; 

app.hybrid=hybrid;
app.iterative=0;
app.tdep = true;
app.wave = true;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.nd   = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu  = 1;                       % Number of components of U
app.itmax = 1;

mesh = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

t = 0;
UDG = l2eprojection(mesh,master,@exactsol,[],0,4);
[Minv] = massinv(master, mesh);

time = 0;
app.time = time;
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
        
    UDG = ehdgrk(master,mesh,app,Minv,UDG);    
    time = time + dt; 
    app.time = time;    
end
% erru = calerror(UDG,mesh,master,@exactsol,itime*dt);

setapplicationpath('FM/wave');

porder = 2;
ngrid  = 9;
nstage =porder+1;
torder =porder+2;
hybrid = 'hdg';

dt = (1/(ngrid-1))/(2*(2*porder+1));
ntime = round(1/dt);

app.tau = 1;
app.kappa = 1;
app.localsolve=1;
app.arg = {1,1};
app.bcm = [1;1;1;1];
app.bcs = [0;0;0;0]; 

app.hybrid=hybrid;
app.iterative=0;
app.tdep = true;
app.wave = true;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.nd   = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu  = 1;                       % Number of components of U
app.itmax = 1;

mesh = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

t = 0;
UDG = l2eprojection(mesh,master,@exactsol,[],0,4);
[~Minv=massinv(master, mesh);

time = 0;
app.time = time;
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
        
    UDG = ehdgrk(master,mesh,app,Minv,UDG);    
    time = time + dt; 
    app.time = time;    
end
erru = calerror(UDG,mesh,master,@exactsol,itime*dt);
setapplicationpath('FM/wave');

porder = 2;
ngrid  = 9;
nstage =porder+1;
torder =porder+2;
hybrid = 'hdg';

dt = (1/(ngrid-1))/(2*(2*porder+1));
ntime = round(1/dt);

app.tau = 1;
app.kappa = 1;
app.localsolve=1;
app.arg = {1,1};
app.bcm = [1;1;1;1];
app.bcs = [0;0;0;0]; 

app.hybrid=hybrid;
app.iterative=0;
app.tdep = true;
app.wave = true;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.nd   = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu  = 1;                       % Number of components of U
app.itmax = 1;

mesh = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

t = 0;
UDG = l2eprojection(mesh,master,@exactsol,[],0,4);
[~Minv=massinv(master, mesh);

time = 0;
app.time = time;
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
        
    UDG = ehdgrk(master,mesh,app,Minv,UDG);    
    time = time + dt; 
    app.time = time;    
end
erru = calerror(UDG,mesh,master,@exactsol,itime*dt);
setapplicationpath('FM/wave');

porder = 2;
ngrid  = 9;
nstage =porder+1;
torder =porder+2;
hybrid = 'hdg';

dt = (1/(ngrid-1))/(2*(2*porder+1));
ntime = round(1/dt);

app.tau = 1;
app.kappa = 1;
app.localsolve=1;
app.arg = {1,1};
app.bcm = [1;1;1;1];
app.bcs = [0;0;0;0]; 

app.hybrid=hybrid;
app.iterative=0;
app.tdep = true;
app.wave = true;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.nd   = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu  = 1;                       % Number of components of U
app.itmax = 1;

mesh = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

t = 0;
UDG = l2eprojection(mesh,master,@exactsol,[],0,4);
[~Minv=massinv(master, mesh);

time = 0;
app.time = time;
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
        
    UDG = ehdgrk(master,mesh,app,Minv,UDG);    
    time = time + dt; 
    app.time = time;    
end
erru = calerror(UDG,mesh,master,@exactsol,itime*dt);



