setapplicationpath('FM/wave');

porder = 2;
% ngrid  = 9;
nstage = 3;
torder = 3;
hybrid = 'hdg';
nref = 2;

dt = 0.0001;
ntime = 100;

k = [1, pi/6];
app.k = k;

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

app.nstage = nstage;
app.torder = torder;

mesh = mkmesh_square_imex(porder);
mesh = refinemesh(mesh,porder,nref);
master = mkmaster(mesh,2*porder);

[master,mesh] = preprocess(master,mesh,hybrid);

t = 0;
UDG = l2eprojection(mesh,master,@exactsol,[],0,4);
Minv=massinv(master, mesh);
mesh.t2t = mkt2t(mesh.t,mesh.elemtype);
% Determine if an element is explicit or implicit
mesh = create_mesh_imex(mesh,dt);

% Find implicit/explicit boundary
mesh = find_imex_bound(mesh);

% Find implicit faces
mesh = create_mesh_imface(mesh);
ind = find(mesh.imex == 1); % Implicit element indicies
impE = length(ind)

dt = 0.0001;
time = 0;
app.dt = dt;
app.torder = torder;
app.time = time;
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
      UDG = hdg_vibration_imex_coupled_UPDATED_v4(master,mesh,app,UDG(:,1:3,:),UDG(:,4,:),time,dt,nstage,torder,Minv);
%     UDG = ehdgrk_test(master,mesh,app,Minv,UDG);    
    time = time + dt; 
    app.time = time;    
end
erru = calerror(UDG,mesh,master,@exactsol,[],itime*dt,[1 2 3 4]);

uex = exactsol(mesh.dgnodes, [], itime*dt);
figure(1); clf; scaplot(mesh,uex(:,4,:)-UDG(:,4,:),[],2,1); axis off;  
figure(2); clf; scaplot(mesh,uex(:,4,:),[],2,1); axis off;
figure(3); clf; scaplot(mesh,UDG(:,4,:),[],2,1); axis off;

% figure(1); clf; scaplot(mesh,UDG(:,4,:),[-1 1],2,1); axis off; colormap jet;
