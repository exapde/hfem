setapplicationpath('FM/wave');

porder = 3;
nstage = 3;
torder = 3;
hybrid = 'hdg';

c2 = 1;
k = [6,0];
% dt = 0.09;
dt = 0.002;
ntime = 2000; % explicit number of steps
%ntime = 20;  % implicit number of steps
tau = 2;

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.uqpk = false;
app.linear = 1;

app.localsolve=1;
app.arg = {c2,k,tau};
app.bcm = [4;4;4;4;3];
app.bcs = [0;0;0;0;0]; 
app.bcd = app.bcm;
app.bcv = app.bcs;

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

% hdg_wave parameters,,
app.kappa = 1;
app.tau = tau;
app.k = k;

% Stage parameter
app.nstage = nstage;
app.torder = torder;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mesh = mkmesh_circleinsquare(porder,1,1);
mesh=mkmesh_triangleinrect(porder);
master = mkmaster(mesh,2*porder);

%mesh1 = mkmesh_circleinsquare(porder+1,1,1);
mesh1=mkmesh_triangleinrect(porder+1);
master1 = mkmaster(mesh1,2*(porder+1));

[master,mesh] = preprocess(master,mesh,hybrid);
[master1,mesh1] = preprocess(master1,mesh1,hybrid);

UDG = initu(mesh,{0;0;0});
UH = inituhat(master,mesh.elcon,UDG,app.ncu);
PDG = initu(mesh,{0});

% Explicit Code Precomputation
[M,Minv]=massinv(master, mesh);
UDG(:,4,:)=PDG; % explicit HDG
mesh.t2t = mkt2t(mesh.t,mesh.elemtype);

% Determine if an element is explicit or implicit
mesh = create_mesh_imex(mesh,dt);

% Find implicit/explicit boundary
mesh = find_imex_bound(mesh);

% Find implicit faces
mesh = create_mesh_imface(mesh);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HDG solver
time=0;
%figure(1);
set(gcf,'color','black');
app.torder = 3;


  
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
    
    % Implicit Code
    UDG = hdg_wave_imex_coupled_final(master,mesh,app,UDG(:,1:3,:),UDG(:,4,:),time,dt,nstage,torder,Minv);   

    %UDG = cat(2,UDGimp,PDG); % Fixes notation differences between explicit and implicit codes
    
    figure(3); clf; scaplot(mesh,UDG(:,4,:),[-1 1],2,1); axis off; colormap jet; 
    
    time = time + dt;
    
    
    
          
%     pause(0.01)
end


