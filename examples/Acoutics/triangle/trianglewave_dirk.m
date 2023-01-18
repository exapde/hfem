setapplicationpath('FM/wave');

porder = 4;
nstage = 3;
torder = 3;
hybrid = 'hdg';

c2 = 1;
k = [6,0];
% dt = 0.1/1e2;
dt = 0.01;
ntime = 20; % explicit number of steps
% ntime = 1000;  % implicit number of steps
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
[M,Minv]=massinv(master, mesh);
% UDG(:,4,:)=PDG; % explicit HDG
mesh.t2t = mkt2t(mesh.t,mesh.elemtype);

% HDG solver
time=0;
%figure(1);
set(gcf,'color','black');
app.torder = 3;
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);

%     [UDG,UH,PDG] = hdg_solve_dirk(master,mesh,app,UDG,UH,PDG,time,dt,nstage,torder);            
    [UDG,UH,PDG] = hdg_wave(master,mesh,app,UDG,PDG,time,dt,nstage,torder);    
    
%     UDGt = cat(2,PDG0(:,:,:,1),UDG(:,2:end,:));
%     Ustar = postprocessnd(master,mesh,master1,mesh1,UDGt);

    figure(1); clf; scaplot(mesh,PDG(:,1,:),[-1 1],2,1); axis off; colormap jet;            
%     figure(2); clf; scaplot(mesh,PDG0(:,:,:,1),[-2.5 2.5],2,1); axis off;              
    %figure(3); clf; scaplot(mesh1,Ustar(:,1,:),[-2.5 2.5],2,1); axis off;        
    
%     app.dt = dt;   
%     app.time = time;   
%     UDG = ehdgrk(master,mesh,app,Minv,UDG);
    time = time + dt; 
%     
%     figure(1); clf; scaplot(mesh,UDG(:,1,:),[-1 1],2,1); axis off; colormap jet;            
%     figure(2); clf; scaplot(mesh,UDG(:,4,:),[-1 1],2,1); axis off; colormap jet;            
end


