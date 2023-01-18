setapplicationpath('FM/rans')

porder = 3;
nstage = 1;
torder = 1;
hybrid = 'hdg';

ntime = 20;
dt=0.5e-5*2.^(0:ntime);
dt=repmat(dt,[2 1]);
dt=dt(:);

gam   = 1.4;
epslm = 0.0;
Minf  = 0.3;                  % Infinity conditions
pinf  = 1/(gam*Minf^2);
alpha = 0;
Re    = 1.85e6;
Pr    = 0.72;

ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1), 0.2/Re];

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.arg = {gam,epslm,Re,Pr,Minf};
app.bcm  = [2,1];  
app.bcs  = [ui; ui];
app.bcd  = [2,1];  % 2: Wall, 1: Far-field
app.bcv  = [ui;ui];

app.uqpk=0;
app.localsolve = 1;
app.hybrid = hybrid;
app.iterative = 0;
app.wave = false;
app.tdep = true;
app.alag = false;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_q = 1;
app.fc_u = 0;
app.fc_p = 0;

app.nd   = 2;
app.ncu  = 3+app.nd;               % Number of components of U
app.nch  = app.ncu;                % Number of componets of UH
app.nq   = app.nch*app.nd;         % Number of componets of Q
app.nc   = app.nch+app.nq;         % Number of componeents of UDG

mesh = mkmesh_naca0012(porder,1,1,1);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
dist = meshdist(mesh,1);
mesh.dgnodes(:,3,:) = dist;

UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),ui(5); 0,0,0,0,0; 0,0,0,0,0});
%UH = inituhat(master,mesh,app,UDG);
UH=inituhat(master,mesh.elcon,UDG,app.ncu);

size(UDG)
time = 0;
for itime = 33:length(dt)    %% WARNING START TO 33
    fprintf('Timestep :  %d,   Time :  %g\n', itime, time+dt(itime));
        
    %function [Un,Hn,Pn] = hdg_solve_dirk(master,mesh,app,UDG,UH,PDG,time,dt,nstage,torder)
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);    
    time = time + dt(itime);     
    
    figure(1); clf; scaplot(mesh,UDG(:,2,:),[],2,0); 
    axis on; axis square; axis tight; colorbar; colormap jet;
    
    figure(2); clf; scaplot(mesh,UDG(:,3,:),[],2,0); 
    axis on; axis square; axis tight; colorbar; colormap jet;
    
    figure(3); clf; scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],2,0); 
    axis on; axis square; axis tight; colorbar; colormap jet;   
    
    figure(4); clf; scaplot(mesh,UDG(:,5,:),[],2,0); 
    axis on; axis square; axis tight; colorbar; colormap jet;
        
    %save naca0012.mat
    fn = ['naca0012' num2str(itime) '.mat'];
    save(fn,'UDG','UH');
end

fprintf('\nSteady State Solve\n');

app.time = [];
app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
save naca0012.mat
% 
% figure(1); clf; scaplot(mesh,UDG(:,1,:),[],2,0); 
% axis on; axis square; axis tight; colorbar;
% 
% figure(2); clf; scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],2,0); 
% axis on; axis square; axis tight; colorbar;    
% 
% figure(3); clf; scaplot(mesh,UDG(:,5,:),[],2,0); 
% axis on; axis square; axis tight; colorbar
% 
% 
% 
