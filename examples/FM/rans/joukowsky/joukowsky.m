setapplicationpath('FM/rans')

porder = 2;
nstage = 1;
torder = 1;
hybrid = 'hdg';

ntime = 20;
dt=1e-4*2.^(0:ntime);
%dt=repmat(dt,[2 1]);
dt=dt(:);

gam   = 1.4;
epslm = 0.0;
Minf  = 0.3;                  % Infinity conditions
pinf  = 1/(gam*Minf^2);
alpha = 0;
Re    = 2.8e6; %1.85e6;
Pr    = 0.71;
tau = 2;

ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1), 0.2/Re];

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.arg = {gam,epslm,Re,Pr,Minf,tau};
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

%mesh = mkmesh_naca0012(porder,1,1,1);
%mesh = mkmesh_joukowsky(porder,25,71,2,5);
mesh = mkmesh_joukowsky(porder,51,121,2,8);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);
dist = meshdist(mesh,1);
mesh.dgnodes(:,3,:) = dist;

%UDG = initu(mesh,{ui(1),ui(2),ui(3),ui(4),ui(5); 0,0,0,0,0; 0,0,0,0,0});
%UH = inituhat(master,mesh,app,UDG);
%UDG = dgprojection(master,mesh,UDGp1,1);
%UH = inituhat(master,mesh.elcon,UDG,app.ncu);

time = 0;
for itime = 11:length(dt)    %% WARNING START TO 33
    fprintf('Timestep :  %d,   Time :  %g\n', itime, time+dt(itime));
        
    %function [Un,Hn,Pn] = hdg_solve_dirk(master,mesh,app,UDG,UH,PDG,time,dt,nstage,torder)
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt(itime),nstage,torder);    
    time = time + dt(itime);     
    
    figure(1); clf; scaplot(mesh,UDG(:,2,:)./UDG(:,1,:),[],2,0); 
    axis on; axis equal; axis tight; colorbar; colormap jet; 
    
    figure(2); clf; scaplot(mesh,UDG(:,3,:),[],2,0); 
    axis on; axis square; axis tight; colorbar; colormap jet;
    
    figure(3); clf; scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],2,0); 
    axis on; axis square; axis tight; colorbar; colormap jet;   
    
    figure(4); clf; scaplot(mesh,UDG(:,5,:),[],2,0); 
    axis on; axis square; axis tight; colorbar; colormap jet;
        
    %save naca0012.mat
    fn = ['joukowsky' num2str(itime) '.mat'];
    save(fn,'UDG','UH');
end

fprintf('\nSteady State Solve\n');

app.time = [];
app.fc_q = 1;
app.fc_u = 0;
app.tdep = false;
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
save joukowsky.mat

for k = 1:4
    if k==1
        load('joukowskyp2Re14.mat')
    elseif k==2
        load('joukowskyp2Re28.mat')
    elseif k==3
        load('joukowskyp2Re56.mat')
    elseif k==4
        load('joukowskyp2Re112.mat')
    end        
    figure(k); clf; scaplot(mesh,UDG(:,2,:)./UDG(:,1,:),[],2,0); 
    axis on; axis equal; axis tight; colorbar; colormap jet; axis([-0.7 1.0 -0.5 0.5]);        
end

return;

for k = 1:4
    if k==1
        load('joukowskyp2Re14.mat')
    elseif k==2
        load('joukowskyp2Re28.mat')
    elseif k==3
        load('joukowskyp2Re56.mat')
    elseif k==4
        load('joukowskyp2Re112.mat')
    end        
    vel = UDG(:,2:3,:);
    vel(:,1,:) = UDG(:,2,:)./UDG(:,1,:);
    vel(:,2,:) = UDG(:,3,:)./UDG(:,1,:);
    [udg,pdg,pm,xm,inde] = getfieldinbox(mesh,vel,[0.49 2.52 -0.5 0.5],[1.5 0],10,1);
    s=579.31;
    pdg = s*pdg;
    pm = s*pm;
    udg = permute(udg,[1 3 2]);
    pdg = permute(pdg,[1 3 2]);
    if k==1
        save ransvelr14p2.mat udg pdg pm;
    elseif k==2
        save ransvelr28p2.mat udg pdg pm;
    elseif k==3
        save ransvelr56p2.mat udg pdg pm;
    elseif k==4
        save ransvelr112p2.mat udg pdg pm;
    end            
end



figure(1); clf;
hold on; 
for k = 1:100
plot([pdg(:,1,k); pdg(1,1,k)],[pdg(:,2,k); pdg(1,2,k)],'-b');
plot(pm(k,1),pm(k,2),'or');
end

% 

a = [14 28 56 112];
b = [0.06 0.12 0.24 0.48];
for i=1:length(a)
    load(['joukowskyp2Re' num2str(a(i)) '.mat'])
    figure(1); clf; scaplot(mesh,b(i)*UDG(:,2,:)./UDG(:,1,:),[],2,0); 
    axis on; axis equal; axis tight; colorbar('FontSize',16);
    axis([-1 1.5 -0.8 0.8]); axis off;
    title(['Re = ' num2str(a(i)*1e5) '  ,  '  'U_\infty = ' num2str(b(i))],'FontSize',16);
    colormap(jet);
    print('-dpng',['ux' num2str(i) '.png']);
    
    figure(2); clf; scaplot(mesh,b(i)*UDG(:,3,:)./UDG(:,1,:),[],2,0); 
    axis on; axis equal; axis tight; colorbar('FontSize',16);
    axis([-1 1.5 -0.8 0.8]); axis off;
    title(['Re = ' num2str(a(i)*1e5) '  ,  '  'U_\infty = ' num2str(b(i))],'FontSize',16);
    colormap(jet);
    print('-dpng',['uy' num2str(i) '.png']);
end


[vx,vy] = potentialflow(mesh.dgnodes(:,1,:),mesh.dgnodes(:,2,:));
figure(1); clf; scaplot(mesh,0.25*vx,[],2,0); 
axis on; axis equal; axis tight; colorbar('FontSize',16);
axis([-1 1.5 -0.8 0.8]); axis off;

figure(1); clf; scaplot(mesh,0.25*vy,[-0.84 0.64],2,0); 
axis on; axis equal; axis tight; colorbar('FontSize',16);
axis([-1 1.5 -0.8 0.8]); axis off;

% 
% figure(2); clf; scaplot(mesh,eulereval(UDG(:,1:app.nch,:),'M',gam),[],2,0); 
% axis on; axis square; axis tight; colorbar;    
% 
% figure(3); clf; scaplot(mesh,UDG(:,5,:),[],2,0); 
% axis on; axis square; axis tight; colorbar
% 

figure(1); clf; scaplot(mesh,UDG(:,2,:)./UDG(:,1,:),[],2,0); 
hold on;
plot(pdg(:,:,1),pdg(:,:,2),'o');
axis on; axis equal; axis tight; colorbar('FontSize',16);
axis([-1 1.5 -0.8 0.8]); axis off;

vel = UDG(:,2:3,:);
vel(:,1,:) = UDG(:,2,:)./UDG(:,1,:);
vel(:,2,:) = UDG(:,3,:)./UDG(:,1,:);
[udg,pdg,pm,xm,inde] = getfieldinbox(mesh,vel,[0.49 2.52 -0.5 0.5],[1.5 0],10,1);
figure(1); clf;
meshplot(mesh); 
hold on;
x = pdg(:,:,1); y = pdg(:,:,2);
plot(x(:),y(:),'ob');
plot(pm(:,1),pm(:,2),'or');

v = getvelatx(udg,pdg,pm,[0.75 0]);

s=579.31;
mesh0=mesh;
mesh0.p=s*mesh.p;
figure(1); clf;
meshplot(mesh0); 
hold on;
x = pdg(:,:,1); y = pdg(:,:,2);
plot(x(:)*s,y(:)*s,'ob');
plot(xa,ya);
axis equal;
axis([260 320 -20 20]);

pdg = s*pdg;
pm = s*pm;
