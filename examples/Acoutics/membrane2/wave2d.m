setapplicationpath('FM/wave');

porder = 4;
ngrid  = 5;
nstage = 2;
torder = 3;
hybrid = 'hdg';

c2 = 1;
k = [0,0];
dt = 0.01;
ntime = 20;

app.localsolve=1;
app.arg = {c2,k};
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

%mesh = mkmesh_circleinsquare(porder,1,1);
mesh = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

UDG = initu(mesh,{0;0;0});
%UH = inituhat(master,mesh.elcon,UDG,app.ncu);
PDG = initu(mesh,{0});

x=mesh.dgnodes(:,1,:);
y=mesh.dgnodes(:,2,:);
u=(1/(sqrt(2)*pi))*sin(pi*x).*sin(pi*y).*sin(sqrt(2)*pi*0);    
v=sin(pi*x).*sin(pi*y).*cos(sqrt(2)*pi*0);
qx=-(1/sqrt(2))*cos(pi*x).*sin(pi*y).*sin(sqrt(2)*pi*0);    
qy=-(1/sqrt(2))*sin(pi*x).*cos(pi*y).*sin(sqrt(2)*pi*0);    
UDG(:,1,:)=v;
UDG(:,2,:)=qx;
UDG(:,3,:)=qy;
UH = inituhat(master,mesh.elcon,UDG,app.ncu);
PDG(:,1,:)=u;

% HDG solver
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
        
    [UDG,UH,~,PDG] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],PDG,(itime-1)*dt,dt,nstage,torder);            
    v=sin(pi*x).*sin(pi*y).*cos(sqrt(2)*pi*itime*dt);
    figure(1); clf; scaplot(mesh,UDG(:,1,:),[],2,1); axis off;         
    figure(2); clf; scaplot(mesh,v-UDG(:,1,:),[],2,1); axis off;  
        
    u=(1/(sqrt(2)*pi))*sin(pi*x).*sin(pi*y).*sin(sqrt(2)*pi*itime*dt);    
    figure(3); clf; scaplot(mesh,PDG(:,1,:),[],2,1); axis off;          
    figure(4); clf; scaplot(mesh,u-PDG(:,1,:),[],2,1); axis off;   
    ev = v-UDG(:,1,:);
    eu = u-PDG(:,1,:);
    [itime max(abs(ev(:))) max(abs(eu(:)))]
    pause
end
