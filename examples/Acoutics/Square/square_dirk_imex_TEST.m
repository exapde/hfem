setapplicationpath('FM/wave');


%%%%%% Loop Code
% erruimp = zeros(4,2,3);
% count = 1;
% 
% for i = 1:2
%     for j = 1:3
%         clear mesh
%         clear UDG
%         clear PDG
%         
%     porder = i;
% nrefset = [2;4;8];
% nref = nrefset(j);
% dtset = [0.06; 0.03; 0.015; 0.015; 0.0075; 0.00375];
% dt = dtset(count);

porder = 2;
nstage = 3;
torder = 3;
hybrid = 'hdg';
nref = 1;


dt = 0.03;
c2 = 1;
k = [1, pi/6];
% dt = 0.12; % dt to create mesh when no nref, p=1
% dt = 0.06; %dt when nref = 2
% dt = 0.03; % dt when nref = 4
% dt = 0.0015;
% dt = 0.0075;
% dt = 0.00325;
% dt = 0.015;
ntime = 50;
tau = 1;
c = sqrt(c2);
alpha = k(1);
phi = k(2);

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';
app.uqpk = false;
app.linear = 1;

app.localsolve=1;
app.arg = {c2,k,tau};
app.arg = {1,1};
app.bcm = [1;1;1;1];
app.bcs = [0;0;0;0]; 
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

% vIC = @(dg) sin(pi*dg(:,1)).*sin(pi*dg(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh=mkmesh_square_imex(porder);
mesh = refinemesh(mesh,porder,nref);
master = mkmaster(mesh,2*porder);

mesh1=mkmesh_square_imex(porder+1);
mesh1 = refinemesh(mesh1,porder,nref);
master1 = mkmaster(mesh1,2*(porder+1));

[master,mesh] = preprocess(master,mesh,hybrid);

[master1,mesh1] = preprocess(master1,mesh1,hybrid);


UDG = l2eprojection(mesh,master,@exactsol,[c alpha phi],0,4);
% UDG = l2eprojection(mesh,master,@exactsol,[],0,4);
PDG(:,1,:)=UDG(:,4,:);
UDG(:,4,:)=[];
 
UH = inituhat(master,mesh.elcon,UDG,app.ncu);



% Explicit Code Precomputation
% [M, Minv]=massinv(master, mesh);
[Minv]=massinv(master, mesh);
UDG(:,4,:)=PDG; % explicit HDG
mesh.t2t = mkt2t(mesh.t,mesh.elemtype);

%%%%%%%%%%%%%%%%% DEBUG CHANGE TIME STEP %%%%%%%%%%%%%%%


% Determine if an element is explicit or implicit
mesh = create_mesh_imex(mesh,dt);

% Find implicit/explicit boundary
mesh = find_imex_bound(mesh);

% Find implicit faces
mesh = create_mesh_imface(mesh);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind = find(mesh.imex == 1); % Implicit element indicies
impE = length(ind)
% HDG solver
time=0;
app.time = time;
%figure(1);
% set(gcf,'color','black');
app.torder = torder;
dt = 0.0001;
app.dt = dt;

UDGexp = UDG;

  
for itime = 1:ntime
    fprintf('Timestep :  %d\n', itime);
    
    
    % Implicit Code
    UDG = hdg_vibration_imex_coupled_UPDATED_v5(master,mesh,app,UDG(:,1:3,:),UDG(:,4,:),time,dt,nstage,torder,Minv); 
    
%     [UDGimp,UH,PDG] = hdg_wave(master,mesh,app,UDG(:,1:3,:),PDG,time,dt,nstage,torder); 
    

    
    % Purely explicit, new bx + g method
%     [UDG2,B,G] = ehdgrk_imex_Ru_newform(master,mesh,app,Minv,UDG);
    
    % Purely explicit, old method
    UDGexp = ehdgrk_test(master,mesh,app,Minv,UDGexp);
%     UDGexp = ehdgrk(master,mesh,app,Minv,UDGexp);
    

    time = time + dt;
    app.time = time; 
          
%     pause(0.01)
end

% ustar = postprocess(master,mesh,master1,mesh1,UDG(:,4,:),UDG(:,2:3,:));
erru = calerror(UDG,mesh,master,@exactsol,[],itime*dt,[1 2 3 4]);
erruexp = calerror(UDGexp,mesh,master,@exactsol,[],itime*dt,[1 2 3 4]);
% erruex2 = calerror(UDGex2,mesh,master,@exactsol,[],itime*dt,[1 2 3 4]);

% erruimp = calerror(cat(2,UDGimp,PDG),mesh,master,@exactsol,[],itime*dt_loop,[1 2 3 4]);
% erruimp(:,i,j) = calerror(UDG,mesh,master,@exactsol,[],itime*dt_loop,[1 2 3 4]);
% erru2 = calerror(UDG2,mesh,master,@exactsol,[],itime*dt,[1 2 3 4]);

uex = exactsol(mesh.dgnodes, [], itime*dt);
figure(1); clf; scaplot(mesh,uex(:,1,:),[],2,1); axis off;
figure(2); clf; scaplot(mesh,UDG(:,1,:),[],2,1); axis off;
% figure(3); clf; scaplot(mesh,uex(:,4,:)-UDGimex(:,4,:),[],2,1); axis off;  
% figure(4); clf; scaplot(mesh,uex(:,4,:)-UDGex(:,4,:),[],2,1); axis off;
% figure(3); clf; scaplot(mesh,uex(:,4,:)-UDG2(:,4,:),[],2,1); axis off;
% figure(3); clf; scaplot(mesh,PDG,[],2,1); axis off; 

% count = count + 1;
%     end
% end