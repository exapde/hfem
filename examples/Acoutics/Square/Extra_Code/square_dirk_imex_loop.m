setapplicationpath('FM/wave');

porder = 1;
nstage = 3;
torder = 3;
hybrid = 'hdg';
nrefall = [2,4,6,8];
dtall = [0.006, 0.003, 0.0015, 0.00075];
erru = zeros(4,5);
r = 4;


c2 = 1;
k = [1, pi/6];
% dt = 0.012; % dt to create mesh when no nref
% dt = 0.006; %dt when nref = 2
% dt = 0.003; % dt when nref = 4
% dt = 0.0015;
% dt = 0.00075
% dt = 0.000325;
% dt = 0.015;
ntime = 2000;
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

for i = 1:r
    dt = dtall(i);
    nref = nrefall(i);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mesh=mkmesh_square_imex(porder);
    mesh = refinemesh(mesh,porder,nref);
    master = mkmaster(mesh,2*porder);

    mesh1=mkmesh_square_imex(porder+1);
    master1 = mkmaster(mesh1,2*(porder+1));

    [master,mesh] = preprocess(master,mesh,hybrid);

    [master1,mesh1] = preprocess(master1,mesh1,hybrid);


    UDG = l2eprojection(mesh,master,@exactsol,[c alpha phi],0,4);
    % UDG = l2eprojection(mesh,master,@exactsol,[],0,4);
    PDG(:,1,:)=UDG(:,4,:);
    UDG(:,4,:)=[];

    UH = inituhat(master,mesh.elcon,UDG,app.ncu);



    % Explicit Code Precomputation
    [M, Minv]=massinv(master, mesh);
    UDG(:,4,:)=PDG; % explicit HDG
    mesh.t2t = mkt2t(mesh.t,mesh.elemtype);

    % Determine if an element is explicit or implicit
    mesh = create_mesh_imex(mesh,dt);

    % Find implicit/explicit boundary
    mesh = find_imex_bound(mesh);

    % Find implicit faces
    mesh = create_mesh_imface(mesh);
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % HDG solver
    time=0;
    app.time = time;
    %figure(1);
    set(gcf,'color','black');
    app.torder = torder;
    dt = 0.0005;



    for itime = 1:ntime
        fprintf('Timestep :  %d\n', itime);


        % Implicit Code
%         UDG = hdg_vibration_imex_coupled(master,mesh,app,UDG(:,1:3,:),UDG(:,4,:),time,dt,nstage,torder,Minv); 
        [UDG,UH,PDG] = hdg_wave(master,mesh,app,UDG(:,1:3,:),PDG,time,dt,nstage,torder); 

    %     app.dt = dt;   
    %     app.time = time; 
    %     UDG = ehdgrk_imex_coupled(master,mesh,app,Minv,UDG);

        %UDG = cat(2,UDGimp,PDG); % Fixes notation differences between explicit and implicit codes

%         figure(1); clf; scaplot(mesh,UDG(:,4,:),[-1 1],2,1); axis off; colormap jet; 

        figure(1); clf; scaplot(mesh,PDG,[-1 1],2,1); axis off; colormap jet;


        time = time + dt;


    %     pause(0.01)
    end
    erru(:,i+1) = calerror(UDG,mesh,master,@exactsol,[],itime*dt,[1 2 3 4]);
end

    uex = exactsol(mesh.dgnodes, [], itime*dt);
    figure(2); clf; scaplot(mesh,uex(:,4,:)-UDG(:,4,:),[],2,1); axis off;    