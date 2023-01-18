%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        CYLINDER SHELL PINCHED BY A NODAL FORCE         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setapplicationpath('SM/nelasuq')
setapplicationpath('SM/NeoHookean_uq')
%setapplicationpath('SM/SVK_uq')

% Elements type
hybrid = 'hdg';
elemtype = 1;
nodetype = 1;
porder = 2;

% Mesh parameters
thickness = 12.7;
alpha = 0.1;
R = 2540.;
L = 254.;
nLon = 9;
nL = 9;
nt = 2;

% Mechanical parameters
Eyoung = 3102.75; nu = 0.3;
lambda = Eyoung*nu/((1+nu)*(1-2*nu));
mu0 = Eyoung/(2*(1+nu));
lambda = lambda/mu0;
mu = 1.;
tau = 2.;

% Application parameters
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.iterative = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {mu,lambda,tau};

% Boundary conditions
bodyforce = [0 0 0];
ui = [0, 0, 0];
app.bcm  = [2,5,5,5,4,5];
app.bcs  = [ui; ui; ui; ui; ui; ui; ui];
app.bcd  = [1];
app.bcv  = [ui];

% Informations relative to the point source (if any)
app.ptsource.globpos = [0 R+0.5*thickness 0];
app.ptsource.dir = [0 -1 0];
app.ptsource.locpos = [0 0.5 0];
app.ptsource.ibpl = [4];
app.ptsource.amplmax = [750./mu0];
app.ptsource.surf = 1;

% Application flags
app.wave = false;
app.tdep = false;
app.alag = false;
app.uqpk = true;

app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;
app.fc_q = 1;
app.fc_u = 1;
app.fc_p = 0;

app.nd   = 3;
app.ncu  = 3;               % Number of components of U
app.nch  = app.ncu;                % Number of componets of UH
app.ncq  = app.nch*app.nd;         % Number of componets of Q
app.nc   = app.ncu+app.ncq;         % Number of componeents of UDG
app.ncp = 0;

app.time = [];
app.dtfc = [];
app.alpha = [];

% Mesh generation
mesh = mkmesh_cylindralshell2(nLon,nt,nL,porder,L,R,thickness,alpha,elemtype,nodetype);
master = mkmaster(mesh,2*(porder));
[master,mesh] = preprocess(master,mesh,hybrid);

% Some Diriclet Nodes are remooved from he resolution
app.nodirnodes = true;
mesh.nfixed = locate_Dirichlet_nodes(mesh, [0. 0. 1.], [R*sin(alpha) R*cos(alpha) 0.]);


% Finding closest faces and nodes for the points loads
if isfield(app,'ptsource'), app = locate_ptsour(mesh, app); end

% Initialization of results sets
res.WB1 = []; res.WB2 = [];
res.WA1 = []; res.WA2 = [];
res.Nit = []; res.Loa = [];
res.Ener = [];

UDG0 = initu(mesh,{0;0;0;0;0;0;0;0;0;0;0;0});
x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
z = mesh.dgnodes(:,3,:);
UDG0(:,1,:) = x;
UDG0(:,2,:) = y;
UDG0(:,3,:) = z;
UDG0(:,4,:) = -1;
UDG0(:,5,:) = 0;
UDG0(:,6,:) = 0;
UDG0(:,7,:) = 0;
UDG0(:,8,:) = -1;
UDG0(:,9,:) = 0;
UDG0(:,10,:) = 0;
UDG0(:,11,:) = 0;
UDG0(:,12,:) = -1;
UH0 = inituhat(master,mesh.elcon,UDG0,app.ncu);

fprintf('\nSteady State Solve\n');

DDUH0 = sparse(size(UH0,1),size(UH0,2));
lambda = 0.; DDlambda = 0; DDlambda0 = 1.;
Nit1 = 20; it = 1; Nref = 1; NitMax = 20;
% Arclength parameters :
Dl = 20.;
psi = Dl/(0.1*app.ptsource.amplmax) ;

while (lambda<1)
    
    if lambda>1
        lambda=1;
    end
    
    fprintf('########################################### \n'); 
    fprintf('Iteration : %d  Load fraction : %f  Load incr : %e\n', ...
        [it, lambda, DDlambda]);
    fprintf('########################################### \n'); 
    
    % Newton HDG solution
    [UDG,UH,DDUH,DDlambda,Nit] = hdg_solve_arclength(master,mesh,app,UDG0,UH0,[],DDUH0,lambda,DDlambda0,psi,Dl,NitMax);
    
    
    if (Nit>NitMax-2)
        Dl = 0.5 * Dl;
        Nref = Nref+1;
        if Nref>12, fprintf('#### NO CONVERGENCE #### '); break; end
    else
        if (Nit<5 && Nit1<5 && Dl<2.)
            Dl = 1.5 * Dl;
        end
        % Updates are accepted
        UDG0 = UDG; UH0 = UH; DDUH0 = DDUH;
        lambda = lambda+DDlambda; DDlambda0 = DDlambda;
        Nit1 = Nit; it = it+1; Nref = 1;
        res.Nit = [res.Nit ; Nit];
        res.Loa = [res.Loa ; lambda];
        
        % Compute Energy
        [Epsi,Ehat,Dmax] = comp_energy(UDG,UH,mesh,master,app,'NeoHookean');
        res.Ener = [res.Ener; [mu0*Epsi,mu0*Ehat,Dmax]];
        
        % Save End points position
        res.WA1 = [res.WA1; UDG(porder+1,1:3,nLon-1)];
        res.WA2 = [res.WA2; UDG((porder+1)^2,1:3,nLon-1)];
        res.WB1 = [res.WB1; UDG(1,1:3,1)];
        res.WB2 = [res.WB2; UDG((porder+1)*porder+1,1:3,1)];
        
    end
    
    % Current Load Plot
    run results/plot_loads.m
    drawnow;
    
end

% original configuration
figure(2); clf; meshplot(mesh,1);  axis equal; axis tight; axis on;

% deformed configuration
mesh1=mesh;
mesh1.dgnodes(:,1:3,:)=UDG(:,1:3,:);
figure(3); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;

% Write outputs for Paraview visu
UDG0(:,1:3,:)   = UDG(:,1:3,:)-mesh.dgnodes(:,1:3,:);
UDG0(:,4:end,:) = UDG(:,4:end,:);
writeINPfile('outputs', mesh1, UDG0);

% symVisu(mesh1,135,25,1);
