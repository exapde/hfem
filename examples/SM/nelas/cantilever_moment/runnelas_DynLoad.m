%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       CANTILEVER LOADED WITH A BENDING MOMENT          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setapplicationpath('SM/nelasuq')
setapplicationpath('SM/NeoHookean_uq')
%setapplicationpath('SM/SVK_uq')

% Elements type
hybrid = 'hdg';
elemtype = 1;
porder = 3;

global lamMin;
global lamMax;

% Mesh parameters
thickness = 0.1;
L = 12.;
l = 1. ;
nL = 9;
nl = 2;
nt = 2;
nodetype = 1;
Sur = thickness*l;

% Mechanical parameters
Eyoung = 1.2e6; nu = 0.;
lambda = Eyoung*nu/((1+nu)*(1-2*nu));
mu0 = Eyoung/(2*(1+nu));
lambda = lambda/mu0;
mu = 1.;
tau = 5.;

% Applied Total Load (traction)
applied_M = 50.*pi/3.;
applied_trac = 12.*applied_M/(thickness^3);
applied_trac = 1.*applied_trac/mu0;
Pmax = [0, applied_trac, 0];

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
app.bcm  = [1,2,3,3,3,3];  
app.bcs  = [ui; Pmax; ui; ui; ui; ui; ui];
app.bcd  = [1];
app.bcv  = [ui];

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
mesh = mkmesh_cube(nL,nl,nt,porder,L,l,thickness,elemtype,nodetype);

master = mkmaster(mesh,2*(porder));
[master,mesh] = preprocess(master,mesh,hybrid);

% Initialization of results sets
res.WB1 = []; res.WB2 = [];
res.WA1 = []; res.WA2 = [];
res.Nit = []; res.Loa = [];
res.Ener= [];
res.lmin{1} = []; res.lmax{1}= [];

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

P = [0. 0. 0.]; DP = 0.1*Pmax;
Nit1 = 20; it = 1; Nref = 1;
NitMax = 32;

while (P(2)<Pmax(2))
    
    % Increasing Load
    P = P + DP;
    if P(2)>Pmax(2)
        P = Pmax;
    end
    
    fprintf('########################################### \n'); 
    fprintf('Iteration : %d  Load fraction : %f  Load incr : %e\n', [it, P(2)/Pmax(2), DP(2)/Pmax(2)]); 
    fprintf('########################################### \n'); 
    
    lamMax = 0; lamMin = 0;
    
    % Newton HDG solution
    app.bcs  = [ui; P; ui; ui; ui; ui; ui];
    [UDG,UH,Nit] = hdg_solve(master,mesh,app,UDG0,UH0,[],NitMax);
    
    res.lmin{end} = [res.lmin{end};lamMin];
    res.lmax{end} = [res.lmax{end};lamMax];
    
    if (Nit>NitMax-1)
        P  = P - DP;
        DP = 0.5 * DP;
        Nref = Nref+1;
        if Nref>8, fprintf('#### NO CONVERGENCE #### '); break; end
    else
        if (Nit<5 && Nit1<5)
            DP = 1.5 * DP;
        end
        UDG0 = UDG; UH0 = UH;
        Nit1 = Nit; it = it+1; Nref = 1;
        res.Nit = [res.Nit ; Nit];
        res.Loa = [res.Loa ; P(2)/Pmax(2)];
        res.lmin{end+1} = [];
        res.lmax{end+1} = [];
        
        % Compute Energy
        [Epsi,Ehat,Dmax] = comp_energy(UDG,UH,mesh,master,app,'NeoHookean');
        res.Ener = [res.Ener; [mu0*Epsi,mu0*Ehat,Dmax]];

        % Save End points position
        switch porder
            case 1
                res.WA1 = [res.WA1; UDG(2,1:3,nL-1)];
                res.WA2 = [res.WA2; UDG(6,1:3,nL-1)];
                res.WB1 = [res.WB1; UDG(4,1:3,nL-1)];
                res.WB2 = [res.WB2; UDG(8,1:3,nL-1)];
            case 2
                res.WA1 = [res.WA1; UDG(3,1:3,nL-1)];
                res.WA2 = [res.WA2; UDG(21,1:3,nL-1)];
                res.WB1 = [res.WB1; UDG(9,1:3,nL-1)];
                res.WB2 = [res.WB2; UDG(27,1:3,nL-1)];
            case 3
                res.WA1 = [res.WA1; UDG(4,1:3,nL-1)];
                res.WA2 = [res.WA2; UDG(52,1:3,nL-1)];
                res.WB1 = [res.WB1; UDG(16,1:3,nL-1)];
                res.WB2 = [res.WB2; UDG(64,1:3,nL-1)];
        end
        
    end
    
end

% original configuration
figure(1); clf; meshplot(mesh,1);  axis equal; axis tight; axis on;

% deformed configuration
mesh1=mesh;
mesh1.dgnodes(:,1:3,:)=UDG(:,1:3,:);
figure(2); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;
% figure(3); scaplot(mesh,UDG(:,1,:),[],2,1); axis equal; axis tight;
% figure(4); scaplot(mesh,UDG(:,2,:),[],2,1); axis equal; axis tight;
% figure(5); clf; plot(res.Loa,res.Ener(:,1),'-s');

lminVec = zeros(length(res.Loa),1);
lmaxVec = zeros(length(res.Loa),1);
for i=1:length(res.Loa)
    lminVec(i) = res.lmin{i}(end);
    lmaxVec(i) = res.lmax{i}(end);
end
plot(res.Loa,lminVec,'-s');
plot(res.Loa,lmaxVec,'-s');
