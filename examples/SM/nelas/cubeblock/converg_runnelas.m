setapplicationpath('SM/nelasuq')
%setapplicationpath('SM/leuq')

hybrid = 'hdg';
elemtype = 1;
nodetype = 0;

mu = 1;
lambda = 1;
tau = 2;
ui = [0, 0, 0];

app.source = 'source';
app.flux = 'flux';
app.fbou = 'fbou';
app.fhat = 'fhat';

app.iterative = 0;
app.hybrid = hybrid;
app.localsolve = 1;
app.arg = {mu,lambda,tau};
app.bcm  = [2,2,1,2,3,3];  
%app.bcm  = [1,3,3,3,3,3];  
app.bcs  = [ui; ui; ui; ui; ui; ui];
app.bcd  = [1];  
app.bcv  = [ui];

app.wave = false;
app.tdep = false;
app.alag = false;
app.uqpk = false;

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

% Parameters for convergence study :
n_ref   = 7;        % Number of refinements
ngrid0  = 2;        % Initial nb of element in each directions
refcoef = sqrt(2.); % Refinement coefficient

% Initialization of errors results
L2errors.errors = zeros(n_ref,app.nc);
L2errors.orders = zeros(n_ref,app.nc);
L2errors.Uerror = zeros(n_ref,2);
L2errors.Qerror = zeros(n_ref,2);
L2errors.h      = zeros(1,n_ref);

for i_ref = 1:n_ref
    
    ngrid = round(ngrid0*refcoef^(i_ref-1))+1;
    fprintf('Number of elements in each direction :  %d\n', ngrid-1);
    porder= 2;
    
    mesh = mkmesh_cube(ngrid,ngrid,ngrid,porder,1,1,1,elemtype,nodetype);
    master = mkmaster(mesh,2*(porder));
    [master,mesh] = preprocess(master,mesh,hybrid);
    
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
    [UDG,UH] = hdg_solve(master,mesh,app,UDG0,UH0,[]);
    
    % Computation of the errors
    err = calerror(UDG,mesh,master,@exactsolution);
    L2errors.h(i_ref)     = 1./(ngrid-1);   % Mesh size
    L2errors.errors(i_ref,:) = err;         % Errors
    % Convergence Orders
    if i_ref== 1 
        L2errors.orders(i_ref,:) = 0.;
    else
        L2errors.orders(i_ref,:) = log(L2errors.errors(i_ref,:)./L2errors.errors(i_ref-1,:))...
                                 /(log(L2errors.h(i_ref)/L2errors.h(i_ref-1)));
    end
end

% Computing error and convergence orders for vector U and tensor Q
L2errors.Uerror(:,1) = sqrt(L2errors.errors(:,1).^2 + L2errors.errors(:,2).^2 + L2errors.errors(:,3).^2);
L2errors.Uerror(2:end,2) = log(L2errors.Uerror(2:end,1)./L2errors.Uerror(1:end-1,1))...
                        ./(log(L2errors.h(2:end)./L2errors.h(1:end-1))');
L2errors.Qerror(:,1) = sqrt(L2errors.errors(:,4).^2 + L2errors.errors(:,5).^2 + L2errors.errors(:,6).^2 ...
                          + L2errors.errors(:,7).^2 + L2errors.errors(:,8).^2 + L2errors.errors(:,9).^2 ...
                          + L2errors.errors(:,10).^2+ L2errors.errors(:,11).^2+ L2errors.errors(:,12).^2);
L2errors.Qerror(2:end,2) = log(L2errors.Qerror(2:end,1)./L2errors.Qerror(1:end-1,1))...
                        ./(log(L2errors.h(2:end)./L2errors.h(1:end-1))');

% original configuration
figure(1); clf; meshplot(mesh,1); axis equal; axis tight; axis on;

% deformed configuration
mesh1=mesh;
mesh1.dgnodes = UDG(:,1:3,:); 
figure(2); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;

% exact deformed configuration
uex = exactsolution(mesh.dgnodes);
mesh1.dgnodes(:,1:3,:)=uex(:,1:3,:);
figure(3); clf; meshplot(mesh1,1); axis equal; axis tight; axis on;

% Convergence plot
figure(4); clf; loglog(1./L2errors.h,L2errors.Uerror(:,1),'-s'); grid on;
title('Convergence for error on displacements u'); xlabel('1/h'); ylabel('||u-u_h||')


