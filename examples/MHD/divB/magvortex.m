clear

setapplicationpath('MHD/mhd')

% ------ Parameters --------- %
hybrid      = 'hdg';
elemtype    = 0;
nodetype    = 0;
porder      = 1;        % HDG order
nstage      = 3;        % n-stage DIRK
torder      = 3;        % t-order DIRK
T           = 10;        % Final time
ngrid       = 9;       % mesh discretization
dtps        = 1/(ngrid-1);     % Time step
ntime       = T/dtps;   % number time iteration
ntime       = ceil(ntime);
gam         = 5/3;      % adiabatic constant depending on the physical properties of the fluid
icase       = 5;        % 1:scalar case, 2:Orszag-Tang vortex, 3:Rotor problem, 4:Smooth Alfven wave, 5:Smooth vortex problem


app.source  = 'source2d';
app.flux    = 'flux2d';
app.fbou    = 'fbou_mhd';
app.fhat    = 'fhat2d';

app.iterative = 0;
app.hybrid  = hybrid;
app.localsolve = 0;
app.bcm     = [];
app.bcs     = [];
app.bcd     = [];
app.bcv     = [];
app.wave    = false;
app.tdep    = true;
app.alag    = false;
app.uqpk    = 0;
app.flg_q   = 0;
app.flg_p   = 0;
app.flg_g   = 0;
app.fc_q    = 0;
app.fc_u    = 1;
app.fc_p    = 0;

app.nd      = 2;            % Dimension
app.ncu     = 5 + app.nd;   % Number of components of U
app.nch     = app.ncu;      % Number of componets of UH
app.nq      = 0;            % Number of componets of Q
app.nc      = app.ncu;      % Number of componeents of UDG
app.nco     = 0;
app.ncq     = 0; 
app.ncp     = 0; 

app.time    = 0;
app.dt      = dtps*ones(1,ntime);
app.torder  = torder;
app.nstage  = nstage;
app.dtcoef = [1.0;1.0;1.0;1.0;1.0;1.0;1.0];
%app.dtcoef = [1.0;1.0;1.0;1.0;1.0;1.0;0.0];
           
app.adjoint         = 0;               % 1 if adjoint problem. 0 otherwise
app.linearproblem   = 0;               % 0 if problem is linear. 1 if nonlinear
app.appname         = 'mhd';
app.linearSolver    = 1;
app.jacobianStep    = 0;
app.orderingStep    = 0;
app.morder          = [porder porder];
app.porder          = [porder porder];
app.nodetype        = nodetype;
app.pgauss          = 2*[porder porder];
app.pgaussR         = 2*[porder porder];
app.overlappinglevel= 1;
app.preconditioner  = 0;
app.quadtype        = [0 0];
check               = 0;

% ----- Mesh depending of the test case ------- %
if icase == 1 || icase == 2
   mesh = mkmesh_square(ngrid,ngrid,porder,0,2*pi,2*pi);
elseif icase == 3
   mesh = mkmesh_square(ngrid,ngrid,porder,0,1,1);
elseif icase == 4 
   mesh = mkmesh_rect(ngrid,ngrid,porder,0,[0,1/cos(pi/6),0,1/sin(pi/6)]);
elseif icase == 5
   mesh = mkmesh_rect(ngrid,ngrid,porder,0,[-5,5,-5,5]);
elseif icase == 6
   mesh = mkmesh_rect(ngrid,ngrid,porder,0,[-0.5,0.5,-0.5,0.5]);
end
master  = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

% Periodic boundary conditions
mesh = periodic(mesh,{1,'p(:,1)',3,'p(:,1)';4,'p(:,2)',2,'p(:,2)'},2);

% ------ Initialize solution ------------- %
%UDG0    = initsol(mesh.dgnodes,icase,gam);
UDG = initL2(master,mesh,icase,gam);
%UH0  = inituhat(master,mesh.elcon,UDG0,app.ncu);
UH  = initL2f(master,mesh,icase,gam);
%[Q,~,~] = getq(master,mesh,UDG0,UH0);
%  ----- Compute Tau and alpha1 ------ %
app.arg{1}   = gam;
[alpha1,tau] = glm2D(mesh,app,UDG,UH,ngrid);
alpha1 = 1.0;
alpha = 2.0;
%alpha = 0.0;
app.arg      = {gam,tau,alpha1,alpha};

app.ncd     = size(mesh.dgnodes,2);
elementtype = elemtype*ones(mesh.ne,1);

% ------- Newton Parameters ---------- %
app.newtoniter = 100;  % def 10
app.newtontol  = 1e-13; % def 1e-7
% % time step looping 
time = 0;

%TableError = zeros(length(ntime), 10);
for itime = 1:length(app.dt)
    dt = app.dt(itime);
    fprintf('Timestep :  %d,   Time :  %g\n', itime, time + dt);
    
%     [alpha1,alpha2,tau] = glm(mesh,app,UDG,UH,ngrid);
%     app.arg{2} =  tau;
            
    [UDG,UH] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],time,dt,nstage,torder);    
    time = time + dt;  

%     % --- Evaluate L2 errors in time ---
%     err = calerror(UDG,mesh,master,@exactsolAlfven);
%     TableError(itime,:) = [time err'];
    
%          figure(1)
%          clf; scaplot(mesh,UDG(:,1,:),[],2,0); 
%          axis on; axis square; axis tight; colorbar; colormap jet;
%          title('density');
%          drawnow
%          
         figure(2)
         clf; scaplot(mesh,UDG(:,2,:),[],2,0); 
         axis on; axis square; axis tight; colorbar; colormap jet;
         title('rho*ux');
         drawnow       
%         
%          figure(3)
%          clf; scaplot(mesh,UDG(:,3,:),[],2,0); 
%          axis on; axis square; axis tight; colorbar; colormap jet;
%          title('rho*uy');
%          drawnow
%          
%          figure(4)
%          clf; scaplot(mesh,UDG(:,4,:),[],2,0); 
%          axis on; axis square; axis tight; colorbar; colormap jet;
%          title('rho*uz');
%          drawnow
%          
%          figure(5)
%          clf; scaplot(mesh,UDG(:,5,:),[],2,0); 
%          axis on; axis square; axis tight; colorbar; colormap jet;
%          title('e');
%          drawnow 
%          
%          figure(6)
%          clf; scaplot(mesh,UDG(:,6,:),[],2,0); 
%          axis on; axis square; axis tight; colorbar; colormap jet;
%          title('Bx');
%          drawnow
%          
%          figure(7)
%          clf; scaplot(mesh,UDG(:,7,:),[],2,0); 
%          axis on; axis square; axis tight; colorbar; colormap jet;
%          title('By');
%          drawnow  
%         
%          figure(8)
%          clf; scaplot(mesh,UDG(:,8,:),[],2,0); 
%          axis on; axis square; axis tight; colorbar; colormap jet;
%          title('Bz');
%          drawnow
%          
% %          figure(9)
% %          clf; scaplot(mesh,UDG(:,9,:),[],2,0); 
% %          axis on; axis square; axis tight; colorbar; colormap jet;
% %          title('psi');
% %          drawnow           
end
err = calerror(UDG,mesh,master,@exactsolvortex);