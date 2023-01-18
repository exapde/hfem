setapplicationpath('FM/condiff');

% nn = [2 4 8 16 32]+1;
% for porder = 1:4
%     for ii=1:length(nn)
%         ngrid = nn(ii);
%         cd1;
%         err(ii,porder) = calerror(UDG(:,1,:),mesh,master,@exactsol1);
%     end
% end

porder = 1;
ngrid  = 2*4+1;
hybrid = 'hdg';

kappa = 0;
c = [1,1]; 
tau = 1/2+1/2;
fhatform = 3;
Aform = 4;

app.source = 'source1';
app.flux = 'flux';
app.fbou = 'fbou1';
app.fhat = 'fhat';
app.adjoint = 0;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve=1;
app.arg = {kappa,c,tau,fhatform,Aform};
app.bcm = [5;2;2;1];
app.bcs = [0;0;0;0]; 
app.bcd = [1;2;2;1];
app.bcv = [0;0;0;0]; 

app.denseblock = 0;
app.tdep = false;
app.wave = false;
app.alag = false;
app.uqpk=0;
app.flg_q = 1;
app.flg_p = 0;
app.flg_g = 0;

app.fc_q = 1;
app.fc_u = 0;
app.fc_p = 0;

app.np = 2;
app.nd = 2;
app.nch  = 1;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = 1;

app.appname = 'cd';
app.time = [];
app.dtfc = [];
app.alpha = [];

mesh   = mkmesh_square(ngrid,ngrid,porder,1);
% mesh.p(:,1) = logdec(mesh.p(:,1),1);
% mesh.p(:,2) = logdec(mesh.p(:,2),1);
% mesh.dgnodes(:,1,:) = logdec(mesh.dgnodes(:,1,:),1);
% mesh.dgnodes(:,2,:) = logdec(mesh.dgnodes(:,2,:),1);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,app);
%return;

UDG = zeros(master.npv,3,mesh.ne);
UH = inituhat(master,mesh.elcon,UDG,app.ncu);
% UDG = initu(mesh,{0;0;0});
% UH = inituhat(master,mesh.elcon,UDG,app.ncu);

% HDG solver
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

figure(1); clf; scaplot(mesh,UDG(:,1,:),[],2); 
axis equal; axis tight; colormap jet;
x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
uex = sin(pi*(x-y));
ind = find(x<=y);
uex(ind) = 0;
figure(2); clf; scaplot(mesh,uex,[],2); 
axis equal; axis tight; colormap jet;



% % HDG postprocessing 
% mesh1 = mkmesh_square(ngrid,ngrid,porder+1);
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG);
% 
% %figure(2); clf; scaplot(mesh1,UDGstar(:,1,:),[],2); axis equal; axis tight;
% 
% [VDG,VH] = hdg_solve_adjoint(master,mesh,UDG,UH,[],app);
% figure(3);clf; scaplot(mesh,VDG(:,1,:),[],2); axis equal; axis tight;
% 

% app.source = 'source';
% app.flux = 'flux';
% app.ubou = 'ldgubou';
% app.fbou = 'ldgfbou';
% app.fhat = 'ldgfhat';
% tol = 1e-8;
% dt = 1e-1*ones(1000,1);
% [UDGA,UHA,normR] = ardm(master,app,mesh,0*UDG(:,1,:),dt/1.5,tol);
% 
% figure(2); clf; scaplot(mesh,UDGA(:,1,:),[],0,1); axis equal; axis tight; colormap jet;




