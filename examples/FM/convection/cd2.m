setapplicationpath('FM/condiff');
% 
% nn = [5 10 20 40 80]+1;
% for porder = 0:3
%     for ii=1:length(nn)
%         ngrid = nn(ii);
%         cd2;
%         err(ii,porder+1) = calerror2(UDG(:,1,:),mesh,master);        
%         [errh(ii,porder+1),errb(ii,porder+1)] = faceerror2(UDG(:,1,:),UH,mesh,master);
%     end
% end

porder = 3;
ngrid  = 11;
hybrid = 'hdg';

kappa = 0;
c = [1,1]; 
tau = 1/2;

app.source = 'source1';
app.flux = 'flux';
app.fbou = 'fbou1';
app.fhat = 'fhat';
app.adjoint = 0;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve=1;
app.arg = {kappa,c,tau};
app.bcm = [7;2;2;1];
app.bcs = [0;0;0;1]; 
app.bcd = [1;2;2;1];
app.bcv = [0;0;0;0]; 

app.upqk=0;
app.denseblock = 0;
app.tdep = false;
app.wave = false;
app.alag = false;
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
master = mkmaster(mesh,max(2,2*porder));
[master,mesh] = preprocess(master,mesh,app);
%return;

%UDG = initu(mesh,{0;0;0});
UDG = zeros(master.npv,3,mesh.ne);
UH = inituhat(master,mesh.elcon,UDG,app.ncu);

% HDG solver
[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);
figure(1); clf; scaplot3(mesh,UDG(:,1,:),[-0.2 1.2],2); 
%scaplot(mesh,UDG(:,1,:),[0 1],1,1,1)
axis equal; axis tight; colormap jet;
axis([0 1 0 1 -0.2 1.2]);
axis normal; colorbar off;
set(gca,'FontSize',16);
set(gca,'xtick',[0:0.2:1]);
set(gca,'ytick',[0:0.2:1]);
set(gca,'ztick',[-0.2:0.2:1.2]);
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
box on;
return;

[y,u,v,w] = uhplot(mesh,master,UH,UDG);
if porder == 0
    y0 = y; u0 = u; v0 = v; w0 = w;
elseif porder==1
    y1 = y; u1 = u; v1 = v; w1 = w;
elseif porder==2
    y2 = y; u2 = u; v2 = v; w2 = w;
elseif porder==3
    y3 = y; u3 = u; v3 = v; w3 = w;    
end
y0 = reshape(y0,40,[]); w0 = reshape(w0,40,[]);
y1 = reshape(y1,40,[]); w1 = reshape(w1,40,[]);
y2 = reshape(y2,40,[]); w2 = reshape(w2,40,[]);
y3 = reshape(y3,40,[]); w3 = reshape(w3,40,[]);
figure(2); clf;
hold on;
for i = 1:size(y0,2)
    plot(y0(:,i),w0(:,i),'-k','LineWidth',1);
end
for i = 1:size(y1,2)
    plot(y1(:,i),w1(:,i),'-b','LineWidth',1);
end
for i = 1:size(y2,2)
    plot(y2(:,i),w2(:,i),'-r','LineWidth',1);
end
for i = 1:size(y3,2)
    plot(y3(:,i),w3(:,i),'-m','LineWidth',1);
end


figure(2); clf;
plot(y0(:),w0(:),'-g',y1(:),w1(:),'-.b',y2(:),w2(:),'--r',y3(:),w3(:),':k','LineWidth',2);
axis equal;
axis([0 sqrt(2) -0.2 1.2]);
set(gca,'xtick',[0:0.2:1.4]);
set(gca,'ytick',[-0.2:0.2:1.2]);
set(gca,'FontSize',16);
legend({'$k=0$','$k=1$','$k=2$','$k=3$'},'interpreter','latex','FontSize',15);

figure(3); clf;
plot(y0(:),w0(:),'-g',y1(:),w1(:),'-.b',y2(:),w2(:),'--r','LineWidth',2);
axis equal;
axis([0 sqrt(2) -0.1 1.1]);
set(gca,'xtick',[0:0.1:1.4]);
set(gca,'ytick',[-0.2:0.2:1.2]);
set(gca,'FontSize',16);
legend({'$k=0$','$k=1$','$k=2$'},'interpreter','latex','FontSize',15);

%axis tight;

% figure(2); clf;
% plot(y,u,'-k',y,v,'-b',y,w,'-r');
% axis equal;
% axis tight;

% figure(3); clf;
% plot(y(1:end-1),w(1:end-1),'-r');
% axis equal;
% axis tight;

x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
uex = 0*x;
for i = 1:mesh.ne
    xm = mean(x(:,1,i));
    ym = mean(y(:,1,i));
    if xm<=ym+0.2
        uex(:,1,i) = 1;
    end
end
figure(2); clf; scaplot3(mesh,uex,[],0); 
axis equal; axis tight; colormap jet;
axis([0 1 0 1 -0.2 1.2]);
axis normal
colorbar
% 
%err = calerror(UDG(:,1,:),mesh,master,@exactsol1);

% figure(2); clf; 
% meshplot(mesh);
% set(gca,'FontSize',16);
% axis equal; axis tight; 


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




