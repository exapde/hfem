setapplicationpath('FM/wave');

%porder = 2;
%ngrid  = 17;
nstage = 3;
torder = porder+2;
hybrid = 'hdg';

c2 = 1;
k = [1,pi/6];
dt = (fac/(ngrid-1));
ntime = T/dt;

app.tau = 1;
app.kappa = c2;
app.localsolve=1;
app.arg = {c2,k};
app.bcm = [1;1;1;1]*5;
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

mesh = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

c = sqrt(c2);
alpha = k(1);
phi = k(2);

clear Vn;
UDG = l2eprojection(mesh,master,@exactsol,[c alpha phi],0,4);
Vn(:,1,:)=UDG(:,1,:);
UDG(:,1,:) = UDG(:,4,:);
UDG(:,4,:)=[];

%UH = inituhat(master,mesh.elcon,UDG,app.ncu);
%VH = inituhat(master,mesh.elcon,Vn,app.ncu);
UH = l2fprojection(mesh,master,@exactsol,[c alpha phi],0,4);
VH = UH(1,:);
UH = UH(4,:);
UH = UH'; VH = VH';

x=mesh.dgnodes(:,1,:);
y=mesh.dgnodes(:,2,:);
% HDG solver
for itime = 1:ntime
    t = dt*itime;
    fprintf('Timestep :  %d\n', itime);
            
    [UDG,UH,Vn,VH] = hdg_wave_ec(master,mesh,app,UDG,Vn,UH,VH,(itime-1)*dt,dt,nstage,torder);
%    [UDG,UH,Vn] = hdg_wave_ec2(master,mesh,app,UDG,Vn,(itime-1)*dt,dt,nstage,torder);
end
Un = UDG(:,1,:);
[QDG,UHAT] = get_quhat(master, mesh, app, Un, t, 4);
[PDG,VHAT] = get_quhat(master, mesh, app, Vn, t, 1);
UHAT = full(UHAT);
Eh = get_Eh(master, mesh, app, Un, Vn, QDG, reshape(UHAT(mesh.elcon),[(porder+1) 3 mesh.ne]));

UDG(:,4,:) = Un;
UDG(:,1,:) = Vn;
e = calerror(UDG,mesh,master,@exactsol,[c k],t,[1 2 3 4]);

mesh1 = mkmesh_square(ngrid,ngrid,porder+1);
master1 = mkmaster(mesh1,2*(porder+1));
[master1,mesh1] = preprocess(master1,mesh1,hybrid);

clear UDGstar VDG;
VDG(:,1,:) = Un; VDG(:,2:3,:) = -QDG; 
UDGstar(:,1,:) = postprocessnd(master,mesh,master1,mesh1,VDG(:,1:3,:)); % u
VDG(:,1,:) = Vn; VDG(:,2:3,:) = -PDG;
UDGstar(:,2,:) = postprocessnd(master,mesh,master1,mesh1,VDG(:,1:3,:)); % v
es = calerror(UDGstar,mesh1,master1,@exactsol,[c k],t,[4 1]);   


% UDG(:,4,:) = UDG(:,1,:);
% UDG(:,1,:) = Vn;
% e = calerror(UDG,mesh,master,@exactsol,[c k],t);
% 
% mesh1 = mkmesh_square(ngrid,ngrid,porder+1);
% master1 = mkmaster(mesh1,2*(porder+1));
% [master1,mesh1] = preprocess(master1,mesh1,hybrid);
% 
% VDG = -UDG; VDG(:,1,:) = UDG(:,4,:);
% clear UDGstar;
% UDGstar(:,1,:) = postprocessnd(master,mesh,master1,mesh1,VDG(:,1:3,:));
% VDG(:,1,:) = Vn(:,1,:); VDG(:,2:3,:) = getq(mesh,master,Vn,reshape(VH(mesh.elcon),[1 (porder+1)*3 mesh.ne])); 
% UDGstar(:,2,:) = postprocessnd(master,mesh,master1,mesh1,VDG(:,1:3,:));
% es = calerror(UDGstar,mesh1,master1,@exactsol1,[c k],t);   

return;


nn = [2 4 8 16 32]+1;
for porder = 1:2
    for ii=1:length(nn)
        ngrid = nn(ii);
        wave_ec;             
        errv(ii,porder) = e(1); errq(ii,porder) = sqrt(e(2)^2+e(3)^2);
        erru(ii,porder) = e(4);
        ersu(ii,porder) = es(1); 
        ersv(ii,porder) = es(2);
        Erh(ii,porder) = Eh(1);
    end
end
cu=log(erru(1:end-1,:)./erru(2:end,:))/log(2);
cv=log(errv(1:end-1,:)./errv(2:end,:))/log(2);
cq=log(errq(1:end-1,:)./errq(2:end,:))/log(2);
su=log(ersu(1:end-1,:)./ersu(2:end,:))/log(2);
sv=log(ersv(1:end-1,:)./ersv(2:end,:))/log(2);

% a=[erru [0 0 0 0; cu]; errq [0 0 0 0; cq]; errs [0 0 0 0; cs]];
% a=a(:,[1 5 2 6 3 7 4 8]);
% a=[[nn-1 nn-1 nn-1]' a];
% 

