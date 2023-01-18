setapplicationpath('FM/wave');

porder = 2;
ngrid  = 9;
nstage =porder+1;
torder =porder+2;
hybrid = 'hdg';

c2 = 1;
k = [1,pi/6];
dt = (0.25/(ngrid-1));
ntime = 1/dt;

app.tau = c2;
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

%mesh = mkmesh_circleinsquare(porder,1,1);
mesh = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,hybrid);

%UDG = initu(mesh,{0;0;0});
%UH = inituhat(master,mesh.elcon,UDG,app.ncu);
PDG = initu(mesh,{0});

c = sqrt(c2);
alpha = k(1);
phi = k(2);
x=mesh.dgnodes(:,1,:);
y=mesh.dgnodes(:,2,:);
% u=(1/(sqrt(2)*pi))*sin(pi*x).*sin(pi*y).*sin(sqrt(2)*pi*0);    
% v=sin(pi*x).*sin(pi*y).*cos(sqrt(2)*pi*0);
% qx=-(1/sqrt(2))*cos(pi*x).*sin(pi*y).*sin(sqrt(2)*pi*0);    
% qy=-(1/sqrt(2))*sin(pi*x).*cos(pi*y).*sin(sqrt(2)*pi*0);    
t = 0;
% c = sqrt(c2);
% u = -tanh(alpha*(x*cos(phi) - c*t + y*sin(phi)) - 1);
% v = -alpha*c*(tanh(alpha*(x*cos(phi) - c*t + y*sin(phi)) - 1).^2 - 1);
% qx = -(alpha*cos(phi)*(tanh(alpha*(x*cos(phi) - c*t + y*sin(phi)) - 1).^2 - 1));
% qy = -(alpha*sin(phi)*(tanh(alpha*(x*cos(phi) - c*t + y*sin(phi)) - 1).^2 - 1));
% UDG(:,1,:)=v;
% UDG(:,2,:)=qx;
% UDG(:,3,:)=qy;
% PDG(:,1,:)=u;

clear PDG;
UDG = l2eprojection(mesh,master,@exactsol,[c alpha phi],0,4);
PDG(:,1,:)=UDG(:,4,:);
UDG(:,4,:)=[];

UH = inituhat(master,mesh.elcon,UDG,app.ncu);


% HDG solver
for itime = 1:ntime
    t = dt*itime;
    fprintf('Timestep :  %d\n', itime);
        
    %[UDG,UH,~,PDG] = hdg_solve_dirk(master,mesh,app,UDG,UH,[],PDG,(itime-1)*dt,dt,nstage,torder);            
    [UDG,UH,PDG] = hdg_wave(master,mesh,app,UDG,PDG,(itime-1)*dt,dt,nstage,torder);
    %v=sin(pi*x).*sin(pi*y).*cos(sqrt(2)*pi*itime*dt);
%     v=-alpha*c*(tanh(alpha*(x*cos(phi) - c*t + y*sin(phi)) - 1).^2 - 1);
%     figure(1); clf; scaplot(mesh,UDG(:,1,:),[],2,1); axis off;         
%     figure(2); clf; scaplot(mesh,v-UDG(:,1,:),[],2,1); axis off;  
%         
%     %u=(1/(sqrt(2)*pi))*sin(pi*x).*sin(pi*y).*sin(sqrt(2)*pi*itime*dt);       
%     u = -tanh(alpha*(x*cos(phi) - c*t + y*sin(phi)) - 1);
%     figure(3); clf; scaplot(mesh,PDG(:,1,:),[],2,1); axis off;          
%     figure(4); clf; scaplot(mesh,u-PDG(:,1,:),[],2,1); axis off;   
%     ev = v-UDG(:,1,:);
%     eu = u-PDG(:,1,:);
%     [itime max(abs(ev(:))) max(abs(eu(:)))]
%     pause
end
UDG(:,4,:) = PDG;
e = calerror(UDG,mesh,master,@exactsol,[c k],t,[1 2 3 4]);
[e(1) sqrt(e(2)^2+e(3)^2) e(4)]

%err = calerror(UDG,mesh,master,func,param,time,ind)
return;

nn = [2 4 8 16 32]+1;
for porder = 1:2
    for ii=1:length(nn)
        ngrid = nn(ii);
        wave;             
        erra(ii,porder) = e(1); errq(ii,porder) = sqrt(e(2)^2+e(3)^2);
        errb(ii,porder) = e(4);
    end
end

%steepness of traveling wave
% alpha = 1;
% %orientation of traveling wave
% phi = pi/6;

% syms alpha phi x y t c;
% zeta = cos(phi)*x+sin(phi)*y;
% u = tanh(1-alpha*(zeta-t*c));
% ut = simplify(diff(u,'t'));
% ux = simplify(diff(u,'x'));
% uy = simplify(diff(u,'y'));
% utt = simplify(diff(ut,'t'));
% uxx = simplify(diff(ux,'x'));
% uyy = simplify(diff(uy,'y'));
% f = simplify(utt-c*c*uxx-c*c*uyy);

syms alpha phi x y t c;
u = cos(x+y-sqrt(2)*t);
ut = simplify(diff(u,'t'));
ux = simplify(diff(u,'x'));
uy = simplify(diff(u,'y'));
utt = simplify(diff(ut,'t'));
uxx = simplify(diff(ux,'x'));
uyy = simplify(diff(uy,'y'));
f = simplify(utt-uxx-uyy);



