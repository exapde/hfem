setapplicationpath('FM/corona');

% See the paper Neimarlija2009.pdf and Medlin1998b.pdf
 
porder = 4;

% physical parameters
K = 1.5e-4; % (m^2/Vs)
D = 0.1;    %  (m^2/s)
e0 = 8.854e-12; % (F/m)

% lengthscale parameters  
a = 0.01; % radius of the wire  (m)
b = 15;   % distance from the wire to the plate (m)
L0 = 1;   % reference length scale (m)  

T0 = 1;   % reference time (s)
T1 = 20;  % time to reach maximum background electric (s)  
T2 = 30;  % final time (s)

% discharge electric strength
Ec = 3.1*(1 + 0.0308*sqrt(1/a))*1e6;
Etmax = 40*1e3; % maximum background electric (V/m)
E0 = Etmax; % reference electric field (V/m)

% nondimensional parameters
Kstar = K*E0*T0/L0;
Dstar = D*T0/L0;
estar = e0/e0;
Ecstar = Ec/E0;
Etstar = Etmax/E0;
T1star = T1/T0;
T2star = T2/T0;

% applied potential at the wire
phia = Etmax*b;  % (V)
phiastar = phia/(L0*E0);

% applied ion density at the wire
%rhoa = 2.39*1e-8; % (C/m^3)
rhoastar = L0*rhoa/(e0*E0);

% stabilization parameter
beta1 = -1/200;
beta2  = -1/200;
tau = 20;
param = {Kstar,Dstar,estar,Ecstar,Etstar,T1star,T2star,beta1,beta2,tau};
ui = [phiastar rhoastar];

hybrid = 'hdg';
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fboucorona';
app.fhat = 'fhat';
app.adjoint = 0;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve=1;
app.arg = param;
app.bcm = [6;8;7;9;1];
app.bcs = [[0 0];[0 0];[0 0];[0 0];ui];
app.bcd = [];
app.bcv = []; 

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
app.nch  = 2;                       % Number of componets of UH
app.nc   = app.nch*(app.nd+1);    % Number of componeents of UDG
app.ncu = 2;

app.appname = 'corona';
app.time = [];
app.dtfc = [];
app.alpha = [];
app.dtcoef = [0,1]';

mesh = mkmesh_circleinrect2(porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,app);

% UDG = initu(mesh,{1;1;0;0;0;0});
% UH = inituhat(master,mesh.elcon,UDG,app.ncu);

wvel = wind(mesh,windvel);
mesh.dgnodes(:,3,:) = wvel/(K*E0);
mesh.dgnodes(:,4,:) = 0;

[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);



Et = Etstar;     
Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
[~,u] = potentialfield(master,mesh,E,-5,1);    
Em = mean(abs(u(:)));  
[Em Ecstar]
abs(Em-Ecstar)/Ecstar
2*pi*a*K*Ec*rhoa*1e6


return;


t = linspace(0,2*pi,1000);
c = cos(t);


save solswind55.mat
tm = [tm rhoa*1e7];
tm(2:end)-tm(1:end-1)
tm(end)+tm(end)-tm(end-1)

w = 0:5:200;
for i=1:length(w)
    fn = ['solswind' num2str(w(i)) '.mat']
    load(fn);
    cur(i) = 2*pi*a*K*Ec*rhoa*1e6;
end

w = 0;
fn = ['solswind' num2str(w) '.mat']
load(fn);
figure(1); clf; scaplot(mesh,(L0*E0)*UDG(:,1,:),[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['ex5_potential_wind' num2str(w)];
print('-dpng',fn);

figure(2); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:),[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet;    
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['ex5_density_wind' num2str(w)];
print('-dpng',fn);

fn = ['solswind' num2str(w) '.mat']
load(fn);
figure(3); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:),[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet;    
axis([-0.25 0.75 -0.3 0.4])
% axis([-0.5 1.5 -0.6 0.8])
% axis([-1 3 -1.2 1.6])
% axis([-2 6 -2.4 3.2])
%set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['ex5_densityzoom_wind' num2str(w)];
print('-dpng',fn);

w = [0 10 20 50 100 200];
for i=1:length(w)
    fn = ['solswind' num2str(w(i)) '.mat']
    load(fn);
    figure(1); clf; scaplot(mesh,(L0*E0)*UDG(:,1,:),[0 8e5],2); 
    xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
    colorbar('FontSize',15); set(gca,'FontSize',18); box off;
    axis equal; axis tight; colormap jet; 
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    fn = ['ex5_potential_wind' num2str(w(i))];
    print('-dpng',fn);

    figure(2); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:),[0 1e-7],2); 
    xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
    colorbar('FontSize',15); set(gca,'FontSize',18); box off;
    axis equal; axis tight; colormap jet;    
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    fn = ['ex5_density_wind' num2str(w(i))];
    print('-dpng',fn);

    Et = Etstar;     
    Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
    figure(3); clf; 
    scaplot(mesh,E0*E,[0 0.5e5],2); 
    axis equal; axis tight; axis on; colormap jet;
    colorbar('FontSize',15);
    set(gca,'FontSize',18);
    xlabel('x (m)','FontSize',20);
    ylabel('y (m)','FontSize',20);
    box on;
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    fn = ['ex5_electric_wind' num2str(w(i))];
    print('-dpng',fn);
end

w = [0 10 20 50 100 200];
for i=1:length(w)
    fn = ['solswind' num2str(w(i)) '.mat']
    load(fn);
    figure(1); clf; scaplot(mesh,(L0*E0)*UDG(:,1,:),[0 6e5],2); 
    xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
    colorbar('FontSize',15); set(gca,'FontSize',18); box on;
    axis equal; axis tight; colormap jet; 
    axis([-0.25 0.75 -0.3 0.4])
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    fn = ['ex5_potentialzoom_wind' num2str(w(i))];
    print('-dpng',fn);

    figure(2); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:),[],2); 
    xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
    colorbar('FontSize',15); set(gca,'FontSize',18); box on;
    axis equal; axis tight; colormap jet;    
    axis([-0.25 0.75 -0.3 0.4])
    set(gca, 'LooseInset', get(gca, 'TightInset'));    
    fn = ['ex5_densityzoom_wind' num2str(w(i))];
    print('-dpng',fn);

    Et = Etstar;     
    Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; 
    figure(3); clf; 
    scaplot(mesh,E0*Ex,[-4e6 4e6],2); 
    axis equal; axis tight; axis on; colormap jet;
    colorbar('FontSize',15);
    set(gca,'FontSize',18);
    xlabel('x (m)','FontSize',20);
    ylabel('y (m)','FontSize',20);
    box on;
    axis([-0.25 0.75 -0.3 0.4])
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    fn = ['ex5_exzoom_wind' num2str(w(i))];
    print('-dpng',fn);
    
    figure(4); clf; 
    scaplot(mesh,E0*Ey,[-4e6 4e6],2); 
    axis equal; axis tight; axis on; colormap jet;
    colorbar('FontSize',15);
    set(gca,'FontSize',18);
    xlabel('x (m)','FontSize',20);
    ylabel('y (m)','FontSize',20);
    box on;
    axis([-0.25 0.75 -0.3 0.4])
    set(gca, 'LooseInset', get(gca, 'TightInset'));
    fn = ['ex5_eyzoom_wind' num2str(w(i))];
    print('-dpng',fn);    
end


figure(1); clf; 
meshplot(mesh,1); 
axis equal; axis tight; axis on; 
set(gca,'FontSize',18);
set(gca, 'LooseInset', get(gca, 'TightInset'));

w = 0:5:200;
figure(1); clf; 
plot(w,cur,'-ob','LineWidth',1.5);
axis tight; axis on; box on;
set(gca,'FontSize',18);
xlabel('Wind velocity (m/s)','FontSize',20);
ylabel('Corona current (\muA/m)','FontSize',20);

load(fn);
Et = Etstar;     
Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
figure(3); clf; scaplot(mesh,E0*E,[0 2e6],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet;    
axis([-0.25 0.75 -0.3 0.4])


