setapplicationpath('FM/corona');

% See the paper Neimarlija2009.pdf and Medlin1998b.pdf

% fn = ['solfwind' num2str(windvel) '.mat']; 
% load(fn);

porder = 6;

% physical parameters
K = 1.5e-4; % (m^2/Vs)
D = 0.0005;    %  (m^2/s)
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
delta_peek = 0.15/100;

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
%rhoa = 0.045775*1.5796e-04; % (C/m^3)
% rhoastar = L0*rhoa/(e0*E0);
% rhoparam = [L0*rhoa/(e0*E0) 0.9835*pi 0.2495];

% stabilization parameter
beta1 = -1/200;
beta2  = -1/200;
tau = max(10,windvel);
param = {Kstar,Dstar,estar,Ecstar,Etstar,T1star,T2star,beta1,beta2,rhoparam,tau};
ui = [phiastar rhoparam(1)];

hybrid = 'hdg';
app.source = 'source';
app.flux = 'flux';
app.fbou = 'fboucorona2';
app.fhat = 'fhat';
app.adjoint = 0;
app.denseblock = 0;
app.hybrid = hybrid;
app.localsolve=1;
app.arg = param;
if windvel>0
    app.bcm = [10;8;8;9;1];
else
    app.bcm = [10;9;8;9;1];
end
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

mesh = mkmesh_circleinrect5(porder);
master = mkmaster(mesh,2*porder);
[master,mesh] = preprocess(master,mesh,app);
% 
% % wvel = wind(mesh,windvel);
% % mesh.dgnodes(:,3,:) = wvel/(K*E0);
% % mesh.dgnodes(:,4,:) = 0;
% 
[u,v] = potentialvelocity(mesh.dgnodes(:,1,:),mesh.dgnodes(:,2,:),windvel,0.01);
mesh.dgnodes(:,3,:) = u/(K*E0);
mesh.dgnodes(:,4,:) = v/(K*E0);
% load('solfwind200delta015.mat')
% mesh.dgnodes(:,5,:) = UDG(:,2,:);
% app.arg{10}(1) = 0;
% app.bcm = [10;8;8;9;1];
% app.bcs = [[0 0];[0 0];[0 0];[0 0];ui*0];
% 
% 
% UDG = initu(mesh,{1;rhoparam(1);0;0;0;0});
% UH = inituhat(master,mesh.elcon,UDG,app.ncu);

[UDG,UH] = hdg_solve(master,mesh,app,UDG,UH,[]);

% 
% figure(1); clf;
% hold on;
% for k = 1:ns
%     plot3(dgnodes(:,1,k),dgnodes(:,2,k),(u(:,1,k)),'-k');
% end
% hold off;
% axis normal;
% axis tight;
% view(3);
% 
% figure(2); clf;
% hold on;
% for k = 1:ns
%     plot3(dgnodes(:,1,k),dgnodes(:,2,k),(u(:,1,k)-Ecstar),'-k');
% end
% hold off;
% axis normal;
% axis tight;
% view(3);
% 
% figure(3); clf;
% hold on;
% for k = 1:ns
%     plot3(dgnodes(:,1,k),dgnodes(:,2,k),abs(u(:,1,k)-Ecstar),'-k');
% end
% hold off;
% axis normal;
% axis tight;
% view(3);

Et = Etstar;     
Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
[dgnodes,u] = potentialfield(master,mesh,E,-5,1);    
[~,rho] = potentialfield(master,mesh,UDG(:,2,:),-5,1);    

x = dgnodes(:,1,:);x=x(:);
y = dgnodes(:,2,:);y=y(:);
f = u(:,1,:);f=f(:);
q = rho(:,1,:);q=q(:);
t = cart2pol(x,y);
i = find(t<0);
t(i) = 2*pi+t(i);
[tj,jj] = sort(t);
[~,ii] = sort(abs(f-(1-delta_peek)*Ecstar));
[~,i1]=sort(t(ii(1:8)));
ii = ii(i1([1 end]));

figure(1);clf;
%fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(f) min(f) max(f) max(f)]*E0/1e6,[0.4 0.4 0.4]*1.5);
%fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[f(ii(1)) f(ii(1)) max(f) max(f)]*E0/1e6,[0.4 0.4 0.4]*1.5);
%hold on;
%fill([0 2*pi 2*pi 0],[f(ii(1)) f(ii(1)) max(f) max(f)]*E0/1e6,[0.4 0.4 0.4]*1.5);
fill([0 2*pi 2*pi 0],[(1-delta_peek)*E0*Ecstar/1e6 (1-delta_peek)*E0*Ecstar/1e6 (1+delta_peek)*E0*Ecstar/1e6 (1+delta_peek)*E0*Ecstar/1e6],[0.4 0.4 0.4]*1.5);
hold on;
plot(tj,E0*f(jj)/1e6,'b-',tj,E0*Ecstar*ones(size(t))/1e6,'--k','LineWidth',1.5);
plot(tj,(1+delta_peek)*E0*Ecstar*ones(size(t))/1e6,'-k','LineWidth',1.5);
plot(tj,(1-delta_peek)*E0*Ecstar*ones(size(t))/1e6,'-k','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));

figure(2);clf;
fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(q) min(q) max(q) max(q)]*e0*E0/L0*1e6,[0.4 0.4 0.4]*1.5);
hold on;
plot(tj,1e6*e0*E0/L0*q(jj),'-b','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Charge density (\muC/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));

c0 = E0*Ecstar/1e6;
c1 = max(E0*f(jj)/1e6);
100*abs(c1-c0)/c0

[t(ii(1)) t(ii(2))]
[1e6*e0*E0/L0*q(ii(1)) 1e6*e0*E0/L0*q(ii(2))]

return;

figure(1);clf;
plot(tj,E0*f(jj)/1e6,'b-','LineWidth',1.5);


t = linspace(0,2*pi,1000);
a = pi; b = 0.2;
z = c*exp(-(t-a).^2/(b^2))+ exp(-(t-2*pi).^2/(b^2));
figure(1); clf; plot(t,z);

% z1 = (1/sqrt(2*pi*b^2))*exp(-(t-a)^2)/(2*b^2);
% z2 = c1*exp(-(t-c2)^2)/(c3^2) = c1*(sqrt(pi*c3^2)/sqrt(pi*c3^2))*exp(-(t-c2)^2)/(2*(c3/sqrt(2))^2);
% area of z2 = c1*c3*sqrt(2*pi)


windv=[160:-10:50];
const=[0.8:-0.05:0.25];
for kk=1:12; 
    windvel = windv(kk); coa = const(kk); 
    [kk windvel coa]
    fwirewind; 
    fn = ['solfwind' num2str(windvel) '.mat']; 
    save(fn,'app', 'UDG', 'UH', 'rhoparam', 'rhoa'); 
end

for i =200:-10:10
    fn = ['solfwind' num2str(i) '.mat']; 
    load(fn);
    fn = ['solfwind' num2str(i*6) '.mat']; 
    save(fn,'app', 'UDG', 'UH', 'rhoparam'); 
end


Et = Etstar;     
Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
[dgnodes,u] = potentialfield(master,mesh,E,-5,1);    
[~,rho] = potentialfield(master,mesh,UDG(:,2,:),-5,1);    

x = dgnodes(:,1,:);x=x(:);
y = dgnodes(:,2,:);y=y(:);
f = u(:,1,:);f=f(:);
q = rho(:,1,:);q=q(:);
t = cart2pol(x,y);
i = find(t<0);
t(i) = 2*pi+t(i);
[tj,jj] = sort(t);
[~,ii] = sort(abs(f-Ecstar));
[~,i1]=sort(t(ii(1:8)));
ii = ii(i1([1 end]));

figure(1);clf;
fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(f) min(f) max(f) max(f)]*E0/1e6,[0.4 0.4 0.4]*1.5);
hold on;
plot(tj,E0*f(jj)/1e6,'b-',tj,E0*Ecstar*ones(size(t))/1e6,'--k','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_electric_onwire_wind' num2str(windvel)];
print('-dpng',fn);

figure(2);clf;
fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(q) min(q) max(q) max(q)]*e0*E0/L0*1e6,[0.4 0.4 0.4]*1.5);
hold on;
plot(tj,1e6*e0*E0/L0*q(jj),'-b','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Charge density (\muC/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_density_onwire_wind' num2str(windvel)];
print('-dpng',fn);





% figure(1);clf;
% plot(t,f,'*',t,Ecstar*ones(size(t)));
% 
% figure(2);clf;
% i1 = find(y<0);
% g = zeros(size(x));
% g(i1) = exp(-abs(1+x(i1)/0.01).^1.0/(0.25^2));
% i2 = find(y>=0);
% g(i2) = exp(-abs(1+x(i2)/0.01).^1.0/(0.25^2));
% h = exp(-(t-0.969*pi).^2/(0.35^2));
% plot(t,g,'ob',t,h,'xr');

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

figure(2); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:),[0 1e-6],2); 
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


