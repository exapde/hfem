% 
% % parameter study
% delta = [0.25 0.5 1 2]/100;
% Einf = [20 22.5 25 27.5 30 32.5 35 37.5 40]*1e3;
% uinf = [200:-10:10 5];
% 
% %wcase = 1;
% 
% if wcase == 1    
%     delta_peek = delta(1);
%     Etmax = Einf(kk);
%     betall = 0*uinf;
%     for jj = 1:length(uinf)
%         windvel = uinf(jj);     
%         if jj>1
%             beta = windvel*(betall(jj-1)/uinf(jj-1));
%         end        
%         iterwire;    
%         betall(jj) = app.arg{10}(1);    
%         fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
%         save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
%     end
% end
% 
% Einf = [25 30 35 40]*1e3;
% for kk=1:4
% uinf = [9 8 7 6];
% delta_peek = 0.25/100;
% Etmax = Einf(kk);
% fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(10) '.mat']
% load(fn);
% beta = app.arg{10}(1)*9/10;
% for jj = 1:length(uinf)
%     windvel = uinf(jj);     
%     iterwire;        
%     fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
%     save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
% end
% end
% 
% if wcase == 2
%     load('sole040delta25uinf20.mat');
%     beta = app.arg{10}(1);
%     windvel = 20;
%     Etmax = 40*1e3;    
%     delta = [0.5 1 2 4]/100; 
%     for jj = 1:length(delta)                
%         delta_peek = delta(jj);
%         iterwire;            
%         fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
%         save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
%     end    
% end
% 
% fn = ['sole0' num2str(25) 'delta' num2str(10) 'uinf' num2str(20) '.mat']
% load(fn);
% beta = app.arg{10}(1);   
% windvel = 20;
% delta_peek = 0.05/100;
% Etmax = 25*1e3;
% iterwire;            
% fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
% save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
% 
% fn = ['sole0' num2str(30) 'delta' num2str(10) 'uinf' num2str(20) '.mat']
% load(fn);
% beta = app.arg{10}(1);   
% windvel = 20;
% delta_peek = 0.05/100;
% Etmax = 30*1e3;
% iterwire;            
% fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
% save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
% 
% fn = ['sole0' num2str(35) 'delta' num2str(10) 'uinf' num2str(20) '.mat']
% load(fn);
% beta = app.arg{10}(1);   
% windvel = 20;
% delta_peek = 0.05/100;
% Etmax = 35*1e3;
% iterwire;            
% fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
% save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
% 
% fn = ['sole0' num2str(40) 'delta' num2str(10) 'uinf' num2str(20) '.mat']
% load(fn);
% beta = app.arg{10}(1);   
% windvel = 20;
% delta_peek = 0.05/100;
% Etmax = 40*1e3;
% iterwire;            
% fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
% save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
% 
% 
% if wcase==3
%     load('sole040delta50uinf20.mat');
%     beta = app.arg{10}(1)/2;
%     windvel = 20;
%     delta_peek = 0.5/100;
%     Einf = [37.5 32.5 27.5 22.5]*1e3;
%     for jj = 1:length(Einf)                
%         Etmax = Einf(jj);
%         iterwire;            
%         fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
%         save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
%     end       
% end
% 
% load('sole040delta50uinf200.mat');
% beta = app.arg{10}(1);
% windvel = 200;
% delta_peek = 0.5/100;
% Einf = [45 50 55 60]*1e3;
% for jj = 1:length(Einf)                    
%     Etmax = Einf(jj);
%     iterwire;            
%     fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
%     save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
% end       
% 
% if wcase==4
%     load('sole032.5delta50uinf20.mat');
%     beta = app.arg{10}(1);
%     windvel = 20;
%     delta_peek = 1/100;
%     Einf = [35 30]*1e3;
%     for jj = 1:length(Einf)                
%         Etmax = Einf(jj);
%         iterwire;            
%         fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
%         save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
%     end       
%     
%     load('sole032.5delta50uinf20.mat');
%     beta = app.arg{10}(1);
%     windvel = 20;
%     delta_peek = 2/100;
%     Einf = [35 30]*1e3;
%     for jj = 1:length(Einf)                
%         Etmax = Einf(jj);
%         iterwire;            
%         fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
%         save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
%     end       
%         
%     load('sole032.5delta50uinf200.mat');
%     beta = app.arg{10}(1);
%     windvel = 200;
%     delta_peek = 1/100;
%     Einf = [35 30]*1e3;
%     for jj = 1:length(Einf)                
%         Etmax = Einf(jj);
%         iterwire;            
%         fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
%         save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
%     end       
%         
%     load('sole032.5delta50uinf200.mat');
%     beta = app.arg{10}(1);
%     windvel = 200;
%     delta_peek = 2/100;
%     Einf = [35 30]*1e3;
%     for jj = 1:length(Einf)                
%         Etmax = Einf(jj);
%         iterwire;            
%         fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
%         save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
%     end       
%             
%     windvel = 200;
%     delta_peek = 0.5/100;
%     Einf = [40 37.5 35 32.5 30 27.5 25 22.5]*1e3;
%     for jj = 8:length(Einf)                
%         Etmax = Einf(jj);
%         iterwire;            
%         fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
%         save(fn,'UDG','UH','app','Etmax','delta_peek','windvel');    
%     end       
% end

load meshmaster.mat;

delta_peek = 0.25/100;
e0 = 8.854e-12;
uinf = [200:-10:10 9:-1:0];
Einf = [25 30 35 40]*1e3;
Iwind = zeros(length(uinf),length(Einf));
Icren = 0*Iwind;
for j = 1:length(Einf)
for i = 1:length(uinf)
    windvel = uinf(i);
    Etmax = Einf(j);
    E0 = Etmax;
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    load(fn);
    Iwind(i,j) = windvel*e0*E0;
    Icren(i,j) = current(master,mesh,UDG,-5,Etmax);
end
end
figure(1);clf;
plot(uinf,1e6*Icren(:,1),'b-',uinf,1e6*Icren(:,2),'r-',...
     uinf,1e6*Icren(:,3),'g-',uinf,1e6*Icren(:,4),'k-',...
     'MarkerSize',8,'LineWidth',1.5);
xlabel('u_\infty (m/s)','FontSize',20); 
ylabel('Corona current per unit length (\muA/m)','FontSize',20);
legend({'$E_\infty = 25 \mbox{ KV/m}$','$E_\infty = 30 \mbox{ KV/m}$','$E_\infty = 35 \mbox{ KV/m}$','$E_\infty = 40 \mbox{ KV/m}$'},'interpreter','latex','location','NorthWest','FontSize',18);
set(gca,'FontSize',18); box on;
axis tight;
axis([0 200 0 70]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_current_vs_velocity2'];
print('-dpng',fn);


delta_peek = 0.5/100;
e0 = 8.854e-12;
uinf = [200];
Einf = [60 55 50 45 40 37.5 35 32.5 30 27.5 25 22.5]*1e3;
Iwind = zeros(length(uinf),length(Einf));
Icren = 0*Iwind;
for j = 1:length(Einf)
for i = 1:length(uinf)
    windvel = uinf(i);
    Etmax = Einf(j);
    E0 = Etmax;
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    load(fn);
    Iwind(i,j) = windvel*e0*E0;
    Icren(i,j) = current(master,mesh,UDG,-5,Etmax);
end
end
figure(1);clf;
plot(Einf/1e3,1e6*Iwind(1,:),'b-o',Einf/1e3,1e6*Icren(1,:),'r-s',...     
     'MarkerSize',8,'LineWidth',1.5);
xlabel('E_\infty (KV/m)','FontSize',20); 
ylabel('Corona current per unit length (\muA/m)','FontSize',20);
legend({'$I = \epsilon_0 E_\infty u_\infty$','$I = \int_0^{2\pi} a \rho \mu E  d \theta$'},'interpreter','latex','location','NorthWest','FontSize',18);
set(gca,'FontSize',18); box on;
axis tight;
%axis([0 200 0 75]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_current_vs_einf_uinf200'];
print('-dpng',fn);

figure(2);clf;
plot(Einf/1e3,1e6*Iwind(2,:),'b-o',Einf/1e3,1e6*Icren(2,:),'r-s',...
     'MarkerSize',8,'LineWidth',1.5);
xlabel('E_\infty (KV/m)','FontSize',20); 
ylabel('Corona current per unit length (\muA/m)','FontSize',20);
legend({'$I = \epsilon_0 E_\infty u_\infty$','$I = \int_0^{2\pi} a \rho \mu E  d \theta$'},'interpreter','latex','location','NorthWest','FontSize',18);
set(gca,'FontSize',18); box on;
axis tight;
%axis([0 200 0 75]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_current_vs_einf_uinf20'];
print('-dpng',fn);


windvel = 20;   
delta = [0.05 0.1 0.25 0.5 1 2 4]/100; 
e0 = 8.854e-12;
Einf = [25 40]*1e3;
Iwind = zeros(length(delta),length(Einf));
Icren = 0*Iwind;
for j = 1:length(Einf)
for i = 1:length(delta)
    delta_peek = delta(i);
    Etmax = Einf(j);
    E0 = Etmax;
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    load(fn);
    Iwind(i,j) = windvel*e0*E0;
    Icren(i,j) = current(master,mesh,UDG,-5,Etmax);
end
end
figure(1);clf;
plot(delta,1e6*Icren(:,1),'b-o',delta,1e6*Icren(:,2),'r-s',...
     'MarkerSize',8,'LineWidth',1.5);
xlabel('\delta_E','FontSize',20); 
ylabel('Corona current per unit length (\muA/m)','FontSize',20);
legend({'$E_\infty = 25 \mbox{ KV/m}$','$E_\infty = 40 \mbox{ KV/m}$'},'interpreter','latex','location','NorthWest','FontSize',18);
set(gca,'FontSize',18); box on;
axis tight;
axis([0 4/100 0 6]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_current_vs_delta'];
print('-dpng',fn);


delta_peek = 0.25/100;
uinf = [5 10 20 50 100 200];
Einf = [25 30 35 40]*1e3;
for j = 1:length(Einf)
for i = 1:length(uinf)
windvel = uinf(i);
Etmax = Einf(j);
E0 = Etmax;
Ecstar = Ec/E0;
Etstar = Etmax/E0;
lb =  E0*Ecstar/1e6;
ub = (1+delta_peek)*lb;
fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
load(fn);
Et = Etstar;     
Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
[dgnodes,u] = normalelectricfield(master,mesh,UDG,UH,-5,0);
f = u(:,end,:);f=abs(f(:));
Ef(:,i) = f; 
end
x = dgnodes(:,1,:);x=x(:);
y = dgnodes(:,2,:);y=y(:);
t = cart2pol(x,y);
ind = t<0;
t(ind) = 2*pi+t(ind);
[tj,jj] = sort(t);

figure(j);clf;
fill([0 2*pi 2*pi 0],[lb lb ub ub],[0.4 0.4 0.4]*1.5); 
hold on;
xpl = plot(tj,E0*Ef(jj,1)/1e6,'b-',...
     tj,E0*Ef(jj,2)/1e6,'r-',...
     tj,E0*Ef(jj,3)/1e6,'g-',...
     tj,E0*Ef(jj,4)/1e6,'k-',...
     tj,E0*Ef(jj,5)/1e6,'m-',...
     tj,E0*Ef(jj,6)/1e6,'y-',...
     'LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis([0 2*pi 3.95 4.07]);
if j == 1
legend(xpl,{'$u_\infty = 5\mbox{ m/s}$','$u_\infty = 10\mbox{ m/s}$','$u_\infty = 20\mbox{ m/s}$',...
        '$u_\infty = 50\mbox{ m/s}$', '$u_\infty = 100\mbox{ m/s}$', '$u_\infty = 200\mbox{ m/s}$'},...
        'interpreter','latex','location','Best','FontSize',18);
end
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_electric_onwire_sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) '.png'];
print('-dpng',fn);
end


delta_peek = 0.05/100;
e0 = 8.854e-12;
lb =  Ec/1e6;
ub = (1+delta_peek)*lb;
uinf = 20;
Einf = [25 30 35 40]*1e3;
for i = 1:length(uinf)
for j = 1:length(Einf)
    windvel = uinf(i);
    Etmax = Einf(j);    
    fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
    load(fn);    
    Et = 1;     
    Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
    [dgnodes,u] = normalelectricfield(master,mesh,UDG,UH,-5,0);
    f = u(:,end,:);f=abs(f(:));
    Ef(:,j) = Etmax*f;     
end
x = dgnodes(:,1,:);x=x(:);
y = dgnodes(:,2,:);y=y(:);
t = cart2pol(x,y);
ind = t<0;
t(ind) = 2*pi+t(ind);
[tj,jj] = sort(t);

figure(i);clf;
fill([0 2*pi 2*pi 0],[lb lb ub ub],[0.4 0.4 0.4]*1.5); 
hold on;
xpl = plot(tj,Ef(jj,1)/1e6,'b-',...
     tj,Ef(jj,2)/1e6,'r-',...
     tj,Ef(jj,3)/1e6,'g-',...
     tj,Ef(jj,4)/1e6,'k-',...
     'LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis([0 2*pi 3.95 4.15]);
if delta_peek == 0.25/100
legend(xpl,{'$E_\infty = 25\mbox{ KV/m}$','$E_\infty = 30\mbox{ KV/m}$',...
            '$E_\infty = 35\mbox{ KV/m}$','$E_\infty = 40\mbox{ KV/m}$'},...
        'interpreter','latex','location','Best','FontSize',18);
end
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_electric_onwire_sol' 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.png'];
print('-dpng',fn);
end

delta_peek = 0.25/100;
Etmax = 25*1e3;
uinf = [5 10 20 50 100 200];
dd = [128 64 32 16 8 4]/4;
for i = 1:length(uinf)
windvel = uinf(i);
fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
load(fn);    
E0 = Etmax;
tm = e0*E0/L0*UDG(:,2,:)*1e6;
figure(i); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:)*1e6,[0 max(tm(:))*0.98],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.5 0.5 -0.5 0.5]/5*dd(i));
set(gca,'xtick',[-0.5:0.25:0.5]/5*dd(i));
set(gca,'ytick',[-0.5:0.25:0.5]/5*dd(i));
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_densityzoom_e0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel)];
print('-dpng',fn);
end

for i = 1:length(uinf)
windvel = uinf(i);
fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
load(fn);    
E0 = Etmax;
figure(i); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:)*1e6,[0 0.1e-6/2]*1e6,2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet;    
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_density_e0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel)];
print('-dpng',fn);
end

for i = 1:length(uinf)
windvel = uinf(i);
fn = ['sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
load(fn);    
E0 = Etmax;
Et = 1;     
Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
figure(i); clf; scaplot(mesh,E0*Ey/1e6,[-1e4 1e4]/1e6,2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet;    
%set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_electric_e0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel)];
print('-dpng',fn);
end


