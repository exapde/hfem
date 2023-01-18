Et = Etstar;     
Ex = UDG(:,3,:); Ey = UDG(:,5,:)+Et; E = sqrt(Ex.^2+Ey.^2);
% [dgnodes,u] = potentialfield(master,mesh,E,-5,1);    
[dgnodes,u] = normalelectricfield(master,mesh,UDG,UH,-5,0);
[~,rho] = potentialfield(master,mesh,UDG(:,2,:),-5,1);    

E0 = Etmax; % reference electric field (V/m)
Ecstar = Ec/E0;
Etstar = Etmax/E0;
lb =  E0*Ecstar/1e6;
ub = (1+delta_peek)*lb;

x = dgnodes(:,1,:);x=x(:);
y = dgnodes(:,2,:);y=y(:);
f = u(:,end,:);f=abs(f(:));
q = rho(:,1,:);q=q(:);
t = cart2pol(x,y);
i = find(t<0);
t(i) = 2*pi+t(i);
[tj,jj] = sort(t);
[~,ii] = sort(abs(f-Ecstar));
[~,i1]=sort(t(ii(1:8)));
ii = ii(i1([1 end]));

figure(1);clf;
%fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(f) min(f) max(f) max(f)]*E0/1e6,[0.4 0.4 0.4]*1.5);
%fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[f(ii(1)) f(ii(1)) max(f) max(f)]*E0/1e6,[0.4 0.4 0.4]*1.5);
%hold on;
%fill([0 2*pi 2*pi 0],[f(ii(1)) f(ii(1)) max(f) max(f)]*E0/1e6,[0.4 0.4 0.4]*1.5);
fill([0 2*pi 2*pi 0],[lb lb ub ub],[0.4 0.4 0.4]*1.5);
hold on;
plot(tj,E0*f(jj)/1e6,'b-','LineWidth',1.5);
% plot(tj,(1+delta_peek)*E0*Ecstar*ones(size(t))/1e6,'-k','LineWidth',1.5);
% plot(tj,(1-delta_peek)*E0*Ecstar*ones(size(t))/1e6,'-k','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
axis([0 2*pi min(E0*f(jj)/1e6)/1.1 max(E0*f(jj)/1e6)*1.05]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['wire_electric_onwire_sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.mat'];
print('-dpng',fn);


figure(2);clf;
%fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(q) min(q) max(q) max(q)]*e0*E0/L0*1e6,[0.4 0.4 0.4]*1.5);
hold on;
plot(tj,1e6*e0*E0/L0*q(jj),'-b','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Charge density (\muC/m)','FontSize',20);
set(gca,'FontSize',18); box on;
%axis([0 2*pi 0 15]);
set(gca, 'LooseInset', get(gca, 'TightInset'));

q = app.arg{10}(1)*(E0*f/1e6-lb)*(Ecstar/lb)/2;
q(q<=0) = 0;
figure(3);clf;
%fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(q) min(q) max(q) max(q)]*e0*E0/L0*1e6,[0.4 0.4 0.4]*1.5);
hold on;
plot(tj,1e6*e0*E0/L0*q(jj),'-b','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Charge density (\muC/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis([0 2*pi 0 15]);
set(gca, 'LooseInset', get(gca, 'TightInset'));

return;

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
[ta,i1]=sort(t(ii(1:8)));
ii = ii(i1([1 end]));

figure(1);clf;
fill([0 2*pi 2*pi 0],[(1-delta_peek)*E0*Ecstar/1e6 (1-delta_peek)*E0*Ecstar/1e6 (1+delta_peek)*E0*Ecstar/1e6 (1+delta_peek)*E0*Ecstar/1e6],[0.4 0.4 0.4]*1.5);
hold on;
plot([t(ii(1)) t(ii(1))],[3.96 4.065],'--k','LineWidth',1.5);
plot([t(ii(2)) t(ii(2))],[3.96 4.065],'--k','LineWidth',1.5);
plot(tj,E0*f(jj)/1e6,'b-','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
axis([0 2*pi 3.96 4.065]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_electric_onwire_wind' num2str(windvel)];
print('-dpng',fn);

figure(2);clf;
%fill([t(ii(1)) t(ii(2)) t(ii(2)) t(ii(1))],[min(q) min(q) max(q) max(q)]*e0*E0/L0*1e6,[0.4 0.4 0.4]*1.5);
hold on;
plot([t(ii(1)) t(ii(1))],[0 dmax],'--k','LineWidth',1.5);
plot([t(ii(2)) t(ii(2))],[0 dmax],'--k','LineWidth',1.5);
plot(tj,1e6*e0*E0/L0*q(jj),'-b','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); ylabel('Charge density (\muC/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
axis([0 2*pi 0 dmax]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_density_onwire_wind' num2str(windvel)];
print('-dpng',fn);

[dgnodes,u] = potentialfield(master,mesh,E,-1,1);    
x = dgnodes(:,1,:);x=x(:);
y = dgnodes(:,2,:);y=y(:);
f = u(:,1,:);f=f(:);
[x,jj] = sort(x);
f = f(jj);
figure(3);clf;
plot(x,E0*f/1e6,'b-','LineWidth',1.5);
xlabel('x (m)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
axis([-100 100 0.005 0.04]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_electric_ground_wind' num2str(windvel)];
print('-dpng',fn);

[dgnodes,u] = potentialfield(master,mesh,E,-2,1);    
x = dgnodes(:,1,:);x=x(:);
y = dgnodes(:,2,:);y=y(:);
f = u(:,1,:);f=f(:);
[y,jj] = sort(y);
f = f(jj);
figure(4);clf;
plot(y,E0*f/1e6,'b-','LineWidth',1.5);
xlabel('y (m)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
set(gca,'FontSize',18); box on;
axis tight;
axis([-15 100 0.005 0.05]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_electric_rightbnd_wind' num2str(windvel)];
print('-dpng',fn);

figure(5); clf; scaplot(mesh,(L0*E0)*UDG(:,1,:)/1e6,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_potential_wind' num2str(windvel)];
print('-dpng',fn);

figure(6); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:)*1e6,[0 0.2e-6]*1e6,2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet;    
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_density_wind' num2str(windvel)];
print('-dpng',fn);

figure(7); clf; scaplot(mesh,E0*E/1e6,[0 0.5e5]/1e6,2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet;    
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_electric_wind' num2str(windvel)];
print('-dpng',fn);

figure(8); clf; scaplot(mesh,(L0*E0)*UDG(:,1,:)/1e6,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_potentialzoom_wind' num2str(windvel)];
print('-dpng',fn);

figure(9); clf; scaplot(mesh,e0*E0/L0*UDG(:,2,:)*1e6,[0 dmax],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.5 0.5 -0.5 0.5]/5*100);
set(gca,'xtick',[-0.5:0.25:0.5]/5*100);
set(gca,'ytick',[-0.5:0.25:0.5]/5*100);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_densityzoom_wind' num2str(windvel)];
print('-dpng',fn);

figure(10); clf; scaplot(mesh,E0*UDG(:,3,:)/1e6,[-4e6 4e6]/1e6,2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_exzoom_wind' num2str(windvel)];
print('-dpng',fn);

figure(11); clf; scaplot(mesh,E0*UDG(:,5,:)/1e6,[-4e6 4e6]/1e6,2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_eyzoom_wind' num2str(windvel)];
print('-dpng',fn);

figure(12); clf; scaplot(mesh,E0*E/1e6,[1e6 4.e6]/1e6,2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_electriczoom_wind' num2str(windvel)];
print('-dpng',fn);


figure(13); clf; scaplot(mesh,(L0*E0)*(UDG(:,1,:)+mesh.dgnodes(:,2,:))/1e6,[0 4.5],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_totalpotential_wind' num2str(windvel)];
print('-dpng',fn);

return;



U = [0 5 10 20 50 100 200];
for i = 1:length(U)
    windvel = U(i);
    fn = ['solfwind' num2str(windvel) 'delta015.mat'];
    load(fn);
    Iwind(i) = windvel*e0*E0;
    c1=rhoparam(1);c2=rhoparam(2);c3=rhoparam(3);
    Idens(i) = a*K*(Ec)*((e0*E0/L0)*c1)*(sqrt(pi)/2)*(c3)*(erf(c2/c3)+ erf((2*pi-c2)/c3));
    Icren(i) = current(master,mesh,UDG,-5);
end

figure(3);clf;
plot(U,1e6*Iwind,'b-o',U,1e6*Icren,'r-s','MarkerSize',8,'LineWidth',1.5);
xlabel('u_\infty (m/s)','FontSize',20); 
ylabel('Corona current per unit length (\mu A/m)','FontSize',20);
legend({'$I_1 = \epsilon_0 E_\infty u_\infty$','$I_2 = \int_0^{2\pi} a \rho \mu |E|  d \theta$'},'interpreter','latex','location','NorthWest','FontSize',18);
set(gca,'FontSize',18); box on;
axis tight;
axis([0 200 0 75]);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire_corona_current'];
print('-dpng',fn);




% figure(1); clf; 
% meshplot(mesh,1); 
% axis equal; axis tight; axis on; 
% set(gca,'FontSize',18);
% set(gca, 'LooseInset', get(gca, 'TightInset'));
% axis([-0.4 0.4 -0.4 0.4]/10);
% fn = ['groundedwire_meshzoom'];
% print('-dpng',fn);
% 
