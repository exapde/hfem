%load('meshmaster.mat');
%load('sole025delta10uinf20.mat');

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
E0 = Etmax; % reference electric field (V/m)
Ecstar = Ec/E0;
lb =  E0*Ecstar/1e6;
ub = (1+delta_peek)*lb;

dEcdr = -(47740/(a^2*(1/a)^(1/2)))/1e6;

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
     
E = sqrt(UDG(:,3,:).^2+(UDG(:,5,:)+ Etstar).^2);
E = E0*E/1e6;
[gradE,~] = getdudx(master, mesh, E);

figure(1); clf; scaplot(mesh,E,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/20);
set(gca,'xtick',[-0.4:0.2:0.4]/20);
set(gca,'ytick',[-0.4:0.2:0.4]/20);
set(gca, 'LooseInset', get(gca, 'TightInset'));

figure(2); clf; scaplot(mesh,(1e6*e0*E0/L0)*UDG(:,2,:),[0 0.5]*(1e6*e0*E0/L0),2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/2);
set(gca,'xtick',[-0.4:0.2:0.4]/2);
set(gca,'ytick',[-0.4:0.2:0.4]/2);
set(gca, 'LooseInset', get(gca, 'TightInset'));

x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
t = cart2pol(x,y);
i = find(t<0);
t(i) = 2*pi+t(i);
dxdr = cos(t);
dydr = sin(t);
dEdr = -(gradE(:,1,:).*dxdr + gradE(:,2,:).*dydr);
figure(3); clf; scaplot(mesh,dEdr,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/20);
set(gca,'xtick',[-0.4:0.2:0.4]/20);
set(gca,'ytick',[-0.4:0.2:0.4]/20);
set(gca, 'LooseInset', get(gca, 'TightInset'));








% [dgnodes,u] = normalfield(master,mesh,gradE,-5);
% [~,v] = normalelectricfield(master,mesh,UDG,UH,-5,0);
% 
% x = dgnodes(:,1,:);x=x(:);
% y = dgnodes(:,2,:);y=y(:);
% f = u(:,end,:);f=f(:);
% g = v(:,end,:);g=abs(g(:));
% t = cart2pol(x,y);
% i = find(t<0);
% t(i) = 2*pi+t(i);
% [tj,jj] = sort(t);

% figure(1);clf;
% fill([0 2*pi 2*pi 0],[lb lb ub ub],[0.4 0.4 0.4]*1.5);
% hold on;
% plot(tj,E0*g(jj)/1e6,'b-','LineWidth',1.5);
% xlabel('\theta (rad)','FontSize',20); ylabel('Electric field intensity (MV/m)','FontSize',20);
% set(gca,'FontSize',18); box on;
% axis tight;
% axis([0 2*pi min(E0*g(jj)/1e6)/1.05 max(E0*g(jj)/1e6)*1.05]);
% set(gca, 'LooseInset', get(gca, 'TightInset'));
% % fn = ['wire_electric_onwire_sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.png'];
% % print('-dpng',fn);
% 
% figure(2);clf;
% plot(tj,f(jj),'b-','LineWidth',1.5);
% xlabel('\theta (rad)','FontSize',20); ylabel('Electric field gradient (MV/m^2)','FontSize',20);
% set(gca,'FontSize',18); box on;
% axis tight;
% axis([0 2*pi min(f(jj))*1.05 max(f(jj))/1.05]);
% set(gca, 'LooseInset', get(gca, 'TightInset'));
% % fn = ['wire_electric_gradient_sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.png'];
% % print('-dpng',fn);
 
figure(1); clf; scaplot(mesh,E,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/20);
set(gca,'xtick',[-0.4:0.2:0.4]/20);
set(gca,'ytick',[-0.4:0.2:0.4]/20);
set(gca, 'LooseInset', get(gca, 'TightInset'));
%fn = ['wire_electric_zoom_sole0' num2str(Etmax/1e3) 'delta' num2str(delta_peek*1e4) 'uinf' num2str(windvel) '.png'];
%print('-dpng',fn);
% 
% figure(4); clf; scaplot(mesh,-gradE(:,1,:),[],2); 
% xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
% colorbar('FontSize',15); set(gca,'FontSize',18); box off;
% axis equal; axis tight; colormap jet; 
% axis([-0.4 0.4 -0.4 0.4]/20);
% set(gca,'xtick',[-0.4:0.2:0.4]/20);
% set(gca,'ytick',[-0.4:0.2:0.4]/20);
% set(gca, 'LooseInset', get(gca, 'TightInset'));
% 
% figure(5); clf; scaplot(mesh,-gradE(:,2,:),[],2); 
% xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
% colorbar('FontSize',15); set(gca,'FontSize',18); box off;
% axis equal; axis tight; colormap jet; 
% axis([-0.4 0.4 -0.4 0.4]/20);
% set(gca,'xtick',[-0.4:0.2:0.4]/20);
% set(gca,'ytick',[-0.4:0.2:0.4]/20);
% set(gca, 'LooseInset', get(gca, 'TightInset'));

