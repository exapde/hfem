
load('solfwind200delta015.mat')
U = E0*UDG(:,1,:)/1e6; Ex = E0*UDG(:,3,:)/1e6; Ey = E0*UDG(:,5,:)/1e6;

load('solfwindappliedpotential.mat')
U1 = E0*UDG(:,1,:)/1e6; E1x = E0*UDG(:,3,:)/1e6; E1y = E0*UDG(:,5,:)/1e6;

load('solfwinddensitysourcewind200.mat');
U2 = E0*UDG(:,1,:)/1e6; E2x = E0*UDG(:,3,:)/1e6; E2y = E0*UDG(:,5,:)/1e6;

% figure(1); clf; scaplot(mesh,U-U1-U2,[],2); 
% xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
% colorbar('FontSize',15); set(gca,'FontSize',18); box off;
% axis equal; axis tight; colormap jet; 
% set(gca, 'LooseInset', get(gca, 'TightInset'));
% 
% figure(2); clf; scaplot(mesh,Ex-E1x-E2x,[],2); 
% xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
% colorbar('FontSize',15); set(gca,'FontSize',18); box off;
% axis equal; axis tight; colormap jet; 
% set(gca, 'LooseInset', get(gca, 'TightInset'));
% 
% figure(3); clf; scaplot(mesh,Ey-E1y-E2y,[],2); 
% xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
% colorbar('FontSize',15); set(gca,'FontSize',18); box off;
% axis equal; axis tight; colormap jet; 
% set(gca, 'LooseInset', get(gca, 'TightInset'));

figure(1); clf; scaplot(mesh,Ex,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_Exzoom_' num2str(windvel)];
print('-dpng',fn);

figure(2); clf; scaplot(mesh,E1x,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_E1xzoom_' num2str(windvel)];
print('-dpng',fn);

figure(3); clf; scaplot(mesh,E2x,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_E2xzoom_' num2str(windvel)];
print('-dpng',fn);

figure(4); clf; scaplot(mesh,Ey,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_Eyzoom_' num2str(windvel)];
print('-dpng',fn);

figure(5); clf; scaplot(mesh,E1y,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_E1yzoom_' num2str(windvel)];
print('-dpng',fn);

figure(6); clf; scaplot(mesh,E2y,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_E2yzoom_' num2str(windvel)];
print('-dpng',fn);

figure(7); clf; scaplot(mesh,U,[13 15]*E0/1e6,2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_Uzoom_' num2str(windvel)];
print('-dpng',fn);

figure(8); clf; scaplot(mesh,U1,[11 15]*E0/1e6,2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_U1zoom_' num2str(windvel)];
print('-dpng',fn);

figure(9); clf; scaplot(mesh,U2,[0 1.5]*E0/1e6,2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
axis([-0.4 0.4 -0.4 0.4]/8);
set(gca,'xtick',[-0.4:0.2:0.4]/8);
set(gca,'ytick',[-0.4:0.2:0.4]/8);
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_U2zoom_' num2str(windvel)];
print('-dpng',fn);

figure(10); clf; scaplot(mesh,U,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_U_' num2str(windvel)];
print('-dpng',fn);

figure(11); clf; scaplot(mesh,U1,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_U1_' num2str(windvel)];
print('-dpng',fn);

figure(12); clf; scaplot(mesh,U2,[],2); 
xlabel('x (m)','FontSize',20); ylabel('y (m)','FontSize',20);
colorbar('FontSize',15); set(gca,'FontSize',18); box off;
axis equal; axis tight; colormap jet; 
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_U2_' num2str(windvel)];
print('-dpng',fn);


[dgnodes,ex] = potentialfield(master,mesh,Ex,-5,1);    
[~,ey] = potentialfield(master,mesh,Ey,-5,1);    
[~,e1x] = potentialfield(master,mesh,E1x,-5,1);    
[~,e1y] = potentialfield(master,mesh,E1y,-5,1);    
[~,e2x] = potentialfield(master,mesh,E2x,-5,1);    
[~,e2y] = potentialfield(master,mesh,E2y,-5,1);    
x = dgnodes(:,1,:);x=x(:);
y = dgnodes(:,2,:);y=y(:);
t = cart2pol(x,y);
i = find(t<0);
t(i) = 2*pi+t(i);
[tj,jj] = sort(t);


figure(13);clf;
hold on;
plot(tj,ex(jj),'-k','LineWidth',1.5);
plot(tj,e1x(jj),'--b','LineWidth',1.5);
plot(tj,e2x(jj),'-.r','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); 
ylabel('x-component of electric field (MV/m)','FontSize',20);
legend({'$E_x$','$E^1_x$','$E^2_x$'},'interpreter','latex','location','Best','FontSize',18);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_Ex_onwire' num2str(windvel)];
print('-dpng',fn);

figure(14);clf;
hold on;
plot(tj,ey(jj),'-k','LineWidth',1.5);
plot(tj,e1y(jj),'--b','LineWidth',1.5);
plot(tj,e2y(jj),'-.r','LineWidth',1.5);
xlabel('\theta (rad)','FontSize',20); 
ylabel('y-component of electric field (MV/m)','FontSize',20);
legend({'$E_y$','$E^1_y$','$E^2_y$'},'interpreter','latex','location','Best','FontSize',18);
set(gca,'FontSize',18); box on;
axis tight;
set(gca, 'LooseInset', get(gca, 'TightInset'));
fn = ['groundedwire1_Ey_onwire' num2str(windvel)];
print('-dpng',fn);











