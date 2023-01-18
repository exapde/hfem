function u = wind(mesh,Uinf)

x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
X = x+200+1;
Y = y+15;
nu = 1.6e-5;

Re = Uinf/nu;
ReX = Re*X;
deltaX = 0.37*X./(ReX.^(1/5));
u = Uinf*(Y./deltaX).^(1/7);
ind = Y>deltaX;
u(ind) = Uinf;

b  = 100;
dw = meshdist(mesh,5);
uw = (1-exp(-b*dw));
u  = Uinf.*uw;

figure(1); clf; 
scaplot(mesh,u,[],2); 
axis equal; axis tight; axis on; colormap jet;
colorbar('FontSize',15);
set(gca,'FontSize',18);
xlabel('x (m)','FontSize',20);
ylabel('y (m)','FontSize',20);
box off;
axis([-0.012 0.012 -0.01 0.01]*10)

% ind = find(abs(x+15)<=0);
% yin = y(ind);
% uin = u(ind);
% [yin,ii] = sort(yin);
% uin = uin(ii);
% 
% figure(2); clf;
% plot(yin(:),uin(:),'o-','LineWidth',1);
% 

% [dgnodes,w] = potentialfield(master,mesh,u,-2,[1]);
% figure(2); clf;
% hold on;
% for k = 1:size(w,3)    
%     plot(dgnodes(:,2,k),w(:,1,k),'-k','LineWidth',1.5);
% end
% hold off;
% axis normal; axis tight; axis on; box on;
% xlabel('y (m)','FontSize',20);
% ylabel('Velocity (m/s)','FontSize',20);
% %set(gca,'xtick',[-5:1:5]);
% %set(gca,'ytick',[0:0.2:1]);
% %axis([-5 5 0 2e5])
% set(gca,'FontSize',18);









