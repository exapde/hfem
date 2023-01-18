
% plot matlab solution
figure(1); clf;
scaplot(mesh,eulereval(UDG(:,1:4,:),'M',1.4),[],3,0); 
hold on;
axis off; axis equal; axis tight;    
colormap('jet');
colorbar('FontSize',14);
axis([-0.2 2 -0.8 0.8]);

% get c++ parallel solution and plot it
%function [UDG,UH] = getsolfrombinaryfile(filename,nproc,npv,nc,npf,nch,hybrid)
[UDGpar,UHpar] = getsolfrombinaryfile(apppar.fileout,apppar.nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);

figure(2); clf;
scaplot(mesh,eulereval(UDGpar(:,1:4,:),'M',1.4),[],3,0); 
hold on;
axis off; axis equal; axis tight;    
colormap('jet');
colorbar('FontSize',14);
axis([-0.2 2 -0.8 0.8]);

% get c++ serial solution and plot it
[UDGser,UHser] = getsolfrombinaryfile(appser.fileout,appser.nproc,master.npv,app.nc,master.npf,app.nch,app.hybrid);

figure(3); clf;
scaplot(mesh,eulereval(UDGser(:,1:4,:),'M',1.4),[],3,0); 
hold on;
axis off; axis equal; axis tight;    
colormap('jet');
colorbar('FontSize',14);
axis([-0.2 2 -0.8 0.8]);


