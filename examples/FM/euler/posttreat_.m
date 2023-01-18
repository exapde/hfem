

timeStepNo = 1402;
nproc = 112;
filename = 'eulerout';
[UDG,UH] = getsolfrombinaryfile(filename,nproc,timeStepNo,master.npv,app.nc,master.npf,app.nch,app.hybrid);

%f = figure('visible','off');
scaplot(mesh,UDG(:,1,:),[0.9,2.1],2,0); 
axis on; axis tight; colormap jet;
%print -dpng KH120.png 
%close(f)
%title('rho');
drawnow