Q = qcriterion(UDG);
gam = 1.4;
P = (gam-1)*(UDG(:,5,:) - 0.5*(UDG(:,2,:).*UDG(:,2,:)./UDG(:,1,:) + UDG(:,3,:).*UDG(:,3,:)./UDG(:,1,:) + UDG(:,4,:).*UDG(:,4,:)./UDG(:,1,:)));

[p0,t0,c0] = surfaceplot(mesh,mesh.elcon2,P,-7,1);
[p1,t1,c1,p2,t2,c2,p3,c3,t3] = isosurfaceplot(mesh, Q, -0.603, P, 1);

f = figure('Visible','off');
figure(1);clf;
set(axes,'DataAspectRatio',[1 1 1]);
hold on; 
patch('vertices',p0,'faces',t0,'cdata',c0(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');    
patch('vertices',p1,'faces',t1,'cdata',c1(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');    
patch('vertices',p2,'faces',t2,'cdata',c2(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');        
patch('vertices',p3,'faces',t3,'cdata',c3(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');                      
hold off;    
axis([-0.01 1.4 -0.05 0.8 0 0.1]);
cameratoolbar('SetCoordSys','y');
camlight('headlight');
lighting gouraud
camorbit(5,52,'data',[0 1 0]);
camzoom(1.2);
colormap jet;
axis off;
print -dpng les690c.png
close(f);


figure(1);clf;
set(axes,'DataAspectRatio',[1 1 1]);
hold on; 
patch('vertices',p0,'faces',t0,'cdata',c0(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');    
patch('vertices',p1,'faces',t1,'cdata',c1(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');    
patch('vertices',p2,'faces',t2,'cdata',c2(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');        
patch('vertices',p3,'faces',t3,'cdata',c3(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');                      
hold off;    
axis([0.3 1.4 0.3 0.7 0 0.1]);
colormap jet;
lighting gouraud
cameratoolbar('SetCoordSys','y');
campos([-2    4    3.5]);
camtarget([0.8500  0.5000  0.0500]);
camzoom(1.2);
camlight('right');
%camorbit(10,80,'data',[0 1 0]);
axis off;



figure(1);clf;
set(axes,'DataAspectRatio',[1 1 1]);
hold on; 
patch('vertices',p0,'faces',t0,'cdata',c0(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');    
patch('vertices',p1,'faces',t1,'cdata',c1(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');    
patch('vertices',p2,'faces',t2,'cdata',c2(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');        
patch('vertices',p3,'faces',t3,'cdata',c3(:), ...
   'facecol','interp','edgec','none','FaceLighting','gouraud');                      
hold off;    
axis([-0.05 1.2 -0.05 0.75 0 0.1]);
colormap jet;
lighting gouraud
cameratoolbar('SetCoordSys','y');
campos([-2 5 5]);
camtarget([0.6 0.25 0.05]);
camzoom(1.5);
camlight('right');
axis off;


for n = 100:10:150
    for i = 1:nproc
        %filename = ['/scratch/DIGASO/LES/Naca651810Re250kAoA54/sol3DLESIEDGP2offDesign_t' num2str(n) '_np' num2str(i-1)];
        filename = ['/scratch/DIGASO/LES/Naca651810Re250kAoA54/out3DLESoffDesignHexP2IEDG_t' num2str(n) '_np' num2str(i-1)];
        fileID = fopen([filename,'.bin'],'r');
        data = fread(fileID,'double');    
        fclose(fileID);

        nei = sum(dmd{i}.elempartpts(1:2));
        nfi = sum(dmd{i}.entpartpts(1:2)); 
        inde = dmd{i}.elempart(1:nei);               
        indf = dmd{i}.entpart(1:nfi);

        npf = size(mesh.perm,1);
        nc = app.nc;    
        nch = app.nch;    
        npe = size(mesh.plocal,1);        
        n1 = npe*nc*nei;        

        UDG(:,:,inde) = reshape(data(1:n1),[npe nc nei]);                                 
    end
    
    func = @(p) all(abs(p(:,3)-0.05)<1e-5);
    figure(1); clf; 
    faceplot(mesh,1,UDG(:,2,:)./UDG(:,1,:),func,[-0.1 1.1]);
    hold on;
    faceplot(mesh2,1,UDG(:,2,:)./UDG(:,1,:),func,[-0.1 1.1]);
    axis equal; axis tight; axis off; 
    colormap jet;  colorbar('location','EastOutSide','FontSize',12);    
    fn = ['uvel/fig'];     
    print('-dpng',sprintf('%s%05d.png',fn,n));
    
    figure(2); clf; 
    faceplot(mesh,1,UDG(:,4,:)./UDG(:,1,:),func,[-0.1 0.1]);
    hold on;
    faceplot(mesh2,1,UDG(:,4,:)./UDG(:,1,:),func,[-0.1 0.1]);
    axis equal; axis tight; axis off; 
    colormap jet;  colorbar('location','EastOutSide','FontSize',12);
    fn = ['wvel/fig'];     
    print('-dpng',sprintf('%s%05d.png',fn,n));
            
    Q = qcriterion(UDG);
    gam = 1.4;
    P = (gam-1)*(UDG(:,5,:) - 0.5*(UDG(:,2,:).*UDG(:,2,:)./UDG(:,1,:) + UDG(:,3,:).*UDG(:,3,:)./UDG(:,1,:) + UDG(:,4,:).*UDG(:,4,:)./UDG(:,1,:)));
    [p0,t0,c0] = surfaceplot(mesh,mesh.elcon2,P,-7,1);
    
    npl = size(Q,1);
    Q = reshape(Q,npl,[]);
    umin = min(Q,[],1);
    umax = max(Q,[],1);    
    isoa = linspace(-1,1,401); isob = isoa;
    for j = 1:length(isoa)
        inde = find((umin<isoa(j)) & (isoa(j) < umax));
        isob(j) = length(inde);
    end
    [~,imax] = max(isob);    
    isoa = linspace(isoa(imax-1),isoa(imax+1),401); isob = isoa;
    for j = 1:length(isoa)
        inde = find((umin<isoa(j)) & (isoa(j) < umax));
        isob(j) = length(inde);
    end
    [~,imax] = max(isob);
    isoval = isoa(imax);    
    
%     npl = size(Q,1);
%     Q = reshape(Q,npl,[]);
%     umin = min(Q,[],1);
%     umax = max(Q,[],1);    
%     isoa = linspace(isoval-0.01,isoval+0.01,401); isob = isoa;
%     for j = 1:length(isoa)
%         inde = find((umin<isoa(j)) & (isoa(j) < umax));
%         isob(j) = length(inde);
%     end
%     [~,imax] = max(isob);
%     isoval = isoa(imax);    
%  
    [n isoval]
    max(isob)
    
    [p1,t1,c1,p2,t2,c2,p3,c3,t3] = isosurfaceplot(mesh, Q, isoval, P, 1);    
    figure(3);clf;
    set(axes,'DataAspectRatio',[1 1 1]);
    hold on; 
    patch('vertices',p0,'faces',t0,'cdata',c0(:), ...
       'facecol','interp','edgec','none','FaceLighting','gouraud');    
    patch('vertices',p1,'faces',t1,'cdata',c1(:), ...
       'facecol','interp','edgec','none','FaceLighting','gouraud');    
    patch('vertices',p2,'faces',t2,'cdata',c2(:), ...
       'facecol','interp','edgec','none','FaceLighting','gouraud');        
    patch('vertices',p3,'faces',t3,'cdata',c3(:), ...
       'facecol','interp','edgec','none','FaceLighting','gouraud');                      
    hold off;    
    axis([-0.05 1.2 -0.05 0.75 0 0.1]);
    colormap jet;
    lighting gouraud
    cameratoolbar('SetCoordSys','y');
    campos([-2 5 5]);
    camtarget([0.6 0.25 0.05]);
    camzoom(1.5);
    camlight('right');
    axis off;
    fn = ['qcri/fig'];     
    print('-dpng',sprintf('%s%05d.png',fn,n));
        
end


N = [180 270 300 330 350 410 550 550 730 740 760 780 910 930 1010 1020 1060 1120 1130 1140]-150;
N = [530 550 650 930 1020]-150;
N = 700-150;
%N = [350 380 550 650 930 1020];
% N = 10:10:150;
% N = 20;
% %for n = 120:10:120
% N = 10:10:990;
for m = 1:length(N)
    n = N(m);
    for i = 1:nproc        
        filename = ['/scratch/DIGASO/LES/Naca651810Re250kAoA54/out3DLESoffDesignHexP2IEDG_t' num2str(n) '_np' num2str(i-1)];
        %filename = ['/scratch/DIGASO/LES/Naca651810Re250kAoA54/sol3DLESIEDGP2offDesign_t' num2str(n) '_np' num2str(i-1)];
        fileID = fopen([filename,'.bin'],'r');
        data = fread(fileID,'double');    
        fclose(fileID);

        nei = sum(dmd{i}.elempartpts(1:2));
        nfi = sum(dmd{i}.entpartpts(1:2)); 
        inde = dmd{i}.elempart(1:nei);               
        indf = dmd{i}.entpart(1:nfi);

        npf = size(mesh.perm,1);
        nc = app.nc;    
        nch = app.nch;    
        npe = size(mesh.plocal,1);        
        n1 = npe*nc*nei;        

        UDG(:,:,inde) = reshape(data(1:n1),[npe nc nei]);                                 
    end
    
%     func = @(p) all(abs(p(:,3)-0.05)<1e-5);
%     figure(1); clf; 
%     faceplot(mesh,1,UDG(:,2,:)./UDG(:,1,:),func,[-0.1 1.1]);
%     hold on;
%     faceplot(mesh2,1,UDG(:,2,:)./UDG(:,1,:),func,[-0.1 1.1]);
%     axis equal; axis tight; axis off; 
%     colormap jet;  colorbar('location','EastOutSide','FontSize',12);    
%     fn = ['uvel/fig'];     
%     print('-dpng',sprintf('%s%05d.png',fn,150+n));
%     
%     figure(2); clf; 
%     faceplot(mesh,1,UDG(:,4,:)./UDG(:,1,:),func,[-0.1 0.1]);
%     hold on;
%     faceplot(mesh2,1,UDG(:,4,:)./UDG(:,1,:),func,[-0.1 0.1]);
%     axis equal; axis tight; axis off; 
%     colormap jet;  colorbar('location','EastOutSide','FontSize',12);
%     fn = ['wvel/fig'];     
%     print('-dpng',sprintf('%s%05d.png',fn,150+n));
            
    Q = qcriterion(UDG);
    gam = 1.4;
    P = (gam-1)*(UDG(:,5,:) - 0.5*(UDG(:,2,:).*UDG(:,2,:)./UDG(:,1,:) + UDG(:,3,:).*UDG(:,3,:)./UDG(:,1,:) + UDG(:,4,:).*UDG(:,4,:)./UDG(:,1,:)));
    [p0,t0,c0] = surfaceplot(mesh,mesh.elcon2,P,-7,1);
    
    npl = size(Q,1);
    Q = reshape(Q,npl,[]);
    umin = min(Q,[],1);
    umax = max(Q,[],1);    
    isoa = linspace(-1,1,401); isob = isoa;
    for j = 1:length(isoa)
        inde = find((umin<isoa(j)) & (isoa(j) < umax));
        isob(j) = length(inde);
    end
    [~,imax] = max(isob);    
    isoa = linspace(isoa(imax-1),isoa(imax+1),401); isob = isoa;
    for j = 1:length(isoa)
        inde = find((umin<isoa(j)) & (isoa(j) < umax));
        isob(j) = length(inde);
    end
    [~,imax] = max(isob);
    isoval = isoa(imax);    
   
    [n isoval]
    max(isob)
    
    [p1,t1,c1,p2,t2,c2,p3,c3,t3] = isosurfaceplot(mesh, Q, isoval, P, 1);    
    figure(3);clf;
    set(axes,'DataAspectRatio',[1 1 1]);
    hold on; 
    patch('vertices',p0,'faces',t0,'cdata',c0(:), ...
       'facecol','interp','edgec','none','FaceLighting','gouraud');    
    patch('vertices',p1,'faces',t1,'cdata',c1(:), ...
       'facecol','interp','edgec','none','FaceLighting','gouraud');    
    patch('vertices',p2,'faces',t2,'cdata',c2(:), ...
       'facecol','interp','edgec','none','FaceLighting','gouraud');        
    patch('vertices',p3,'faces',t3,'cdata',c3(:), ...
       'facecol','interp','edgec','none','FaceLighting','gouraud');                      
    hold off;    
    axis([-0.05 1.2 -0.05 0.85 0 0.1]);
    set(gca,'clim',[97.8 99.05]);
    colormap jet;    
    lighting gouraud
    cameratoolbar('SetCoordSys','y');
    campos([-2 5 5]);
    camtarget([0.65 0.325 0.05]);
    camzoom(1.5);
    camlight('right');    
    axis off;
    fn = ['qcri2/fig'];     
    print('-dpng',sprintf('%s%05d.png',fn,150+n));
        
end

