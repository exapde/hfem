
ISOVAL = 0;

% Generate mesh for unit cube
m=41;n=m;o=m;
[X,Y,Z]=ndgrid((0:m-1)/(m-1),(0:n-1)/(n-1),(0:o-1)/(o-1));

MASK = zeros(m,n,o);
MASK(:) = sqrt(sum([X(:)-0.5 Y(:)-0.5 Z(:)-0.5].^2,2)) - 0.3 - interp3([0:1/10:1],[0:1/10:1],[0:1/10:1],0.1*rand(11,11,11),X(:),Y(:),Z(:),'cubic');

DATA = zeros(m,n,o);
DATA(:) = interp3([0:1/10:1],[0:1/10:1],[0:1/10:1],rand(11,11,11),X(:),Y(:),Z(:),'cubic');

[F,V,C]=isosurface(X,Y,Z,MASK,ISOVAL,DATA);
figure(1);clf;
hold on; axis square; axis([0 1 0 1 0 1]); view(3); camlight, lighting gouraud
patch('Faces',F,'Vertices',V,'EdgeColor','none','FaceColor','interp','CData',C);

[F,V,C]=MarchingCubes(X,Y,Z,MASK,ISOVAL,DATA);
figure(2);clf;
hold on; axis square; axis([0 1 0 1 0 1]); view(3); camlight, lighting gouraud
patch('Faces',F,'Vertices',V,'EdgeColor','none','FaceColor','interp','CData',C);

p=[X(:),Y(:),Z(:)];
t = [1 2 m+2 m+1 m*n+1 m*n+2 m*n+m+2 m*n+m+1];
t=kron(t,ones(o-1,1))+kron(ones(size(t)),(0:o-2)'*(m*n));        
t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');            
xdg=reshape(p(t',1),[size(t,2) size(t,1)]);
ydg=reshape(p(t',2),[size(t,2) size(t,1)]);
zdg=reshape(p(t',3),[size(t,2) size(t,1)]);
MASK=MASK(:);DATA=DATA(:);
udg=reshape(MASK(t',1),[size(t,2) size(t,1)]);
cdg=reshape(DATA(t',1),[size(t,2) size(t,1)]);

umin = min(udg,[],1);
umax = max(udg,[],1);

% hexes contain the isosurface
inde = find((umin<ISOVAL) & (ISOVAL < umax));
ne   = length(inde);
udg = udg(:,inde);
cdg = cdg(:,inde);
xdg = xdg(:,inde);
ydg = ydg(:,inde);
zdg = zdg(:,inde);

% cube to tets
c2t = [1 2 4 6; 1 5 4 6; 5 8 4 6; 2 3 4 6; 7 3 4 6; 7 8 4 6]';

udg = reshape(udg(c2t,:),[4 6*ne]);
cdg = reshape(cdg(c2t,:),[4 6*ne]);
xdg = reshape(xdg(c2t,:),[4 6*ne]);
ydg = reshape(ydg(c2t,:),[4 6*ne]);
zdg = reshape(zdg(c2t,:),[4 6*ne]);

marchingtet(xdg, ydg, zdg, udg, ISOVAL, cdg);
figure(1),figure(2),figure(3);

% [X Y Z] = meshgrid([0:1/100:1],[0:1/100:1],[0:1/100:1]);
% % GENERATE RANDOM DATA
% DATA = zeros(101,101,101);
% DATA(:) = interp3([0:1/10:1],[0:1/10:1],[0:1/10:1],rand(11,11,11),X(:),Y(:),Z(:),'cubic');
% 
% % GENERATE A RANDOM MASK
% MASK = zeros(101,101,101);
% MASK(:) = sqrt(sum([X(:)-0.5 Y(:)-0.5 Z(:)-0.5].^2,2)) - 0.3 - interp3([0:1/10:1],[0:1/10:1],[0:1/10:1],0.1*rand(11,11,11),X(:),Y(:),Z(:),'cubic');
% %
% % ACTUAL PROBLEM
% % ==============
% % EXTRACT THE MASK SURFACE
% SURF = isosurface(X,Y,Z,MASK,0);
% % INTERPOLATE DATA ON MASK SURFACE
% DATA_SURF = interp3(X,Y,Z,DATA,SURF.vertices(:,1),SURF.vertices(:,2),SURF.vertices(:,3));
% % PLOT THE MASK SURFACE AND DATA
% 
% figure(1);clf;
% hold on; axis square; axis([0 1 0 1 0 1]); view(3); camlight
% patch('Faces',SURF.faces,'Vertices',SURF.vertices,'EdgeColor','none','FaceColor','interp','FaceVertexCData',DATA_SURF);
% 
% [F,V,C]=MarchingCubes(X,Y,Z,MASK,0,DATA);
% figure(2);clf;
% hold on; axis square; axis([0 1 0 1 0 1]); view(3); camlight
% patch('Faces',F,'Vertices',V,'EdgeColor','none','FaceColor','interp','CData',C);
% figure(1),figure(2)

% let us suppose you have the following matrices:
% 
% x(nx,ny,nz): x coordinates of your grid
% 
% y(nx,ny,nz): y coordinates of your grid
% 
% z(nx,ny,nz): z coordinates of your grid
% 
% dudx(nx,ny,nz): x-wise derivative of x-wise velocity component
% 
% dudy(nx,ny,nz): y-wise derivative of x-wise velocity component
% .
% .
% .
% dwdz(nx,ny,nz): z-wise derivative of z-wise velocity component
% 
% This is, more or less, how i would put it down in Matlab the task of plotting Q isosurfaces (tipically used for flow visualization in DNS/LES):
% 
% %MATLAB CODE
% 
% iso_q=100; %Pick your number here
% 
% %Definition of Q
% q=-0.5*(dudx.^2+dvdy.^2+dwdz.^2)-dudy.*dvdx-dudz.*dwdx-dvdz.*dwdy;
% 
% %Plotting a Q isosurface, Q=iso_q
% figure()
% p=patch(isosurface(x,y,z,q,iso_q));
% set(p,'FaceColor','red','EdgeColor','none');
% daspect([1,1,1])
% axis tight
% ax = -1; ay = 1; az = 1;
% view([ax,ay,az]);
% camroll(240)
% camlight
% lighting gouraud