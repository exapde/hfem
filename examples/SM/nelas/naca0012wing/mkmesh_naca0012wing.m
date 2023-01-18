function mesh = mkmesh_naca0012wing(porder,n2d,d,zz)

% Example: mesh = mkmesh_naca0012wing(4,20,0.005,linspace(0,4,20)); meshplot(mesh,1);

thick = 12;
th = (pi:-pi/200:pi/2)';
xt = (cos(th)+1)*1.0089304129;  
xt = xt(end:-1:1);
yt=naca(xt,thick);  
xb = flipud(xt);   
yb=-naca(xb,thick);
yt(1) = 0; yb(end) = 0;
xf =[xt; xb(2:end)];
yf =[yt; yb(2:end)];
xf(end) = xf(1);
yf(end) = yf(1);

chord = max(xf);
xf = xf/chord;
yf = yf/chord;
mesh2d = mkmesh_naca0012foil( xf, yf, porder, n2d, d);

bndexpr = {'all(p(:,3)<min(p0(:,3))+1e-6)','all(p(:,3)>max(p0(:,3))-1e-6)', ...
           'all(p(:,2)<1e-6)','all(p(:,2)>-1e-6)'};     
%bndexpr = {'all(p(:,3)<min(p0(:,3))+1e-6)','all(p(:,3)>max(p0(:,3))-1e-6)','true'};            
mesh = mkmesh_3dextrudemesh(mesh2d, zz, bndexpr);

pmin = min(mesh.p);
pmax = max(mesh.p);

% map p and dgnodes from the actual geometry to the unit cube
p = mesh.p;
dgnodes = mesh.dgnodes;
for i=1:3
    p(:,i) = (p(:,i)-pmin(i))/(pmax(i)-pmin(i));
    dgnodes(:,i,:) = (dgnodes(:,i,:)-pmin(i))/(pmax(i)-pmin(i));
end

% define the mapping
r = [0,1,0,1,1.2,1.5,1.2,1.5];
s = [0,0,1,1,0.25,0.25,0.75,0.75];
t = [0,0,0,0,1,1,1,1];
map = [r(:) s(:) t(:)];

% compute the mapping
p = mapp(p,map);
x = dgnodes(:,1,:); y = dgnodes(:,2,:); z = dgnodes(:,3,:);
p1 = mapp([x(:) y(:) z(:)],map);

% map p and dgnodes from the unit cube back to the actual geometry
for i=1:3
    p(:,i) = pmin(i) + p(:,i)*(pmax(i)-pmin(i));
    dgnodes(:,i,:) = reshape(p1(:,i), size(mesh.dgnodes(:,i,:)));
    dgnodes(:,i,:) = pmin(i) + dgnodes(:,i,:)*(pmax(i)-pmin(i));
end

mesh.p = p;
mesh.dgnodes = dgnodes;


% xzmap = [[0,1,1.2,1+0.5]' [0,0,1,1]'];
% 
% pxz = mesh.p(:,[1 3]);
% pxz(:,1) = (pxz(:,1)-pmin(1))/(pmax(1)-pmin(1));
% pxz(:,2) = (pxz(:,2)-pmin(3))/(pmax(3)-pmin(3));
% pxz = mapp(pxz,xzmap);
% pxz(:,1) = pmin(1) + (pmax(1)-pmin(1))*pxz(:,1);
% pxz(:,2) = pmin(3) + (pmax(3)-pmin(3))*pxz(:,2);
% mesh.p(:,1) = pxz(:,1);
% mesh.p(:,3) = pxz(:,2);
% 
% x = mesh.dgnodes(:,1,:);
% z = mesh.dgnodes(:,3,:);
% pxz = [x(:) z(:)];
% pxz(:,1) = (pxz(:,1)-pmin(1))/(pmax(1)-pmin(1));
% pxz(:,2) = (pxz(:,2)-pmin(3))/(pmax(3)-pmin(3));
% pxz = mapp(pxz,xzmap);
% pxz(:,1) = pmin(1) + (pmax(1)-pmin(1))*pxz(:,1);
% pxz(:,2) = pmin(3) + (pmax(3)-pmin(3))*pxz(:,2);
% mesh.dgnodes(:,1,:) = reshape(pxz(:,1), size(mesh.dgnodes(:,1,:)));
% mesh.dgnodes(:,3,:) = reshape(pxz(:,2), size(mesh.dgnodes(:,1,:)));








