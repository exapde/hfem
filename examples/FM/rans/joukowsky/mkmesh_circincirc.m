function mesh = mkmesh_circincirc(porder,m,n,r1,r2,scale1,scale2)
% Exampe: mesh = mkmesh_circincirc(2,21,21,1,10);

nodetype = 1;
elemtype = 1;
[p,t] = squaremesh(m,n,1,elemtype);
ind = p(:,1)<=0.5;
p(ind,1) = loginc(p(ind,1),scale1);
ind = p(:,1)>=0.5;
p(ind,1) = logdec(p(ind,1),scale1);
p(:,2) = loginc(p(:,2),scale2);
figure(5);clf;simpplot(p,t); axis on;

bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<1e-3)'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

p = mesh.p;
pnew(:,1) = (r1+(r2-r1)*p(:,2)).*sin(2*pi*p(:,1)+pi/2);
pnew(:,2) = (r1+(r2-r1)*p(:,2)).*cos(2*pi*p(:,1)+pi/2);
figure(6);clf;simpplot(pnew,t); axis on;
mesh.p = pnew;
[mesh.p,mesh.t]=fixmesh(mesh.p,mesh.t);
%figure(6);clf;simpplot(mesh.p,mesh.t); axis on;

p = mesh.dgnodes;   
pnew = zeros(size(p));
pnew(:,1,:) = (r1+(r2-r1)*p(:,2,:)).*sin(2*pi*p(:,1,:)+pi/2);
pnew(:,2,:) = (r1+(r2-r1)*p(:,2,:)).*cos(2*pi*p(:,1,:)+pi/2);
mesh.dgnodes = pnew;
mesh.fcurved = true(size(mesh.f,1),1);
mesh.tcurved = true(size(mesh.t,1),1);

bndexpr = {'all(sqrt(sum(p.^2,2))<1.15)','true'};     
[mesh.f,mesh.t2f] = mkt2f(mesh.t,elemtype);
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);
mesh.nf = size(mesh.f,1);
