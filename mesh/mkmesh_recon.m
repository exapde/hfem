function mesh = mkmesh_recon(h,porder)

if nargin<1, h=1/8; end
if nargin<2, porder=1; end

%L = 1.5;
Lx = 25.6;
Ly = 12.8;
rect1 = [-Lx/2 -Ly/2;  Lx/2 -Ly/2; Lx/2 -2.0; -Lx/2 -2.0];
rect2 = [-Lx/2  -2.0;  Lx/2 -2.0;   Lx/2 2.0; -Lx/2  2.0];
rect3 = [-Lx/2   2.0;  Lx/2  2.0;  Lx/2 Ly/2; -Lx/2 Ly/2];

[p2,t2]=polymesh({rect2},[1],[1,0],[h/8,1.3]);
figure(2);clf;simpplot(p2,t2); axis on;

ind = find(p2(:,2)==-2.0);
pv1 = [rect1(1:3,:); p2(ind(end:-1:3),:); rect1(4,:)];
[p1,t1]=polymesh({pv1},[1],[1,0],[h,1.3]);
figure(1);clf;simpplot(p1,t1); axis on;

ind = find(p2(:,2)==2.0);
pv3 = [rect3(1,:); p2(ind(end:-1:3),:); rect3(2:4,:)];
[p3,t3]=polymesh({pv3},[1],[1,0],[h,1.3]);
figure(3);clf;simpplot(p3,t3); axis on;

[p12,t12] = connectmesh(p1,t1,p2,t2);
[p,t] = connectmesh(p12,t12,p3,t3);
figure(4);clf;simpplot(p,t); axis on;

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     

mesh = mkmesh(p,t,porder,bndexpr,0,1);
