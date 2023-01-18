function [vx,vy,phi,psi] = potentialflow(x,y)

alpha = 0;
scale = 1/4.0552;
Gamma = 1;
Lambda= 0.1;
delta = 0.1;
k1=1.05;
k2=20;

a  = k1*Gamma*sqrt((1+Lambda)^2 + delta^2); % Inner radius 
b  = k2*Gamma*sqrt((1+Lambda)^2 + delta^2); % Outer radius
c  = Gamma*(-Lambda + 1i*delta);  % origin
cx = real(c); cy = imag(c); 

[x,y,dXdx,dYdx,dXdy,dYdy] = airfoil2unitcircle(x,y,a,cx,cy,scale,alpha);

s = cos(-alpha) + 1i*sin(-alpha);
z = x + 1i*y;
g = s*(z + 1./z);
psi = imag(g);
%phi = real(g);
r2 = x.^2+y.^2;
p = x + x./r2; % potential
q = y - y./r2; % streamline
px = 1 + 1./r2 - (2*x.*x)./(r2.^2);
py = -(2*x.*y)./(r2.^2) ;
qx = (2*x.*y)./(r2.^2) ;
qy = 1 - 1./r2 + (2*y.*y)./(r2.^2);

% u1 = -Ea*x + (Ea*a^2)*x./r2;
% Ex = Ea - (Ea*a^2)./r2 + (2*Ea*a^2*x.*x)./(r2.^2);
% Ey = (2*Ea*a^2*x.*y)./(r2.^2);  

phi = real(s)*p - imag(s)*q;
phix = real(s)*px - imag(s)*qx;
phiy = real(s)*py - imag(s)*qy;

%vx = phix; vy = phiy;
vx = phix.*dXdx + phix.*dYdx;
vy = phiy.*dXdy + phiy.*dYdy;


% function [X,Y,dXdx,dYdx,dXdy,dYdy] = airfoil2unitcircle(x,y,alpha,a,cx,cy,scale)
% 
% w = x + y*1i;
% w = exp(1i*alpha)*w;
% c = cx + cy*1i;
% 
% s = w/(2*scale) - sqrt((w/(2*scale)).^2 - 1);
% Z = (s-c)/a;
% X = real(Z);
% Y = imag(Z);
% dZdw = (1/(2*scale) - w./(4*scale^2*(w.^2/(4*scale^2) - 1).^(1/2)))/a;
% dXdx = real(dZdw);
% dYdx = imag(dZdw);
% dXdy = -dYdx;
% dYdy = dXdx;
% 
% ind = find((X).^2+(Y).^2<1-1e-10);
% w = w(ind);
% s = w/(2*scale) + sqrt((w/(2*scale)).^2 - 1);
% Z = (s-c)/a;
% X(ind) = real(Z);
% Y(ind) = imag(Z);
% dZdw = (1/(2*scale) + w./(4*scale^2*(w.^2/(4*scale^2) - 1).^(1/2)))/a;
% dXdx(ind) = real(dZdw);
% dYdx(ind) = imag(dZdw);
% dXdy(ind) = -dYdx(ind);
% dYdy(ind) = dXdx(ind);
% 
return;


scale = 1e3/7;
Gamma = 1;
Lambda= 0.1;
delta = 0.1;
k1=1.05;
k2=10;

a  = k1*Gamma*sqrt((1+Lambda)^2 + delta^2); % Inner radius 
b  = k2*Gamma*sqrt((1+Lambda)^2 + delta^2); % Outer radius
c  = Gamma*(-Lambda + 1i*delta);  % origin

% unit circle mesh
mesh0 = mkmesh_circincirc(2,101,101,a/a,b/a,0.25,0.25);

alpha = 0;
mesh2 = mesh0;
[mesh2.p(:,1), mesh2.p(:,2)] = unitcircle2airfoil(mesh0.p(:,1),mesh0.p(:,2),a,real(c),imag(c),scale,-alpha);
[mesh2.dgnodes(:,1,:),mesh2.dgnodes(:,2,:)] = unitcircle2airfoil(mesh0.dgnodes(:,1,:),mesh0.dgnodes(:,2,:),a,real(c),imag(c),scale,-alpha);

mesh3 = mesh0;
[mesh3.p(:,1), mesh3.p(:,2)] = airfoil2unitcircle(mesh2.p(:,1),mesh2.p(:,2),a,real(c),imag(c),scale,alpha);
[mesh3.dgnodes(:,1,:),mesh3.dgnodes(:,2,:)] = airfoil2unitcircle(mesh2.dgnodes(:,1,:),mesh2.dgnodes(:,2,:),a,real(c),imag(c),scale,alpha);
e = mesh3.dgnodes-mesh0.dgnodes; max(abs(e(:)))

cx = real(c);
cy = imag(c);
[vx,vy,phi,psi] = potentialflow(mesh2.dgnodes(:,1,:),mesh2.dgnodes(:,2,:),alpha,a,cx,cy,scale);

%[ex,ey] = groundedjoukowskiuniformfield(Ri,Einf(end),mesh2.dgnodes(:,1,:),mesh2.dgnodes(:,2,:),Cx,Cy,Ro,Vapp,Scale); 

figure(1); clf; scaplot(mesh2,160*vx,[],1); axis on; colormap jet; 
figure(2); clf; scaplot(mesh2,160*vy,[],1); axis on; colormap jet; 
figure(3); clf; scaplot(mesh2,phi,[],1); axis on; colormap jet; 
figure(4); clf; scaplot(mesh2,psi,[],1); axis on; colormap jet; 
figure(5); clf; scaplot(mesh2,psi,[],1); 
hold on;
quiver(mesh2.dgnodes(:,1,:),mesh2.dgnodes(:,2,:),vx,vy);
axis on; colormap jet; 

figure(1); clf; scaplot(mesh2,ex,[0 0.4],1); axis off; colormap jet; set(gca,'FontSize',18);
figure(2); clf; scaplot(mesh2,ey,[],1); axis off; colormap jet; set(gca,'FontSize',18);

contour(squeeze(mesh2.dgnodes(:,1,:)),squeeze(mesh2.dgnodes(:,2,:)),squeeze(psi),100);
