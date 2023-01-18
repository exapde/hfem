function [mesh,mesh0] = mkmesh_joukowsky(porder,nr,nt,scale1,scale2)

scale = 1/4.0552;
Gamma = 1;
Lambda= 0.1;
delta = 0.1;
k1=1.05;
k2=20;

a  = k1*Gamma*sqrt((1+Lambda)^2 + delta^2); % Inner radius 
b  = k2*Gamma*sqrt((1+Lambda)^2 + delta^2); % Outer radius
c  = Gamma*(-Lambda + 1i*delta);  % origin

% unit circle mesh
mesh0 = mkmesh_circincirc(porder,nt,nr,a/a,b/a,scale1,scale2);

mesh = mesh0;
[mesh.p(:,1), mesh.p(:,2)] = unitcircle2airfoil(mesh0.p(:,1),mesh0.p(:,2),a,real(c),imag(c),scale);
[mesh.dgnodes(:,1,:),mesh.dgnodes(:,2,:)] = unitcircle2airfoil(mesh0.dgnodes(:,1,:),mesh0.dgnodes(:,2,:),a,real(c),imag(c),scale);
 
