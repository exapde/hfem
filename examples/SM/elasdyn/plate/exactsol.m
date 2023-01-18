function uex = exactsol(p,t)

x = p(:,1);
y = p(:,2);
z = p(:,3);

A = 0.4;

uex = zeros(size(x,1),12);

% Velocities
uex(:,1) = 0;
uex(:,2) = 0;
uex(:,3) = A*pi*cos(pi*t)*sin(pi*x).*sin(pi*y);
% Deformation Gradient
uex(:,4) = -1;
uex(:,5) = 0;
uex(:,6) = 0;
uex(:,7) = 0;
uex(:,8) = -1;
uex(:,9) = 0;
uex(:,10)= -A*pi*sin(pi*t)*cos(pi*x).*sin(pi*y);
uex(:,11)= -A*pi*sin(pi*t)*cos(pi*y).*sin(pi*x);
uex(:,12)= -1;

return;

syms A x y z t rho mu lambda;

nd = 3;

% displacement field
u1 = 0;
u2 = 0;
u3 = A*sin(pi*x)*sin(pi*y)*sin(pi*t);

% position field
p1 = u1+x;
p2 = u2+y;
p3 = u3+z;

% velocity field
v1 = diff(p1,'t');
v2 = diff(p2,'t');
v3 = diff(p3,'t');

% acceleration field
a1 = diff(v1,'t');
a2 = diff(v2,'t');
a3 = diff(v3,'t');

% deformation gradient
F11 = diff(p1,'x');
F12 = diff(p1,'y');
F13 = diff(p1,'z');
F21 = diff(p2,'x');
F22 = diff(p2,'y');
F23 = diff(p2,'z');
F31 = diff(p3,'x');
F32 = diff(p3,'y');
F33 = diff(p3,'z');
F = [F11 F12 F13; F21 F22 F23; F31 F32 F33];

 % identity matrix
I = eye(nd); 

% First Piola Kirchhoff stress :
% Linear Elasticity
e = 0.5*(F + F.' - 2*I); % small strain tensor
P = 2.*mu*e + lambda*trace(e)*I;  
% Saint Venant-Kirchhoff
E = 0.5*(F.'* F - I);
P = F * (2.*mu*E + lambda*trace(E)*I);

% body force
b1 = rho*a1 - diff(P(1,1),'x') - diff(P(1,2),'y') - diff(P(1,3),'z');
b2 = rho*a2 - diff(P(2,1),'x') - diff(P(2,2),'y') - diff(P(2,3),'z');
b3 = rho*a3 - diff(P(3,1),'x') - diff(P(3,2),'y') - diff(P(3,3),'z');



