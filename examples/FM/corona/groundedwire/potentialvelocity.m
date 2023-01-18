function [u,v] = potentialvelocity(x,y,uinf,a)

r = sqrt(x.^2 + y.^2);
%phi = (r+a^2./r).*(x./r);
u = (a^2*r.^2 + r.^4 - 2*a^2*x.^2)./(r.^4);
v = -(2*a^2*x.*y)./(r.^4);
u = uinf*u; v = uinf*v;

