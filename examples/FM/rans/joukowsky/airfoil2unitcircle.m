function [X,Y,dXdx,dYdx,dXdy,dYdy] = airfoil2unitcircle(x,y,a,cx,cy,scale,alpha)

% [X,Y,dXdx,dYdx,dXdy,dYdy] = inverseJoukowsky(x,y,a,cx,cy,scale);
% X = X-cx;
% Y = Y-cy;
% X = X/a;
% Y = Y/a;
% dXdx = dXdx/a;
% dYdx = dYdx/a;
% dXdy = dXdy/a;
% dYdy = dYdy/a;

if nargin<7
    alpha = 0;
end

w = x + y*1i;
w = exp(1i*alpha)*w;
c = cx + cy*1i;

s = w/(2*scale) - sqrt((w/(2*scale)).^2 - 1);
Z = (s-c)/a;
X = real(Z);
Y = imag(Z);
dZdw = (1/(2*scale) - w./(4*scale^2*(w.^2/(4*scale^2) - 1).^(1/2)))/a;
dZdw = exp(1i*alpha)*dZdw;
dXdx = real(dZdw);
dYdx = imag(dZdw);
dXdy = -dYdx;
dYdy = dXdx;

ind = find((X).^2+(Y).^2<1-1e-10);
w = w(ind);
s = w/(2*scale) + sqrt((w/(2*scale)).^2 - 1);
Z = (s-c)/a;
X(ind) = real(Z);
Y(ind) = imag(Z);
dZdw = (1/(2*scale) + w./(4*scale^2*(w.^2/(4*scale^2) - 1).^(1/2)))/a;
dZdw = exp(1i*alpha)*dZdw;
dXdx(ind) = real(dZdw);
dYdx(ind) = imag(dZdw);
dXdy(ind) = -dYdx(ind);
dYdy(ind) = dXdx(ind);

return;

syms a c scale w s Z 
s = w/(2*scale) - sqrt((w/(2*scale)).^2 - 1);
Z = (s-c)/a;



