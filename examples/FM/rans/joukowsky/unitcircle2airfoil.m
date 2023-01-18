function [x,y] = unitcircle2airfoil(X,Y,a,cx,cy,scale,alpha)

% X = X*a+cx;
% Y = Y*a+cy;
% [x,y] = Joukowsky(X,Y,scale);

if nargin<7
    alpha = 0;
end
    
Z = X + Y*1i;            % unit circle
s = a*Z + (cx + cy*1i);  % circle with radius a and center (cx, cy)
w = scale*(s + 1./s);     % scaled Joukowski airfoil
w = exp(1i*alpha)*w;
x = real(w);
y = imag(w);

