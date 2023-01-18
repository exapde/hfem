function uex = exactsolution(p)

mu=1;
kappa=1;
beta=1;

x = p(:,1,:);
y = p(:,2,:);

u = x + 0.5*y.^3 + 0.5*sin(0.5*pi*y);
v = beta*y;

F11 = -ones(size(x));
F12 = -(pi*cos((pi*y)/2))/4 - (3*y.^2)/2;
F21 = -zeros(size(x));
F22 = -beta*ones(size(x));

P11 = ((beta*(beta - 1))/kappa)*ones(size(x));
P12 = -mu*((pi*cos((pi*y)/2))/4 + (3*y.^2)/2);
P21 = -(mu*((pi*cos((pi*y)/2))/4 + (3*y.^2)/2))/beta + ((beta - 1)*((pi*cos((pi*y)/2))/4 + (3*y.^2)/2))/kappa;
P22 = (beta*mu - mu/beta + (beta - 1)/kappa)*ones(size(x));
  
uex(:,1,:) = u;
uex(:,2,:) = v;
uex(:,3,:) = F11;
uex(:,4,:) = F21;
uex(:,5,:) = F12;
uex(:,6,:) = F22;
uex(:,7,:) = P11;
uex(:,8,:) = P21;
uex(:,9,:) = P12;
uex(:,10,:) = P22;

% syms x y mu kappa beta
% 
% u = x + 0.5*y.^3 + 0.5*sin(0.5*pi*y);
% %u = x + 0.5*y.^2;
% v = beta*y;
% 
% F11 = diff(u,'x');
% F12 = diff(u,'y');
% 
% F21 = diff(v,'x');
% F22 = diff(v,'y');
% 
% F = [F11 F12; F21 F22];
% 
% detF = (F11*F22 - F12*F21);
% 
% Pi = [F11.*mu - (F22.*mu)./(detF), F12.*mu + (F21.*mu)./(detF);...
%       F21.*mu + (F12.*mu)./(detF), F22.*mu - (F11.*mu)./(detF)];
%   
% bi = [diff(Pi(1,1),'x') + diff(Pi(1,2),'y'); diff(Pi(2,1),'x') + diff(Pi(2,2),'y')];
% 
% pres  = (detF-1)/kappa;
% % px = diff(p,'x');
% % px = diff(p,'x');
% 
% Pv = [F22.*pres, -F21.*pres; - F12.*pres, F11.*pres]; 
% 
% bv = [diff(Pv(1,1),'x') + diff(Pv(1,2),'y'); diff(Pv(2,1),'x') + diff(Pv(2,2),'y')];
% 
% P  = Pi + Pv;
% 
% b  = [diff(P(1,1),'x') + diff(P(1,2),'y'); diff(P(2,1),'x') + diff(P(2,2),'y')];
% 
% 
% 
% 
