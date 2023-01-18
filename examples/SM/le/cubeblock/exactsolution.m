function uex = exactsolution(p)

mu=1;
kappa=1;
beta=1;

x = p(:,1,:);
y = p(:,2,:);
z = p(:,3,:);

u = x + 0.5*y.^3 + 0.5*sin(0.5*pi*y);
v = beta*y;
w = z; 

F11 = -ones(size(x));
F12 = -(pi*cos((pi*y)/2))/4 - (3*y.^2)/2;
F13 = -zeros(size(x));
F21 = -zeros(size(x));
F22 = -beta*ones(size(x));
F23 = -zeros(size(x));
F31 = -zeros(size(x));
F32 = -zeros(size(x));
F33 = -ones(size(x));

P11 = ((beta*(beta - 1))/kappa)*ones(size(x));
P12 = -mu*((pi*cos((pi*y)/2))/4 + (3*y.^2)/2);
P13 = 0;
P21 = -(mu*((pi*cos((pi*y)/2))/4 + (3*y.^2)/2))/beta + ((beta - 1)*((pi*cos((pi*y)/2))/4 + (3*y.^2)/2))/kappa;
P22 = (beta*mu - mu/beta + (beta - 1)/kappa)*ones(size(x));
P23 = 0;
P31 = 0;
P32 = 0;
P33 = 0;

uex(:,1,:) = u;
uex(:,2,:) = v;
uex(:,3,:) = w;

uex(:,4,:) = F11;
uex(:,5,:) = F21;
uex(:,6,:) = F31;
uex(:,7,:) = F12;
uex(:,8,:) = F22;
uex(:,9,:) = F32;
uex(:,10,:) = F13;
uex(:,11,:) = F23;
uex(:,12,:) = F33;

uex(:,13,:) = P11;
uex(:,14,:) = P21;
uex(:,15,:) = P31;
uex(:,16,:) = P12;
uex(:,17,:) = P22;
uex(:,18,:) = P32;
uex(:,19,:) = P13;
uex(:,20,:) = P23;
uex(:,21,:) = P33;

return;

syms x y z mu kappa

u = x + 0.5*y.^3 + 0.5*sin(0.5*pi*y);
v = y;
w = z;

F11 = diff(u,'x');
F12 = diff(u,'y');
F13 = diff(u,'z');
F21 = diff(v,'x');
F22 = diff(v,'y');
F23 = diff(v,'z');
F31 = diff(w,'x');
F32 = diff(w,'y');
F33 = diff(w,'z');

% F = [F11 F12 F13 F21 F22 F23 F31 F32 F33].';
% JF=det(F);
% CF=F.'*F;
% W  = (mu/2)*(trace(CF)+C-3) - mu*log(JF) + (kappa/2)*(JF-1)^2;
% P = jacobian(W,F);

F = [F11 F12 F13; F21 F22 F23; F31 F32 F33];
JF = det(F);
P11 = F11*mu - kappa*(F22*F33 - F23*F32)*(1-JF) - (mu*(F22*F33 - F23*F32))/JF;
P21 = F21*mu + kappa*(F12*F33 - F13*F32)*(1-JF) + (mu*(F12*F33 - F13*F32))/JF;
P31 = F31*mu - kappa*(F12*F23 - F13*F22)*(1-JF) - (mu*(F12*F23 - F13*F22))/JF;
P12 = F12*mu + kappa*(F21*F33 - F23*F31)*(1-JF) + (mu*(F21*F33 - F23*F31))/JF;
P22 = F22*mu - kappa*(F11*F33 - F13*F31)*(1-JF) - (mu*(F11*F33 - F13*F31))/JF;
P32 = F32*mu + kappa*(F11*F23 - F13*F21)*(1-JF) + (mu*(F11*F23 - F13*F21))/JF;
P13 = F13*mu - kappa*(F21*F32 - F22*F31)*(1-JF) - (mu*(F21*F32 - F22*F31))/JF;
P23 = F23*mu + kappa*(F11*F32 - F12*F31)*(1-JF) + (mu*(F11*F32 - F12*F31))/JF;
P33 = F33*mu - kappa*(F11*F22 - F12*F21)*(1-JF) - (mu*(F11*F22 - F12*F21))/JF;
 
F = [F11 F12 F13 F21 F22 F23 F31 F32 F33].';
P = [P11 P12 P13 P21 P22 P23 P31 P32 P33].';

b  = [diff(P11,'x') + diff(P12,'y') + diff(P13,'z'); ...
      diff(P21,'x') + diff(P22,'y') + diff(P23,'z');...
      diff(P31,'x') + diff(P32,'y') + diff(P33,'z')];

  



