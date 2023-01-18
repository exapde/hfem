function uex = exactsol1(p, param, time)

x = p(:,1,:);
y = p(:,2,:);

u=(1/(sqrt(2)*pi))*sin(pi*x).*sin(pi*y).*sin(sqrt(2)*pi*time);    
v=sin(pi*x).*sin(pi*y).*cos(sqrt(2)*pi*time);

uex(:,1,:) = u;
uex(:,2,:) = v;


