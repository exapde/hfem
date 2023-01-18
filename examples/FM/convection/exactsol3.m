function uex = exactsol3(p)

x = p(:,1,:);
y = p(:,2,:);
uex = sin(pi*x).*sin(pi*y);
uex(:,2,:) = pi*sin(pi*(x + y));
uex(:,3,:) = pi*sin(pi*(x + y));


