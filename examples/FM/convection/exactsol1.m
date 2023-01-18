function uex = exactsol1(p)

x = p(:,1,:);
y = p(:,2,:);
uex = sin(pi*(x-y));
ind = x<=y;
uex(ind) = 0;

