function uex = exactsol2(p)

x = p(:,1,:);
y = p(:,2,:);
uex = 0*x;
ind = find(x<=y+0.2);
uex(ind) = 1;

