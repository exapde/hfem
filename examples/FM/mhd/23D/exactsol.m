function UDG = exactsol(p,t)

p1 = size(p,1);
p3 = size(p,3);
UDG = zeros(p1,9,p3);
x = p(:,1,:);
y = p(:,2,:);
 
UDG(:,1,:) = 2 + sin(x + y - 2*t); 
UDG(:,2,:) = ones(size(UDG(:,2,:))); 
UDG(:,3,:) = ones(size(UDG(:,2,:)));
UDG(:,5,:) = 5*ones(size(UDG(:,2,:)));

