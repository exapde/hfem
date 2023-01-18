function UDG = exactsolAlfven(p)

%   udg(:,1,:) :        rho
%   udg(:,2,:) :        rho*ux
%   udg(:,3,:) :        rho*uy
%   udg(:,4,:) :        rho*uz
%   udg(:,5,:) :        rho*e...E
%   udg(:,6,:) :        Bx
%   udg(:,7,:) :        By
%   udg(:,8,:) :        Bz
%   udg(:,9,:) :        phi

p1 = size(p,1);
p3 = size(p,3);
UDG = zeros(p1,9,p3);
x = p(:,1,:);
y = p(:,2,:);
theta = pi/6;
gam = 5/3;

xx = x*cos(theta) + y*sin(theta);
c = 0.1*sin(2*pi*xx);
UDG(:,1,:) = ones(size(UDG(:,1,:)));
UDG(:,2,:) = -c*sin(theta);
UDG(:,3,:) = c*cos(theta);
UDG(:,4,:) = 0.1*cos(2*pi*xx);
UDG(:,7,:) = UDG(:,1,:)*sin(theta)+c*cos(theta);
UDG(:,6,:) = UDG(:,1,:)*cos(theta)-c*sin(theta);
UDG(:,8,:) = UDG(:,4,:);
b = UDG(:,6,:).*UDG(:,6,:) + UDG(:,7,:).*UDG(:,7,:) + UDG(:,8,:).*UDG(:,8,:);
v = UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:) + UDG(:,4,:).*UDG(:,4,:);
UDG(:,5,:) = 0.1/(gam-1) + v*0.5*UDG(:,1,:) + b*0.5;
