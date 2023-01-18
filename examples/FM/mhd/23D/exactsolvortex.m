function UDG = exactsolvortex(p,t)

%   UDG(:,1,:) :        rho
%   UDG(:,2,:) :        rho*ux
%   UDG(:,3,:) :        rho*uy
%   UDG(:,4,:) :        rho*uz
%   UDG(:,5,:) :        rho*e
%   UDG(:,6,:) :        Bx
%   UDG(:,7,:) :        By
%   UDG(:,8,:) :        Bz
%   UDG(:,9,:) :        phi

p1 = size(p,1);
p3 = size(p,3);
UDG = zeros(p1,9,p3);
x = p(:,1,:)+t*ones(size(p(:,1,:)));
y = p(:,2,:)+t*ones(size(p(:,2,:)));
gam = 5/3;

r = sqrt(x.*x + y.*y);
dx = -exp(0.5*(1-r.*r))*0.5.*y/pi;
dy = exp(0.5*(1-r.*r))*0.5.*x/pi;
dp = (1-r.*r).*exp(1-r.*r)/(8*pi^2) - exp(1-r.*r)/(8*pi^2);
UDG(:,1,:) = ones(size(UDG(:,1,:)));           
UDG(:,2,:) = UDG(:,1,:) + dx;                
UDG(:,3,:) = UDG(:,1,:) + dy;
UDG(:,6,:) = dx; 
UDG(:,7,:) = dy;
p = UDG(:,1,:) + dp;
b = UDG(:,6,:).*UDG(:,6,:) + UDG(:,7,:).*UDG(:,7,:) + UDG(:,8,:).*UDG(:,8,:);
v = UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:) + UDG(:,4,:).*UDG(:,4,:);
UDG(:,5,:) = p/(gam-1) + v*0.5 + b*0.5;
