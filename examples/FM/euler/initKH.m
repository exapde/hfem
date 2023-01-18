function UDG = initKH(p,param,time)

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
UDG = zeros(p1,4,p3);
y = p(:,2,:);
x = p(:,1,:);
gam = 1.4;
p = 2.5;
lambda = 1/6;
A = 0.025;
for i=1:p3
   for j=1:p1
      if (abs(y(j,i))>0.25)
         UDG(j,1,i) = 1;
         UDG(j,2,i) = -0.5;
         UDG(j,3,i) = 0.0;
      else
         UDG(j,1,i) = 2;
         UDG(j,2,i) = 1;
         UDG(j,3,i) = 0.0;
      end
   end
end
for i=1:p3
   for j=1:p1
      if (abs(y(j,i)-0.25)<0.025)
         UDG(j,3,i) = UDG(j,3,i)+UDG(j,1,i)*(A*sin(-2*pi*(x(j,i)+0.5)/lambda));
      elseif (abs(y(j,i)+0.25)<0.025)
         UDG(j,3,i) = UDG(j,3,i)+UDG(j,1,i)*(A*sin(2*pi*(x(j,i)+0.5)/lambda));
      end
   end
end

v = UDG(:,2,:).*UDG(:,2,:) + UDG(:,3,:).*UDG(:,3,:);
UDG(:,4,:) = p/(gam-1)*ones(size(UDG(:,1,:))) + v*0.5.*UDG(:,1,:);