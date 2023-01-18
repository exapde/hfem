function uex = exactsol(p, param, time)

x = p(:,1,:);
y = p(:,2,:);

u=(1/(sqrt(2)*pi))*sin(pi*x).*sin(pi*y).*sin(sqrt(2)*pi*time);    
v=sin(pi*x).*sin(pi*y).*cos(sqrt(2)*pi*time);
qx=(1/sqrt(2))*cos(pi*x).*sin(pi*y).*sin(sqrt(2)*pi*time);    
qy=(1/sqrt(2))*sin(pi*x).*cos(pi*y).*sin(sqrt(2)*pi*time);    

uex(:,1,:) = v;
uex(:,2,:) = qx;
uex(:,3,:) = qy;
uex(:,4,:) = u;

% x = p(:,1,:);
% y = p(:,2,:);
% c = param(1);
% alpha = param(2);
% phi = param(3);
% % 
% % U(:,4,:) = -tanh(alpha*(x*cos(phi) - c*t + y*sin(phi)) - 1); % u
% % U(:,1,:) = -alpha*c*(tanh(alpha*(x*cos(phi) - c*t + y*sin(phi)) - 1).^2 - 1);
% % U(:,2,:) = (alpha*cos(phi)*(tanh(alpha*(x*cos(phi) - c*t + y*sin(phi)) - 1).^2 - 1));
% % U(:,3,:) = (alpha*sin(phi)*(tanh(alpha*(x*cos(phi) - c*t + y*sin(phi)) - 1).^2 - 1));
% 
% U(:,4,:) = cos(x+y-sqrt(2)*t); % u
% U(:,1,:) = sqrt(2)*sin(x+y-sqrt(2)*t);
% U(:,2,:) = -sin(x+y-sqrt(2)*t);
% U(:,3,:) = -sin(x+y-sqrt(2)*t);
