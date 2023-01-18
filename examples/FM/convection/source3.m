function [sr,sr_udg] = source3(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

x = p(:,1);
y = p(:,2);
sr = -udg(:,1)+sin(pi*x).*sin(pi*y)+pi*cos(pi*x).*sin(pi*y)+pi*cos(pi*y).*sin(pi*x);
sr_udg = zeros(ng,nch,nc); 
sr_udg(:,:,1) = -1;


% if nc==3
%     xg = p(:,1);
%     yg = p(:,2);
%     %sr = 2*pi^2*sin(pi*xg).*sin(pi*yg);
% elseif nc==4
%     xg = p(:,1);
%     yg = p(:,2);
%     zg = p(:,3);
%     sr = 3*pi^2*sin(pi*xg).*sin(pi*yg).*sin(pi*zg);
% end
