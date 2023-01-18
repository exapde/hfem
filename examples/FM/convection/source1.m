function [sr,sr_udg] = source1(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

sr = zeros(ng,nch);

sr_udg = zeros(ng,nch,nc); 

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
