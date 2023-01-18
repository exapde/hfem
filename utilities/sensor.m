function [f] = sensor(UDG)

gam  = 1.4;
gam1 = gam-1;

% r   = UDG(:,1,:);
% ru  = UDG(:,2,:);
% rv  = UDG(:,3,:);
% rE  = UDG(:,4,:);
% 
% rx  = UDG(:,5,:);
% rux = UDG(:,6,:);
% 
% ry  = UDG(:,9,:);    
% rvy = UDG(:,11,:);    
% 
% r1  = 1./r;
% u   = ru.*r1;
% v   = rv.*r1;
% 
% ux  = (rux - rx.*u).*r1;
% vy  = (rvy - ry.*v).*r1;
% % uy  = (ruy - ry.*u).*r1;
% % vx  = (rvx - rx.*v).*r1;
% 
% af  = (u.*u+v.*v);
% p   = gam1*(rE - 0.5*r.*af);
% c2  = gam* p.*r1;
% c   = sqrt(c2);
% 
% f   = -(ux+vy)./(c); % divergence sensor
% %g   = (uy-vx)./(sqrt(af)+c); % vorticity sensor


r    = UDG(:,1,:);
ru   = UDG(:,2,:);
rv   = UDG(:,3,:);
rE   = UDG(:,4,:);

rx   = UDG(:,5,:);
rux  = UDG(:,6,:);

ry   = UDG(:,9,:);    
rvy  = UDG(:,11,:);    

r1   = 1./r;
u    = ru.*r1;
v    = rv.*r1;

ux  = (rux - rx.*u).*r1;
vy  = (rvy - ry.*v).*r1;

af    = 0.5*(u.*u+v.*v);
Duv   = ux+vy; 
p     = gam1*(rE -r.*af);
c2    = gam* p.*r1;

f = -Duv./sqrt(c2);
