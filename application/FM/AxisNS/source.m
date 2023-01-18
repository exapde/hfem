function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);

if nc==12
    nch = 4;
elseif nc==20
    nch = 5;
end

gam  = param{1};
gam1 = gam - 1.0;
Re   = param{3};
Re1  = 1/Re;

sr = zeros(ng,nch);
sr_udg = zeros(ng,nch,nc); 

xr   = p(:,2);
r    = udg(:,1);
ru   = udg(:,2);
rv   = udg(:,3);
rE   = udg(:,4);

rx   = udg(:,5);
rux  = udg(:,6);

ry   = udg(:,9);
rvy  = udg(:,11);

r1   = 1./r;
uv   = ru.*r1;
vv   = rv.*r1;
af   = 0.5*(uv.*uv+vv.*vv);
pr    = gam1*(rE-r.*af);

ux  = (rux - rx.*uv).*r1;
vy  = (rvy - ry.*vv).*r1;

% u_r  = -uv.*r1;
% u_ru =  r1;
v_r  = -vv.*r1;
v_rv =  r1;

ux_r  = (2*rx.*uv-rux).*(r1.^2);
ux_ru = -rx.*r1.^2;
ux_rx = -uv.*r1;
ux_rux = r1;
vy_r  = (2*ry.*vv-rvy).*(r1.^2);
vy_rv = -ry.*r1.^2;
vy_ry = -vv.*r1;
vy_rvy = r1;

sr(:,3) = pr-Re1*(4.0/3.0)*vv./xr+Re1*(2.0/3.0)*(ux+vy);

sr_udg(:,3,1) = gam1*af-Re1*(4.0/3.0)*v_r./xr+Re1*(2.0/3.0)*(ux_r+vy_r);
sr_udg(:,3,2) = -gam1*uv+Re1*(2.0/3.0)*(ux_ru);
sr_udg(:,3,3) = -gam1*vv-Re1*(4.0/3.0)*v_rv./xr+Re1*(2.0/3.0)*(vy_rv);
sr_udg(:,3,4) = gam1;

sr_udg(:,3,5) = Re1*(2.0/3.0)*(ux_rx);
sr_udg(:,3,6) = Re1*(2.0/3.0)*(ux_rux);

sr_udg(:,3,9) = Re1*(2.0/3.0)*(vy_ry);
sr_udg(:,3,11) = Re1*(2.0/3.0)*(vy_rvy);

% sr = 0*sr;
% sr_udg = 0*sr_udg;


% xr   = p(:,1);
% r    = udg(:,1);
% ru   = udg(:,2);
% rv   = udg(:,3);
% rE   = udg(:,4);
% 
% rx   = udg(:,5);
% rux  = udg(:,6);
% 
% ry   = udg(:,9);
% rvy  = udg(:,11);
% 
% r1   = 1./r;
% uv   = ru.*r1;
% vv   = rv.*r1;
% af   = 0.5*(uv.*uv+vv.*vv);
% pr    = gam1*(rE-r.*af);
% 
% ux  = (rux - rx.*uv).*r1;
% vy  = (rvy - ry.*vv).*r1;
% 
% u_r  = -uv.*r1;
% u_ru =  r1;
% % v_r  = -vv.*r1;
% % v_rv =  r1;
% 
% ux_r  = (2*rx.*uv-rux).*(r1.^2);
% ux_ru = -rx.*r1.^2;
% ux_rx = -uv.*r1;
% ux_rux = r1;
% vy_r  = (2*ry.*vv-rvy).*(r1.^2);
% vy_rv = -ry.*r1.^2;
% vy_ry = -vv.*r1;
% vy_rvy = r1;
% 
% sr(:,2) = pr-Re1*(4.0/3.0)*uv./xr+Re1*(2.0/3.0)*(ux+vy);
% 
% sr_udg(:,2,1) = gam1*af-Re1*(4.0/3.0)*u_r./xr+Re1*(2.0/3.0)*(ux_r+vy_r);
% sr_udg(:,2,2) = -gam1*uv-Re1*(4.0/3.0)*u_ru./xr+Re1*(2.0/3.0)*(ux_ru);
% sr_udg(:,2,3) = -gam1*vv+Re1*(2.0/3.0)*(vy_rv);
% sr_udg(:,2,4) = gam1;
% 
% sr_udg(:,2,5) = Re1*(2.0/3.0)*(ux_rx);
% sr_udg(:,2,6) = Re1*(2.0/3.0)*(ux_rux);
% 
% sr_udg(:,2,9) = Re1*(2.0/3.0)*(vy_ry);
% sr_udg(:,2,11) = Re1*(2.0/3.0)*(vy_rvy);
% 
% 
% 
