function [alpha1,eiga] = glm(mesh,app,UDG,UH,ngrid)

% r = UDG(:,1,:);
% vx = UDG(:,2,:)./UDG(:,1,:);
% vy = UDG(:,2,:)./UDG(:,1,:);
% re = UDG(:,4,:);
% bx = UDG(:,5,:);
% by = UDG(:,6,:);

r = UH(1,:);
vx = UH(2,:)./UH(1,:);
vy = UH(3,:)./UH(1,:);
vz = UH(4,:)./UH(1,:);
re = UH(5,:);
bx = UH(6,:);
by = UH(7,:);
bz = UH(7,:);
b = bx.*bx + by.*by + bz.*bz;
gam = app.arg{1};
q = 0.5*(vx.*vx + vy.*vy + vz.*vz);
p = (gam - 1)*(re - r.*q - b*0.5);

a    = sqrt(gam*p./r); % sound speed
ca   = b./sqrt(r);   % Alfven speed
cax  = bx./sqrt(r);
cay  = by./sqrt(r);
caz  = bz./sqrt(r);

% cf = g(B,v); 
cfx = a.*a +  b./r + sqrt((a.*a +  b./r).^2 - 4*a.*a.*bx.*bx./r);
cfy = a.*a +  b./r + sqrt((a.*a +  b./r).^2 - 4*a.*a.*by.*by./r);
cfz = a.*a +  b./r + sqrt((a.*a +  b./r).^2 - 4*a.*a.*bz.*bz./r);
cfx = sqrt(cfx*0.5);
cfy = sqrt(cfy*0.5);
cfz = sqrt(cfz*0.5);
vxcfx = abs(vx) + cfx;
vycfy = abs(vy) + cfy;
vzcfz = abs(vz) + cfz;

alpha1 = max(max(max([vxcfx,vycfy,vzcfz])));
% eigax = max(max(max(vxcfx+cax)));
% eigay = max(max(max(vycfy+cay)));
% eiga = eigax + eigay;
%eiga = max(max(max([vxcfx+cax,vycfy+cay])));
eiga = max(max(max(vxcfx + cax + vycfy + cay + vzcfz + caz)));

%%%% Dedner & al.
%alpha2 = sqrt(0.18 * alpha1); 

%%%% Mignone & al.
% coefdiffadv = 1; % [0,1]
% dh = min(max(mesh.p(:,1))/(ngrid-1),max(mesh.p(:,2))/(ngrid-1));
% alpha2 = sqrt(dh * alpha1 / coefdiffadv);