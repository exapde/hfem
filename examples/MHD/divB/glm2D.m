function [alpha1,eiga] = glm2D(mesh,app,UDG,UH,ngrid)

% r = UDG(:,1,:);
% vx = UDG(:,2,:)./UDG(:,1,:);
% vy = UDG(:,3,:)./UDG(:,1,:);
% re = UDG(:,4,:);
% bx = UDG(:,5,:);
% by = UDG(:,6,:);

r = UH(1,:);
vx = UH(2,:)./UH(1,:);
vy = UH(3,:)./UH(1,:);
re = UH(4,:);
bx = UH(5,:);
by = UH(6,:);
psi = UH(7,:);
b = bx.*bx + by.*by;
gam = app.arg{1};
q = 0.5*(vx.*vx + vy.*vy);
p = (gam - 1)*(re - r.*q - b*0.5-0.5*psi.^2);

a    = sqrt(gam*p./r); % sound speed
ca   = b./sqrt(r);   % Alfven speed
cax  = bx./sqrt(r);
cay  = by./sqrt(r);

% cf = g(B,v); 
cfx = a.*a +  b./r + sqrt((a.*a +  b./r).^2 - 4*a.*a.*bx.*bx./r); 
cfy = a.*a +  b./r + sqrt((a.*a +  b./r).^2 - 4*a.*a.*by.*by./r);
cfx = sqrt(cfx*0.5);
cfy = sqrt(cfy*0.5);
vxcfx = abs(vx) + cfx;
vycfy = abs(vy) + cfy;
alpha1 = max(max(max([vxcfx,vycfy])));
%alpha1 = 0.1;
%M = cat(2,vxcfx(:),vycfy(:),vzcfz(:));
%alpha1 = reshape(max(M,[],2),size(vxcfx));
% eigax = max(max(max(vxcfx+cax)));
% eigay = max(max(max(vycfy+cay)));
% eiga = eigax + eigay;
%eiga = max(max(max([vxcfx+cax,vycfy+cay])));
eiga = max(max(max(vxcfx + cax + vycfy + cay))); %why this choice of tau?  should be lambda_max if LF
%eiga = alpha1;
%%%% Dedner & al.
%alpha2 = sqrt(0.18 * alpha1); 

%%%% Mignone & al.
% coefdiffadv = 1; % [0,1]
% dh = min(max(mesh.p(:,1))/(ngrid-1),max(mesh.p(:,2))/(ngrid-1));
% alpha2 = sqrt(dh * alpha1 / coefdiffadv);