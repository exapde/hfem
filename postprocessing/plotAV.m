
function epsilon = plotAV(mesh,app,UDG)

k_h = 1.5;
alpha = 1.0e4;
beta = 0.01;
porder = mesh.porder;
gam = app.param(1);

if mesh.nd == 2
    h = mesh.dgnodes(:,3,:);
    wallDistance = mesh.dgnodes(:,4,:);

    r = UDG(:,1,:);
    ru = UDG(:,2,:);
    rv = UDG(:,3,:);
    rE = UDG(:,4,:);
    rx = UDG(:,5,:);
    rux = UDG(:,6,:);
    rvx = UDG(:,7,:);
    rEx = UDG(:,8,:);
    ry = UDG(:,9,:);
    ruy = UDG(:,10,:);
    rvy = UDG(:,11,:);
    rEy = UDG(:,12,:);
    
    u = ru ./ r;
    v = rv ./ r;
    velMag = sqrt(u.^2+v.^2);
    
    ux = (rux - u.*rx ) ./ r;
    uy = (ruy - u.*ry ) ./ r;
    vx = (rvx - v.*rx ) ./ r;
    vy = (rvy - v.*ry ) ./ r;
    div_v = ux + vy;

    p = (gam-1)*(rE-r.*velMag.^2/2);
    rH = rE + p;
    H = rH./r;

    c_cr = sqrt((2*(gam-1)*H)/(gam+1));            % Speed of sound at critical conditions.
    x = ((k_h*h/porder).*div_v)./c_cr;
    f = x - beta;
    epsilon = (k_h*h/porder).*(velMag+sqrt(gam*max(1.0e-16,p)./r)).*f;
    scaplot(mesh,epsilon);
    
elseif mesh.nd == 3
    error('AV visualization in 3D not implemented yet.');
    
end
