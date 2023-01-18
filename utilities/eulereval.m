
function sca = eulereval(u,str,gam,mach)
%EULERVAL Calculates derived quantities for the Euler equation variables.
%   SCA=EULERVAL(U,STR,GAM)
%
%      UP(npl,4,nt):   np plus states
%      STR:            String used to specify requested quantity
%                      - STR: 'r' Density
%                      - STR: 'u' u_x velocity
%                      - STR: 'v' u_y velocity
%                      - STR: 'p' Pressure
%                      - STR: 'p0' Total pressure
%                      - STR: 'M' Mach number
%                      - STR; 's' Entropy
%                      - STR; 'T0' Total temperature
%      GAM:            Value of Gamma
%      SCA(npl,4,nt):  Scalar field requested by STR 
%

if nargin<3; gam = 1.4; end

gam1 = gam-1;

if strcmp(str,'r')
    sca = u(:,1,:);
elseif strcmp(str,'p')
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    sca = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
elseif strcmp(str,'p0')
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    u2 = sqrt(uv.^2+vv.^2);
    p = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    M = u2./sqrt(gam*max(0,p./u(:,1,:)));
    p0 = p .* (1+0.5*gam1*M.^2).^(gam/gam1);
    sca = p0;
elseif strcmp(str,'c')
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p = max(0,(gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv)));
    sca = sqrt(gam*p./u(:,1,:));
elseif strcmp(str,'M')
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    u2 = sqrt(uv.^2+vv.^2);
    p = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    sca = u2./sqrt(gam*max(0,p./u(:,1,:)));
elseif strcmp(str,'s')
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    sca = p./(u(:,1,:).^gam);
elseif strcmp(str,'u')
    sca = u(:,2,:)./u(:,1,:);
elseif strcmp(str,'v')
    sca = u(:,3,:)./u(:,1,:);
elseif strcmp(str,'velMag')
    uu = u(:,2,:)./u(:,1,:);
    uv = u(:,3,:)./u(:,1,:);
    sca = sqrt(uu.^2+uv.^2);
elseif strcmp(str,'c2')
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    sca = gam*(gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv))./u(:,1,:);
elseif strcmp(str,'divV')
    r   = u(:,1,:);
    rx  = u(:,5,:);
    rux = u(:,6,:); 
    ry  = u(:,9,:);
    rvy = u(:,11,:);
    
    r1   = 1./r;
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);

    ux  = (rux - rx.*uv).*r1;
    vy  = (rvy - ry.*vv).*r1;
    
    sca = - (ux + vy);
elseif strcmp(str,'ss')
    uv  = u(:,2,:)./u(:,1,:);
    vv  = u(:,3,:)./u(:,1,:);
    c2  = gam*(gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv))./u(:,1,:);

    r   = u(:,1,:);
    rx  = u(:,5,:);
    rux = u(:,6,:); 
    ry  = u(:,9,:);
    rvy = u(:,11,:);

    r1   = 1./r;

    ux  = (rux - rx.*uv).*r1;
    vy  = (rvy - ry.*vv).*r1;

    sca = (ux+vy)./sqrt(c2);
elseif strcmp(str,'vort')
    r   = u(:,1,:);
    uv  = u(:,2,:)./u(:,1,:);
    vv  = u(:,3,:)./u(:,1,:);
    rx  = u(:,5,:);
    rvx = u(:,7,:); 
    ry  = u(:,9,:);
    ruy = u(:,10,:);
    
    vx  = (rvx - rx.*vv)./r;
    uy  = (ruy - ry.*uv)./r;
    
    sca = uy - vx;
elseif strcmp(str,'t')
    r  = u(:,1,:);
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p  = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    sca = (gam*mach^2)*p./r;
elseif strcmp(str,'T')
    r  = u(:,1,:);
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p  = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    T = p ./((gam-1)*r);
    sca = T;
elseif strcmp(str,'gradT_norm')
    r  = u(:,1,:);
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    rE = u(:,4,:);
    rx  = u(:,5,:);
    rux = u(:,6,:);
    rvx = u(:,7,:);
    rEx = u(:,8,:);
    ry  = u(:,9,:);
    ruy = u(:,10,:);
    rvy = u(:,11,:);
    rEy = u(:,12,:);
    
    q = 0.5*(uv.*uv+vv.*vv);
    p = (gam-1)*(rE-r.*q);
    
    ux = (rux - uv.*rx)./r;
    vx = (rvx - vv.*rx)./r;
    qx  = uv.*ux + vv.*vx;
    px  = (gam-1)*(rEx - rx.*q - r.*qx);
    Tx  = (px.*r - p.*rx)./((gam-1)*r.^2);

    uy = (ruy - uv.*ry)./r;
    vy = (rvy - vv.*ry)./r;
    qy  = uv.*uy + vv.*vy;
    py  = (gam-1)*(rEy - ry.*q - r.*qy);
    Ty  = (py.*r - p.*ry)./((gam-1)*r.^2);
    
    sca = sqrt(Tx.^2+Ty.^2);
elseif strcmp(str,'T0')
    r  = u(:,1,:);
    ru  = u(:,2,:);
    rv  = u(:,3,:);
    uv = ru./r;
    vv = rv./r;
    rE  = u(:,4,:);
    
    p  = (gam-1)*(rE - 0.5*(ru.*uv + rv.*vv));
    T = p ./((gam-1)*r);
    T0 = T + 0.5 * (uv.*uv + vv.*vv) / gam;
    
    sca = T0;
%     sca = rE./r-0.5*((gam-1)/gam)*(ru.^2+rv.^2)./r;
elseif strcmp(str,'S_fro')
    r   = u(:,1,:);
    uv  = u(:,2,:)./u(:,1,:);
    vv  = u(:,3,:)./u(:,1,:);
    rx  = u(:,5,:);
    rux = u(:,6,:); 
    rvx = u(:,7,:); 
    ry  = u(:,9,:);
    ruy = u(:,10,:);
    rvy = u(:,11,:);
    
    ux  = (rux - rx.*uv)./r;
    vx  = (rvx - rx.*vv)./r;
    uy  = (ruy - ry.*uv)./r;
    vy  = (rvy - ry.*vv)./r;

    S_xx = ux;
    S_xy = 0.5*(uy+vx);
    S_yx = S_xy;
    S_yy = vy;
    S_fro = sqrt(S_xx.^2 + S_xy.^2 + S_yx.^2 + S_yy.^2);
    
    sca = S_fro;
elseif strcmp(str,'shear')
    r   = u(:,1,:);
    uv  = u(:,2,:)./u(:,1,:);
    vv  = u(:,3,:)./u(:,1,:);
    rx  = u(:,5,:);
    rvx = u(:,7,:); 
    ry  = u(:,9,:);
    ruy = u(:,10,:);
    
    vx  = (rvx - rx.*vv)./r;
    uy  = (ruy - ry.*uv)./r;
    
    sca = abs(uy+vx);
elseif strcmp(str,'H')
    uv  = u(:,2,:)./u(:,1,:);
    vv  = u(:,3,:)./u(:,1,:);
    E   = u(:,4,:)./u(:,1,:);
    
    q = 0.5*(uv.*uv+vv.*vv);
    H = gam*E - (gam-1)*q;
    
    sca = H;
else
    error('unknonw case');
end
