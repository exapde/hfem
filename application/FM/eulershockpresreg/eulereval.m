function sca = eulereval(u,str,gam)
%EULERVAL Calculates derived quantities for the Euler equation variables.
%   SCA=EULERVAL(U,STR,GAM)
%
%      UP(npl,4,nt):   np plus states
%      STR:            String used to specify requested quantity
%                      - STR: 'r' Density
%                      - STR: 'u' u_x velocity
%                      - STR: 'v' u_y velocity
%                      - STR: 'p' Density
%                      - STR: 'M' Density
%                      - STR; 's' Entropy
%      GAM:            Value of Gamma
%      SCA(npl,4,nt):  Scalar field requested by STR 
%
if str == 'r'
    sca = u(:,1,:);
elseif str == 'p'
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    sca = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
elseif str == 'c'
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    sca = sqrt(gam*p./u(:,1,:));
elseif str == 'Jp'
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    c = sqrt(gam*p./u(:,1,:));
    sca = u(:,2,:) + 2*c/(gam-1);
elseif str == 'Jm'
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    c = sqrt(gam*p./u(:,1,:));
    sca = u(:,2,:) - 2*c/(gam-1);
elseif str == 'M'
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    u2 = sqrt(uv.^2+vv.^2);
    p = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    sca = u2./sqrt(gam*p./u(:,1,:));
elseif str == 's'
    uv = u(:,2,:)./u(:,1,:);
    vv = u(:,3,:)./u(:,1,:);
    p = (gam-1)*(u(:,4,:) - 0.5*(u(:,2,:).*uv + u(:,3,:).*vv));
    sca = p./(u(:,1,:).^gam);
elseif str == 'u'
    sca = u(:,2,:)./u(:,1,:);
elseif str == 'v'
    sca = u(:,3,:)./u(:,1,:);
else
    error('unknonw case');
end