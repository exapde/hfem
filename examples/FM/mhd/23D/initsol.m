function udg = initsol(p,icase,gam)

%   udg(:,1,:) :        rho
%   udg(:,2,:) :        rho*ux
%   udg(:,3,:) :        rho*uy
%   udg(:,4,:) :        rho*uz
%   udg(:,5,:) :        rho*e
%   udg(:,6,:) :        Bx
%   udg(:,7,:) :        By
%   udg(:,8,:) :        By
%   udg(:,9,:) :        phi

p1 = size(p,1);
p3 = size(p,3);
udg = zeros(p1,9,p3);

x = p(:,1,:);
y = p(:,2,:);

switch (icase)
    case 1 % scalar case
        udg(:,1,:) = 2 + sin(x + y);           
        udg(:,2,:) = udg(:,1,:);                
        udg(:,3,:) = udg(:,1,:); 
        udg(:,5,:) = 5*ones(size(udg(:,1,:))) + udg(:,1,:); 
    case 2 % Orszag-Tang vortex
        udg(:,1,:) = gam^2 * ones(size(udg(:,1,:))); 
        ux = -sin(y);
        udg(:,2,:) = ux.*udg(:,1,:); 
        uy = sin(x);
        udg(:,3,:) = uy.*udg(:,1,:);
        udg(:,6,:) = -sin(y);
        udg(:,7,:) = sin(2*x);
        b = udg(:,6,:).*udg(:,6,:) + udg(:,7,:).*udg(:,7,:) + udg(:,8,:).*udg(:,8,:);
        v = udg(:,2,:).*udg(:,2,:) + udg(:,3,:).*udg(:,3,:) + udg(:,4,:).*udg(:,4,:);
        udg(:,5,:) = gam/(gam-1) + v*0.5 + b*0.5;
    case 3 % Rotor problem
        r0 = 0.1;
        r1 = 0.115;
        udg(:,5,:) = 5/(sqrt(4*pi));
        for ie = 1:p3
            for ip = 1:p1
                r = sqrt((x-0.5*ones(size(x))).^2 + (y-0.5*ones(size(x))).^2);
                if r(ip,1,ie) < r0
                    udg(ip,1,ie) = 10;
                    ux = -(y(ip,1,ie) - 0.5)/r0;
                    udg(ip,2,ie) = udg(ip,1,ie)*ux;
                    uy = (x(ip,1,ie) - 0.5)/r0;
                    udg(ip,3,ie) = udg(ip,1,ie)*uy;
                    normu = norm([ux uy])^2;
                    normb = norm([udg(ip,5,ie) udg(ip,6,ie)])^2;
                    udg(ip,4,ie) = 5/(gam-1) + udg(ip,1,ie)*normu/2 + normb/2;
                elseif r0 < r(ip,1,ie) && r(ip,1,ie) < r1
                    lam = (r1-r(ip,1,ie))/(r1-r0);
                    udg(ip,1,ie) = 1 + 9*lam;
                    ux = -lam*(y(ip,1,ie)-0.5)/r(ip,1,ie);
                    udg(ip,2,ie) = udg(ip,1,ie)*ux;
                    uy = lam*(x(ip,1,ie)-0.5)/r(ip,1,ie);
                    udg(ip,3,ie) = udg(ip,1,ie)*uy;
                    normu = norm([ux uy])^2;
                    normb = norm([udg(ip,5,ie) udg(ip,6,ie)])^2;
                    udg(ip,4,ie) = 5/(gam-1) + udg(ip,1,ie)*normu/2 + normb/2;
                else
                    udg(ip,1,ie) = 1;
                    ux = udg(ip,1,ie)*0;
                    uy = udg(ip,1,ie)*0;
                    normu = norm([ux uy])^2;
                    normb = norm([udg(ip,5,ie) udg(ip,6,ie)])^2;
                    udg(ip,4,ie) = 5/(gam-1) + udg(ip,1,ie)*normu/2 + normb/2;
                end
            end
        end
    case 4 % Smooth Alfven wave
        theta = pi/6;
        xx = x*cos(theta) + y*sin(theta);
        c = 0.1*sin(2*pi*xx);
        udg(:,1,:) = ones(size(udg(:,1,:)));
        udg(:,2,:) = -c*sin(theta);
        udg(:,3,:) = c*cos(theta);
        udg(:,4,:) = 0.1*cos(2*pi*xx);
        udg(:,7,:) = udg(:,1,:)*sin(theta)+c*cos(theta);
        udg(:,6,:) = udg(:,1,:)*cos(theta)-c*sin(theta);
        udg(:,8,:) = udg(:,4,:);
        b = udg(:,6,:).*udg(:,6,:) + udg(:,7,:).*udg(:,7,:) + udg(:,8,:).*udg(:,8,:);
        v = udg(:,2,:).*udg(:,2,:) + udg(:,3,:).*udg(:,3,:) + udg(:,4,:).*udg(:,4,:);
        udg(:,5,:) = 0.1/(gam-1) + v*0.5*udg(:,1,:) + b*0.5;
        
    case 5 % Smooth vortex problem
        r = sqrt(x.*x + y.*y);
        dx = -exp(0.5*(1-r.*r))*0.5.*y/pi;
        dy = exp(0.5*(1-r.*r))*0.5.*x/pi;
        dp = (1-r.*r).*exp(1-r.*r)/(8*pi^2) - exp(1-r.*r)/(8*pi^2);
        udg(:,1,:) = ones(size(udg(:,1,:)));           
        udg(:,2,:) = udg(:,1,:) + dx;                
        udg(:,3,:) = udg(:,1,:) + dy;
        udg(:,6,:) = dx; 
        udg(:,7,:) = dy;
        p = udg(:,1,:) + dp;
        b = udg(:,6,:).*udg(:,6,:) + udg(:,7,:).*udg(:,7,:) + udg(:,8,:).*udg(:,8,:);
        v = udg(:,2,:).*udg(:,2,:) + udg(:,3,:).*udg(:,3,:) + udg(:,4,:).*udg(:,4,:);
        udg(:,5,:) = p/(gam-1) + v*0.5 + b*0.5;
end