function p0 = segmentintersection(p1, v1,  p2, v2)

% x1 + s*a1 = x2 + t*a2
% y1 + s*b1 = y2 + t*b2
% b2*(x1 + s*a1) = b2*(x2 + t*a2)
% a2*(y1 + s*b1) = a2*(y2 + t*b2)
% b2*x1 - a2*y1 + s*(a1*b2-a2*b1) = b2*x2 - a2*y2
%  s*(a1*b2-a2*b1) = b2*x2 - a2*y2 - b2*x1 + a2*y1

nd = length(p1);
if nd == 2
    a1 = v1(1); b1 = v1(2);
    a2 = v2(1); b2 = v2(2);    
    x1 = p1(1); y1 = p1(2);
    x2 = p2(1); y2 = p2(2);
    if abs(a1*b2 - a2*b1)>1e-12
        s = (x2*b2 - x1*b2 + y1*a2 - y2*a2)/(a1*b2 - a2*b1);
        if abs(a2)>abs(b2)
            t = (x1+s*a1-x2)/a2;
        else
            t = (y1+s*b1-y2)/b2;
        end
        p0 = [s;t];
    else
        p0 = [];
    end
%     [x y]
%     figure(1); clf;
%     plot(x,y,'or','MarkerSize',8);
%     hold on;    
%     c = 100;
%     pa = p1(:) + c*v1(:);
%     pb = p1(:) - c*v1(:);    
%     plot([pa(1),pb(1)],[pa(2),pb(2)],'-b');
%     pa = p2(:) + c*v2(:);
%     pb = p2(:) - c*v2(:);
%     plot([pa(1),pb(1)],[pa(2),pb(2)],'-b');
%     pause
else
    a1 = v1(1); b1 = v1(2); c1 = v1(3);
    a2 = v2(1); b2 = v2(2); c2 = v2(3);   
    x1 = p1(1); y1 = p1(2); z1 = p1(3);
    x2 = p2(1); y2 = p2(2); z2 = p2(3);
    if abs(a1*b2 - a2*b1)>1e-12
        y = (x2*b2*b1 - x1*b2*b1 + y1*a1*b2 - y2*a2*b1)/(a1*b2 - a2*b1);
        if abs(b1)>abs(b2)
            x = x1 + (y-y1)*a1/b1;                        
        else
            x = x2 + (y-y2)*a2/b2;            
        end        
        if abs(b1)>abs(a1)
            za = z1 + (y-y1)*c1/b1;
        else
            za = z1 + (x-x1)*c1/a1;
        end
        if abs(b2)>abs(a2)            
            zb = z2 + (y-y2)*c2/b2;
        else
            zb = z2 + (x-x2)*c2/a2;
        end
        if abs(za-zb)<1e-12
            z = 0.5*(za+zb);
            p0 = [x y z];
        else
            p0 = [];
        end        
    elseif abs(c1*b2 - c2*b1)>1e-12
        y = (z2*b2*b1 - z1*b2*b1 + y1*c1*b2 - y2*c2*b1)/(c1*b2 - c2*b1);
        if abs(b1)>abs(b2)            
            z = z1 + (y-y1)*c1/b1;
        else            
            z = z2 + (y-y2)*c2/b2;
        end
        if abs(b1)>abs(c1)
            xa = x1 + (y-y1)*a1/b1;
        else
            xa = x1 + (z-z1)*a1/c1;
        end
        if abs(b2)>abs(c2)            
            xb = x2 + (y-y2)*a2/b2;
        else
            xb = x2 + (z-z2)*a2/c2;
        end        
        if abs(xa-xb)<1e-12
            x = 0.5*(xa+xb);
            p0 = [x y z];
        else
            p0 = [];
        end        
    elseif abs(c1*a2 - c2*a1)>1e-12
        x = (z2*a2*a1 - z1*a2*a1 + x1*c1*a2 - x2*c2*a1)/(c1*a2 - c2*a1);
        if abs(a1)>abs(a2)            
            z = z1 + (x-x1)*c1/a1;
        else            
            z = z2 + (x-x2)*c2/a2;
        end
        if abs(a1)>abs(c1)
            ya = y1 + (x-x1)*b1/a1;
        else
            ya = y1 + (z-z1)*b1/c1;
        end
        if abs(a2)>abs(c2)            
            yb = y2 + (x-x2)*b2/a2;
        else
            yb = y2 + (z-z2)*b2/c2;
        end        
        if abs(ya-yb)<1e-12
            y = 0.5*(ya+yb);
            p0 = [x y z];
        else
            p0 = [];
        end         
    else
        p0 = [];
    end    
end


return;

nd = 2;
p1 = rand(1,nd);
v1 = rand(1,nd);
p2 = rand(1,nd);
v2 = rand(1,nd);
p0 = lineintersection(p1, v1,  p2, v2);

