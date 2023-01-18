function err = calerror2(UDG,mesh,master)

if mesh.porder==0
    err = 0;
    ne = mesh.ne; 
    jac = 1/ne;
    for i = 1:ne
        p1 = mesh.p(mesh.t(i,:),:);
        pm = mean(p1,1);    
        if abs(pm(1)-pm(2)-0.2)<1e-10    
            err = err + 0.5*jac*((UDG(1,1,i)-1).^2+UDG(1,1,i).^2);
%             if abs(p1(1,1)-p1(1,2)-0.2)<1e-10
%                 pn = mean(p1(2:3,:),1);
%                 pa = [p1(1,:); pn; p1(2,:)];
%                 pb = [p1(1,:); pn; p1(3,:)];
%             elseif abs(p1(2,1)-p1(2,2)-0.2)<1e-10
%                 pn = mean(p1([1 3],:),1);
%                 pa = [p1(2,:); pn; p1(1,:)];
%                 pb = [p1(2,:); pn; p1(3,:)];
%             elseif abs(p1(3,1)-p1(3,2)-0.2)<1e-10
%                 pn = mean(p1([1 2],:),1);            
%                 pa = [p1(3,:); pn; p1(1,:)];
%                 pb = [p1(3,:); pn; p1(2,:)];
%             else
%                 error('something wrong');
%             end
%             pc = mean(pa,1);   
%             if abs(pc(1)-pc(2)-0.2)<1e-10        
%                 error('something wrong');
%             end
%             if pc(1)<pc(2)+0.2
%                 uex = 1;                
%                 err = err + 0.5*jac*(UDG(1,1,i)-uex).^2;                 
%             else
%                 uex = 0;                
%                 err = err + 0.5*jac*(UDG(1,1,i)-uex).^2;                 
%             end
%             pc = mean(pb,1);   
%             if abs(pc(1)-pc(2)-0.2)<1e-10        
%                 error('something wrong');
%             end
%             if pc(1)<pc(2)+0.2
%                 uex = 1;                
%                 err = err + 0.5*jac*(UDG(1,1,i)-uex).^2;                 
%             else
%                 uex = 0;                
%                 err = err + 0.5*jac*(UDG(1,1,i)-uex).^2;                 
%             end
        elseif pm(1)<pm(2)+0.2
            uex = 1;                
            err = err + jac*(UDG(1,1,i)-uex).^2;                 
        else
            uex = 0;                
            err = err + jac*(UDG(1,1,i)-uex).^2;                 
        end       
    end
    err = sqrt(err);
    return;
end

[npv, nc, ne] = size(UDG);
ngv = master.ngv;
nd  = master.nd;

shapvt    = squeeze(master.shapvt(:,:,1));
dshapvt   = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);

plocal = mesh.plocal;
xi  = plocal(:,1);
eta = plocal(:,2);    
philocal(:,1) = 1 - xi - eta;
philocal(:,2) = xi;
philocal(:,3) = eta;

err = zeros(nc,1);
for i = 1:ne
    dg = mesh.dgnodes(:,:,i);
    
    % compute the Jacobian matrix at Gauss points: dx/dxi
%     Jg = dshapvt*dg(:,1:nd);
%     Jg = reshape(Jg,[ngv nd nd]);        
%     jab = volgeom(Jg);    
    jac = 2/ne;    
                
    p1 = mesh.p(mesh.t(i,:),:);
    x1 = p1(1,1); y1 = p1(1,2);
    x2 = p1(2,1); y2 = p1(2,2);
    x3 = p1(3,1); y3 = p1(3,2);
    d1 = (y3-y1)*(x2-x1) - (x3-x1)*(y2-y1);
    pm = mean(p1,1);    
    if abs(pm(1)-pm(2)-0.2)<1e-10        
        if abs(p1(1,1)-p1(1,2)-0.2)<1e-10
            pn = mean(p1(2:3,:),1);
            pa = [p1(1,:); pn; p1(2,:)];
            pb = [p1(1,:); pn; p1(3,:)];            
        elseif abs(p1(2,1)-p1(2,2)-0.2)<1e-10
            pn = mean(p1([1 3],:),1);
            pa = [p1(2,:); pn; p1(1,:)];
            pb = [p1(2,:); pn; p1(3,:)];
        elseif abs(p1(3,1)-p1(3,2)-0.2)<1e-10
            pn = mean(p1([1 2],:),1);       
            pa = [p1(3,:); pn; p1(1,:)];
            pb = [p1(3,:); pn; p1(2,:)];
        else
            error('something wrong');
        end
        % dgnodes on the element a
        dga = philocal*pa;        
        xia = ((y3-y1)*(dga(:,1)-x1) - (x3-x1)*(dga(:,2)-y1))/d1;
        eta = -((y2-y1)*(dga(:,1)-x1) - (x2-x1)*(dga(:,2)-y1))/d1;
        sha = mkshape(mesh.porder,master.plocvl,[xia eta],mesh.elemtype);
        ua = sha(:,:,1)'*UDG(:,:,i);
        udga = shapvt*ua;
        
        % dgnodes on the element b
        dgb = philocal*pb;
        xib = ((y3-y1)*(dgb(:,1)-x1) - (x3-x1)*(dgb(:,2)-y1))/d1;
        etb = -((y2-y1)*(dgb(:,1)-x1) - (x2-x1)*(dgb(:,2)-y1))/d1;
        shb = mkshape(mesh.porder,master.plocvl,[xib etb],mesh.elemtype);         
        ub = shb(:,:,1)'*UDG(:,:,i);
        udgb = shapvt*ub;
        
        pc = mean(pa,1);   
        if abs(pc(1)-pc(2)-0.2)<1e-10        
            error('something wrong');
        end
        if pc(1)<pc(2)+0.2
            uea = 1;                            
        else
            uea = 0;                            
        end
        
        pc = mean(pb,1);   
        if abs(pc(1)-pc(2)-0.2)<1e-10        
            error('something wrong');
        end
        if pc(1)<pc(2)+0.2
            ueb = 1;                            
        else
            ueb = 0;                            
        end        
        
        if uea==ueb
            error('something wrong');
        else                     
            err = err + 0.5*(master.gwvl.*jac)'*(udga-uea).^2;  
            err = err + 0.5*(master.gwvl.*jac)'*(udgb-ueb).^2;  
        end
                    
    elseif pm(1)<pm(2)+0.2
        uex = 1;
        udgg = shapvt*UDG(:,:,i);      
        err = err + (master.gwvl.*jac)'*(udgg-uex).^2;          
    else
        uex = 0;
        udgg = shapvt*UDG(:,:,i);    
        err = err + (master.gwvl.*jac)'*(udgg-uex).^2;          
    end        
end
err  = sqrt(err);

% function [jac] = volgeom(Jg)
% 
% nd  = size(Jg,2);
% switch nd
%     case 1
%         jac = Jg;        
%     case 2
%         jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);        
%     case 3
%         jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
%               Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
%               Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);                    
%     otherwise
%         error('Dimension is not implemented');
% end
% 
% 
% 
% 
