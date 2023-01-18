function err = calerror_IM_reference(UDG,mesh,master,UDGimp,meshAll,time)

[npv nc ne] = size(UDG);
ngv = master.ngv;
nd  = master.nd;

shapvt    = squeeze(master.shapvt(:,:,1));
dshapvt   = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);

submesh_index = find_submesh_index(mesh,meshAll);

err = zeros(nc,1);
for i = 1:ne
    dg = mesh.dgnodes(:,:,i);
    si = submesh_index(i);
    
    % compute the Jacobian matrix at Gauss points: dx/dxi
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[ngv nd nd]);        
    jac = volgeom(Jg);
                
%     pg = shapvt*dg;    
%     if nargin<5
%         udge = feval(func,pg);
%     else
%         udge = feval(func,pg,time);
%     end
    
    udgg = shapvt*UDG(:,:,i);
    udge = shapvt*UDGimp(:,:,si);
    for j = 1:nc
        err(j) = err(j) + (master.gwvl.*jac)'*(udgg(:,j)-udge(:,j)).^2;  
    end    
end
err  = sqrt(err);

function [jac] = volgeom(Jg)

nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;        
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);        
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);                    
    otherwise
        error('Dimension is not implemented');
end




