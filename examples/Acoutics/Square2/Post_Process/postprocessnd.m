function UDGstar = postprocessnd(master,mesh,master1,mesh1,UDG)

%HDG_POSTPROCESS postprocesses the HDG solution to obtain a better solution.
%   [ustarh]=hdg_postprocess(mesh,master,uh,qh,uhath)
%
%      MASTER:       Master structure of porder
%      MESH:         Mesh structure of porder
%      MASTER1:      Master structure of porder+1
%      MESH1:        Mesh structure of porder+1
%      UH:           Approximate scalar variable
%      QH:           Approximate flux
%      USTARH:       Postprocessed scalar variable


nd    = master.nd;
nc    = size(UDG,2);
nch   = nc/(nd+1);
npv   = master.npv;
ngv   = master.ngv;
npv1  = master1.npv;
ngv1  = master1.ngv;
ne    = size(mesh1.dgnodes,3);

%shapvt    = master.shapvt(:,:,1); 
dshapvt   = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);
dshapvt1  = reshape(permute(master1.shapvt(:,:,2:nd+1),[1 3 2]),[ngv1*nd npv1]);

shapt = mkshape(master.porder,master.plocvl,master1.gpvl,mesh.elemtype);
shap  = squeeze(shapt(:,:,1));

UDGstar = zeros(npv1,nch,ne);
for i=1:ne
    dg  = mesh.dgnodes(:,:,i);
    dg1 = mesh1.dgnodes(:,:,i);        

    % compute the Jacobian matrix at Gauss points: dx/dxi
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[ngv nd nd]);        
    [jac,Xx] = volgeom(Jg);
    
    Jg1 = dshapvt1*dg1(:,1:nd);
    Jg1 = reshape(Jg1,[ngv1 nd nd]);        
    [jac1,Xx1] = volgeom(Jg1);        
        
    F  = master.shapvl(:,:,1)*(master.gwvl.*jac);
    F1 = master1.shapvl(:,:,1)*(master1.gwvl.*jac1);
    
% Xx = -inv(Jg)*jac: minus the inverse of the Jacobian matrix times 
%                          the determinant of the Jacobian matrix 
       
    dshapdx  = zeros(npv,ngv,nd);
    dshapdx1 = zeros(npv1,ngv1,nd);
    for ii=1:nd                
        dshapdx(:,:,ii)  = bsxfun(@times,master.shapvl(:,:,2),Xx(:,ii,1)');
        dshapdx1(:,:,ii) = bsxfun(@times,master1.shapvl(:,:,2),Xx1(:,ii,1)');        
        for jj=2:nd
            dshapdx(:,:,ii)  = dshapdx(:,:,ii) + bsxfun(@times,master.shapvl(:,:,1+jj),Xx(:,ii,jj)');
            dshapdx1(:,:,ii) = dshapdx1(:,:,ii) + bsxfun(@times,master1.shapvl(:,:,1+jj),Xx1(:,ii,jj)');
        end            
    end
    
    K1 = dshapdx1(:,:,1)*diag(master1.gwvl./jac1)*(dshapdx1(:,:,1)');
    C1 = dshapdx1(:,:,1)*diag(master1.gwvl)*shap';
    L1 = C1*UDG(:,nch+1:2*nch,i);
    for ii=2:nd
        K1 = K1 + dshapdx1(:,:,ii)*diag(master1.gwvl./jac1)*(dshapdx1(:,:,ii)');
        C1 = dshapdx1(:,:,ii)*diag(master1.gwvl)*shap';
        L1 = L1 + C1*UDG(:,ii*nch+1:(ii+1)*nch,i);
    end    

    K1(end,:) = F1';    
    L1(end,:) = F'*UDG(:,1:nch,i);
    UDGstar(:,1:nch,i) = K1\L1;    
end


function [jac,Xx] = volgeom(Jg)

ngv = size(Jg,1); 
nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;
        Xx = -ones(ngv,1);
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);
        Xx(:,1,1) = -Jg(:,2,2);
        Xx(:,2,1) = Jg(:,2,1);
        Xx(:,1,2) = Jg(:,1,2);
        Xx(:,2,2) = -Jg(:,1,1);
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);            
        Xx(:,1,1) = Jg(:,2,3).*Jg(:,3,2) - Jg(:,2,2).*Jg(:,3,3);
        Xx(:,2,1) = Jg(:,2,1).*Jg(:,3,3) - Jg(:,2,3).*Jg(:,3,1);
        Xx(:,3,1) = Jg(:,2,2).*Jg(:,3,1) - Jg(:,2,1).*Jg(:,3,2);
        Xx(:,1,2) = Jg(:,1,2).*Jg(:,3,3) - Jg(:,1,3).*Jg(:,3,2);
        Xx(:,2,2) = Jg(:,1,3).*Jg(:,3,1) - Jg(:,1,1).*Jg(:,3,3);
        Xx(:,3,2) = Jg(:,1,1).*Jg(:,3,2) - Jg(:,1,2).*Jg(:,3,1);
        Xx(:,1,3) = Jg(:,1,3).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,3);
        Xx(:,2,3) = Jg(:,1,1).*Jg(:,2,3) - Jg(:,1,3).*Jg(:,2,1);
        Xx(:,3,3) = Jg(:,1,2).*Jg(:,2,1) - Jg(:,1,1).*Jg(:,2,2);
    otherwise
        error('Dimension is not implemented');
end




