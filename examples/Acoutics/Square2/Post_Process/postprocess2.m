function ustar = postprocess(master,mesh,master1,mesh1,u,q)
% POSTPROCESS postprocesses the HDG solution to obtain a better solution
%
%      ustar = postprocess(mesh,master,uh,qh)
%
%      MASTER:       Master structure of porder
%      MESH:         Mesh structure of porder
%      MASTER1:      Master structure of porder+1
%      MESH1:        Mesh structure of porder+1
%      U:            Approximate scalar variable
%      Q:            Approximate flux
%      USTAR:        Postprocessed scalar variable


ne    = size(mesh.dgnodes,3);
npl1  = size(mesh1.dgnodes,1);
ustar = zeros(npl1,1,ne);




master.shap = permute(master.shapvl,[1 3 2]);
shapxi = squeeze(master.shap(:,2,:))*diag(master.gwvl);
shapet = squeeze(master.shap(:,3,:))*diag(master.gwvl); 
shap = permute(master.shap,[1 3 2]);

master1.shap1 = permute(master1.shapvl,[1 3 2]);
shapxi1 = squeeze(master1.shap1(:,2,:))*diag(master1.gwvl);
shapet1 = squeeze(master1.shap1(:,3,:))*diag(master1.gwvl); 
shapxi2 = squeeze(master1.shap1(:,2,:));
shapet2 = squeeze(master1.shap1(:,3,:)); 
shap1 = permute(master1.shap1,[1 3 2]);


% shapt = shape2d(master.porder,master.plocal,master1.gpts);
shapt = mkshape(master.porder,master.plocvl,master1.gpvl,mesh.elemtype);
% shap0 = squeeze(shapt(:,1,:));
shap0  = squeeze(shapt(:,:,1));

for i=1:ne
    dg  = mesh.dgnodes(:,:,i);    
    dg1 = mesh1.dgnodes(:,:,i);    

    xxi = (shap(:,:,2))'*(mesh.dgnodes(:,1,i));
    xet = (shap(:,:,3))'*(mesh.dgnodes(:,1,i));
    yxi = (shap(:,:,2))'*(mesh.dgnodes(:,2,i));
    yet = (shap(:,:,3))'*(mesh.dgnodes(:,2,i));
    jac = xxi.*yet - xet.*yxi;
    
    xxi1 = (shap1(:,:,2))'*(mesh1.dgnodes(:,1,i));
    xet1 = (shap1(:,:,3))'*(mesh1.dgnodes(:,1,i));
    yxi1 = (shap1(:,:,2))'*(mesh1.dgnodes(:,2,i));
    yet1 = (shap1(:,:,3))'*(mesh1.dgnodes(:,2,i));
    jac1 = xxi1.*yet1 - xet1.*yxi1;
    
    shapx1 = - shapxi1*diag(yet1) + shapet1*diag(yxi1);
    shapy1 =   shapxi1*diag(xet1) - shapet1*diag(xxi1);
    shapx2 = - shapxi2*diag(yet1) + shapet2*diag(yxi1);
    shapy2 =   shapxi2*diag(xet1) - shapet2*diag(xxi1);

    F  = shap(:,:,1)*(master.gwvl.*jac);
    F1 = shap1(:,:,1)*(master1.gwvl.*jac1);
    K1 = (shapx2*diag(1./jac1)*shapx1'+shapy2*diag(1./jac1)*shapy1');
    Cx = shapx1*shap0';
    Cy = shapy1*shap0';

    K1(end,:) = F1';
    L = (Cx*q(:,1,i)+Cy*q(:,2,i));
    L(end,:) = F'*u(:,1,i);
    ustar(:,1,i) = K1\L;    
end