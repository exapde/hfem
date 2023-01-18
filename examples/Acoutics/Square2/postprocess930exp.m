function vstar = postprocess930exp(master,mesh,master1,mesh1,v,vh)
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
npv = size(mesh.dgnodes,1);
nps = size(master.plocfc,1);

elcon = elconnectivities(mesh,size(mesh.f,1));

vstar = zeros(npl1,1,ne);
perm = master.perm;
sh1d = master.shapfc;

master.shap = permute(master.shapvl,[1 3 2]);
master1.shap = permute(master1.shapvl,[1 3 2]);

shap    = squeeze(master.shap(:,1,:));
shapxi = squeeze(master.shap(:,2,:))*diag(master.gwvl);
shapet = squeeze(master.shap(:,3,:))*diag(master.gwvl); 
shap1   = squeeze(master1.shap(:,1,:));
shapxi1 = squeeze(master1.shap(:,2,:))*diag(master1.gwvl);
shapet1 = squeeze(master1.shap(:,3,:))*diag(master1.gwvl); 
shapxi2 = squeeze(master1.shap(:,2,:));
shapet2 = squeeze(master1.shap(:,3,:)); 

shapt = mkshape(master.porder,master.plocvl,master1.gpvl,mesh.elemtype);
shap0  = squeeze(shapt(:,:,1));


for i=1:ne
    dg  = mesh.dgnodes(:,:,i);    
    dg1 = mesh1.dgnodes(:,:,i);  
    A_K = zeros(2*npv,2*npv);
    B_K = zeros(2*npv,npv);
    C_K = zeros(2*npv,3*nps);

    xxi = squeeze(master.shap(:,2,:))'*squeeze(dg(:,1));
    xet = squeeze(master.shap(:,3,:))'*squeeze(dg(:,1));
    yxi = squeeze(master.shap(:,2,:))'*squeeze(dg(:,2));
    yet = squeeze(master.shap(:,3,:))'*squeeze(dg(:,2));
    jac = xxi.*yet - xet.*yxi;
    
    shapx =   shapxi*diag(yet) - shapet*diag(yxi);
    shapy = - shapxi*diag(xet) + shapet*diag(xxi);
    Mass  = shap(:,:,1)*diag(master.gwvl.*jac)*shap(:,:,1)';
    
     % form A_K
    A_K(1:npv,1:npv) = Mass;
    A_K((npv+1):2*npv,(npv+1):2*npv) = Mass; 
    % form B_K
    B_K(1:npv,1:npv) = shapx*shap(:,:,1)';
    B_K((npv+1):2*npv,1:npv) = shapy*shap(:,:,1)';
    
    
    for j = 1:3 % loop over each face of the element i
        
        I = perm(:,j,1);
        J = ((j-1)*nps+1):j*nps;           
        xxi = (sh1d(:,:,2))'*(mesh.dgnodes(I,1,i)); 
        yxi = (sh1d(:,:,2))'*(mesh.dgnodes(I,2,i));  
        dsdxi = sqrt(xxi.^2+yxi.^2);
        nl = [yxi./dsdxi,-xxi./dsdxi];
        dws = master.gwfc.*dsdxi;
     
     % form C_K
        C_K(I,J) = C_K(I,J) + sh1d(:,:,1)*diag(dws.*nl(:,1))*sh1d(:,:,1)';
        C_K(npv+I,J) = C_K(npv+I,J) + sh1d(:,:,1)*diag(dws.*nl(:,2))*sh1d(:,:,1)';
    end
    
    % Solve for p
%     imap=elcon(:,i);
    p = A_K\(-B_K*v(:,1,i)+C_K*vh(:,i));
    p = reshape(p,[npv,2]);
    
    
    xxi1 = squeeze(master1.shap(:,2,:))'*squeeze(dg1(:,1));
    xet1 = squeeze(master1.shap(:,3,:))'*squeeze(dg1(:,1));
    yxi1 = squeeze(master1.shap(:,2,:))'*squeeze(dg1(:,2));
    yet1 = squeeze(master1.shap(:,3,:))'*squeeze(dg1(:,2));
    jac1 = xxi1.*yet1 - xet1.*yxi1;
    shapx1 =   shapxi1*diag(yet1) - shapet1*diag(yxi1);
    shapy1 = - shapxi1*diag(xet1) + shapet1*diag(xxi1);
    shapx2 =   shapxi2*diag(yet1) - shapet2*diag(yxi1);
    shapy2 = - shapxi2*diag(xet1) + shapet2*diag(xxi1);

    F  = shap*(master.gwvl.*jac);
    F1 = shap1*(master1.gwvl.*jac1);
    K1 = (shapx2*diag(1./jac1)*shapx1'+shapy2*diag(1./jac1)*shapy1');
    Cx = shapx1*shap0';
    Cy = shapy1*shap0';

    K1(end,:) = F1';
    L = (Cx*p(:,1)+Cy*p(:,2));
    L(end,:) = F'*v(:,1,i);
    vstar(:,1,i) = K1\L;    
end


