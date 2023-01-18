function Eh = get_Eh(master, mesh, app, UDG, VDG, QDG, UHAT)
%
%      MASTER:       Master structure
%      MESH:         Mesh structure

tau = app.tau;
nt  = size(mesh.t,1);
perm = master.perm;
shap = permute(master.shapvl,[1 2 3]);
sh1d = permute(master.shapfc,[1 2 3]);

Eh = 0;
for i = 1:nt % loop over each element
    xxi = (shap(:,:,2))'*(mesh.dgnodes(:,1,i));
    xet = (shap(:,:,3))'*(mesh.dgnodes(:,1,i));
    yxi = (shap(:,:,2))'*(mesh.dgnodes(:,2,i));
    yet = (shap(:,:,3))'*(mesh.dgnodes(:,2,i));
    jac = xxi.*yet - xet.*yxi;
    Mass  = shap(:,:,1)*diag(master.gwvl.*jac)*shap(:,:,1)';
    
    uhat = UHAT(:,:,i);
    uh = UDG(:,1,i);
    vh = VDG(:,1,i);
    qx = QDG(:,1,i);
    qy = QDG(:,2,i);        
    Eh = Eh + vh'*Mass*vh + qx'*Mass*qx + qy'*Mass*qy;
    
    for j = 1:3 % loop over each face of the element i
        I = perm(:,j,1);                        
        xxi = (sh1d(:,:,2))'*(mesh.dgnodes(I,1,i)); 
        yxi = (sh1d(:,:,2))'*(mesh.dgnodes(I,2,i));  
        dsdxi = sqrt(xxi.^2+yxi.^2);        
        dws = master.gwfc.*dsdxi;                                        
        tmp = sh1d(:,:,1)*diag(tau*dws)*sh1d(:,:,1)';
                
        eh = uh(I,1) - uhat(:,j);
        Eh = Eh + eh'*tmp*eh;               
    end        
end
Eh = 0.5*Eh;

