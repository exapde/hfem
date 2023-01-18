function [divdg,divdu,err] = divergence(UDG,UH,mesh,master,c)

nt  = size(mesh.t,1);
%nf  = size(mesh.f,1);
npv = size(mesh.dgnodes,1);
nps = size(master.plocfc,1);
ne = mesh.ne;

%shap = squeeze(master.shap(:,1,:));
shapxi = squeeze(master.shapvl(:,:,2))*diag(master.gwvl);
shapet = squeeze(master.shapvl(:,:,3))*diag(master.gwvl); 
%sh1d = squeeze(master.sh1d(:,1,:));
perm = master.perm;

% allocate memory for the elemental matrices
divdg = zeros(npv,1,ne);
divdu = zeros(npv,1,ne);

shap = permute(master.shapvl,[1 2 3]);
sh1d = permute(master.shapfc,[1 2 3]);

err = [0 0];
for i = 1:nt % loop over each element
    xxi = (shap(:,:,2))'*(mesh.dgnodes(:,1,i));
    xet = (shap(:,:,3))'*(mesh.dgnodes(:,1,i));
    yxi = (shap(:,:,2))'*(mesh.dgnodes(:,2,i));
    yet = (shap(:,:,3))'*(mesh.dgnodes(:,2,i));
    jac = xxi.*yet - xet.*yxi;
    shapx =   shapxi*diag(yet) - shapet*diag(yxi);
    shapy = - shapxi*diag(xet) + shapet*diag(xxi);
    Mass  = shap(:,:,1)*diag(master.gwvl.*jac)*shap(:,:,1)';
    
    Bx = shapx*shap(:,:,1)';
    By = shapy*shap(:,:,1)';    
    
    % form F_K
    xg = (shap(:,:,1))'*mesh.dgnodes(:,:,i);
    fg = pi*sin(pi*(xg(:,1) + xg(:,2)));
    F = shap(:,:,1)*(master.gwvl.*jac.*fg);
    divpu = Mass\F;
    
    % allocate memory for the elemental matrices
    Cx = zeros(npv,3*nps);    
    Cy = zeros(npv,3*nps);    
    EJ = zeros(npv,3*nps);    
    Dx = zeros(npv,npv);    
    Dy = zeros(npv,npv);    
    EI = zeros(npv,npv);    
    for j = 1:3 % loop over each face of the element i
        I = perm(:,j,1);
        J = ((j-1)*nps+1):j*nps;
                
        xxi = (sh1d(:,:,2))'*(mesh.dgnodes(I,1,i)); 
        yxi = (sh1d(:,:,2))'*(mesh.dgnodes(I,2,i));  
        dsdxi = sqrt(xxi.^2+yxi.^2);
        nl = [yxi./dsdxi,-xxi./dsdxi];
        dws = master.gwfc.*dsdxi;                
                
        % form C_K
        Cx(I,J) = Cx(I,J) + sh1d(:,:,1)*diag(dws.*nl(:,1))*sh1d(:,:,1)';
        Cy(I,J) = Cy(I,J) + sh1d(:,:,1)*diag(dws.*nl(:,2))*sh1d(:,:,1)';
        
        Dx(I,I) = Dx(I,I) + sh1d(:,:,1)*diag(dws.*nl(:,1))*sh1d(:,:,1)';
        Dy(I,I) = Dy(I,I) + sh1d(:,:,1)*diag(dws.*nl(:,2))*sh1d(:,:,1)';
        
        An = nl(:,1)*c(1)+nl(:,2).*c(2);
        tau = 0.5*(abs(An))+1e-10;
        EI(I,I) = EI(I,I) + sh1d(:,:,1)*diag(dws.*tau)*sh1d(:,:,1)';
        EJ(I,J) = EJ(I,J) + sh1d(:,:,1)*diag(dws.*tau)*sh1d(:,:,1)';
    end        
    
%     tau = 0.5*(abs(An))+1e-10;
% %tau = 1/2*ones(ng,1);
% fh = An.*uh + kappa*reshape(mapContractK(nl,q,[],2,1,3,2,1),[ng nch]) + tau.*(u-uh);

    u = UDG(:,1,i);
    imap = mesh.elcon(:,i);
    uh = UH(1,imap)';
    divdg(:,1,i) = Mass\(-(Bx*c(1)+By*c(2))*u + (c(1)*Cx+c(2)*Cy)*uh);    
    divdg(:,1,i) = divdg(:,1,i) + Mass\(EI*u -EJ*uh);
    
    divdu(:,1,i) = Mass\(c(1)*(Bx'*u) + c(2)*(By'*u));
         
    %divdg(:,1,i) = divdu(:,1,i) + Mass\(-(Dx*c(1)+Dy*c(2))*u + (c(1)*Cx+c(2)*Cy)*uh);
    
    eu = divpu-divdg(:,1,i);
    err(1) = err(1) + eu'*Mass*eu;
    eu = divpu-divdu(:,1,i);
    err(2) = err(2) + eu'*Mass*eu;    
end
err = sqrt(err);





