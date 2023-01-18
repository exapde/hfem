function [QDG,UH] = get_quhat(master, mesh, app, UDG, time, indsol)
%
%      MASTER:       Master structure
%      MESH:         Mesh structure

kappa = app.kappa;
tau = app.tau;

nt  = size(mesh.t,1);
%nf  = size(mesh.f,1);
npv = size(mesh.dgnodes,1);
nps = size(master.plocfc,1);

%shap = squeeze(master.shap(:,1,:));
shapxi = squeeze(master.shapvl(:,:,2))*diag(master.gwvl);
shapet = squeeze(master.shapvl(:,:,3))*diag(master.gwvl); 
%sh1d = squeeze(master.sh1d(:,1,:));
perm = master.perm;

% allocate memory for the elemental matrices
A_K = zeros(2*npv,2*npv);
B_K = zeros(2*npv,npv);

% allocate memory 
P = zeros(2*npv,3*nps,nt);
L = zeros(2*npv,nt);
H = zeros(3*nps,3*nps,nt);
R = zeros(3*nps,nt);

shap = permute(master.shapvl,[1 2 3]);
sh1d = permute(master.shapfc,[1 2 3]);
for i = 1:nt % loop over each element
    xxi = (shap(:,:,2))'*(mesh.dgnodes(:,1,i));
    xet = (shap(:,:,3))'*(mesh.dgnodes(:,1,i));
    yxi = (shap(:,:,2))'*(mesh.dgnodes(:,2,i));
    yet = (shap(:,:,3))'*(mesh.dgnodes(:,2,i));
    jac = xxi.*yet - xet.*yxi;
    shapx =   shapxi*diag(yet) - shapet*diag(yxi);
    shapy = - shapxi*diag(xet) + shapet*diag(xxi);
    Mass  = shap(:,:,1)*diag(master.gwvl.*jac)*shap(:,:,1)';
        
    % form A_K
    A_K(1:npv,1:npv) = (1/kappa)*Mass;
    A_K((npv+1):2*npv,(npv+1):2*npv) = (1/kappa)*Mass;
    
    % form B_K
    B_K(1:npv,1:npv) = shapx*shap(:,:,1)';
    B_K((npv+1):2*npv,1:npv) = shapy*shap(:,:,1)';    
        
    % allocate memory for the elemental matrices
    C_K = zeros(2*npv,3*nps);    
    E_K = zeros(npv,3*nps);
    M_K = zeros(3*nps,3*nps);    
    G_K = zeros(3*nps,1);    
    N_K = zeros(3*nps,2*npv);
    P_K = zeros(3*nps,npv);
    for j = 1:3 % loop over each face of the element i
        I = perm(:,j,1);
        J = ((j-1)*nps+1):j*nps;
                
        xxi = (sh1d(:,:,2))'*(mesh.dgnodes(I,1,i)); 
        yxi = (sh1d(:,:,2))'*(mesh.dgnodes(I,2,i));  
        dsdxi = sqrt(xxi.^2+yxi.^2);
        nl = [yxi./dsdxi,-xxi./dsdxi];
        dws = master.gwfc.*dsdxi;                
                        
        tmp = sh1d(:,:,1)*diag(tau*dws)*sh1d(:,:,1)';
        tmpx = sh1d(:,:,1)*diag(dws.*nl(:,1))*sh1d(:,:,1)';
        tmpy = sh1d(:,:,1)*diag(dws.*nl(:,2))*sh1d(:,:,1)';
        
        % form C_K
        C_K(I,J) = C_K(I,J) + tmpx;
        C_K(npv+I,J) = C_K(npv+I,J) + tmpy;
        
        % form E_K
        E_K(I,J) = E_K(I,J) + tmp;

        % form M_K
        M_K(J,J) = M_K(J,J) + tmp;     
        
        bf = mesh.f(abs(mesh.t2f(i,j)),end);
        if bf<0                        
            pg = (sh1d(:,:,1))'*(mesh.dgnodes(I,:,i));              
            fg = exactsol(pg, [app.arg{1} app.arg{2}], time);
            G_K(J,1) = G_K(J,1) + tau*sh1d(:,:,1)*(dws.*fg(:,indsol));                                        
        else
            % form N_K
            N_K(J,I) = N_K(J,I) + tmpx;
            N_K(J,npv+I) = N_K(J,npv+I) + tmpy;
            
            % form P_K
            P_K(J,I) = P_K(J,I) + tmp;
        end                        
    end        
    
    % these are needed to compute u and q
    P(:,:,i) = A_K\C_K;
    L(:,i) = A_K\(-B_K*UDG(:,1,i));

    % form elemental stiffness matrix 
    H_K = M_K + N_K*P(:,:,i);
            
    % form elemental residual vector 
    R_K = G_K + P_K*UDG(:,1,i) - N_K*L(:,i);    
        
    H(:,:,i) = H_K;
    R(:,i)   = R_K; 
end
        
% compute the index mapping
elcon = mesh.elcon;
il = zeros(3*nps,3*nps,nt);
jl = zeros(3*nps,3*nps,nt);
for i=1:nt
    con = repmat((elcon(:,i)'-1),1,1)+ones(1,3*nps);
    con = reshape(con,3*nps,1);
    il(:,:,i) = repmat(con ,1,3*nps);
    jl(:,:,i) = repmat(con',3*nps,1);        
end

H = sparse(reshape(il,(nps*3)^2*nt,1),reshape(jl,(nps*3)^2*nt,1),reshape(H,(nps*3)^2*nt,1));
R = sparse(reshape(il(:,1,:),(nps*3)*nt,1),ones((nps*3)*nt,1),reshape(R,(nps*3)*nt,1));                    

% solve for vhat
UH = H\R;

% compute q
QDG = zeros(npv,2,nt);
for i = 1:nt
    imap = elcon(:,i);
    QV = L(:,i) + P(:,:,i)*UH(imap);
    QDG(:,:,i) = reshape(QV(1:2*npv),[npv 2]);    
end
