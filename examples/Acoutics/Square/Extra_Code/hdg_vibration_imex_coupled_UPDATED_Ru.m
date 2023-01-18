function UDG = hdg_vibration_imex_coupled_UPDATED(master,mesh,app,VDGn,Un,tn,dt,nstage,torder,Minv)

% dirk coefficients
[~,b,c,D] = dirkcoeff(nstage,torder);

% Finds implicit elements
ind = find(mesh.imex == 1); % Implicit element indicies

% Finds imex boundary explicit elements
E_ind = find(mesh.bound_ele(:,1) > 0);
E_ele = mesh.bound_ele(E_ind,1);
expE = length(E_ind);

ind = sort(cat(1,ind,E_ele));
% ind2 = setdiff(ind,mesh.bound_ele(:,3)); % Implicit elements not on the imex boundary
impE = length(ind); % Number of implicit elements

[npv,nc,ne] = size(VDGn);
ncu = 1;

UDGall = zeros(npv,nc+ncu,ne,nstage);
UDGall(:,:,:,1) = cat(2, VDGn, Un); % add on column for explicit calc
UDG = cat(2, VDGn, Un);
% UDGni = UDG;
% Uall = zeros(npv,ncu,ne,nstage);

% Define time
 app.dt = dt;   
 app.time = tn;  

for istage = 1:nstage % for each dirk stage    
    fprintf('DIRK stage :  %d\n', istage);
   
    % intermidiate time
    tni = tn+dt*c(istage);
    
    % Explicit Code
    UDGni = ehdgrk_imex_coupled(master,mesh,app,Minv,UDGall(:,:,:,istage));
%     UDGni = ehdgrk_imex_Ru_newform(master,mesh,app,Minv,UDGall(:,:,:,istage));
    % Break up code for implicit solve
    VDGni = UDGni(:,1:3,:);
    Uni = UDGni(:,4,:);
    
    % update the timestep factor
    fc = D(istage,istage)/(dt);        
    
    % update the source term
    SDG = fc*VDGn(:,:,ind);
    for jstage = 1:istage-1
        SDG = SDG - (D(istage,jstage)/dt)*(UDGall(:,1:nc,ind,jstage)-VDGn(:,:,ind));
    end    
    
    % compute (v,q,vhat)
    [VDGni,VHni] = hdg_solve(master,mesh,app,SDG,fc,tni, VDGni, Minv, Uni);
    
    % update the intermediate states for (v, q)
    UDGall(:,1:nc,ind,istage) = VDGni;  % uses ind to place implicit sol. in global datastructure
    
    SH = UDGall(:,1,ind,istage) + fc*Un(:,:,ind);
    for jstage = 1:istage-1
        SH = SH - (D(istage,jstage)/dt)*(UDGall(:,end,ind,jstage)-Un(:,:,ind));
    end        
    
    % update intermediate states for u
    Uni = (dt/D(istage,istage))*SH;
    UDGall(:,end,ind,istage) = Uni;
    
end

% Recover explicit solution
[UDG] = ehdgrk_imex_coupled(master,mesh,app,Minv,UDG);

% VHnew = VHni; % Test to assign

VDGnew=VDGn;
Unew = Un;

% VDGnew=VDGni;
% Unew = Uni;

% Unew=UDG(:,end,:);

for i=1:nstage
   for j=1:nstage     
      VDGnew(:,:,ind)=VDGnew(:,:,ind)+b(i)*D(i,j).*((UDGall(:,1:nc,ind,j)-VDGn(:,:,ind)));
      Unew(:,:,ind)=Unew(:,:,ind)+b(i)*D(i,j).*((UDGall(:,end,ind,j)-Un(:,:,ind)));
   end
end

UDG(:,1:nc,ind) = VDGnew(:,:,ind);
UDG(:,end,ind) = Unew(:,:,ind);
a = 2;



function [VDG,VH] = hdg_solve(master, mesh, app, SDG, fc, time, UDG, Minv, Uni)
%
%      MASTER:       Master structure
%      MESH:         Mesh structure

kappa = app.kappa;
tau = app.tau;
kx  = app.k(1);
ky  = app.k(2);
kk  = sqrt(kx^2+ky^2);

% Finds implicit elements
ind = find(mesh.imex == 1); % Implicit element indicies
impE = length(ind); % Number of implicit elements

% Finds implicit/explicit boundary
indf = find(mesh.fimex == 1); % Implicit element indicies
faceimex = length(indf); % Number of implicit elements

% Finds imex boundary explicit elements
E_ind = find(mesh.bound_ele(:,1) > 0);
E_ele = mesh.bound_ele(E_ind,1);
expE = length(E_ind);

% nt = impE;

nt  = size(mesh.t,1);
nf  = size(mesh.f,1);
npv = size(mesh.dgnodes,1);
nps = size(master.plocfc,1);

%shap = squeeze(master.shap(:,1,:));
master.shap = permute(master.shapvl,[1 3 2]);
shapxi = squeeze(master.shap(:,2,:))*diag(master.gwvl);
shapet = squeeze(master.shap(:,3,:))*diag(master.gwvl); 
%sh1d = squeeze(master.sh1d(:,1,:));
perm = master.perm;

% allocate memory for the elemental matrices
A_K = zeros(2*npv,2*npv);
B_K = zeros(2*npv,npv);

% Different from uncoupled code
shap = permute(master.shap,[1 3 2]);
% sh1d = permute(master.sh1d,[1 3 2]);
sh1d = master.shapfc;


% allocate memory 
P = zeros(3*npv,3*nps,impE + expE);
L = zeros(3*npv,impE + expE);
H = zeros(3*nps,3*nps,impE + expE);
R = zeros(3*nps,impE + expE);

for count = 1:impE
    i = ind(count); % loop over each implicit element
    
    xxi = (shap(:,:,2))'*(mesh.dgnodes(:,1,i));
    xet = (shap(:,:,3))'*(mesh.dgnodes(:,1,i));
    yxi = (shap(:,:,2))'*(mesh.dgnodes(:,2,i));
    yet = (shap(:,:,3))'*(mesh.dgnodes(:,2,i));
    jac = xxi.*yet - xet.*yxi;
    shapx =   shapxi*diag(yet) - shapet*diag(yxi);
    shapy = - shapxi*diag(xet) + shapet*diag(xxi);
    Mass  = shap(:,:,1)*diag(master.gwvl.*jac)*shap(:,:,1)';
        
    % form A_K
    A_K(1:npv,1:npv) = (fc/kappa)*Mass;
    A_K((npv+1):2*npv,(npv+1):2*npv) = (fc/kappa)*Mass;
    
    % form B_K
    B_K(1:npv,1:npv) = shapx*shap(:,:,1)';
    B_K((npv+1):2*npv,1:npv) = shapy*shap(:,:,1)';    
    
    % form F_K
    xg = (shap(:,:,1))'*mesh.dgnodes(:,:,i);
    fg = 0*xg(:,1);
    F_K = shap(:,:,1)*(master.gwvl.*jac.*fg);
    
    % allocate memory for the elemental matrices
    C_K = zeros(2*npv,3*nps);
    D_K = fc*Mass;
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
                
        % form C_K
        C_K(I,J) = C_K(I,J) + sh1d(:,:,1)*diag(dws.*nl(:,1))*sh1d(:,:,1)';
        C_K(npv+I,J) = C_K(npv+I,J) + sh1d(:,:,1)*diag(dws.*nl(:,2))*sh1d(:,:,1)';

        tmp = sh1d(:,:,1)*diag(tau*dws)*sh1d(:,:,1)';
        
        % form D_K
        D_K(I,I) = D_K(I,I) + tmp;

        % form E_K
        E_K(I,J) = E_K(I,J) + tmp;

        % form M_K
        M_K(J,J) = M_K(J,J) + tmp;     
        
%         % absorbing condition
%         if mesh.f(abs(mesh.t2f(i,j)),end)==-1 || mesh.f(abs(mesh.t2f(i,j)),end)==-2 || mesh.f(abs(mesh.t2f(i,j)),end)==-3 || mesh.f(abs(mesh.t2f(i,j)),end)==-4
%             M_K(J,J) = M_K(J,J) + tmp/tau;
%         end
        
%         % reflection condition
%         if mesh.f(abs(mesh.t2f(i,j)),end)==-5
%             pg = (sh1d(:,:,1))'*(mesh.dgnodes(I,:,i));             
%             u1 = (kx)*cos(kx*pg(:,1)+ky*pg(:,2) - kk*time); 
%             u2 = (ky)*cos(kx*pg(:,1)+ky*pg(:,2) - kk*time); 
%             gn = u1.*nl(:,1) + u2.*nl(:,2);            
%             G_K(J,1) = G_K(J,1) - sh1d(:,:,1)*(dws.*gn);
%         end

        if mesh.f(abs(mesh.t2f(i,j)),end)==-1 || mesh.f(abs(mesh.t2f(i,j)),end)==-2 || mesh.f(abs(mesh.t2f(i,j)),end)==-3 || mesh.f(abs(mesh.t2f(i,j)),end)==-4                    
            pg = (sh1d(:,:,1))'*(mesh.dgnodes(I,:,i));              
            fg = exactsol(pg, [app.arg{1} app.arg{2}], time);
            G_K(J,1) = G_K(J,1) + tau*sh1d(:,:,1)*(dws.*fg(:,1));                                        
        else
            % form N_K
            N_K(J,I) = N_K(J,I) + sh1d(:,:,1)*diag(dws.*nl(:,1))*sh1d(:,:,1)';
            N_K(J,npv+I) = N_K(J,npv+I) + sh1d(:,:,1)*diag(dws.*nl(:,2))*sh1d(:,:,1)';
            
            % form P_K
            P_K(J,I) = P_K(J,I) + tmp;
        end      
        
%         % IMEX boundary
%         if mesh.fimex(abs(mesh.t2f(i,j)))== 1
%             M_K(J,J) = M_K(J,J) + 2*tmp;
%             en = mesh.t2t(i,j); % neighboring explicit element             
%             ug = (sh1d(:,:,1))'*(UDG(I,:,en)); % solution on the explicit element                   
%             gn = ug(:,2).*nl(:,1) + ug(:,3).*nl(:,2) + tau*ug(:,1); % q_E dot n + tau_v_E            
%             G_K(J,1) = G_K(J,1) + sh1d(:,:,1)*(dws.*gn);            
%         end
        
    end        
    
    % these are needed to compute u and q
    P(:,:,count) = [A_K B_K; -B_K' D_K]\[C_K; E_K];
    L(:,count) = [A_K B_K; -B_K' D_K]\[Mass*SDG(:,2,count)/kappa; Mass*SDG(:,3,count)/kappa; F_K + Mass*SDG(:,1,count)];
   
    % form elemental stiffness matrix 
    H_K = M_K + [N_K -P_K]*P(:,:,count);
            
    % form elemental residual vector 
    R_K = G_K-[N_K -P_K]*L(:,count);    
    
    H(:,:,count) = H_K;
    R(:,count)   = R_K; 
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP OVER BOUNDARY EXPLICIT ELEMENTS


% Find adjacent elements
E_neighbors = mesh.t2t(E_ele,:);

for count = 1:expE
    i = E_ele(count); % loop over each explicit boundary element
    [UDG,B,G] = ehdgrk_imex_Ru_newform(master,mesh,app,Minv,cat(2,UDG,Uni));
%     [B,G] = ehdgrk_imex_coupled_UPDATED(master,mesh,app,Minv,UDG, i,3);
%     if E_neighbors(count,:) > 0
%         imp_exp = mesh.imex(E_neighbors(count,:)); % determines if neighboring element is implicit or explicit
    
    xxi = (shap(:,:,2))'*(mesh.dgnodes(:,1,i));
    xet = (shap(:,:,3))'*(mesh.dgnodes(:,1,i));
    yxi = (shap(:,:,2))'*(mesh.dgnodes(:,2,i));
    yet = (shap(:,:,3))'*(mesh.dgnodes(:,2,i));
    jac = xxi.*yet - xet.*yxi;
    shapx =   shapxi*diag(yet) - shapet*diag(yxi);
    shapy = - shapxi*diag(xet) + shapet*diag(xxi);
    Mass  = shap(:,:,1)*diag(master.gwvl.*jac)*shap(:,:,1)';
    
    % form A_K
    A_K(1:npv,1:npv) = (fc/kappa)*Mass;
    A_K((npv+1):2*npv,(npv+1):2*npv) = (fc/kappa)*Mass;
    
    % form B_K
    B_K(1:npv,1:npv) = shapx*shap(:,:,1)';
    B_K((npv+1):2*npv,1:npv) = shapy*shap(:,:,1)';
    
    % form F_K
    xg = (shap(:,:,1))'*mesh.dgnodes(:,:,i);
    fg = 0*xg(:,1);
    F_K = shap(:,:,1)*(master.gwvl.*jac.*fg);
    
    % allocate memory for the elemental matrices
    C_K = zeros(2*npv,3*nps);
    D_K = fc*Mass;
    E_K = zeros(npv,3*nps);
    M_K = zeros(3*nps,3*nps);
    G_K = zeros(3*nps,1);
    N_K = zeros(3*nps,2*npv);
    P_K = zeros(3*nps,npv);
    for j = 1:3 % loop over each face of the element i
        
        if E_neighbors(count,j) > 0
            imp_exp = mesh.imex(E_neighbors(count,j)); % determines if neighboring element is implicit or explicit
        else
            imp_exp = -1;  % account for element on domain boundary
        end
        
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
        
        tmp = sh1d(:,:,1)*diag(tau*dws)*sh1d(:,:,1)';
        
        % form D_K
        D_K(I,I) = D_K(I,I) + tmp;
        
        % form E_K
        E_K(I,J) = E_K(I,J) + tmp;
        
        % form M_K
        M_K(J,J) = M_K(J,J) + tmp;
        
        
%         if imp_exp == 0   % Condition to determine if apply explicit solution
%             neighbor_local_face = find(mesh.t2t(E_neighbors(count,j),:) == i);
%             pg = (sh1d(:,:,1))'*(mesh.dgnodes(I,:,i));
%             fg = UDG(I,neighbor_local_face,E_neighbors(count,j)); % explicit solution
%             % fg is obtained by using mkfhat
%             G_K(J,1) = G_K(J,1) + tau*sh1d(:,:,1)*(dws.*fg(:,1));            
%         if mesh.f(abs(mesh.t2f(i,j)),end)==-1 || mesh.f(abs(mesh.t2f(i,j)),end)==-2 || mesh.f(abs(mesh.t2f(i,j)),end)==-3 || mesh.f(abs(mesh.t2f(i,j)),end)==-4                    
%            pg = (sh1d(:,:,1))'*(mesh.dgnodes(I,:,i));              
%            fg = exactsol(pg, [app.arg{1} app.arg{2}], time);
%            G_K(J,1) = G_K(J,1) + tau*sh1d(:,:,1)*(dws.*fg(:,1));   
%         else
%             % form N_K
%             N_K(J,I) = N_K(J,I) + sh1d(:,:,1)*diag(dws.*nl(:,1))*sh1d(:,:,1)';
%             N_K(J,npv+I) = N_K(J,npv+I) + sh1d(:,:,1)*diag(dws.*nl(:,2))*sh1d(:,:,1)';
%             
%             % form P_K
% %             if imp_exp ~= 0
% %             P_K(J,I) = P_K(J,I) + tmp;
% %             end
%         end
        
        % IMEX boundary
        if mesh.fimex(abs(mesh.t2f(i,j)))== 1
            %M_K(J,J) = M_K(J,J) + 2*tmp;
            en = mesh.t2t(i,j); % neighboring explicit element
            ug = (sh1d(:,:,1))'*(UDG(I,:,en)); % solution on the explicit element
            gn = ug(:,2).*nl(:,1) + ug(:,3).*nl(:,2) + tau*ug(:,1); % q_E dot n + tau_v_E
            G_K(J,1) = G_K(J,1) + sh1d(:,:,1)*(dws.*gn);
        end  
    end
    
    count2 = count + impE;
    
    % these are needed to compute u and q
%     P(:,:,count2) = [A_K B_K; -B_K' D_K]\[C_K; E_K];
%     L(:,count2) = [A_K B_K; -B_K' D_K]\[Mass*SDG(:,2,count2)/kappa; Mass*SDG(:,3,count2)/kappa; F_K + Mass*SDG(:,1,count2)];
    Minv_Exp = Minv(:,:,i);
    ZeroBlk = zeros(size(Minv_Exp));
    Minv_Blk = [Minv_Exp, ZeroBlk, ZeroBlk; ZeroBlk, Minv_Exp, ZeroBlk; ZeroBlk, ZeroBlk, Minv_Exp];
    P(:,:,count2) = -Minv_Blk*B;
    L(:,count2) = Minv_Blk*G;
   
    % form elemental stiffness matrix 
    H_K = M_K + [N_K -P_K]*P(:,:,count2);
            
    % form elemental residual vector 
    R_K = G_K-[N_K -P_K]*L(:,count2);    
    
    H(:,:,count2) = H_K;
    R(:,count2)   = R_K;
end
    


impE = impE + expE;
ind = sort(cat(1,ind,E_ele));

% compute the index mapping
% elcon = elconnectivities_orig(nps,mesh.t2f);
% change to compute over only implicit elementsfaceimex
[elcon,ndof,edg,hybridn] = elconnectivities_imex_UPDATED(mesh,nf,E_ele);
il = zeros(3*nps,3*nps,impE);
jl = zeros(3*nps,3*nps,impE);
for count = 1:impE
    i = ind(count); % loop over each implicit element
    
%     con = repmat((elcon(:,i)'-1),1,1)+ones(1,3*nps);
    con = repmat((elcon(:,count)'-1),1,1)+ones(1,3*nps);
    con = reshape(con,3*nps,1);
    il(:,:,count) = repmat(con ,1,3*nps);
    jl(:,:,count) = repmat(con',3*nps,1);        
end

H = sparse(reshape(il,(nps*3)^2*impE,1),reshape(jl,(nps*3)^2*impE,1),reshape(H,(nps*3)^2*impE,1));
R = sparse(reshape(il(:,1,:),(nps*3)*impE,1),ones((nps*3)*impE,1),reshape(R,(nps*3)*impE,1));                    

% solve for vhat
VH = H\R;



% compute v and q
VDG = zeros(npv,3,impE);
for count = 1:impE
    i = ind(count); 
    imap = elcon(:,count);
    QV = L(:,count) + P(:,:,count)*VH(imap);
    VDG(:,2:3,count) = reshape(QV(1:2*npv),[npv 2]);
    VDG(:,1,count) = QV(2*npv+1:end);    
end

function [A,b,c,D] = dirkcoeff(q,p)
% q - number of stages
% p - order of accuracy

if q == 1 && p == 1
    A = 1;
    b = 1;
    c = 1;
elseif q == 1 && p == 2 % implcit midpoint   
    A = 0.5;
    b = 1;
    c = 0.5;    
elseif q == 2 && p == 2
    A = [1-0.5*sqrt(2), 0;
           0.5*sqrt(2), 1-0.5*sqrt(2)];
    b = [  0.5*sqrt(2), 1-0.5*sqrt(2)];
    c = [1-0.5*sqrt(2), 1];
elseif q == 2 && p == 3    
    A = [0.5+0.5/sqrt(3), 0;
              -1/sqrt(3), 0.5 + 0.5/sqrt(3)];
    b = [            0.5, 0.5];
    c = [0.5+0.5/sqrt(3), 0.5-0.5/sqrt(3)];    
elseif q == 3 && p == 3
    a1 = 0.4358665215;
    t1 = (1+a1)/2;
    b1 = -(6*a1^2-16*a1+1)/4;
    b2 = (6*a1^2-20*a1+5)/4;    
    A = [a1   ,  0, 0;
         t1-a1, a1, 0;
         b1   , b2, a1];
    b = [b1, b2, a1];
    c = [a1, t1, 1];    
elseif q == 3 && p == 4    
    a1 = 2*cos(pi/18)/sqrt(3);
    A = [ 0.5*(1+a1),            0,          0;
             -0.5*a1,   0.5*(1+a1),          0;
                1+a1,    -(1+2*a1), 0.5*(1+a1)];
    b = [ 1/(6*a1^2), 1-1/(3*a1^2), 1/(6*a1^2)];
    c = [ 0.5*(1+a1),          0.5, 0.5*(1-a1)];    
else
    error('Invalid (q,p) combination');
end

D = inv(A);


