function [UDGnew,UHnew,Vnew,VHnew] = hdg_wave_ec(master,mesh,app,UDGn,Vn,UHn,VHn,tn,dt,nstage,torder)

% dirk coefficients
[~,b,c,D] = dirkcoeff2(nstage,torder);

[npv,nc,ne] = size(UDGn);
ncu = 1;

% intermediate solutions
UDGall = zeros(npv,nc,ne,nstage);
Vall = zeros(npv,ncu,ne,nstage);
UHall = zeros((mesh.porder+1)*mesh.nf,ncu,nstage);
VHall = UHall;

for istage = 1:nstage % for each dirk stage    
    fprintf('DIRK stage :  %d\n', istage);
    
    % intermidiate time
    tni = tn+dt*c(istage);
    
    % update the timestep factor
    fc = D(istage,istage)/(dt);            
        
    % update the source term
    SDG = fc*(fc*UDGn(:,1,:)+Vn);
    for jstage = 1:istage-1
        SDG = SDG - (D(istage,jstage)/dt)*(Vall(:,:,:,jstage)-Vn);
        SDG = SDG - fc*(D(istage,jstage)/dt)*(UDGall(:,1,:,jstage)-UDGn(:,1,:));
    end    
    
    % compute (u,q,uhat)
    [UDGni,UHni] = hdg_solve(master,mesh,app,SDG,fc,tni);
    
    % update the intermediate states for (u, q)
    UDGall(:,:,:,istage) = UDGni;    
    UHall(:,:,istage) = UHni;        
    
    % update intermediate states for v
    Vall(:,:,:,istage) = fc*(UDGall(:,1,:,istage) - UDGn(:,1,:));
    VHall(:,:,istage) = fc*(UHall(:,:,istage) - UHn(:,1));
    for jstage = 1:istage-1
        Vall(:,:,:,istage) = Vall(:,:,:,istage) + (D(istage,jstage)/dt)*(UDGall(:,1,:,jstage)-UDGn(:,1,:));
        VHall(:,:,istage) = VHall(:,:,istage) + (D(istage,jstage)/dt)*(UHall(:,:,jstage)-UHn);
    end            
end

% Solution at the new timestep n+1
% UHnew = UHni;
% UDGnew = UDGall(:,:,:,end);
% Vnew = Vall(:,:,:,end);
% VHnew = VHall(:,:,end);

UHnew=UHn;
UDGnew=UDGn;
Vnew=Vn;
VHnew=VHn;
for i=1:nstage
   for j=1:nstage     
      UDGnew=UDGnew+b(i)*D(i,j).*(UDGall(:,:,:,j)-UDGn(:,:,:));
      UHnew=UHnew+b(i)*D(i,j).*(UHall(:,:,j)-UHn);
      Vnew=Vnew+b(i)*D(i,j).*(Vall(:,1,:,j)-Vn);
      VHnew=VHnew+b(i)*D(i,j).*(VHall(:,:,j)-VHn);
   end
end

function [UDG,UH] = hdg_solve(master, mesh, app, SDG, fc, time)
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
P = zeros(3*npv,3*nps,nt);
L = zeros(3*npv,nt);
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
    
    % form F_K
    xg = (shap(:,:,1))'*mesh.dgnodes(:,:,i);
    fg = 0*xg(:,1);
    F_K = shap(:,:,1)*(master.gwvl.*jac.*fg);
    
    % allocate memory for the elemental matrices
    C_K = zeros(2*npv,3*nps);
    D_K = (fc*fc)*Mass;
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
        
        % form D_K
        D_K(I,I) = D_K(I,I) + tmp;

        % form E_K
        E_K(I,J) = E_K(I,J) + tmp;

        % form M_K
        M_K(J,J) = M_K(J,J) + tmp;     
        
        bf = mesh.f(abs(mesh.t2f(i,j)),end);
        if bf<0                        
            pg = (sh1d(:,:,1))'*(mesh.dgnodes(I,:,i));              
            fg = exactsol(pg, [app.arg{1} app.arg{2}], time);
            G_K(J,1) = G_K(J,1) + tau*sh1d(:,:,1)*(dws.*fg(:,4));                                        
        else
            % form N_K
            N_K(J,I) = N_K(J,I) + tmpx;
            N_K(J,npv+I) = N_K(J,npv+I) + tmpy;
            
            % form P_K
            P_K(J,I) = P_K(J,I) + tmp;
        end                        
    end        
    
    % these are needed to compute u and q
    P(:,:,i) = [A_K B_K; -B_K' D_K]\[C_K; E_K];
    L(:,i) = [A_K B_K; -B_K' D_K]\[0*F_K; 0*F_K; F_K + Mass*SDG(:,1,i)];
   
    % form elemental stiffness matrix 
    H_K = M_K + [N_K -P_K]*P(:,:,i);
            
    % form elemental residual vector 
    R_K = G_K-[N_K -P_K]*L(:,i);    
    
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

% compute v and q
UDG = zeros(npv,3,nt);
for i = 1:nt
    imap = elcon(:,i);
    QV = L(:,i) + P(:,:,i)*UH(imap);
    UDG(:,2:3,i) = reshape(QV(1:2*npv),[npv 2]);
    UDG(:,1,i) = QV(2*npv+1:end);    
end

function [A,b,c,D] = dirkcoeff(q,p)
% q - number of stages
% p - order of accuracy

if q == 1 && p == 1
    A = 1;
    b = 1;
    c = 1;
elseif q == 1 && p == 2
    % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
    A = 0.5;
    b = 1;
    c = 0.5;    
elseif q == 2 && p == 2
    A = [1-0.5*sqrt(2), 0;
           0.5*sqrt(2), 1-0.5*sqrt(2)];
    b = [  0.5*sqrt(2), 1-0.5*sqrt(2)];
    c = [1-0.5*sqrt(2), 1];
elseif q == 2 && p == 3
    % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
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
    % fprintf('Warning: Innacutate for non-wave formuations - use (2,2) or (3,3)');
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


function [A,b,c,D] = dirkcoeff2(q,p)
% q - number of stages
% p - order of accuracy
% symplectic dirk schemes
if  q == 1 && p == 2 % implicit midpoint  
    A = 0.5;
    b = 1;
    c = 0.5;   
elseif q == 2 && p == 2
    A = [1/4, 0;
           0.5, 1/4];
    b = [  0.5, 0.5];
    c = [1/4, 3/4];  
elseif q == 3 && p == 3
    a = 1.351207191959658; 
    A = [1/2*a   ,  0, 0;
         a, 1/2*a, 0;
         a   , a, 1/2-a];
    b = [a, a, 1-2*a];
    c = [a/2, 3/2*a, 1/2+a];   
elseif q == 3 && p == 4   
    a = 1.351207191959658;
    A = [ 0.5*a,            0,          0;
             a,   0.5-a,          0;
                a,    1-2*a, 0.5*a];
    b = [ a, 1-2*a, a];
    c = [ 0.5*a,          0.5, 1-1/2*a];
elseif q == 4 && p == 4
    b1 = -2.70309412;
    b2 = -0.53652708;
    b3 = 2.37893931;
    b4 = 1.8606818856;
    A = [b1/2, 0, 0, 0;
         b1, b2/2, 0, 0;
         b1, b2, b3/2, 0;
         b1, b2, b3, b4/2];
     b = [b1, b2, b3, b4];
     c = [b1/2, b1+b2/2, b1+b2+b3/2, b1+b2+b3+b4/2];
elseif q == 6 && p == 5
    b1 = 0.5080048194000274;
    b2 = 1.360107162294827;
    b3 = 2.0192933591817224;
    b4 = 0.5685658926458251;
    b5 = -1.4598520495864393;
    b6 = -1.9961191839359627;
    A = [b1/2, 0, 0, 0, 0, 0;
         b1, b2/2, 0, 0, 0, 0;
         b1, b2, b3/2, 0, 0, 0;
         b1, b2, b3, b4/2, 0, 0;
         b1, b2, b3, b4, b5/2, 0;
         b1, b2, b3, b4, b5, b6/2];
     b = [b1, b2, b3, b4, b5, b6];
     c = [b1/2, b1+b2/2, b1+b2+b3/2, b1+b2+b3+b4/2, b1+b2+b3+b4+b5/2, b1+b2+b3+b4+b5+b6/2];
elseif q == 7 && p == 5
    b1 = -1.9119528632302611;
    b2 = 1.7330487520694027;
    b3 = 1.7978839093830516;
    b4 = 1.2522794817510428;
    b5 = 2.87965671639611;
    b6 = -2.8389631331390848;
    b7=b1;
    A = [b1/2, 0, 0, 0, 0, 0, 0;
         b1, b2/2, 0, 0, 0, 0, 0;
         b1, b2, b3/2, 0, 0, 0, 0;
         b1, b2, b3, b4/2, 0, 0, 0;
         b1, b2, b3, b4, b5/2, 0, 0;
         b1, b2, b3, b4, b5, b6/2, 0;
         b1, b2, b3, b4, b5, b6, b7/2];
     b = [b1, b2, b3, b4, b5, b6, b7];
     c = [b1/2, b1+b2/2, b1+b2+b3/2, b1+b2+b3+b4/2, b1+b2+b3+b4+b5/2, b1+b2+b3+b4+b5+b6/2, b1+b2+b3+b4+b5+b6+b7/2];
elseif q == 7 && p == 6
    b1 = 7.8451361047755652e-01;
    b2 = 2.3557321335935860e-01;
    b3 = -1.1776799841788705e+00;
    b4 = 1.3151863206839107e+00;
    b5 = b3;
    b6 = b2;
    b7 = b1;
    A = [b1/2, 0, 0, 0, 0, 0, 0;
         b1, b2/2, 0, 0, 0, 0, 0;
         b1, b2, b3/2, 0, 0, 0, 0;
         b1, b2, b3, b4/2, 0, 0, 0;
         b1, b2, b3, b4, b5/2, 0, 0;
         b1, b2, b3, b4, b5, b6/2, 0;
         b1, b2, b3, b4, b5, b6, b7/2];
     b = [b1, b2, b3, b4, b5, b6, b7];
     c = [b1/2, b1+b2/2, b1+b2+b3/2, b1+b2+b3+b4/2, b1+b2+b3+b4+b5/2, b1+b2+b3+b4+b5+b6/2, b1+b2+b3+b4+b5+b6+b7/2];
elseif q == 9 && p == 6
    b1 = 2.44398640327406;
    b2 = -2.46929010453909;
    b3 = 0.28158632623993;
    b4 = 0.50745789725108;
    b5 = 1.17888214306555;
    b6 = -2.31558614555863;
    b7 = 2.35136242638295;
    b8 = -1.24653876689005;
    b9 = 0.26813982077420;
    A = [b1/2, 0, 0, 0, 0, 0, 0, 0, 0;
         b1, b2/2, 0, 0, 0, 0, 0, 0, 0;
         b1, b2, b3/2, 0, 0, 0, 0, 0, 0;
         b1, b2, b3, b4/2, 0, 0, 0, 0, 0;
         b1, b2, b3, b4, b5/2, 0, 0, 0, 0;
         b1, b2, b3, b4, b5, b6/2, 0, 0, 0;
         b1, b2, b3, b4, b5, b6, b7/2, 0, 0;
         b1, b2, b3, b4, b5, b6, b7, b8/2, 0;
         b1, b2, b3, b4, b5, b6, b7, b8, b9/2];
     b = [b1, b2, b3, b4, b5, b6, b7, b8, b9];
     c = [b1/2, b1+b2/2, b1+b2+b3/2, b1+b2+b3+b4/2, b1+b2+b3+b4+b5/2, b1+b2+b3+b4+b5+b6/2, b1+b2+b3+b4+b5+b6+b7/2, ...
         b1+b2+b3+b4+b5+b6+b7+b8/2, b1+b2+b3+b4+b5+b6+b7+b8+b9/2];
else
    error('Invalid (q,p) combination');
end

D = inv(A);
