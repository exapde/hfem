function [UHATo,B,G,D] = ehdgrk_imex_coupled_UPDATED_v2(master,mesh,app,Minv,UDG,i,j)
% i: is explicit element
% j: is the stage

tau = app.tau;
dt = app.dt;
nstage = app.nstage;
torder = app.torder;
D = ERKcoeff(nstage,torder);

[npv,ncu,ne] = size(UDG); 
% UDGS = zeros(npv,ncu+1,nstage+1);
UDGS = zeros(npv,ncu+1,1);
% UDGS(:,:,1) = UDG(:,:,i); % element i

iu = 1:3;  %%% Temporarily hard coding in 3 instead of ncu-1
% iu = 1:(ncu);
UDGS(:,iu,1) = UDG(:,iu,i); % element i %%%%% CHANGED TO UDG(:,IU,I)

ind = find(mesh.imex == 0); % Explicit element indicies
expE = length(ind); % Number of explicit elements

dtj = dt/D(j,j);

%%%%%%%%%%%%%%%%%%% New computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
      
% UHAT = mkfhat(app,mesh,UDGS(:,iu,1),tau); % compute uhat at the implicit time   
[npf,nfe] = size(mesh.perm);

 UHAT = mkfhat(app,mesh,UDG(:,1:3,:),tau);
 UHAT = reshape(UHAT,[npf,nfe,ne]);
 UHATo = UHAT;
%  UHAT = UHAT(:,i);

% Compute QHAT
time = app.time;
kx  = app.k(1);
ky  = app.k(2);
kk  = sqrt(kx^2+ky^2);

shapfc = mkshape(mesh.porder,mesh.plocfc,mesh.plocfc,mesh.elemtype);
dshapft  = shapfc(:,:,2)';

nd = mesh.nd;
[npf,nfe] = size(mesh.perm);
[~,~,ne] = size(UDG(:,iu,:));
perm = mesh.perm;
t2t = mesh.t2t; %mkt2t(mesh.t,mesh.elemtype);

QHAT = zeros(npf,nfe,1);
U1 = zeros(npf,nfe,1);
e1 = i;    
    for j = 1:nfe % for each face of element i        
        e2 = t2t(i,j);
        u1 = UDG(perm(:,j),iu,e1);        % trying e2
        pb = mesh.dgnodes(perm(:,j),:,e1);
        dpg = reshape(dshapft*pb,[npf nd]);     
        jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nl   = [dpg(:,2)./jac,-dpg(:,1)./jac];
        
           
        qn = sum(u1(:,2:end).*nl,2);
        QHAT(:,j,1) = qn - tau*(u1(:,1) - UHAT(:,j,e1));
        U1(:,j,1) = qn - tau*(u1(:,1));
    end
    QHAT = reshape(QHAT,[npf*nfe,1]);
    U1 = reshape(U1,[npf*nfe,1]);
    
    % Ru = B*UHAT + G
    % need to compute B and G
    
[Ru,B,G,B2,G2] = hdgres_imex_bound(app,master,mesh,UDGS(:,iu,1),tau,UHAT,QHAT,U1); 
% [Ru,B,G,B2,G2] = hdgres_imex_bound(app,master,mesh,UDGS(:,:,1),tau,UHAT,QHAT,U1); 

%%%%%%%%%%%%%%%%DEBUGGING
% Bcoupled = B(:,:,1)

UHAT = reshape(UHAT(:,:,i),[npf*nfe,1]);
debug = 1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% UDG1 = UDGS(:,:,1);
% for k=1:j-1
%    UDG1 = UDG1 - (D(j,k)/D(j,j))*(UDGS(:,:,k+1) - UDGS(:,:,1));
% end    
%         
% %  UDGS(:,iu,j+1) = UDG1(:,iu,1) - (dtj*Minv(:,:,i))*(B*UHAT + G);
%  UDGS(:,iu,j+1) = UDG1(:,iu,1) - (dtj*Minv(:,:,i))*(Ru);
%  UDGS(:,end,j+1) = UDG1(:,end) + dtj*UDGS(:,1,j);   
% %  UDGS(:,iu,j+1) = (dtj*Minv(:,:,i)*B2)*UHAT - dtj*Minv(:,:,i)*G2 + UDG1(:,iu);
% 
%     
% if expE == 0
%     UDGS(:,:,:,4) = UDG;
% end
% 
% UDG = UDGS(:,:,:,end);
% a = 2;


function D = ERKcoeff(q,p)
% q - number of stages
% p - order of accuracy

if q == 1 && p == 1    
    A = 0;
    b = 1; 
    c = 1;
elseif q == 3 && p == 2
    alpha = 1 - (sqrt(2)/2);
    delta = -(2*sqrt(2))/3;
    
    A = [0,     0,       0;
         alpha, 0,       0;
         delta, 1-delta, 0];
    b = [0, 1-delta, 0]; 
    c = [0, alpha,    1];

elseif q == 3 && p == 3
    alpha = (3 + sqrt(3))/6;
    
    A = [0,       0,           0;
         alpha,   0,           0;
         alpha-1, 2*(1-alpha), 0];
    b = [0, 0.5,      0.5]; 
    c = [0, alpha,    1-alpha];    
    
elseif q == 4 && p == 3    
    A = [0,            0,            0,            0;
         0.4358665215, 0,            0,            0;
         0.3212788860, 0.3966543747, 0,            0;
         -0.105858296, 0.5529291479, 0.5529291479, 0];
     
     b = [0, 1.208496649, -0.644363171, 0.4358665215];
     c = [0, 0.4358665215, 0.7179332608, 1];
elseif q == 4 && p == 4    
    A = [0,            0,            0,            0;
         1/2, 0,            0,            0;
         0, 1/2, 0,            0;
         0, 0, 1, 0];
     
     b = [1/6, 1/3, 1/3, 1/6];
     c = [0, 1/2, 1/2, 1];     
else
    error('Invalid (q,p) combination');
end

Dinv = A(2:end,1:end);
Dinv = cat(1,Dinv,b);
D = inv(Dinv);




function [Ru,B,G,B2,G2] = hdgres_imex_bound(app,master,mesh,UDG,tau,UHAT,QHAT,U1)


nfe = size(mesh.perm,2);
ne = 1;
npv = master.npv;
ngv = master.ngv;
nd = mesh.nd;
npf = master.npf;
ngf = master.ngf;


dshapvt(:,:,1) = master.shapvl(:,:,2)';
dshapvt(:,:,2) = master.shapvl(:,:,3)';
dshapvt = reshape(permute(dshapvt,[1 3 2]),[ngv*nd npv]);
%shapft  = master.shapfc(:,:,1)';
dshapft = master.shapfc(:,:,2)';
perm = master.perm;
%bcm  = mesh.bcm;

shapvgdotshapvl  = zeros(npv*npv,ngv,nd+1);      
for d=1:nd+1    
    shapvg = master.shapvl(:,:,d)*diag(master.gwvl);    
    for ii=1:npv
        for jj = 1:npv
            shapvgdotshapvl((ii-1)*npv+jj,:,d) = shapvg(jj,:).*master.shapvl(ii,:,1);                    
        end
    end            
end
shapfgdotshapfc  = zeros(npf*npf,ngf,nd);      
for d=1:nd    
    shapfg = master.shapfc(:,:,d)*diag(master.gwfc);
    for ii=1:npf
        for jj = 1:npf
            shapfgdotshapfc((ii-1)*npf+jj,:,d) = shapfg(jj,:).*master.shapfc(ii,:,1);                    
        end
    end            
end
% 

Ru = zeros(npv,3,ne);
nblk(1,:) = 1:100:ne;
nblk(2,:) = [nblk(1,2:end) ne];
for n = 1:size(nblk,2)
    e1 = nblk(1,n);
    e2 = nblk(2,n);                         
    ns = (e2-e1+1);

    pn = reshape(mesh.dgnodes(:,1:nd,e1:e2),[npv nd*ns]);    
    Jg = reshape(dshapvt*pn,[ngv nd nd ns]);  
    Jg = permute(Jg,[1 4 2 3]);
    jac = Jg(:,:,1,1).*Jg(:,:,2,2) - Jg(:,:,1,2).*Jg(:,:,2,1);
    xxi = Jg(:,:,1,1)./jac;
    xet = Jg(:,:,2,1)./jac;
    yxi = Jg(:,:,1,2)./jac;
    yet = Jg(:,:,2,2)./jac;
        
    EE = zeros(npv,npf*nfe,ns);
    CX = zeros(npv,npf*nfe,ns);
    CY = zeros(npv,npf*nfe,ns);
    
%     tm = shapvgdotshapvl(:,:,1)*reshape(jac,[ngv ns]);
%     AA = reshape(tm,[npv npv ns]);

    tm = shapvgdotshapvl(:,:,2)*reshape(jac.*yet,[ngv ns])-...
         shapvgdotshapvl(:,:,3)*reshape(jac.*yxi,[ngv ns]);
    BX = reshape(tm,[npv npv ns]);

    tm = -shapvgdotshapvl(:,:,2)*reshape(jac.*xet,[ngv ns])+...
          shapvgdotshapvl(:,:,3)*reshape(jac.*xxi,[ngv ns]);
    BY = reshape(tm,[npv npv ns]);
    for jf=1:nfe
        I = perm(:,jf);        
        J = ((jf-1)*npf+1):jf*npf; 
        pb = pn(I,:);
        %pg = reshape(shapft*pb,[ngf nd ns]);            
        dpg = reshape(dshapft*pb,[ngf nd ns]);        
        detf = reshape(sqrt(dpg(:,1,:).^2+dpg(:,2,:).^2),[ngf ns]);
        nx = reshape(dpg(:,2,:),[ngf ns])./detf;
        ny = reshape(-dpg(:,1,:),[ngf ns])./detf;                    
        
        tm = shapfgdotshapfc(:,:,1)*reshape(detf,[ngf ns]);
        EE(I,J,:) = EE(I,J,:) + reshape(tm,[npf npf ns]);  
        tm = shapfgdotshapfc(:,:,1)*reshape(detf.*nx,[ngf ns]);
        CX(I,J,:) = CX(I,J,:) + reshape(tm,[npf npf ns]); 
        tm = shapfgdotshapfc(:,:,1)*reshape(detf.*ny,[ngf ns]);
        CY(I,J,:) = CY(I,J,:) + reshape(tm,[npf npf ns]);        
    end
    
    % B Matrix 
    B1 = -tau*EE;
    B2 = -CX;
    B3 = -CY;
    
    B = [B1;B2;B3];
%     B = [B3;B2;B1];
    B2 = [B1, B2, B3];
    
    % G Matrix
    G1 = BX*UDG(:,2) + BY*UDG(:,3) - EE*U1;
    G2 = BX*UDG(:,1);
    G3 = BY*UDG(:,1);
    
    G = [G1;G2;G3];
%     G = [G3;G2;G1];
    G2 = [G1, G2, G3];
    UHAT = reshape(UHAT(:,:,e1),[npf*nfe,1]);
%     for i = 1:ns
%         ii = e1 + i - 1; 
        Ru(:,1) = BX(:,:)*UDG(:,2) + BY(:,:)*UDG(:,3) - EE(:,:)*QHAT;
        Ru(:,2) = BX(:,:)*UDG(:,1) - CX(:,:)*UHAT;
        Ru(:,3) = BY(:,:)*UDG(:,1) - CY(:,:)*UHAT;
%         %UDG(:,:,ii) = UDG(:,:,ii) - dt*AA(:,:,i)\Ru;                
%     end    
end

