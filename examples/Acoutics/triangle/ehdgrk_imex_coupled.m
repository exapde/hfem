function UDG = ehdgrk_imex_coupled(master,mesh,app,Minv,UDG)

tau = app.tau;
dt = app.dt;
nstage = app.nstage;
torder = app.torder;
D = ERKcoeff(nstage,torder);

[npv,ncu,ne] = size(UDG); 
UDGS = zeros(npv,ncu,ne,nstage+1);
UDGS(:,:,:,1) = UDG;
iu = 1:(ncu-1);

ind = find(mesh.imex == 0); % Explicit element indicies
expE = length(ind); % Number of explicit elements

for j = 1:nstage
    dtj = dt/D(j,j);
    Ru = hdgres(app,master,mesh,UDGS(:,iu,:,j),tau);
    
    UDG1 = UDGS(:,:,:,1);
    for k=1:j-1
       UDG1 = UDG1 - (D(j,k)/D(j,j))*(UDGS(:,:,:,k+1) - UDGS(:,:,:,1));
    end    
    
    for count = 1:expE
        i = ind(count);
        UDGS(:,end,i,j+1) = UDG1(:,end,i) + dtj*UDGS(:,1,i,j);
        UDGS(:,iu,i,j+1) = UDG1(:,iu,i) - dtj*Minv(:,:,i)*Ru(:,:,i); 
    end
    
    if expE == 0
        UDGS(:,:,:,4) = UDG;
    end
    
%     for i = 1:ne
%         if mesh.imex(i) == 1
%             UDGS(:,end,i,j+1) = UDG1(:,end,i) + dtj*UDGS(:,1,i,j);
%             UDGS(:,iu,i,j+1) = UDG1(:,iu,i) - dtj*Minv(:,:,i)*Ru(:,:,i);  
%         end
%     end
end

UDG = UDGS(:,:,:,end);
a = 2;


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

function Ru = hdgres(app,master,mesh,UDG,tau)

[UHAT,QHAT] = mkfhat(app,mesh,UDG,tau);

nfe = size(mesh.perm,2);
ne = mesh.ne;
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
    
    for i = 1:ns
        ii = e1 + i - 1; 
        Ru(:,1,ii) = BX(:,:,i)*UDG(:,2,ii) + BY(:,:,i)*UDG(:,3,ii) - EE(:,:,i)*QHAT(:,ii);
        Ru(:,2,ii) = BX(:,:,i)*UDG(:,1,ii) - CX(:,:,i)*UHAT(:,ii);
        Ru(:,3,ii) = BY(:,:,i)*UDG(:,1,ii) - CY(:,:,i)*UHAT(:,ii);
        %UDG(:,:,ii) = UDG(:,:,ii) - dt*AA(:,:,i)\Ru;                
    end    
end

