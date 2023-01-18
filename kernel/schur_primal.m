function [DUDG, DUDG_DUH, AE, FE] = schur_primal(M, C, E, L, Q, BD, F, GK, H, Rq, Ru, Rh)        

warning('OFF');

if isempty(M)==1
    %[DUDG, DUDG_DUH, AE, FE] = schur_u2(BD, F, GK, H, Ru, Rh);        
    if isempty(L)==1
        [DUDG, DUDG_DUH, AE, FE] = schur_u2(BD, F, GK, H, Ru, Rh);        
    else
        [DUDG, DUDG_DUH, AE, FE] = schur_u3(BD, F, L, GK, H, Ru, Rh); 
    end
elseif isempty(Q)==1    
    [DUDG, DUDG_DUH, AE, FE] = schur_uq2(M, C, E, BD, F, GK, H, Rq, Ru, Rh);        
else
    [DUDG, DUDG_DUH, AE, FE] = schur_upq(M, C, E, L, Q, BD, F, GK, H, Rq, Ru, Rh);
end

function  [DUDG, DUDG_DUH, AE, FE] = schur_u2(BD, F, GK, H, Ru, Rh)        
% schur complement for the following system
% [D F] [du] = [Ru]
% [K H] [dh] = [Rh]

% BD: npv*ncu*npv*nc*ne
% F:  npv*ncu*ncu*ndf*ne
% GK: ncu*ndf*npv*nc*ne
% H : ncu*ndf*ncu*ndf*ne
% Ru: npv*ncu*ne
% Rh: ncu*ndf*ne

npv = size(BD,1);
ncu = size(BD,2);
nc  = size(BD,4);
ne  = size(BD,5);
ndf = size(H,2);

DUDG     = zeros(npv,nc,ne);
DUDG_DUH = zeros(npv,nc,ncu,ndf,ne);
FE  = zeros(ncu,ndf,ne);
AE  = zeros(ncu,ndf,ncu,ndf,ne);
for i = 1:ne    
    SD    = reshape(BD(:,:,:,1:ncu,i),npv*ncu,npv*ncu);    
    
    du   = SD\(reshape(Ru(:,:,i),[npv*ncu 1]));       
    DUDG(:,:,i) = reshape(du,[npv ncu]);
    
    dlu = SD\(-reshape(F(:,:,:,:,i),[npv*ncu ncu*ndf]));        
    DUDG_DUH(:,:,:,:,i) = reshape(dlu,npv,ncu,ncu,ndf);
    
    tmp = reshape(GK(:,:,:,:,i),[ncu*ndf npv*nc])*reshape(DUDG(:,:,i),[npv*nc 1]);        
    FE(:,:,i) = Rh(:,:,i) - reshape(tmp,[ncu ndf]);

    tmp = reshape(GK(:,:,:,:,i),[ncu*ndf npv*nc])*reshape(DUDG_DUH(:,:,:,:,i),[npv*nc ncu*ndf]);
    AE(:,:,:,:,i) = H(:,:,:,:,i) + reshape(tmp,[ncu ndf ncu ndf]);      
end

function  [DUDG, DUDG_DUH, AE, FE] = schur_u3(BD, F, L, GK, H, Ru, Rh)        
% schur complement for the following system
% [D F] [du] = [Ru]
% [K H] [dh] = [Rh]

% BD: npv*ncu*npv*nc*ne
% F:  npv*ncu*ncu*ndf*ne
% GK: ncu*ndf*npv*nc*ne
% H : ncu*ndf*ncu*ndf*ne
% Ru: npv*ncu*ne
% Rh: ncu*ndf*ne

npv = size(BD,1);
ncu = size(BD,2);
nc  = size(BD,4);
ne  = size(BD,5);
ndf = size(H,2);

DUDG     = zeros(npv,nc,ne);
DUDG_DUH = zeros(npv,nc,ncu,ndf,ne);
FE  = zeros(ncu,ndf,ne);
AE  = zeros(ncu,ndf,ncu,ndf,ne);
for i = 1:ne    
    SD = reshape(BD(:,:,:,1:ncu,i),npv*ncu,npv*ncu);    
    
    if (i==ne)
        [size(SD) rank(SD) 0*cond(SD)]
    end

    LD = L(:,i)';        
    SD(end+1,(end-npv+1):end) = LD;
    
    du   = SD\[reshape(Ru(:,:,i),[npv*ncu 1]); 0];       
    DUDG(:,:,i) = reshape(du,[npv ncu]);
    
    dlu = SD\[-reshape(F(:,:,:,:,i),[npv*ncu ncu*ndf]); zeros(1,ncu*ndf)];        
    DUDG_DUH(:,:,:,:,i) = reshape(dlu,npv,ncu,ncu,ndf);
    
%     LD = L(:,i)';        
%     SD(end,(end-npv+1):end) = LD;
    
%     du   = SD\reshape(Ru(:,:,i),[npv*ncu 1]);       
%     DUDG(:,:,i) = reshape(du,[npv ncu]);
%     
%     dlu = SD\(-reshape(F(:,:,:,:,i),[npv*ncu ncu*ndf]));        
%     DUDG_DUH(:,:,:,:,i) = reshape(dlu,npv,ncu,ncu,ndf);

    if (i==ne)
        [size(SD) rank(SD) 0*cond(SD)]
    end
    
    tmp = reshape(GK(:,:,:,:,i),[ncu*ndf npv*nc])*reshape(DUDG(:,:,i),[npv*nc 1]);        
    FE(:,:,i) = Rh(:,:,i) - reshape(tmp,[ncu ndf]);
    
    tmp = reshape(GK(:,:,:,:,i),[ncu*ndf npv*nc])*reshape(DUDG_DUH(:,:,:,:,i),[npv*nc ncu*ndf]);
    AE(:,:,:,:,i) = H(:,:,:,:,i) + reshape(tmp,[ncu ndf ncu ndf]);      
end

function [DUDG, DUDG_DUH, AE, FE] = schur_uq2(M, C, E, BD, F, GK, H, Rq, Ru, Rh)        
% schur complement for the following system
% [M -C  E] [dq] = [Rq]
% [B  D  F] [du] = [Ru]
% [G  K  H] [dh] = [Rh]

% M:  npv*npv*ne
% C:  npv*npv*nd*ne
% E:  npv*ndf*nd*ne
% BD: npv*ncu*npv*nc*ne
% F:  npv*ncu*ncu*ndf*ne
% GK: ncu*ndf*npv*nc*ne
% H : ncu*ndf*ncu*ndf*ne
% Rq: npv*ncu*nd*ne
% Ru: npv*ncu*ne
% Rh: ncu*ndf*ne

npv = size(M,1);
ne  = size(M,3);
ncu = size(BD,2);
nc  = size(BD,4);
ndf = size(E,2);
nd  = size(E,3);

MiC = zeros(npv,npv,nd);
MiE = zeros(npv,ndf,nd);
MiCdu = zeros(npv, ncu, nd);
MiCdlu = zeros(npv,ncu,nd,ncu,ndf);
MiEdlu = zeros(npv,ncu,nd,ncu,ndf);
DUDG     = zeros(npv,nc,ne);
DUDG_DUH = zeros(npv,nc,ncu,ndf,ne);
FE  = zeros(ncu,ndf,ne);
AE  = zeros(ncu,ndf,ncu,ndf,ne);
for i = 1:ne
    Mi = inv(M(:,:,i));
    MiRq = Mi*reshape(Rq(:,:,:,i),[npv ncu*nd]);
    for d = 1:nd
        MiC(:,:,d)  = Mi*C(:,:,d,i);
        MiE(:,:,d)  = Mi*E(:,:,d,i);
    end    
    
    D    = reshape(BD(:,:,:,1:ncu,i),npv,ncu,npv,ncu);
    B    = reshape(BD(:,:,:,ncu+1:end,i),[npv,ncu,npv,ncu,nd]);           
    BMiR = reshape(B,[npv*ncu npv*ncu*nd])*reshape(MiRq,[npv*ncu*nd 1]);        
    BMiC = mapContractK(B,MiC,[1 2 4],[3 5],[],[1 3],2,[]);   
    BMiC = permute(BMiC,[1 2 4 3]);
    BMiE = mapContractK(B,MiE,[1 2 4],[3 5],[],[1 3],2,[]);   

    %reshape(D,[npv*ncu,npv*ncu])
    %reshape(BMiC,[npv*ncu,npv*ncu])
    %reshape(GK(:,:,:,1:ncu,i),[ncu*ndf npv*ncu]);
%     G = reshape(GK(:,:,:,ncu+1:end,i),[ncu,ndf,npv,ncu,nd]);    
%     GMiC = mapContractK(G,MiC,[1 2 4],[3 5],[],[1 3],2,[]);   
%     GMiC = permute(GMiC,[1 2 4 3]);
%     reshape(GMiC,[ndf*ncu,npv*ncu])
%     max(abs(G(:)))
%     pause
        
    SD = reshape(D+BMiC,[npv*ncu,npv*ncu]);    
    du   = SD\(reshape(Ru(:,:,i),[npv*ncu 1]) - BMiR(:));       
    for d=1:nd
        MiCdu(:,:,d) = MiC(:,:,d)*reshape(du,[npv ncu]);            
    end        
    dq   = MiRq + reshape(MiCdu,[npv ncu*nd]);
    DUDG(:,:,i) = [reshape(du,[npv ncu]), reshape(dq,[npv ncu*nd])];  
    
    dlu = SD\(reshape(BMiE,[npv*ncu ncu*ndf])-reshape(F(:,:,:,:,i),[npv*ncu ncu*ndf]));     
%    dlu    
%     fileID = fopen(['D.bin'],'r');
%     Dt = fread(fileID,npv*ncu*npv*ncu,'double');
%     max(abs(SD(:)-Dt(:)))
%     fileID = fopen(['F.bin'],'r');
%     Ft = fread(fileID,npv*ncu*ncu*ndf,'double');
%     tmp = reshape(F(:,:,:,:,i),[npv*ncu ncu*ndf])-reshape(BMiE,[npv*ncu ncu*ndf]);    
%     max(abs(tmp(:)-Ft(:)))
%     fileID = fopen(['H.bin'],'r');
%     Ht = fread(fileID,ndf*ncu*ncu*ndf,'double');
%     G = reshape(GK(:,:,:,ncu+1:end,i),[ncu,ndf,npv,ncu,nd]);    
%     GMiE = mapContractK(G,MiE,[1 2 4],[3 5],[],[1 3],2,[]);   
%     tmp = reshape(H(:,:,:,:,i),[ndf*ncu ncu*ndf]) - reshape(GMiE,[ndf*ncu ncu*ndf]);    
%     max(abs(tmp(:)-Ht(:)))
%     %reshape(H(:,:,:,:,i),[ndf*ncu ncu*ndf])    
%     pause

    for j=1:nd
        MiCdlu(:,:,j,:,:) = reshape(MiC(:,:,j)*reshape(dlu,[npv ncu*ncu*ndf]),[npv ncu 1 ncu ndf]);            
        MiEdlu(:,:,j,:,:) = bsxfun(@times,reshape(MiE(:,:,j),[npv 1 1 ndf]),reshape(eye(ncu),[1 ncu ncu 1]));    
    end                                
    dlq = reshape(MiCdlu,[npv ncu*nd ncu ndf]) - reshape(MiEdlu,[npv ncu*nd ncu ndf]);    
    DUDG_DUH(:,:,:,:,i) = cat(2,reshape(dlu,npv,ncu,ncu,ndf), reshape(dlq,npv,ncu*nd,ncu,ndf));                
    
    tmp = reshape(GK(:,:,:,:,i),[ncu*ndf npv*nc])*reshape(DUDG(:,:,i),[npv*nc 1]);        
    FE(:,:,i) = Rh(:,:,i) - reshape(tmp,[ncu ndf]);

    tmp = reshape(GK(:,:,:,:,i),[ncu*ndf npv*nc])*reshape(DUDG_DUH(:,:,:,:,i),[npv*nc ncu*ndf]);
    AE(:,:,:,:,i) = H(:,:,:,:,i) + reshape(tmp,[ncu ndf ncu ndf]);     
%     reshape(AE(:,:,:,:,i),[ncu*ndf ncu*ndf])
%     pause
end


