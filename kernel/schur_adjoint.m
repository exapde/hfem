function [DUDG, DUDG_DUH, AE, FE] = schur_adjoint(M, C, E, L, Q, BD, F, GK, H, Jq, Ju, Jh)        

if isempty(M)==1
    [DUDG, DUDG_DUH, AE, FE] = schur_u_adjoint(BD, F, GK, H, Ju, Jh);        
elseif isempty(Q)==1
%    [DUDG2, DUDG2_DUH, AE2, FE2] = schur_uq_adjoint2(M, C, E, BD, F, GK, H, Jq, Ju, Jh);
    [DUDG, DUDG_DUH, AE, FE] = schur_uq_adjoint(M, C, E, BD, F, GK, H, Jq, Ju, Jh);
%     max(abs(DUDG2(:)-DUDG(:)))
%     max(abs(DUDG2_DUH(:)-DUDG_DUH(:)))
%     max(abs(AE2(:)-AE(:)))
%     max(abs(FE2(:)-FE(:)))
%     pause
else
    [DUDG, DUDG_DUH, AE, FE] = schur_upq_adjoint(M, C, E, L, Q, BD, F, GK, H, Jq, Ju, Jh);
end

function  [DUDG, DUDG_DUH, AE, FE] = schur_u_adjoint2(BD, F, GK, H, Ju, Jh)        
% schur complement for the following system
% [D^T K^T] [du] = [Ju]
% [F^T H^T] [dh] = [Jh]

% BD: npv*ncu*npv*nc*ne
% F:  npv*ncu*ncu*ndf*ne
% GK: ncu*ndf*npv*nc*ne
% H : ncu*ndf*ncu*ndf*ne
% Ju: npv*ncu*ne
% Jh: ncu*ndf*ne

%[~, ~, ~, ~, SDE] = schur_uq(M, C, E, BD, F, GK, H, Jq, Ju, Jh);

npv = size(BD,1);
ncu = size(BD,2);
nc  = size(BD,4);
ne  = size(BD,5);
ndf = size(H,2);

% Taking transpose
BD = permute(BD,[3 4 1 2 5]);
GK = permute(GK,[3 4 1 2 5]);
F  = permute(F,[3 4 1 2 5]);
H  = permute(H,[3 4 1 2 5]);

DUDG     = zeros(npv,nc,ne);
DUDG_DUH = zeros(npv,nc,ncu,ndf,ne);
FE  = zeros(ncu,ndf,ne);
AE  = zeros(ncu,ndf,ncu,ndf,ne);    
for i = 1:ne    
    Dt = reshape(BD(:,1:ncu,:,:,i),npv,ncu,npv,ncu);     
    Kt = reshape(GK(:,1:ncu,:,:,i),npv,ncu,ncu,ndf);       
    
    SDt  = reshape(Dt(:),[npv*ncu,npv*ncu]);     
    dv   = SDt\(reshape(Ju(:,:,i),[npv*ncu 1]));            
    DUDG(:,:,i) = reshape(dv,[npv ncu]);
    
    dlv = SDt\(-reshape(Kt,[npv*ncu ncu*ndf]));     
    DUDG_DUH(:,:,:,:,i) = reshape(dlv,npv,ncu,ncu,ndf);
    
    EF = reshape(F(:,:,:,:,i),[ncu*ndf npv*ncu]);

    tmp = EF*reshape(DUDG(:,:,i),[npv*nc 1]);        
    FE(:,:,i) = Jh(:,:,i) - reshape(tmp,[ncu ndf]);

    tmp = EF*reshape(DUDG_DUH(:,:,:,:,i),[npv*nc ncu*ndf]);
    AE(:,:,:,:,i) = H(:,:,:,:,i) + reshape(tmp,[ncu ndf ncu ndf]);    
end


function  [DUDG, DUDG_DUH, AE, FE] = schur_uq_adjoint2(M, C, E, BD, F, GK, H, Jq, Ju, Jh)        
% schur complement for the following system
% [ M^T B^T G^T] [dq] = [Jq]
% [-C^T D^T K^T] [du] = [Ju]
% [ E^T F^T H^T] [dh] = [Jh]

% M:  npv*npv*ne
% C:  npv*npv*nd*ne
% E:  npv*ndf*nd*ne
% BD: npv*ncu*npv*nc*ne
% F:  npv*ncu*ncu*ndf*ne
% GK: ncu*ndf*npv*nc*ne
% H : ncu*ndf*ncu*ndf*ne
% Jq: npv*ncu*nd*ne
% Ju: npv*ncu*ne
% Jh: ncu*ndf*ne


npv = size(M,1);
ne  = size(M,3);
ncu = size(BD,2);
nc  = size(BD,4);
ndf = size(E,2);
nd  = size(E,3);

%Jq  = reshape(Jq, [npv ncu nd ne]);
    
% Taking transpose
C  = permute(C,[2 1 3 4]);
E  = permute(E,[2 1 3 4]);
BD = permute(BD,[3 4 1 2 5]);
GK = permute(GK,[3 4 1 2 5]);
F  = permute(F,[3 4 1 2 5]);
H  = permute(H,[3 4 1 2 5]);

DUDG     = zeros(npv,nc,ne);
DUDG_DUH = zeros(npv,nc,ncu,ndf,ne);
FE  = zeros(ncu,ndf,ne);
AE  = zeros(ncu,ndf,ncu,ndf,ne);
Et = zeros(ncu,ndf,npv,ncu,nd); 
for i = 1:ne
    Mi = inv(M(:,:,i));    
    Dt = reshape(BD(:,1:ncu,:,:,i),npv,ncu,npv,ncu);
    Bt = reshape(BD(:,ncu+1:end,:,:,i),[npv,ncu,nd,npv,ncu]);  
    Kt = reshape(GK(:,1:ncu,:,:,i),npv,ncu,ncu,ndf);
    Gt = reshape(GK(:,ncu+1:end,:,:,i),[npv,ncu,nd,ncu,ndf]);     

    MiJq = reshape(Mi*reshape(Jq(:,:,i),[npv ncu*nd]),[npv ncu nd]);
    MiBt = reshape(Mi*reshape(Bt,[npv ncu*nd*npv*ncu]),[npv ncu nd npv ncu]);    
    MiGt = reshape(Mi*reshape(Gt,[npv ncu*nd*ncu*ndf]),[npv ncu nd ncu ndf]);                
    CMiJ = C(:,:,1,i)*MiJq(:,:,1);
    CMiBt = C(:,:,1,i)*reshape(MiBt(:,:,1,:,:),[npv ncu*npv*ncu]);
    CMiGt = C(:,:,1,i)*reshape(MiGt(:,:,1,:,:),[npv ncu*ncu*ndf]);          
    for d = 2:nd
        CMiJ = CMiJ + C(:,:,d,i)*MiJq(:,:,d);
        CMiBt = CMiBt + C(:,:,d,i)*reshape(MiBt(:,:,d,:,:),[npv ncu*npv*ncu]);
        CMiGt = CMiGt + C(:,:,d,i)*reshape(MiGt(:,:,d,:,:),[npv ncu*ncu*ndf]);        
    end
        
    SDt  = reshape(Dt(:)+CMiBt(:),[npv*ncu,npv*ncu]);     
    dv   = SDt\(reshape(Ju(:,:,i),[npv*ncu 1]) + CMiJ(:));        
    dq   = MiJq - reshape(reshape(MiBt,[npv*ncu*nd npv*ncu])*dv,[npv ncu nd]);
    DUDG(:,:,i) = [reshape(dv,[npv ncu]), reshape(dq,[npv ncu*nd])];      
    
    dlv = SDt\(-reshape(CMiGt,[npv*ncu ncu*ndf])-reshape(Kt,[npv*ncu ncu*ndf])); 
    MiBdlv = reshape(MiBt,[npv*ncu*nd npv*ncu])*reshape(dlv,[npv*ncu ncu*ndf]);    
    dlq = -reshape(MiGt,[npv ncu*nd ncu ndf]) - reshape(MiBdlv,[npv ncu*nd ncu ndf]);    
    DUDG_DUH(:,:,:,:,i) = cat(2,reshape(dlv,npv,ncu,ncu,ndf), reshape(dlq,npv,ncu*nd,ncu,ndf));                    
    
    for d = 1:nd
        Et(:,:,:,:,d) = bsxfun(@times,reshape(E(:,:,d,i),[1 ndf npv 1]),reshape(eye(ncu),[ncu 1 1 ncu]));   
    end        
    EF = cat(2,reshape(F(:,:,:,:,i),[ncu*ndf npv*ncu]),reshape(Et,[ncu*ndf npv*ncu*nd]));

    tmp = EF*reshape(DUDG(:,:,i),[npv*nc 1]);        
    FE(:,:,i) = Jh(:,:,i) - reshape(tmp,[ncu ndf]);

    tmp = EF*reshape(DUDG_DUH(:,:,:,:,i),[npv*nc ncu*ndf]);
    AE(:,:,:,:,i) = H(:,:,:,:,i) + reshape(tmp,[ncu ndf ncu ndf]);    
end






