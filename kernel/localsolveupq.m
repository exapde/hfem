function [DUDG, DUDG_duh] = localsolveupq(BD1, Ru1, Rlu1, MiC1, MiE1,...
    MiL1, MiQ1, UDG, UH, SH, fc_q, fc_p, alag, alpha)

npv = size(UDG,1);
nc  = size(UDG,2);
[ndf nch ne] = size(UH);
ncu = size(Ru1,2);
nd = size(MiC1,2);

if alag==0
    alpha = 0;
    DUDG_dp = [];
else
    DUDG_dp  = zeros(npv,nc,npv,ne);
end

Rlu1 = permute(Rlu1,[1 2 4 3 5]);
DUDG     = zeros(npv,nc,ne);
DUDG_duh = zeros(npv,nc,ndf,nch,ne);
for k=1:1:ne % For each element     
    BD  = BD1(:,:,:,:,k);    
    Ru  = Ru1(:,:,k);                 
    Ru  = reshape(Ru,npv*ncu,1);
    Rlu = Rlu1(:,:,:,:,k);       
    MiC = reshape(MiC1(:,:,:,k)/fc_q, [npv*nd npv]);
    MiE = reshape(MiE1(:,:,:,k)/fc_q, [npv*nd ndf]);
    MiQ = reshape(MiQ1(:,:,:,k)/(fc_p+alpha), [npv npv*nch*nd]);
    MiL = reshape(MiL1(:,k)/(fc_p+alpha), [npv 1]);
        
    % obtain the current solution and numerical trace    
    u    = UDG(:,1:nch,k);
    p    = UDG(:,nch+1,k);
    q    = UDG(:,nch+2:end,k);
    uh   = UH(:,:,k);     
    sh   = SH(:,:,k);                
    
    MiCu = permute(reshape(MiC*u,[npv nd nch]), [1 3 2]);
    MiEu = permute(reshape(MiE*uh,[npv nd nch]), [1 3 2]);
    MiR  = sh(:,nch+2:end)/fc_q + reshape(MiEu,[npv nch*nd]) - reshape(MiCu,[npv nch*nd]) - q;
        
    % obtain B and D
    D    =  reshape(BD(:,:,:,1:nch),npv,nch,npv,nch);
    B    = -reshape(BD(:,:,:,nch+2:end),[npv,nch,npv,nch,nd]);       
    
    BMiR =  reshape(B,[npv*ncu npv*nch*nd])*reshape(MiR,[npv*nch*nd 1]);    
        
    B    =  reshape(permute(B,[1 2 4 3 5]),[npv*nch*nch npv*nd]);    
    BMiC =  reshape(B*MiC,[npv nch nch npv]);
    BMiC =  permute(BMiC,[1 2 4 3]);
    BMiE =  reshape(B*MiE,[npv nch nch ndf]);
    BMiE =  permute(BMiE,[1 2 4 3]);
        
    Bp    =  -reshape(BD(:,:,:,nch+1),[npv*nch,npv]);               
    MiRp  =  -fc_p*p/(fc_p+alpha) + MiL + reshape(MiQ,[npv npv*nch*nd])*MiR(:);
    BMiRp =  Bp*MiRp;    
    
    MiQt  =  MiQ;
    MiQ   =  permute(reshape(MiQ,[npv npv nch nd]),[1 3 2 4]);
    
    MiQC  =  reshape(reshape(MiQ,[npv*nch npv*nd])*MiC,[npv nch*npv]);            
    BMiQC =  reshape(Bp*MiQC,[npv nch nch npv]);
    BMiQC =  permute(BMiQC,[1 2 4 3]);
%     MiQC
%     pause
    
    MiQE  =  reshape(reshape(MiQ,[npv*nch npv*nd])*MiE,[npv nch*ndf]);        
    BMiQE =  reshape(Bp*MiQE,[npv nch nch ndf]);
    BMiQE =  permute(BMiQE,[1 2 4 3]);
        
    % the reduced local matrix
    SD = reshape(BMiQC+BMiC+D,[npv*nch,npv*nch]);        
%     reshape(D,[npv*nch,npv*nch])       
%     pause
    
    % solve the local system with the residual on the RHS
    du   = SD\(Ru(:) + BMiR(:) + BMiRp(:));
    MiCdu= permute(reshape(MiC*reshape(du,[npv nch]),[npv nd nch]), [1 3 2]);
    dq   = MiR - reshape(MiCdu,[npv nch*nd]);
    dp   = -fc_p*p/(fc_p+alpha) + MiL + reshape(MiQt,[npv npv*nch*nd])*dq(:);
    dudg = [reshape(du,[npv nch]), reshape(dp,[npv 1]), reshape(dq,[npv nch*nd])];
                
    % solve the local system with the numerical trace on the RHS            
    dlu = SD\(reshape(Rlu,[npv*ncu ndf*nch])+reshape(BMiE + BMiQE,[npv*nch ndf*nch]));        
    
    
    MiCdlu = permute(reshape(MiC*reshape(dlu,[npv nch*ndf*nch]),[npv nd nch ndf nch]), [1 3 2 4 5]);
    AiE = bsxfun(@times,reshape(MiE,[npv 1 nd ndf 1]),reshape(eye(nch),[1 nch 1 1 nch]));    
    dlq = reshape(AiE,[npv nch*nd ndf nch]) - reshape(MiCdlu,[npv nch*nd ndf nch]);    
    dlp = reshape(MiQt,[npv npv*nch*nd])*reshape(dlq,[npv*nch*nd ndf*nch]);
    dudg_duh = cat(2,reshape(dlu,npv,nch,ndf,nch), reshape(dlp,npv,1,ndf,nch));
    dudg_duh = cat(2,dudg_duh, reshape(dlq,npv,nch*nd,ndf,nch));        
                
    % obtain local solutions            
    DUDG(:,:,k) = dudg;
    DUDG_duh(:,:,:,:,k) = dudg_duh;                      
            
%     MiCx = reshape(MiC1(:,1,:,k)/fc_q, [npv npv]);
%     MiCy = reshape(MiC1(:,2,:,k)/fc_q, [npv npv]);    
%     AiC = [MiCx/(fc_p+alpha) MiCy/(fc_p+alpha); MiCx 0*MiCy; 0*MiCy MiCx; MiCy 0*MiCy; 0*MiCx MiCy];    
%     D    =  reshape(BD(:,:,:,1:nch),npv*nch,npv*nch);
%     B    = -reshape(BD(:,:,:,nch+1:end),[npv*nch,npv*(nd*nch+1)]);       
%     SD2 = B*AiC+D;
%     
%     max(abs(SD(:)-SD2(:)))
%     SD2
%     B    
%     pause    
% 
%     MiCx = reshape(MiC1(:,1,:,k)/fc_q, [npv npv]);
%     MiCy = reshape(MiC1(:,2,:,k)/fc_q, [npv npv]);    
%     MiCz = reshape(MiC1(:,3,:,k)/fc_q, [npv npv]);    
%     AiC = [MiCx/(fc_p+alpha) MiCy/(fc_p+alpha) MiCz/(fc_p+alpha);...
%            MiCx 0*MiCy 0*MiCz;...
%            0*MiCy MiCx 0*MiCz;...
%            0*MiCy 0*MiCx MiCx;...
%            MiCy 0*MiCy 0*MiCz;...
%            0*MiCy MiCy 0*MiCz;...
%            0*MiCy 0*MiCx MiCy;...
%            MiCz 0*MiCy 0*MiCz;...
%            0*MiCy MiCz 0*MiCz;...
%            0*MiCy 0*MiCx MiCz;];    
%     D    =  reshape(BD(:,:,:,1:nch),npv*nch,npv*nch);
%     B    = -reshape(BD(:,:,:,nch+1:end),[npv*nch,npv*(nd*nch+1)]);       
%     SD2 = B*AiC+D;
%     
%     max(abs(SD(:)-SD2(:)))
%     pause    
end
DUDG_duh = permute(DUDG_duh,[1 2 4 3 5]);
