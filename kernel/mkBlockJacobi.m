function [Hg, Rg] = mkBlockJacobi(mesh,Hg,Rg)
% Apply the block Jacobi preconditioner to the system Hg x = Rg

nrows = mesh.cbsr_nrows;
rowpts= mesh.cbsr_rowpts;
bsz   = size(Rg,1);
Id    = eye(bsz);
for i = 1:nrows
    k = (rowpts(i)+1):rowpts(i+1);  
    Hi = inv(Hg(:,:,k(1)));
    Hg(:,:,k(1)) = Id;
    Rg(:,i) = Hi*Rg(:,i);
    for j = 2:length(k)
        Hg(:,:,k(j)) = Hi*Hg(:,:,k(j));
    end    
end
