function Xg = matvec(mesh,Hg,Rg)

nrows = mesh.cbsr_nrows;
rowpts    = mesh.cbsr_rowpts;
colind    = mesh.cbsr_colind;
bsz   = size(Rg,1);
Xg    = 0*Rg;
for i = 1:nrows
    k = (rowpts(i)+1):rowpts(i+1);
    j = colind(k);    
    Xg(:,i) = reshape(Hg(:,:,k),[bsz bsz*length(k)])*reshape(Rg(:,j),[bsz*length(k) 1]);
end





