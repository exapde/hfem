function [rp, ci, ai, ncols] = sparse_to_csr(A)
% SPARSE_TO_CSR Convert a sparse matrix into compressed row storage arrays
% 
% [rp ci ai] = sparse_to_csr(A) returns the row pointer (rp), column index
% (ci) and value index (ai) arrays of a compressed sparse representation of
% the matrix A.

[nzi, nzj, nzv] = find(A);
nrows=max(nzi);    
ncols=max(nzj);    
nz = length(nzi);
ci = zeros(nz,1); 
ai = zeros(nz,1); 
rp = zeros(nrows+1,1);

if max(nzi) ~= unique(nzi) | max(nzj) ~= unique(nzj)
    error('Warning: There are ghost entities in the mesh.');
end

for i=1:nz
    rp(nzi(i)+1)=rp(nzi(i)+1)+1;
end
rp=cumsum(rp);
for i=1:nz
    ai(rp(nzi(i))+1)=nzv(i);
    ci(rp(nzi(i))+1)=nzj(i);
    rp(nzi(i))=rp(nzi(i))+1;
end
for i=nrows:-1:1
    rp(i+1)=rp(i);
end
rp(1)=0;
