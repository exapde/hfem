function [cg2dg, dg2cg, ncf, isEDGelement] = getcg2dg(perm,elcon)

perm = perm(:);
ndf = length(perm);
cg2dg = zeros(ndf,1);
dg2cg = zeros(ndf,1);

ncf = 0;
for i=1:ndf
    elemNode = perm(i);

    alreadyExistentCGnode = 0;
    for j = 1:ncf
        if (cg2dg(j) == elemNode)
            alreadyExistentCGnode = 1;
            break;
        end
    end

    if (alreadyExistentCGnode == 0) 
        cgNode = ncf+1;
        cg2dg(cgNode) = elemNode;    
        ncf = ncf+1;
    elseif (alreadyExistentCGnode == 1) 
        cgNode = j;
    end
    dg2cg(i) = cgNode;    
end

[cg2dg dg2cg perm]

% Convert "cg2elemNode" to "cg2dg":
for i=1:ncf
    elemNode = cg2dg(i);
    for j=1:ndf
        if (perm(j) == elemNode) 
            dgNode = j;
            break;
        end
    end    
    cg2dg(i) = dgNode;
end
[cg2dg dg2cg perm]

size(elcon)
ne = numel(elcon)/ndf;    
isEDGelement = zeros(ne,1);
for i=1:ne
    isEDGelement(i) = 1;    
    for j=1:ndf        
        if (elcon(ndf*(i-1)+cg2dg(dg2cg(j))) ~= elcon(ndf*(i-1)+j))            
            isEDGelement(i) = 0;
            break;
        end
    end
end   

%         for (i = 0; i < ndf; i++) {
%             elemNode = mesh.perm[i];
% 
%             alreadyExistentCGnode = 0;
%             for (j = 0; j < ncf; j++) {
%                 if (mesh.cg2dg[j] == elemNode) {
%                     alreadyExistentCGnode = 1;
%                     break;
%                 }
%             }
% 
%             if (alreadyExistentCGnode == 0) {
%                 cgNode = ncf;
%                 mesh.cg2dg[cgNode] = elemNode;      // Temporarily, cg2dg stores cg2elemNode
%                 ncf++;
%             }
%             else if (alreadyExistentCGnode == 1) {
%                 cgNode = j;
%             }
%             mesh.dg2cg[i] = cgNode;
%         }
%         
%         // Convert "cg2elemNode" to "cg2dg":
%         for (i = 0; i < ncf; i++) {
%             elemNode = mesh.cg2dg[i];
%             for (j = 0; j < ndf; j++) {
%                 if (mesh.perm[j] == elemNode) {
%                     dgNode = j;
%                     break;
%                 }
%             }
%             if (j == ndf) {
%                 printf("Error when computing mesh.cg2dg and mesh.dg2cg arrays\n");
%                 exit(1);
%             }
%             mesh.cg2dg[i] = dgNode;
%         }
%         mesh.ncf = ncf;            // TODO: Include in binary files!!        
