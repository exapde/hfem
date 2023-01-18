function [ meshIM, meshEX ] = split_IMEX_mesh( mesh )
%SPLIT_IMEX_MESH Summary of this function goes here
%   Detailed explanation goes here

% The boundary expression is modified to create a new label for the
% faces that will be on the IMEX interface
bndexpr = mesh.bndexpr;
bndexpr{end+1} = 'true';

% Gets all the elements of the implicit mesh
ind = find(mesh.imex==1);
tIM = mesh.t(ind,:);
% Gets all the elements of the explicit mesh
ind = find(mesh.imex==0);
tEX = mesh.t(ind,:);
% Remove Unused Nodes
[pIM,tIM] = fixmesh(mesh.p,tIM);
[pEX,tEX] = fixmesh(mesh.p,tEX);

% Make full mesh Structures
meshIM = mkmesh(pIM,tIM,mesh.porder,bndexpr,mesh.elemtype,mesh.nodetype);
meshEX = mkmesh(pEX,tEX,mesh.porder,bndexpr,mesh.elemtype,mesh.nodetype);

% Create array of IMEX faces ordered the same way for both meshes
imexB = mesh.imexB; % Finds correct number for the IMEX boundary
fimex = find(meshIM.f(:,end)==-imexB);
meshIM.fimex = fimex;
centerIM = 0.5*(meshIM.p(meshIM.f(fimex,1),:) + meshIM.p(meshIM.f(fimex,2),:));
fimex = find(meshEX.f(:,end)==-imexB);
meshEX.fimex = fimex;
centerEX = 0.5*(meshEX.p(meshEX.f(fimex,1),:) + meshEX.p(meshEX.f(fimex,2),:));
% Sorting meshEX.fimex
for i=1:length(fimex)
    distCent = centerEX - ones(length(fimex),1)*centerIM(i,:);
    distCent = sum(distCent.^2,2);
    [~,ind] = min(distCent);
    meshEX.fimex(i) = fimex(ind);
end

end

