
% % 2D slice of 3D solution with MATLAB:
% velmin = -1;
% velmax = 1;
% vel = UDG(:,2,:)./UDG(:,1,:);
% func = @(p) all(abs(mesh.p(:,3)-0.05)<1e-5);
% figure(1); clf;
% faceplot(mesh,1,vel,func,[velmin velmax]);
% axis equal; axis tight; axis off; 
% colormap jet;  colorbar('location','EastOutSide','FontSize',12);

function writeINPfile(fileName, mesh, UDG)

% INPUTS:
% UDG: [npv, numFields, numElements]

UDG = real(UDG);
ncu = size(UDG,2) / (mesh.nd+1);

projectionFlag = 0;     % Convert from DG to CG nodes by taking the average of overlapping nodes
[mesh_p1,UCG_p1] = getP1meshAndUCG(mesh, UDG, projectionFlag);

nump1Nodes = mesh_p1.np;
nump1Elements = mesh_p1.ne;
numDimensions = mesh_p1.nd;
elementType = mesh_p1.elemtype;
numSolutionFields = size(UCG_p1,2);
defaultMaterialId = 0;
elementTypeString = getElementTypeString(numDimensions, elementType);

file = fopen([fileName,'.inp'],'w');

% Save header:
toWrite = num2str([nump1Nodes, nump1Elements, numSolutionFields, 0, 0]);
% toWrite = num2str([nump1Nodes, nump1Elements, 2, 0, 0]);
writeStringMatrix(file, toWrite);

% Save nodes:
if mesh.nd == 2; p_aux = [mesh_p1.p, zeros(nump1Nodes,1)];
elseif mesh.nd == 3; p_aux = mesh_p1.p; end
toWrite = num2str([(1:nump1Nodes)', p_aux]);
writeStringMatrix(file, toWrite);

% Save elements:
% for i=1:nump1Elements
%     toWrite = [num2str([i, defaultMaterialId]), '   ' elementTypeString '   ', num2str(mesh_p1.t(i,:))];
%     writeStringMatrix(file, toWrite);
% end
toWrite = [num2str((1:nump1Elements)'), repmat('   ',[nump1Elements,1]), num2str(repmat(defaultMaterialId,[nump1Elements,1])), repmat(['   ' elementTypeString '   '],[nump1Elements,1]), num2str(mesh_p1.t(:,:))];
writeStringMatrix(file, toWrite);

% Save nodal fields:
solutionMatrix = UCG_p1;


% toWrite = ['2   1'];
% writeStringMatrix(file, toWrite);
% toWrite = 'pointField, None';        % None refers to the units of the variable
% writeStringMatrix(file, toWrite);
% toWrite = num2str([(1:nump1Nodes)', r]);
% writeStringMatrix(file, toWrite);
% toWrite = 'pressure, None';        % None refers to the units of the variable
% writeStringMatrix(file, toWrite);
% toWrite = num2str([(1:nump1Nodes)', p]);
% writeStringMatrix(file, toWrite);

toWrite = ['1   ', num2str(size(solutionMatrix,2))];
writeStringMatrix(file, toWrite);
toWrite = 'pointField, None';        % None refers to the units of the variable
writeStringMatrix(file, toWrite);
toWrite = num2str([(1:nump1Nodes)', solutionMatrix]);
writeStringMatrix(file, toWrite);

fclose(file);

end


function writeStringMatrix(file, toWrite)

numLines = size(toWrite,1);

for i=1:numLines
    str = toWrite(i,:);
    fprintf(file, '%s\n', str);
end

end


function [mesh_p1,UCG_p1] = getP1meshAndUCG(mesh, UDG, projectionFlag)

% UDG: [npv, numSolutionFields, numElements]
% UCG_p1: [mesh_p1.np, numSolutionFields]

% NOTE:
% - This function should work for periodic boundaries

if nargin < 3; projectionFlag = 0; end

numDimensions = mesh.nd;
npv = size(mesh.dgnodes,1);
numElements = size(mesh.dgnodes,3);
numSolutionFields = size(UDG,2);

mesh_p1.p = [];
mesh_p1.t = [];
elementPoint_To_p_p1 = (1:npv*numElements);

[~, tlocal,~,~,~,~,~] = masternodes(mesh.porder, mesh.nd, mesh.elemtype, mesh.nodetype);

tmp = mesh.dgnodes(:, 1:numDimensions, :);
tmp = permute(tmp,[1,3,2]);
mesh_p1.p = reshape(tmp,[npv*numElements,numDimensions]);

mesh_p1.t = kron(ones(numElements,1), tlocal) + kron(npv*(0:numElements-1)',ones(size(tlocal)));

% Remove duplicated nodes:
snap = max(max(mesh_p1.p,[],1)-min(mesh_p1.p,[],1),[],2)*1024*eps;
[~,ix,jx] = unique(round(mesh_p1.p/snap)*snap,'rows');
mesh_p1.p = mesh_p1.p(ix,:);
elementPoint_To_p_p1 = jx(elementPoint_To_p_p1);
mesh_p1.t = jx(mesh_p1.t);
if size(mesh_p1.t,2) == 1, mesh_p1.t = mesh_p1.t'; end  % This lines ensures the function works for one element

% Remove nodes that are not contained in t:
[pix,ix,jx] = unique(mesh_p1.t);
mesh_p1.t = reshape(jx,size(mesh_p1.t));
mesh_p1.p = mesh_p1.p(pix,:);
aux = zeros(npv*numElements,1);
aux(pix) = 1:length(pix);
elementPoint_To_p_p1 = aux(elementPoint_To_p_p1);                   % npv * numElements -> p_p1 index

mesh_p1.np = size(mesh_p1.p,1);
mesh_p1.ne = size(mesh_p1.t,1);
mesh_p1.npv = size(mesh_p1.t,2);
mesh_p1.nd = numDimensions;
mesh_p1.porder = 1;
mesh_p1.elemtype = mesh.elemtype;
mesh_p1.nodetype = 0;                   % Since p=1, all node distributions are equivalent

tmp = permute(mesh_p1.t,[2,1]);
tmp = reshape(mesh_p1.p(tmp(:),:), [mesh_p1.npv,mesh_p1.ne,mesh.nd]);
mesh_p1.dgnodes = permute(tmp,[1,3,2]);


% Create UCG for p=1 mesh:
if projectionFlag == 0
    % Take average from DG to CG nodes
    
    UCG_p1 = zeros(mesh_p1.np,numSolutionFields);
    numDGnodesInCGnode = zeros(mesh_p1.np,1);

    UDG_aux = permute(UDG,[1,3,2]);
    UDG_aux = reshape(UDG_aux,[],numSolutionFields);
    for i=1:npv*numElements
        p_index = elementPoint_To_p_p1(i);
        UCG_p1(p_index,:) = UCG_p1(p_index,:) + UDG_aux(i,:);
        numDGnodesInCGnode(p_index) = numDGnodesInCGnode(p_index) + 1;
    end
    for i=1:mesh_p1.np; UCG_p1(i,:) = UCG_p1(i,:) / numDGnodesInCGnode(i); end
    
elseif projectionFlag == 1
    % Galerkin projection
    error('Galerkin projection from DG to CG solution not implemented yet.');
    
else
    error('Projection method from DG to CG solution not implemented.');
end

end


function elementTypeString = getElementTypeString(numDimensions, elementType)

if numDimensions == 1
    elementTypeString = 'line';
elseif numDimensions == 2 && elementType == 0
    elementTypeString = 'tri';
elseif numDimensions == 2 && elementType == 1
    elementTypeString = 'quad';
elseif numDimensions == 3 && elementType == 0
    elementTypeString = 'tet';
elseif numDimensions == 3 && elementType == 1
	elementTypeString = 'hex';
else
    error('Type of element not available.');
end

end
