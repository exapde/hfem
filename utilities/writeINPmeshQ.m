function writeINPmeshQ( fileName, master, mesh)
%WRITEINPMESHQ Summary of this function goes here
%   Detailed explanation goes here

defaultMaterialId = 0;
elementTypeString = getElementTypeString(mesh.nd, mesh.elemtype);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           WRITE DEFORMED MESH              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file = fopen([fileName,'_meshQual.inp'],'w');

% Write header
toWrite = num2str([mesh.np, mesh.ne, 0, mesh.ne, 0]); 
writeStringMatrix(file, toWrite);

% Compute final positions of vertices (via simple averaging at vertices)
vertpos = zeros(size(mesh.p));
[vertoccur,~]=hist(mesh.t(:),unique(mesh.t(:)));

for ii=1:mesh.ne
    vertpos(mesh.t(ii,:),:) = vertpos(mesh.t(ii,:),:) + squeeze(mesh.dgnodes(mesh.permnode,1:3,ii));
end
vertpos(:,1) = vertpos(:,1) ./ vertoccur(:);
vertpos(:,2) = vertpos(:,2) ./ vertoccur(:);
vertpos(:,3) = vertpos(:,3) ./ vertoccur(:);

% Save nodes Coordinates
toWrite = num2str([(1:mesh.np)', vertpos(:,:)]);
%toWrite = num2str([(1:mesh.np)', mesh.p(:,:)]); % <-- Write undef vertexs
writeStringMatrix(file, toWrite);

% Save Element Connectivities
toWrite = [num2str((1:mesh.ne)'), repmat('   ',[mesh.ne,1]), num2str(repmat(defaultMaterialId,[mesh.ne,1])), repmat(['   ' elementTypeString '   '],[mesh.ne,1]), num2str(mesh.t(:,:))];
writeStringMatrix(file, toWrite);

% Getting the mesh quality
mesh_qual = quality_mesh(mesh, master);
toWrite = ['1   ', num2str(1)];
writeStringMatrix(file, toWrite);
toWrite = 'Mesh Distortion, None';        % None refers to the units of the variable
writeStringMatrix(file, toWrite);
toWrite = num2str([(1:mesh.ne)', mesh_qual]);
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

