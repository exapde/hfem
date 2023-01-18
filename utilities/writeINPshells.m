function writeINPshells(fileName, pos, elem, elem2sol, elemfield)
%WRITEINPSHELLS : Write a INP file format for Paraview
%   Detailed explanation goes here


% Number of points and of elements
np = size(pos,1);
ne = size(elem,1);

% Get element type
if size(elem,2)==3
    elementTypeString = 'tri';
elseif size(elem,2)==4
    elementTypeString = 'quad';
elseif size(elem,2)==8
	elementTypeString = 'hex'; 
else
    error('Type of element not available.');
end


file = fopen([fileName,'.inp'],'w');

% Write header
toWrite = num2str([np, ne, 0, ne, 0]); 
writeStringMatrix(file, toWrite);

% Save nodes Coordinates
toWrite = num2str([(1:np)', pos]);
writeStringMatrix(file, toWrite);

% Save Element Connectivities
toWrite = [num2str((1:ne)'), repmat('   ',[ne,1]), num2str(elem2sol), repmat(['   ' elementTypeString '   '],[ne,1]), num2str(elem)];
writeStringMatrix(file, toWrite);

% Writing the quads thickness
toWrite = ['1   ', num2str(1)];
writeStringMatrix(file, toWrite);
toWrite = 'Thickness, None';        % None refers to the units of the variable
writeStringMatrix(file, toWrite);
toWrite = num2str([(1:ne)', elemfield]);
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


