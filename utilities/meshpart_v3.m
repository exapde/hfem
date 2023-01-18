
function [elementPartition, entityPartition, ent2ent, ent2entStart, ...
    elemWeightInProcessor, entWeightInProcessor] = meshpart_v3(mesh,numProc)

% The entity graph is partitioned using METIS. Each entity has different 
% weight based on the number of neighboring entities (this determines the
% cost of matrix-vector product and preconditioner solves)
% The element partition is determined afterwards. Each element has 
% different weight depending on whether all its faces are EDG or not.

% current directory
cdir = pwd;

% move to directory that contains METIS programs
if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end
ps=strcat(pwd,sslash);
is=find(ps==sslash);
up=0;
while ~(strcmp(ps(is(end-up-1)+1:is(end-up)-1),'hdgv1.0') || strcmp(ps(is(end-up-1)+1:is(end-up)-1),'HDGv1.0'))
    up = up+1;
end
cd(strcat(ps(1:is(end-up)),'metis'));

ent2ent = mesh.cbsr_colind;
ent2entStart = mesh.cbsr_rowpts;
if ent2entStart(1) == 0; ent2entStart = ent2entStart + 1; end

% numEntities = length(ent2entStart)-1;
if strcmp(mesh.hybrid,'hdg')
    numTotalEntities = max(mesh.t2f(:));
    numRealEntities = length(unique(mesh.t2f(:)));
elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg') || strcmp(mesh.hybrid,'hedg')
    numTotalEntities = max(mesh.elcon(:));
    numRealEntities = length(unique(mesh.elcon(:)));
end

if numRealEntities ~= numTotalEntities
    error('There are ghost (missing) entities in the mesh.');
end

numEdges = (length(ent2ent) - numTotalEntities) / 2;

if floor(numTotalEntities) ~= ceil(numTotalEntities) || floor(numEdges) ~= ceil(numEdges)
    error('Error No. 2 in meshpart_v2');
end

% Generate a temporary file to be used in METIS
% % % % fid = fopen('temp.txt','w');
% % % % fprintf(fid,'%d %d 010\n', floor(numTotalEntities), floor(numEdges));
% % % % fclose(fid);
% % % % for i=1:numTotalEntities
% % % %     neighEntities = ent2entStart(i+1) - ent2entStart(i) - 1;
% % % %     toWrite = [neighEntities; ent2ent(ent2entStart(i)+1:ent2entStart(i+1)-1)]';      % +1 to avoid self-connectivity
% % % %     dlmwrite('temp.txt', toWrite, '-append', 'delimiter', ' ','precision',10);
% % % % end
fid = fopen('temp.txt','w');
fprintf(fid,'%d %d 010\n', floor(numTotalEntities), floor(numEdges));
maxLenStr = 100000;
str = char(zeros(1,maxLenStr));
lenStr = 0;
for i=1:numTotalEntities
    neighEntities = ent2entStart(i+1) - ent2entStart(i) - 1;
    newStr = [num2str([neighEntities; ent2ent(ent2entStart(i)+1:ent2entStart(i+1)-1)]'),'\n'];
    lenNewStr = length(newStr);
    str(lenStr+1:lenStr+lenNewStr) = newStr;
    lenStr = lenStr + lenNewStr;
    if lenStr > maxLenStr
        fprintf(fid,str(1:lenStr));
        lenStr = 0;
    end
end
fprintf(fid,str);
fclose(fid);

% call gpmetis
str = ['!./gpmetis temp.txt ' num2str(numProc)];
eval(str);

% get node partitioning data
str = ['temp.txt.part.' num2str(numProc)];
entityPartition = textread(str,'%d');

% remove files
delete('temp.txt');
str = ['temp.txt.part.' num2str(numProc)];
delete(str);

% move back to current directory
cd(cdir);

% get element partitioning
[elementPartition, elemWeightInProcessor] = computeElementPartition(mesh,entityPartition,numProc);

% Compute weight of entities in processor:
entWeightInProcessor = zeros(numProc,1);
for i=1:numProc
    entitiesInProcessor = find(entityPartition == i-1);
    entWeightInProcessor(i) = sum(ent2entStart(entitiesInProcessor+1) - ent2entStart(entitiesInProcessor)) - length(entitiesInProcessor);
end

end


function [elementPartition, elemWeightInProcessor] = computeElementPartition(mesh,entityPartition,numProc)

if strcmp(mesh.hybrid,'hdg')
    [elementPartition, elemWeightInProcessor] = computeElementPartition_HDG(mesh,entityPartition,numProc);
elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg') || strcmp(mesh.hybrid,'hedg')
    [elementPartition, elemWeightInProcessor] = computeElementPartition_EDG(mesh,entityPartition,numProc);
else
    error('Hybrid flag has invalid value.');
end

end


function [elementPartition, elemWeightInProcessor] = computeElementPartition_HDG(mesh,entityPartition,numProc)

numElements = mesh.ne;

entityPartition = entityPartition + 1;

elemWeights = ones(numElements,1);      % For HDG, all elements have the same weight

for i=1:numProc
    % list of faces in subdomain i
    facesInProcessor{i} = find(entityPartition == i); 
end

elementsAlreadyAllocated = [];
% loop over each processor
for i=1:numProc
    % find all interface faces
    f2f = mesh.f2f(:,facesInProcessor{i});
    i1 = find(f2f(:)>0);
    in = ones(size(f2f));
    in(i1) = ismember(f2f(i1),facesInProcessor{i});
    [~,j1] = find(in==0);
    j1 = unique(j1);
    interfaceFaces{i} = facesInProcessor{i}(j1);
    
    % find interior faces
    interiorFaces{i} = setdiff(facesInProcessor{i},interfaceFaces{i});
    
    % find all elements that are connected to interior faces
    interiorElements{i} = mesh.f(interiorFaces{i},end-1:end);
    i1 = find(interiorElements{i}(:)>0);
    interiorElements{i} = interiorElements{i}(i1);
    interiorElements{i} = unique(interiorElements{i}(:));
    
    elementsInProcessor{i} = interiorElements{i}(:);
    elemWeightInProcessor(i) = sum(elemWeights(interiorElements{i}(:)));
    elementsAlreadyAllocated = [elementsAlreadyAllocated; interiorElements{i}(:)];
    
    % find all neighboring faces not in the processor i
    f2f = f2f(:,j1);
    i1 = find(f2f(:)>0);
    in = ones(size(f2f));
    in(i1) = ismember(f2f(i1),facesInProcessor{i});
    neighboringFaces{i} = unique(f2f(in==0));  % faces do not belong to processor i
    
    interfAndNeighFaces{i} = [interfaceFaces{i}(:); neighboringFaces{i}(:)];
    if length(interfAndNeighFaces{i}) ~= length(unique(interfAndNeighFaces{i}))
        error('Error No. 1 in computeElementPartition_HDG.m');
    end
    
    % find all other required elements (elements connected to interface faces)
    otherRequiredElements{i} = mesh.f(interfaceFaces{i},end-1:end);
    i1 = find(otherRequiredElements{i}(:)>0);
    otherRequiredElements{i} = otherRequiredElements{i}(i1);
    otherRequiredElements{i} = unique(otherRequiredElements{i}(:));
    otherRequiredElements{i} = setdiff(otherRequiredElements{i},elementsInProcessor{i});
end
elementsToBeAllocated = setdiff(1:numElements,elementsAlreadyAllocated(:)');

% Allocate elements in between two processors
for i=elementsToBeAllocated
    % Find processors that neighbor the element
    neighboringProcessors = [];
    elemWeightInNeighboringProcessors = [];

    for j=1:numProc
        if any(i == otherRequiredElements{j})
            neighboringProcessors = [neighboringProcessors; j];
            elemWeightInNeighboringProcessors = [elemWeightInNeighboringProcessors; elemWeightInProcessor(j)];
        end
    end

    % Assign element to neighboring processor that currently has least
    % weight
    [~,proc_tmp] = min(elemWeightInNeighboringProcessors);
    elementsInProcessor{neighboringProcessors(proc_tmp)} = [elementsInProcessor{neighboringProcessors(proc_tmp)}; i];
    elemWeightInProcessor(neighboringProcessors(proc_tmp)) = elemWeightInProcessor(neighboringProcessors(proc_tmp)) + elemWeights(i);
    elementsAlreadyAllocated = [elementsAlreadyAllocated; i];
end
if ~isempty(setdiff(1:numElements,elementsAlreadyAllocated(:)')); error('Some element was not allocated to any processor.'); end

% Compute elementPartition from elementsInProcessor{}
elementPartition = zeros(numElements,1);
for i=1:numProc
    for j=1:length(elementsInProcessor{i})
        elementPartition(elementsInProcessor{i}(j)) = i;
    end
end

if min(elementPartition) <= 0
    error('Error No. 2 in computeElementPartition_HDG.m');
end

elementPartition = elementPartition - 1;

end


function [elementPartition, elemWeightInProcessor] = computeElementPartition_EDG(mesh,entityPartition,numProc)
% This function assumes all elements are of the same type (e.g. triangles).
% This is due to the way the EDG or non-EDG nature of an element is
% computed (i.e. through mesh.perm)

numElements = mesh.ne;

entityPartition = entityPartition + 1;

for i=1:numProc
    % list of nodes in subdomain i
    nodesInProcessor{i} = find(entityPartition == i); 
end

% Compute weight of each element:
numEDGnodesPerElem = length(unique(mesh.perm(:)));
EDGweight = 1;      % Weight of elements with EDG faces only
nonEDGweight = 2;   % Weight of elements with at least one non-EDG face
EDGelements = [];
for elem=1:numElements
    if length(unique(mesh.elcon(:,elem))) == numEDGnodesPerElem
        EDGelements = [EDGelements, elem];
    elseif length(unique(mesh.elcon(:,elem))) < numEDGnodesPerElem
        error('The number of traced nodes detected in some element was smaller than the EDG case. This is not possible.');
    end
end
nonEDGelements = setdiff(1:numElements,EDGelements);
elemWeights = zeros(numElements,1);
if length(EDGelements) > 1; elemWeights(EDGelements) = EDGweight; end
if length(nonEDGelements) > 1; elemWeights(nonEDGelements) = nonEDGweight; end
if min(elemWeights(:)) == 0; error('Some elements was not detected as EDG or non-EDG.'); end


elementsAlreadyAllocated = [];
% loop over each processor
for i=1:numProc
    % find all interface edgnodes
    numEntInProc = length(nodesInProcessor{i});
    edgnumber = zeros(numEntInProc,1); 
    for j = 1:numEntInProc
        nj = nodesInProcessor{i}(j);
        rj = (mesh.cbsr_rowpts(nj)+1):mesh.cbsr_rowpts(nj+1);
        edg1 = mesh.cbsr_colind(rj); % list of neirghboring edg nodes
        if all(ismember(edg1,nodesInProcessor{i}))
            % edgnode indent{i}(j) is fully inside the subdomain i (all
            % neighboring edg nodes are in the subdomain i)
            edgnumber(j) = 2;
        else
            % edgnode indent{i}(j) is on the interface between two subdomains
            edgnumber(j) = 1;            
        end
    end
    interfaceNodes = nodesInProcessor{i}(edgnumber==1);
    
    % find interior edgnodes
    interiorNodes = setdiff(nodesInProcessor{i},interfaceNodes);
    
    % find all elements that are connected to interior edgnodes
    interiorElements{i} = [];
    for j = 1:length(interiorNodes)
        nj = interiorNodes(j);
        rj = (mesh.cbsr_rowent2elem(nj)+1):mesh.cbsr_rowent2elem(nj+1);
        newElem = mesh.cbsr_colent2elem(rj);
        interiorElements{i} = [interiorElements{i}; newElem(:)];
    end
    interiorElements{i} = unique(interiorElements{i});
    
    elementsInProcessor{i} = interiorElements{i}(:);
    elemWeightInProcessor(i) = sum(elemWeights(interiorElements{i}(:)));
    elementsAlreadyAllocated = [elementsAlreadyAllocated; interiorElements{i}(:)];
    
    % find all neighboring edgnodes not in the processor i
    neighboringNodes = [];
    for j = 1:length(interfaceNodes)
        nj = interfaceNodes(j);
        rj = (mesh.cbsr_rowpts(nj)+1):mesh.cbsr_rowpts(nj+1);
        edg1 = mesh.cbsr_colind(rj); % list of neirghboring edg nodes
        edg1 = unique(edg1(:));
        in = ismember(edg1,nodesInProcessor{i});
        neighboringNodes = [neighboringNodes; edg1(in==0)];
    end
    neighboringNodes = unique(neighboringNodes);        % list of all edg nodes that are not in the subdomain but are connected to edg nodes in the subdomain
    
    interfAndNeighNodes = [interfaceNodes(:); neighboringNodes(:)];
    if length(interfAndNeighNodes) ~= length(unique(interfAndNeighNodes))
        error('Error No. 1 in computeElementPartition_EDG.m');
    end
    
    % find all other required elements (elements connected to interface edgnodes)
    otherRequiredElements{i} = [];
    for j = 1:length(interfaceNodes)
        nj = interfaceNodes(j);
        rj = (mesh.cbsr_rowent2elem(nj)+1):mesh.cbsr_rowent2elem(nj+1);
        newElem = mesh.cbsr_colent2elem(rj);
        otherRequiredElements{i} = [otherRequiredElements{i}; newElem(:)];
    end
    otherRequiredElements{i} = unique(otherRequiredElements{i});
    otherRequiredElements{i} = setdiff(otherRequiredElements{i},elementsInProcessor{i});
end
elementsToBeAllocated = setdiff(1:numElements,elementsAlreadyAllocated(:)');

% Allocate elements in between two processors
for i=elementsToBeAllocated
    % Find processors that neighbor the element
    neighboringProcessors = [];
    elemWeightInNeighboringProcessors = [];

    for j=1:numProc
        if any(i == otherRequiredElements{j})
            neighboringProcessors = [neighboringProcessors; j];
            elemWeightInNeighboringProcessors = [elemWeightInNeighboringProcessors; elemWeightInProcessor(j)];
        end
    end

    % Assign element to neighboring processor that currently has least
    % weight
    [~,proc_tmp] = min(elemWeightInNeighboringProcessors);
    elementsInProcessor{neighboringProcessors(proc_tmp)} = [elementsInProcessor{neighboringProcessors(proc_tmp)}; i];
    elemWeightInProcessor(neighboringProcessors(proc_tmp)) = elemWeightInProcessor(neighboringProcessors(proc_tmp)) + elemWeights(i);
    elementsAlreadyAllocated = [elementsAlreadyAllocated; i];
end
if ~isempty(setdiff(1:numElements,elementsAlreadyAllocated(:)')); error('Some element was not allocated to any processor.'); end

% Compute elementPartition from elementsInProcessor{}
elementPartition = zeros(numElements,1);
for i=1:numProc
    for j=1:length(elementsInProcessor{i})
        elementPartition(elementsInProcessor{i}(j)) = i;
    end
end

if min(elementPartition) <= 0
    error('Error No. 2 in computeElementPartition_EDG.m');
end

elementPartition = elementPartition - 1;

end
