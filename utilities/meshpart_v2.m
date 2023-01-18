
function [elementPartition,entityPartition, ent2ent, ent2entStart] = meshpart_v2(mesh,numProc,ent2entWeight,p0,t0)

% METIS computes entity partition by minimizing the edge-weight cutting.
% The element partition is determined afterwards based on random allocation
% among candidates. Random allocation is fine since no element communication is required in the code.

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

% % Compute ent2ent and ent2entStart
% [ent2ent, ent2entStart] = computeEnt2ent(mesh);
% ent2entStart = ent2entStart + 1;        % Make indices start at 1 (instead of 0, as in C)

ent2ent = mesh.cbsr_colind;
ent2entStart = mesh.cbsr_rowpts;
if ent2entStart(1) == 0; ent2entStart = ent2entStart + 1; end

if nargin < 3; ent2entWeight = ones(length(ent2ent),1); end

if min(ent2entWeight) <= 0; error('Non-positive edge weight was detected.'); end
if length(ent2ent) ~= length(ent2entWeight)
    error('Error No. 1 in meshpart_v2. The number of weights in ent2entWeights does not match the number of entity-to-entity connectivities in the mesh');
end

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

% Create unique ent2entWeight based on contributions from both blocks
ent2entWeightNew = zeros(length(ent2entWeight),1);
for i=1:numTotalEntities
    ent2entWeightNew(ent2entStart(i) + 0) = ent2entWeightNew(ent2entStart(i) + 0) + ent2entWeight(ent2entStart(i) + 0);     % Self-block
    for j=1:(ent2entStart(i+1) - ent2entStart(i) - 1)       % Cross-blocks
        adjEntity = ent2ent(ent2entStart(i) + j);
        ent2entWeightNew(ent2entStart(i) + j) = ent2entWeightNew(ent2entStart(i) + j) + ent2entWeight(ent2entStart(i) + j);
        
        detected = 0;
        for k=1:(ent2entStart(adjEntity+1) - ent2entStart(adjEntity) - 1)
            if ent2ent(ent2entStart(adjEntity) + k) == i
                ent2entWeightNew(ent2entStart(adjEntity) + k) = ent2entWeightNew(ent2entStart(adjEntity) + k) + ent2entWeight(ent2entStart(i) + j);
                detected = 1;
                break;
            end
        end
        if detected == 0; error('No dual-direction connectivity was detected.'); end
    end
end
ent2entWeight = ent2entWeightNew;

% Set zero blocks to a small value
bin = ent2entWeight == 0;
zeroBlocks = find(bin);
if ~isempty(zeroBlocks)
    nonZeroBlocks = find(~bin);
    minNonZeroBlockNorm = min(ent2entWeight(nonZeroBlocks));
    ent2entWeight(zeroBlocks) = minNonZeroBlockNorm / 10;
end

% Convert ent2entWeight to positive integer array
ent2entWeight = ceil(ent2entWeight/min(ent2entWeight));

% Generate a temporary file to be used in METIS
fid = fopen('temp.txt','w');
fprintf(fid,'%d %d 001\n', floor(numTotalEntities), floor(numEdges));
fclose(fid);
% maxNeighEntities = max(ent2entStart(2:end)-ent2entStart(1:end-1))-1;
% toWriteChar = [];
% entireToWrite = zeros(2*maxNeighEntities, numTotalEntities);
for i=1:numTotalEntities
%     neighEntities = ent2entStart(i+1) - ent2entStart(i) - 1;
%     numSpaces = 2*(maxNeighEntities - neighEntities);
    toWriteOdd = ent2ent(ent2entStart(i)+1:ent2entStart(i+1)-1);      % +1 to avoid self-connectivity
    toWriteEven = ent2entWeight(ent2entStart(i)+1:ent2entStart(i+1)-1);
    toWrite = zeros(1,2*length(toWriteOdd));
    toWrite(1:2:2*length(toWriteOdd)) = toWriteOdd;
    toWrite(2:2:2*length(toWriteEven)) = toWriteEven;
%     entireToWrite(1:2:2*length(toWriteOdd),i) = toWriteOdd;
%     entireToWrite(2:2:2*length(toWriteEven),i) = toWriteEven;
%     spaces = repmat(' ',1,numSpaces);
%     toWriteChar = [toWriteChar; char(toWrite), spaces];
    dlmwrite('temp.txt', toWrite, '-append', 'delimiter', ' ','precision',10);
end
% toMakeSpace = find(entireToWrite == num2str(0));
% entireToWrite(toMakeSpace) = ' ';

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
elementPartition = computeElementPartition(mesh,entityPartition,numProc);

% plot the partitioned mesh 
if nargin>3
    p = p0;
    t = t0;
    bcol = [1 1 0;... % yellow 
            1 0 1;... % magneta
            0 1 1;... % cyan
            1 0 0;... % red
            0 1 0;... % green
            0 0 1;... % blue
            1,0.4,0.6;...
            0.4,0.6,1;...
           ];
    figure(2); clf;
    hold on;        
    for i=0:numProc-1
        ind = find(elementPartition==i);
        ti = t(ind,:);
        simpplot(p,ti,[],bcol(i+1,:));                       
    end
%     for it=1:size(t,1)
%         pmid=mean(p(t(it,:),:),1);
%         txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
%         text(pmid(1),pmid(2),num2str(it),txtpars{:});
%     end

    hold off;
    axis equal;      
    axis tight;
    axis on;  
end

end


% function [ent2ent, ent2entStart] = computeEnt2ent(mesh)
% 
% if strcmp(mesh.hybrid,'hdg')
%     
%     ne    = size(mesh.t2f,1);
%     nfe   = size(mesh.t2f,2); 
%     il = zeros(nfe,nfe,ne);
%     jl = zeros(nfe,nfe,ne);
%     for i=1:ne    
%         con = mesh.t2f(i,:);
%         il(:,:,i) = repmat(con' , 1, nfe);
%         jl(:,:,i) = repmat(con , nfe, 1);        
%     end
%     il(:,:,1);
%     jl(:,:,1);
%     sl = ones(nfe,nfe,ne);
% 
% elseif strcmp(mesh.hybrid,'edg')
%     
%     [nn,ne] = size(mesh.elcon);
%     %ne = mesh.ne;
%     il = zeros(nn,nn,ne);
%     jl = zeros(nn,nn,ne);
%     for i=1:ne
%         con = mesh.elcon(:,i);
%         il(:,:,i) = repmat(con ,1,nn);
%         jl(:,:,i) = repmat(con',nn,1);
%     end    
%     sl = ones(nn,nn,ne);
%     
% end
% [rp, cj] = sparse_to_csr(sparse(il(:),jl(:),sl(:)));
% 
% ct = cj;
% for i = 1:length(rp)-1    
%     k = rp(i)+1:rp(i+1);    
%     ct(k) = [i; sort(setdiff(ct(k),i))]; % the first index is the block itself   
% end
% 
% if length(ct) ~= length(cj)
%     error('Matrix must have nonzero diagonal entries.');
% else
%     cj = ct;
% end
% 
% ent2entStart = rp;
% ent2ent = cj;
% 
% end


function elementPartition = computeElementPartition(mesh,entityPartition,numProc)

if strcmp(mesh.hybrid,'hdg')
    elementPartition = computeElementPartition_HDG(mesh,entityPartition,numProc);
elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg') || strcmp(mesh.hybrid,'hedg')
    elementPartition = computeElementPartition_EDG(mesh,entityPartition,numProc);
else
    error('Hybrid flag has invalid value.');
end

end


function elementPartition = computeElementPartition_HDG(mesh,entityPartition,numProc)

numElements = mesh.ne;
maxElementsPerProcessor = ceil(numElements / numProc);

entityPartition = entityPartition + 1;

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
    numElementsInNeighboringProcessors = [];

    for j=1:numProc
        if any(i == otherRequiredElements{j})
            neighboringProcessors = [neighboringProcessors; j];
            numElementsInNeighboringProcessors = [numElementsInNeighboringProcessors; length(elementsInProcessor{j})];
        end
    end

    % Assign element to neighboring processor that currently has least
    % elements
    [~,proc_tmp] = min(numElementsInNeighboringProcessors);
    elementsInProcessor{neighboringProcessors(proc_tmp)} = [elementsInProcessor{neighboringProcessors(proc_tmp)}; i];
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


function elementPartition = computeElementPartition_EDG(mesh,entityPartition,numProc)

numElements = mesh.ne;
maxElementsPerProcessor = ceil(numElements / numProc);

entityPartition = entityPartition + 1;

for i=1:numProc
    % list of nodes in subdomain i
    nodesInProcessor{i} = find(entityPartition == i); 
end

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
    numElementsInNeighboringProcessors = [];

    for j=1:numProc
        if any(i == otherRequiredElements{j})
            neighboringProcessors = [neighboringProcessors; j];
            numElementsInNeighboringProcessors = [numElementsInNeighboringProcessors; length(elementsInProcessor{j})];
        end
    end

    % Assign element to neighboring processor that currently has least
    % elements
    [~,proc_tmp] = min(numElementsInNeighboringProcessors);
    elementsInProcessor{neighboringProcessors(proc_tmp)} = [elementsInProcessor{neighboringProcessors(proc_tmp)}; i];
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
