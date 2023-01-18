function [ptCell, tCell, pp, pb] = mkCornerCell(p, pp, tp, pb, tb, ib, corner, sh2th)
% MKCORNERCELL : creates the hexa Cell at the junction of different beams
%
% SYNTAX:  [ptCell, tCell, pp, pb] = mkCornerCell(p, pp, tp, pb, tb, ib, corner, tol)
%
% INPUTS:
%    p       - Surface nodes positions
%    pp      - Volumetric nodes positions from shell extrusions
%    tp      - Connectivities for the volumetric shells
%    pb      - Volumetric nodes positions from beam extrusions
%    tb      - Connectivities for the volumetric beams
%    ib      - Sets of intersecting shells
%    corner  - Current corner structure
%    sh2th   - Shell to thickness information
%
% OUTPUTS:
%    ptCell  - Nodes positions of the new corner Cell
%    tCell   - Nodes connectivity for the corner Hexa
%    pp      - Volumetric shell nodes positions (modified)
%    pb      - Volumetric beams nodes positions (modified)
%

nd = 3;
% Cell center position
pc = p(corner.node,:);
% Cell element connectivities
tCell = [1:8];
% Number of neighbouring beams
nb = length(corner.beams);


% Get the outgoing normals of the neighbouring beams the the corner
% These normals are obtained from the edges connected to the corner
normals = zeros(nb,nd);
EdNodes = zeros(nb,1);
for i=1:nb
    % Get the neighbouring edge node
    auxVert    = corner.edges(i,:);
    EdNodes(i) = auxVert(auxVert~=corner.node);
    % Build the normal
    normals(i,:) = pc - p(EdNodes(i),:);
    normals(i,:) = normals(i,:)/norm(normals(i,:));
end

% Set the vertical beam as the first beam
[~,ind] = sort(abs(normals(:,3)),'descend');
corner.beams = corner.beams(ind);
normals(:,:) = normals(ind,:);
corner.edges = corner.edges(ind,:);
EdNodes = EdNodes(ind);
ib1 = corner.beams(1);


% List of neighbouring shells
shlist = [];
for i=1:nb
    shlist = [shlist; ib{corner.beams(i)}];
end
shlist = unique(shlist);
nsh = length(shlist);
% Creation of shell to beam connectivities
sh2beam = zeros(nsh,2);
for i=1:nb
    ibe = corner.beams(i);
    for j=1:length(ib{ibe})
        ind = find(shlist==ib{ibe}(j));
        if sh2beam(ind,1)==0
            sh2beam(ind,1) = i;
        else
            sh2beam(ind,2) = i;
        end
    end
end

% Creating the nodes redirection matrices for shells and beams
nShellTarget = zeros(nsh,2);
nShellReloc  = zeros(nsh,2);
nBeamTarget  = zeros(nb, 4);
nBeamReloc   = zeros(nb, 4);

% Get the collapsed nodes for the neighbouring beams
for i=1:nb
    ibe = corner.beams(i);
    % Get the distances from center beam hexes to corner
    HexCentDist = elemCenterDist(pb{ibe},tb{ibe},pc);
    % Get the nearest hex to corner
    [~,ind] = min(HexCentDist);
    % Get the 4 nearest hex nodes to corner
    indNode = findNearestNodes(pb{ibe}(tb{ibe}(ind,:),:),pc,4);
    auxNode = tb{ibe}(ind,indNode);
    % Update table of relocated beam nodes
    nBeamReloc(i,:) = auxNode;
    % Reorder nodes turning around the normal in a direct way
    indNode = checkOrderingQuad(pb{ibe}(auxNode,:),normals(i,:));
    nBeamReloc(i,:) = nBeamReloc(i,indNode);
    if i==1
        % Get Target Nodes (simple reordering for 1st beam)
        nBeamTarget(1,:) = [1:4];
        % Init positions of the Cell corner as a flat hex
        ptCell = zeros(8,3);
        ptCell(1:4,:) = pb{ibe}(nBeamReloc(1,:),:);
        ptCell(5:8,:) = pb{ibe}(nBeamReloc(1,:),:);
    else
        % Find beam target nodes
        nBeamTarget(i,:) = findBeamNodesTarget(ptCell,pb{ibe},normals(i,:),normals(1,:),nBeamReloc(i,:),p(EdNodes(i),:));
    end
end

% Get the collapsed nodes for the neighbouring shells
k=0; avThickness=0;
for i=1:nsh
    ish = shlist(i);
    % Get the distances from center shell hexes to corner
    HexCentDist = elemCenterDist(pp{ish},tp{ish},pc);
    % Get the nearest hex to corner
    [~,ind] = min(HexCentDist);
    % Get the 2 nearest hex nodes to corner
    indNode = findNearestNodes(pp{ish}(tp{ish}(ind,:),:),pc,2);
    auxNode = tp{ish}(ind,indNode);
    % Update table of relocated shell nodes
    nShellReloc(i,:) = auxNode;
    
    if (ismember(1,sh2beam(i,:)))
        ind = sh2beam(i,:)~=1;
        ibe = sh2beam(i,ind);
        ind = nBeamTarget(ibe,:)<5;
        indNode = nBeamTarget(ibe,ind);
        if dot(ptCell(indNode(2),:)-ptCell(indNode(1),:),pp{ish}(nShellReloc(i,2),:)-pp{ish}(nShellReloc(i,1),:))<0
            nShellTarget(i,1) = indNode(2);
            nShellTarget(i,2) = indNode(1);
        else
            nShellTarget(i,:) = indNode;
        end
    else % Shells not connected with 1st beam
        % Get common target nodes shared by the 2 adjacent beams
        indNode = intersect(nBeamTarget(sh2beam(i,1),:),nBeamTarget(sh2beam(i,2),:));
        if length(indNode)<2
            fprintf('problem intersection\n');
            nShellTarget(i,:) = hack_shouldbeprism(ib1,normals(1,3));
        elseif (indNode(2)-indNode(1))*(pp{ish}(nShellReloc(i,2),3)-pp{ish}(nShellReloc(i,1),3))*normals(1,3)<0
            nShellTarget(i,1) = indNode(2);
            nShellTarget(i,2) = indNode(1);
        else
            nShellTarget(i,:) = indNode;
        end  
        % Sum thicknesses for averaging
        avThickness = avThickness + sh2th(ish);
        k = k+1;
    end
end

% Average Thickness
avThickness = avThickness/k;

% Move the horizontal faces of the corner along the 1st beam normal
ptCell(1:4,:) = ptCell(1:4,:) - 0.5*avThickness*ones(4,1)*normals(1,:);
ptCell(5:8,:) = ptCell(5:8,:) + 0.5*avThickness*ones(4,1)*normals(1,:);


% Assign shells and beams collapsed nodes to the new positions
for i=1:nb
    pb{corner.beams(i)}(nBeamReloc(i,:),:) = ptCell(nBeamTarget(i,:),:);
end
for i=1:nsh
    pp{shlist(i)}(nShellReloc(i,:),:) = ptCell(nShellTarget(i,:),:);
end

end




function [ind] = checkOrderingQuad(p,normal)
% This function checks the ordering of the nodes relatively to the normal
% direction in order to create a quad that is untangled and with right
% nodes ordering

% Selection of point n2
for n2=2:4
    if n2==2
        nA = 3; nB = 4;
    elseif n2==3
        nA = 4; nB = 2;
    else %(n2==4)
        nA = 2; nB = 3;
    end
    Vtest = p(n2,:) - p(1,:);
    Vaux  = cross(normal,Vtest);
    auxA = dot(Vaux,p(nA,:)-p(1,:));
    auxB = dot(Vaux,p(nB,:)-p(1,:));
    if (auxA>0 && auxB>0)
        break
    end
end
% Selection of points n3 and n4
Vaux = p(nB,:) - p(nA,:);
V12  = p(n2,:) - p(1, :);
if (dot(Vaux,V12)>0)
    n3 = nB; n4 = nA;
else
    n3 = nA; n4 = nB;
end
ind = [1,n2,n3,n4];

end


function [indTarget] = findBeamNodesTarget(pCell,pb,n,n1,BeamReloc,pEd)
% Find the local beam taget nodes using lines intersections
dists = zeros(4,1);
n2D = n(1:2)/norm(n(1:2));
for i=1:4
    Vedge = pCell(mod(i,4)+1,1:2)-pCell(i,1:2);
    Vedge = Vedge/norm(Vedge);
    paux  = lineintersection(pEd(1:2), n2D,  pCell(i,1:2), Vedge);
    if isempty(paux)
        dists(i) = 1000.;
    else
        % Check if the intersection is in the segment
        if dot(paux'-pCell(mod(i,4)+1,1:2),paux'-pCell(i,1:2))<0.
            dists(i) = norm(paux'-pEd(1:2));
        else % in that case, the intersection is outside the segment
            dists(i) = 1000.;
        end
    end
end
[~,nEdge] = min(dists);

% Check / fix quad ordering : 
indNode = checkOrderingQuad(pb(BeamReloc,:),n);
if ~isequal(indNode,[1:4])
    error('Error. IndNode should be [1 2 3 4] %s.');
end

% Find the two lowest nodes in local n1-repair
Zp = n1(3)*pb(BeamReloc(indNode),3); 

switch(nEdge)
    case 1
        if(max(Zp(1),Zp(2))<min(Zp(3),Zp(4)))
            indTarget = [2 1 5 6];
        elseif(max(Zp(1),Zp(4))<min(Zp(2),Zp(3)))
            indTarget = [1 5 6 2];
        elseif(max(Zp(3),Zp(4))<min(Zp(2),Zp(1)))
            indTarget = [5 6 2 1];
        elseif(max(Zp(3),Zp(2))<min(Zp(4),Zp(1)))
            indTarget = [6 2 1 5]; 
        end
    case 2
        if(max(Zp(1),Zp(2))<min(Zp(3),Zp(4)))
            indTarget = [3 2 6 7];
        elseif(max(Zp(1),Zp(4))<min(Zp(2),Zp(3)))
            indTarget = [2 6 7 3];
        elseif(max(Zp(3),Zp(4))<min(Zp(2),Zp(1)))
            indTarget = [6 7 3 2];
        elseif(max(Zp(3),Zp(2))<min(Zp(4),Zp(1)))
            indTarget = [7 3 2 6]; 
        end
    case 3
        if(max(Zp(1),Zp(2))<min(Zp(3),Zp(4)))
            indTarget = [4 3 7 8];
        elseif(max(Zp(1),Zp(4))<min(Zp(2),Zp(3)))
            indTarget = [3 7 8 4];
        elseif(max(Zp(3),Zp(4))<min(Zp(2),Zp(1)))
            indTarget = [7 8 4 3];
        elseif(max(Zp(3),Zp(2))<min(Zp(4),Zp(1)))
            indTarget = [8 4 3 7]; 
        end
    case 4
       if(max(Zp(1),Zp(2))<min(Zp(3),Zp(4)))
            indTarget = [1 4 8 5];
        elseif(max(Zp(1),Zp(4))<min(Zp(2),Zp(3)))
            indTarget = [4 8 5 1];
        elseif(max(Zp(3),Zp(4))<min(Zp(2),Zp(1)))
            indTarget = [8 5 1 4];
        elseif(max(Zp(3),Zp(2))<min(Zp(4),Zp(1)))
            indTarget = [5 1 4 8]; 
       end
end
end


function [nearNodes] = findNearestNodes(p,pc,nn)
% This function finds the nn nearest nodes from point pc
dist = (p(:,1) - pc(1)).^2 + (p(:,2) - pc(2)).^2 + (p(:,3) - pc(3)).^2;
[~,ind] = sort(dist);
nearNodes = ind(1:nn);
end


function [centerDist] = elemCenterDist(p,t,pc)

centerDist = zeros(size(t,1),1);
for i=1:size(t,1)
    CenterHex = mean(p(t(i,:),:));
    centerDist(i) = norm(CenterHex-pc);
end
end

function [ind] = hack_shouldbeprism(ib1,n3)

if ib1==1
    if n3>0
        ind = [1,5];
    else
        ind = [5,1];
    end
elseif ib1==2
    if n3>0
        ind = [5,1];
    else
        ind = [1,5];
    end
elseif ib1==3
    if n3>0
        ind = [1,5];
    else
        ind = [5,1];
    end
elseif ib1==4
    if n3>0
        ind = [2,6];
    else
        ind = [4,8];
    end
else
    error('Case not taken into account...');
end

end
