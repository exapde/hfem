
function [p3,t3,dg3] = QtoT( p2, q, thickness, nlayers, dg2, porder)

if nargin < 5 && nargout > 2; error('Six (6) input arguments are required and less than five (5) were provided'); end

if nargin > 4 && nargout > 2
    plc1d = masternodes(porder,1,1,1);
    np1d = length(plc1d);
    np2d = size(dg2,1);
end

dz = thickness/nlayers;

lnd = true(size(p2,1),1);
lnd(q) = 0;
np = size(p2,1)-sum(lnd);

p2(lnd,1) = Inf;
[~,I] = sort(p2(:,1));
J(I) = 1:size(p2,1);
q = J(q);       % Updated version of q without non-existent points

p2 = p2(I(1:np),:);     % Updated version of p without non-existent points. Points are ordered based on x-coordinate
t2 = reshape([q(:,1:3)'; q(:,1)'; q(:,3:4)'],3,2*size(q,1))';       % t for 2D tri mesh

% Re-order element nodes (node with smallest y-coordinate is ordered first.
% The rest are ordered counterclockwise)
ne = size(t2,1);
[~,I] = min(t2,[],2);
[~,II] = min(q,[],2);
if nargin > 4 && nargout > 2
    triNodePosition = zeros(ne,3);
    for i = 1:ne/2
        dg2Indices = reshape(1:np2d,[np1d,np1d]);
        rot90Angle = 1-II(i);
        dg2Indices = rot90(dg2Indices,rot90Angle);
        dg2(:,:,i) = dg2(dg2Indices(:),:,i);
        if II(i) == 1
            triNodePosition(2*i-1:2*i,:) = [1, 2, 3; 1, 3, 4];
%             triNodePosition(2*i-1,1) = 1;
%             triNodePosition(2*i-1,2) = 2;
%             triNodePosition(2*i-1,3) = 3;
%             triNodePosition(2*i,1) = 1;
%             triNodePosition(2*i,2) = 3;
%             triNodePosition(2*i,3) = 4;
        elseif II(i) == 2
            triNodePosition(2*i-1:2*i,:) = [4, 1, 2; 4, 2, 3];
%             triNodePosition(2*i-1,1) = 4;
%             triNodePosition(2*i-1,2) = 1;
%             triNodePosition(2*i-1,3) = 2;
%             triNodePosition(2*i,1) = 4;
%             triNodePosition(2*i,2) = 2;
%             triNodePosition(2*i,3) = 3;
        elseif II(i) == 3
            triNodePosition(2*i-1:2*i,:) = [3, 4, 1; 3, 1, 2];
%             triNodePosition(2*i-1,1) = 3;
%             triNodePosition(2*i-1,2) = 4;
%             triNodePosition(2*i-1,3) = 1;
%             triNodePosition(2*i,1) = 3;
%             triNodePosition(2*i,2) = 1;
%             triNodePosition(2*i,3) = 2;
        elseif II(i) == 4
            triNodePosition(2*i-1:2*i,:) = [2, 3, 4; 2, 4, 1];
%             triNodePosition(2*i-1,1) = 2;
%             triNodePosition(2*i-1,2) = 3;
%             triNodePosition(2*i-1,3) = 4;
%             triNodePosition(2*i,1) = 2;
%             triNodePosition(2*i,2) = 4;
%             triNodePosition(2*i,3) = 1;
        end
    end
    triNodePosition_tmp = [triNodePosition,triNodePosition];
end

t = [t2,t2];
for i = 1:ne
    t2(i,:) = [t(i,I(i)),t(i,I(i)+1),t(i,I(i)+2)];
    if nargin > 4 && nargout > 2
        triNodePosition(i,:) = [triNodePosition_tmp(i,I(i)), triNodePosition_tmp(i,I(i)+1), triNodePosition_tmp(i,I(i)+2)];
    end
end

% Create 3D mesh
p3 = zeros(np*(nlayers+1),3);
zc = ones(np,1)*(0:dz:thickness);
p3 = [repmat(p2,nlayers+1,1), zc(:)];

% First layer
t3 = zeros(3*ne*nlayers,4);

if nargin > 4 && nargout > 2
    pm = dz*repmat(plc1d',[np2d 1]);
    pm = pm(:);

    [plocalTetra,~,~,~,~,~,~] = masternodes(porder,3,0,0);
    [plocalHexa,~,~,~,~,~,~] = masternodes(porder,3,1,0);

    dg3 = zeros(size(plocalTetra,1),3,3*ne*nlayers);
end

for i = 1:ne
    iQuad = floor((i-1)/2)+1;
    in = t2(i,:);
    jn = in + np;
    t3(3*(i-1)+1,:) = [jn(1), jn(3), jn(2), in(1)];
    if in(2) < in(3)
        t3(3*(i-1)+2,:) = [in(1), in(2), in(3), jn(3)];
        t3(3*(i-1)+3,:) = [in(1), in(2), jn(3), jn(2)];
    else
        t3(3*(i-1)+2,:) = [in(1), in(2), in(3), jn(2)];
        t3(3*(i-1)+3,:) = [in(3), in(1), jn(2), jn(3)];
    end
    if nargin > 4 && nargout > 2
        ind(1,:) = [triNodePosition(i,1)+4, triNodePosition(i,3)+4, triNodePosition(i,2)+4, triNodePosition(i,1)];
        if in(2) < in(3)
            ind(2,:) = [triNodePosition(i,1), triNodePosition(i,2), triNodePosition(i,3), triNodePosition(i,3)+4];
            ind(3,:) = [triNodePosition(i,1), triNodePosition(i,2), triNodePosition(i,3)+4, triNodePosition(i,2)+4];
        else
            ind(2,:) = [triNodePosition(i,1), triNodePosition(i,2), triNodePosition(i,3), triNodePosition(i,2)+4];
            ind(3,:) = [triNodePosition(i,3), triNodePosition(i,1), triNodePosition(i,2)+4, triNodePosition(i,3)+4];
        end
        mapping = getMapping(ind,plocalTetra,plocalHexa);
        dg3dHexa = [repmat(dg2(:,1:2,iQuad),[np1d 1]) pm(:)];
        dg3(:,:,3*(i-1)+1) = dg3dHexa(mapping(:,1),:);
        dg3(:,:,3*(i-1)+2) = dg3dHexa(mapping(:,2),:);
        dg3(:,:,3*(i-1)+3) = dg3dHexa(mapping(:,3),:);
    end
end

for i = 2:nlayers
    t3(3*ne*(i-1)+1:3*ne*i,:) = t3(3*ne*(i-2)+1:3*ne*(i-1),:) + np;
    if nargin > 4 && nargout > 2
        dg3(:,1:3,3*ne*(i-1)+1:3*ne*i) = [dg3(:,1:2,3*ne*(i-2)+1:3*ne*(i-1)) dg3(:,3,3*ne*(i-2)+1:3*ne*(i-1))+dz];
    end
end

end

function mapping = getMapping(ind,plocalTetra,plocalHexa)

p(1,:) = [0, 0, 0];
p(2,:) = [1, 0, 0];
p(3,:) = [1, 1, 0];
p(4,:) = [0, 1, 0];
p(5,:) = [0, 0, 1];
p(6,:) = [1, 0, 1];
p(7,:) = [1, 1, 1];
p(8,:) = [0, 1, 1];

pMasterTetras(:,:,1) = [p(ind(1,1),:);p(ind(1,2),:);p(ind(1,3),:);p(ind(1,4),:)];
pMasterTetras(:,:,2) = [p(ind(2,1),:);p(ind(2,2),:);p(ind(2,3),:);p(ind(2,4),:)];
pMasterTetras(:,:,3) = [p(ind(3,1),:);p(ind(3,2),:);p(ind(3,3),:);p(ind(3,4),:)];

mapping = zeros(size(plocalTetra,1),3);
for j=1:3
    changeBasisMatrix = [(pMasterTetras(2,:,j)-pMasterTetras(1,:,j))' , (pMasterTetras(3,:,j)-pMasterTetras(1,:,j))', (pMasterTetras(4,:,j)-pMasterTetras(1,:,j))'];
    plocalTetra_j = changeBasisMatrix * plocalTetra' + repmat(reshape(pMasterTetras(1,:,j),[],1),[1,size(plocalTetra,1)]);
    plocalTetra_j = plocalTetra_j';
    for i=1:size(plocalTetra,1)
        aux = find(sqrt(sum((repmat(plocalTetra_j(i,:),[size(plocalHexa,1),1]) - plocalHexa).^2, 2)) < 1.0e-8);
        if length(aux) ~= 1; error('Node from master tetrahedra to master hexahedra was not mapped properly.'); end
        mapping(i,j) = aux;
    end
end

end
