
function bouDist = getDistance2Boundary(mesh, bouIndex)

npv = size(mesh.dgnodes,1);
npf = size(mesh.perm,1);
nd = mesh.nd;
ne = mesh.ne;

% Find faces and elements near boundary:
facesOnBou = find(mesh.f(:,end) == - bouIndex);
[elemNearBou, localIndex] = find(ismember(mesh.t2f,facesOnBou));

% Find DG nodes on boundary of interest:
DGnodesOnBou = zeros(npf*length(elemNearBou), nd);
numDGnodesOnBou = npf*length(elemNearBou);
for i=1:length(elemNearBou)
    DGnodesOnBou((i-1)*npf+1:i*npf,1:nd) = mesh.dgnodes(mesh.perm(:,localIndex(i)), 1:nd, elemNearBou(i));
end

% Get smallest distance to the boundary:
dgNodes_tmp = reshape(permute(mesh.dgnodes(:,1:nd,:),[1,3,2]), [npv*ne, nd]);
bouDist = 1e10*ones(npv*ne, 1);
for i=1:numDGnodesOnBou
    dg2pOnBou = dgNodes_tmp - repmat(DGnodesOnBou(i,:),[npv*ne,1]);
    dist2pOnBou = sum(dg2pOnBou.^2,2);
    ind = dist2pOnBou < bouDist;
    bouDist(ind) = dist2pOnBou(ind);
end
bouDist = sqrt(reshape(bouDist,[npv,ne]));

end
