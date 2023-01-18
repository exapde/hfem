function [ nfixed ] = locate_Dirichlet_nodes(mesh, Vaxis, Paxis)
%LOCATE_DIRICHLET_NODES Locate the nodes that are on a rotation axis
%      MESH:                    Mesh structure
%
%      VAXIS(3):                Vector direction of the axis of rotation
%      PAXIS(3):                Point belonging to the rotation axis

tol = 1.e-6; nfixed = [];

% Norm of the vector Vaxis
Vnorm = norm(Vaxis);

% Loop on the HDG interfaces nodes
for i=1:size(mesh.edgnodes,1)
    % Vector from P to the current node
    PM = mesh.edgnodes(i,:) - Paxis;
    % Distance between the current node and the axis
    d  = norm(cross(PM,Vaxis))/Vnorm;
    if d<tol
        nfixed = [nfixed i];
    end
end




