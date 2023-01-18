function [ptbeam,hexa,eb,pp] = mkbeam_new(p, q2p, pp, sh2q, sh2th, ind, tol)
% MKBEAM_NEW : get the points of beam line, normals, tangents, centers
%
% SYNTAX:  [ptbeam,normal,tangent,shellcenter] = mkbeam_new(p, q2p, pp, sh2q, sh2th, ind, tol)
%
% INPUTS:
%    p       - Surface nodes positions
%    q2p     - Quad to Nodes connectivity
%    pp      - Volumetric nodes positions (from shell extrusion)
%    sh2q    - Shell to Quads connectivity
%    sh2th   - Shell to Thickness information
%    ind     - Indexes of intersecting shells
%    tol     - Tolerance for node collapse
%
% OUTPUTS:
%    ptbeam      - Nodes positions of the new 3D beam
%    hexa        - Hexa to Nodes connectivity
%    eb          - edges informations
%    pp          - Volumetric nodes positions (modified)
%

% Get the 'beamline' for the current given intersection, including nodes
% positions on the beamline, normals and tangent vectors at these nodes.
[beamline,blinedir,normal,tangent] = getbeamline_new(p, q2p, sh2q, ind, tol);

[np, nd, ns] = size(normal);
n = size(beamline,1);
joinp = zeros(nd, 4, np);
joint = zeros(np, 4);
edge = zeros(nd, 2, ns, np);
for i = 1:n
    % Rotate such that the vector given by the beamline is the third direction
    % Rotation matrix to the local coords (normal,tangent,blinedir)
    rotmat = [normal(i,:,1); tangent(i,:,1); blinedir(i,:)];
    % Rotation of vectors to the local coordinates
    loctan = rotmat*squeeze(tangent(i,:,:));
    % Creation of the joint in the local plane
    [pi, joint(i,:), ei] = mkjoint([0;0], loctan(1:2,:), sh2th(ind));
    % Rotation of joint's positions back the global coordinates
    joinploc(1:2,:) = pi;
    joinploc(3,:)   = 0;
    joinp(:,:,i) = rotmat'*joinploc;
    joinp(1,:,i) = joinp(1,:,i) + beamline(i,1);
    joinp(2,:,i) = joinp(2,:,i) + beamline(i,2);
    joinp(3,:,i) = joinp(3,:,i) + beamline(i,3);
    % Connectivities
    joint(i,:) = joint(i,:) + (i-1)*4;
    % Edges (to correct the plates position nodes)
    edgetmp(1:2,:,:) = ei;
    edgetmp(3,:,:)   = 0 ;
    edgetmp = reshape(edgetmp,[nd, 2*ns]);
    edgetmp = rotmat'*edgetmp;
    edgetmp = reshape(edgetmp,[nd, 2, ns]);
    edge(1,:,:,i) = edgetmp(1,:,:) + beamline(i,1);
    edge(2,:,:,i) = edgetmp(2,:,:) + beamline(i,2);
    edge(3,:,:,i) = edgetmp(3,:,:) + beamline(i,3);
end

ptbeam = reshape(permute(joinp, [2 3 1]),[4*np nd]);
hexa = zeros(np-1,8);
for i = 1:np-1
    hexa(i,:) = [joint(i,:) joint(i+1,:)];
end
eb = reshape(permute(edge, [2 4 1 3]),[2*np nd ns]);

for i = 1:ns
    pp{ind(i)} = fixp(pp{ind(i)}, eb(:,:,i));
end   


end



