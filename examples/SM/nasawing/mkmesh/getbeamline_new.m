function [beamline,blinedir,normal,tangent,shellcenter] = getbeamline_new(p, q2p, sh2q, ind, tol)
% GETBEAMLINE_NEW : get the points of beam line, normals, tangents, centers
%
% SYNTAX:  [beamline,normal,tangent,shellcenter] = getbeamline_new(p, q2p, sh2q, ind, tol)
%
% INPUTS:
%    p       - Nodes positions
%    q2p     - Quad to Nodes connectivity
%    sh2q    - Shell to Quads connectivity
%    ind     - Indexes of intersecting shells
%    tol     - Tolerance for node collapse
%
% OUTPUTS:
%    beamline    - positions of points on the beam line
%    blinedir    - unitary vectors in the beamline direction
%    normal      - outgoing normals of shells at beamline points
%    tangent     - tangeant of shells at beamline points
%    shellcenter - coordinates of intersecting shell centers
%
% REMARKS: the vectors (beamline,normal,tangent) are such that they form an
%          orthogonal direct basis.
%

n = length(ind);
if n==2
    % ts1 and ts2 are the points of the 2 intersecting shells
    ts1 = q2p(sh2q{ind(1)},:);
    ts2 = q2p(sh2q{ind(2)},:);
    [f1, ~, n1, n2, c1, c2] = intersectshells(p,ts1,p,ts2,tol);    
    normal = n1;
    normal(:,:,2) = n2;    
    shellcenter = [c1; c2];
elseif n==3
    % ts1, ts2 and ts3 are the points of the 3 intersecting shells
    ts1 = q2p(sh2q{ind(1)},:);
    ts2 = q2p(sh2q{ind(2)},:);
    ts3 = q2p(sh2q{ind(3)},:);
    [f1, ~, n1, n3, c1, c3] = intersectshells(p,ts1,p,ts3,tol);
    [~, ~, n2, ~, c2, ~] = intersectshells(p,ts2,p,ts3,tol);    
    normal = n1;
    normal(:,:,2) = n2;
    normal(:,:,3) = n3;
    shellcenter = [c1; c2; c3];
elseif n==4
    % ts1, ts2, ts3 and ts4 are the points of the 4 intersecting shells
    ts1 = q2p(sh2q{ind(1)},:);
    ts2 = q2p(sh2q{ind(2)},:);
    ts3 = q2p(sh2q{ind(3)},:);
    ts4 = q2p(sh2q{ind(4)},:);
    [f1, ~, n1, n4, c1, c4] = intersectshells(p,ts1,p,ts4,tol);
    [~, ~, n2, n3, c2, c3] = intersectshells(p,ts2,p,ts3,tol);    
    normal = n1;
    normal(:,:,2) = n2;
    normal(:,:,3) = n3;    
    normal(:,:,4) = n4;    
    shellcenter = [c1; c2; c3; c4];    
end

% Get the positions of the beam nodes
beamline = p([f1(:,1); f1(end,2)],:);
blinedir = zeros(size(beamline));

% Get the tangent at each beam node
m = size(beamline,1);
tangent = normal;
for i = 1:n
    for j = 1:m
        % Vector vbline giving local direction of beamline
        if j==m; j2=j-1; else; j2=j; end
        vbline = beamline(j2,:) - beamline(j2+1,:);
        vbline = vbline / norm(vbline);
        auxtan = cross(vbline,normal(j,:,i));
        blinedir(j,:) = vbline;
        vcelin = beamline(j,:)-shellcenter(i,:);
        % Check orientation of tangent : tangent has to point inwards
        if (auxtan*vcelin' < 0)
            tangent(j,:,i) =  auxtan;
        else
            tangent(j,:,i) = -auxtan;
            normal(j,:,i)  = -normal(j,:,i);
        end
    end
end

end




