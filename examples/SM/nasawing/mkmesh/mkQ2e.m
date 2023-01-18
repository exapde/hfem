function [ edge2p, edge2quad, quad2edge ] = mkQ2e( quads )
%MKQ2E Summary of this function goes here
%   Detailed explanation goes here

% Number of points and Quads
np = max(quads(:));
ne = size(quads,1);
nv = size(quads,2);

edge2p   = zeros(ne*nv,2);
eind     = zeros(ne*nv,1);

for i=1:ne
    for j=1:nv
        np1 = quads(i,j);
        np2 = quads(i,mod(j,nv)+1);
        edge2p((i-1)*nv+j,1:2) = [np1 np2];
        % Vector of indices to singularize edges
        % assigning an unique ID to the edges
        eind((i-1)*nv+j) = min(np1,np2)*np+max(np1,np2);
    end
end

% Eliminate multiple edges
[~, ia, ic] = unique(eind);
edge2p = edge2p(ia,:);

% Nb of edges :
ned = size(edge2p,1);

% Creation connectivities Quad2Edges
quad2edge = 1:ne*nv;
quad2edge = quad2edge(ic);
quad2edge = reshape(quad2edge,[nv ne])';

% Creation connectivities Edges2Quad
%nConMax = 5; % edge2quad = sparse(ne,nConMax);
edge2quad = cell(ned,1);
for i=1:ne
    for j=1:nv
        ie = quad2edge(i,j);
        edge2quad{ie} = [edge2quad{ie}; i];
    end
end


