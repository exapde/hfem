function [ corner ] = mkCorners( quads,edge2p,e2b,np )
%mkCorners Summary of this function goes here
%   Detailed explanation goes here

% Number of and Quads, vertices, edges
ne = size(quads,1);
nv = size(quads,2);
ned= size(edge2p,1);

% Init the "point" structures
p2edge = cell(np,1);

% Creating p2edge structure
for i=1:ned
    p2edge{edge2p(i,1)} = [p2edge{edge2p(i,1)}; i];
    p2edge{edge2p(i,2)} = [p2edge{edge2p(i,2)}; i];
end

% Creation of the Corner structure. By definition, a "corner" is a vertex
% whose adjacent edges belong to different beams.
corner = []; j=1;
for i=1:np
    beams = e2b(p2edge{i});
    ind = beams==0 ;
    beams(ind) = [];
    p2edge{i}(ind) = [];
    if length(unique(beams))>1
        corner{j}.beams = beams;
        corner{j}.edges = edge2p(p2edge{i},:);
        corner{j}.node  = i;
        j = j+1;
    end
end

end

