function [rowent2elem,colent2elem,rowent2ent,colent2ent,ent,ent2ind] = elcon2entconmpi(elcon,elem)
% (rowent2elem,colent2elem) entity-to-element connectivity
% (rowent2ent,colent2ent) entity-to-entity connectivity

% IMPORTANT: For matvecImplementation == 1 in our HDG code to work properly
% with the RAS preconditioner, the entities in each row of "colent2elem"/"colent2ent" need to be
% ordered in increasing order.

[np,ne] = size(elcon);
ent = unique(elcon(:));
ndof = length(ent);
entmax = ent(end); % largest entity

% entity-to-index mapping
ent2ind = zeros(entmax,1);
ent2ind(ent) = (1:ndof);

% element-to-index mapping
elem2ind = zeros(max(elem),1);
elem2ind(elem) = (1:ne);

% ent'
% ent2ind'
% elem'
% elem2ind'
% pause
% 
 % store number of neighboring elements for each entity
rowent2elem = zeros(ndof,1);
for i = 1:ne  % for each element i
    elc = elcon(:,i);     % entities on element i  
    k   = unique(elc(:)); % remove duplicate entities on element i  
    ind = ent2ind(k);   % get entity indices 
    rowent2elem(ind) = rowent2elem(ind) + 1;
end
rowent2elem=[0; cumsum(rowent2elem)];

% rowent2elem'
% pause

 % store neighboring-element indices for each entity
colent2elem = zeros(rowent2elem(end),1); 
inc = ones(ndof,1);
for i = 1:ne
    elc = elcon(:,i);   % entities on element i  
    k = unique(elc(:)); % remove duplicate entities on element i  
    ind = ent2ind(k);   % get entity indices 
    % rowent2elem(ind): pointer to the list of entities k
    colent2elem(rowent2elem(ind)+inc(ind)) = elem(i);
    inc(ind) = inc(ind) + 1;
end
if min(colent2elem(:)) < 1; error('Entity-to-element connectivity went wrong.'); end

% colent2elem'
% pause

me = ceil(mean(rowent2elem(2:end)-rowent2elem(1:end-1)));
 % store number of neighboring entities for each entity
rowent2ent = zeros(ndof,1);
 % store neighboring-entity indices for each entity
colent2ent = zeros(ndof*np*me,1);
m = 1;
for i = 1:ndof   
    ei = ent(i);  % entity ei
    j = (rowent2elem(i)+1):rowent2elem(i+1);
    e = colent2elem(j);       % elements neighboring the entity ei
    n = elcon(:,elem2ind(e)); % entities neighboring entity ei
    k = unique(n(:));   % remove duplicate entities neighboring entity ei        
    l = length(k);    
    if ~ismember(ei,k) || (l==1) 
        error('Entity-to-entity connectivity went wrong.'); 
    end    
    rowent2ent(i) = l;        
    colent2ent(m:(m+l-1)) = [ei; k(k~=ei)];
    m = m + l;
end
rowent2ent=[0; cumsum(rowent2ent)];
in = colent2ent==0;
colent2ent(in) = [];

% rowent2ent'
% colent2ent'
% pause

% [length(rowent2elem) length(colent2elem) length(rowent2ent) length(colent2ent) length(ent2ind)]
% pause