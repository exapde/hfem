function [re,ce,rp,cj] = elcon2entcon(elcon)

% IMPORTANT: For matvecImplementation == 1 in our HDG code to work properly
% with the RAS preconditioner, the entities in each row of "ce"/"cj" need to be
% ordered in increasing order.

[np,ne] = size(elcon);
ndof = max(elcon(:));
ndof2 = length(unique(elcon(:)));
if ndof ~= ndof2; error('It looks like there are ghost entities in the mesh.'); end
if min(elcon(:)) < 1; error('elcon <= 0 was detected.'); end

re = zeros(ndof,1); % store number of neighboring elements for each entity
for i = 1:ne
    elc = elcon(:,i);
    elc = unique(elc(:));
    re(elc) = re(elc) + 1;
end
re_aux = re;
re=[0; cumsum(re)];

ce = zeros(re(end),1);  % store neighboring-element indices for each entity
in = ones(ndof,1);
for i = 1:ne
    k = elcon(:,i);   % entities on element i  
    k = unique(k(:)); % remove duplicate entities on element i  
    % re(k): pointer to the list of entities k
    ce(re(k)+in(k)) = i;
    in(k) = in(k) + 1;
end
if min(ce(:)) < 1; error('Entity-to-element connectivity went wrong.'); end
if any(in-1-re_aux); error('Entity-to-element connectivity went wrong.'); end

me = ceil(mean(re(2:end)-re(1:end-1)));
rp = zeros(ndof,1); % store number of neighboring entities for each entity
cj = zeros(ndof*np*me,1); % store neighboring-entity indices for each entity
m = 1;
for i = 1:ndof    
    j = (re(i)+1):re(i+1);
    e = ce(j);      % elements neighboring the entity i
    n = elcon(:,e);     % entities neighboring entity i
    k = unique(n(:));   % entities neighboring entity i
    l = length(k);
    rp(i) = l;    
    if ~ismember(i,k); error('Entity-to-entity connectivity went wrong.'); end
    cj(m:(m+l-1)) = [i; k((k~=i))];
    m = m + l;
end
rp=[0; cumsum(rp)];
in = cj==0;
cj(in) = [];
