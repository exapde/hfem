function intent = mkintent(facecon,nin,face2cpu,my_rank)

ent = facecon(:,1:nin);
ent = unique(ent(:));
ent = ent(ent>0);
nent = length(ent);

intent = zeros(nent,1);
%ent2cpu = intent;
for i = 1:nent
    enti = ent(i);
    [~,ci]=find(facecon == enti);
    cpus = face2cpu(ci);
    a = unique(cpus);
    b = 0*a;
    for j=1:length(b)
        b(j) = length(find(cpus==a(j)));
    end    
    m = b==max(b);
    nbcpu = a(m);    
    icpu = nbcpu(1); % this cpus own enti   
    %ent2cpu(i) = icpu;
    if icpu == my_rank
        intent(i) = 1;
    end        
end
intent = ent(intent==1);


function intent = mkintent2(facecon,extintface,nin,face2cpu,my_rank)

ent = facecon(:,1:nin);
ent = unique(ent(:));
ent = ent(ent>0);
nent = length(ent);

[rowent2face,colent2face,ent2ind] = facecon2entcon(facecon,extintface);

intent = zeros(nent,1);
for i = 1:nent
    enti = ent(i);
    ij = ent2ind(enti);
    rj = (rowent2face(ij)+1):rowent2face(ij+1);
    faces = colent2face(rj);
    
%     enti 
%     faces
%     [ri,ci]=find(facecon == enti);
%     extintface(ci)
%     pause
    
    cpus = face2cpu(faces);    
    a = unique(cpus);
    b = 0*a;
    for j=1:length(b)
        b(j) = length(find(cpus==a(j)));
    end    
    m = find(b==max(b));
    nbcpu = a(m);
    icpu = nbcpu(1); % this cpus own enti   
    if icpu == my_rank
        intent(i) = enti;
    end    
end
intent = intent(intent>0);

function [rowent2elem,colent2elem,ent2ind] = facecon2entcon(elcon,elem)

[~,ne] = size(elcon);
ent = unique(elcon(:));
ent = ent(ent>0);
ndof = length(ent);
entmax = ent(end); % largest entity

% entity-to-index mapping
ent2ind = zeros(entmax,1);
ent2ind(ent) = (1:ndof);

 % store number of neighboring elements for each entity
rowent2elem = zeros(ndof,1);
for i = 1:ne  % for each element i
    elc = elcon(:,i);     % entities on element i  
    k   = unique(elc(:)); % remove duplicate entities on element i  
    k   = k(k>0);
    ind = ent2ind(k);   % get entity indices 
    rowent2elem(ind) = rowent2elem(ind) + 1;
end
rowent2elem=[0; cumsum(rowent2elem)];

 % store neighboring-element indices for each entity
colent2elem = zeros(rowent2elem(end),1); 
inc = ones(ndof,1);
for i = 1:ne
    elc = elcon(:,i);   % entities on element i  
    k = unique(elc(:)); % remove duplicate entities on element i  
    k   = k(k>0);
    ind = ent2ind(k);   % get entity indices 
    % rowent2elem(ind): pointer to the list of entities k
    colent2elem(rowent2elem(ind)+inc(ind)) = elem(i);
    inc(ind) = inc(ind) + 1;
end
