function mesh = computeExactMDFordering(mesh, Hg)

nrows = mesh.cbsr_nrows;
nblks = mesh.cbsr_nblks;
maxcols = mesh.maxBlocksPerRow;
DeltaC = zeros(maxcols,maxcols);
rowpts    = mesh.cbsr_rowpts;
colind    = mesh.cbsr_colind;
w     = zeros(nrows,1);
%bsz  = size(Hg,1);

for k = 1:nrows
    n = (rowpts(k)+1):rowpts(k+1); % pointers to row k 
    c = colind(n);             % connectivity indices for row k    
    DeltaC = 0*DeltaC;
    for i = 2:length(n)
        ri  = c(i);              % row ri is a neighbor of row k  
        m = (rowpts(ri)+1):rowpts(ri+1); % pointers to row ri         
        a = find(colind(m)==k);      % find the position of row k in the neighboring list of row ri        
        ik = rowpts(ri)+a;  % find the connectivity index ik
        for j = 2:length(n)
            rj  = c(j);
            if j~=i                    % for each neigboring row rj that is different from row ri                
                kj = rowpts(k)+j; % find the connectivity index kj                                                             
                if length(find(colind(m)==rj)) == 0
                    DeltaC(i,j) = norm(squeeze(Hg(:,:,ik))*(squeeze(Hg(:,:,rowpts(k)+1))\squeeze(Hg(:,:,kj))),'fro');
                end
            end
        end
    end
    w(k) = norm(DeltaC(:));    
end

bign = 1e100;
p = zeros(nrows,1);
for l = 1:nrows
    [~,p(l)] = min(w);  
    w(p(l)) = bign;                % do not choose pl again
    rpl = (rowpts(p(l))+1):rowpts(p(l)+1); % pointers to row pl
    cpl = colind(rpl);                 % connectivity indices for row pl
    for e = 2:length(rpl)          % for each neighbor of row pl
        k = cpl(e);                % row k is a neighbor to row pl  
        % if row k is not numbered, then recompute the weight for row k
        if w(k) < bign                         
            n = (rowpts(k)+1):rowpts(k+1); % pointers to row k 
            c = colind(n);         % connectivity indices for row k
            nk = find(w(c)<bign);  % find neighbors that are not numbered
            ck = c(nk);            % unnumbered connectivity indices for row k
            DeltaC = 0*DeltaC;
            for i = 2:length(ck)
                ri  = ck(i);             % neighboring row ri to row k  
                m = (rowpts(ri)+1):rowpts(ri+1); % pointers to row ri         
                a = find(colind(m)==k);      % find the position of row k in the neighboring list of row ri                                
                ik = rowpts(ri)+a;  % find the connectivity index ik
                for j = 2:length(ck)
                    rj  = ck(j);
                    if i~=j              % for each neigboring row rj that is different from row ri                                
                        kj = rowpts(k)+nk(j); 
                        if length(find(colind(m)==rj)) == 0
                            DeltaC(i,j) = norm(squeeze(Hg(:,:,ik))*(squeeze(Hg(:,:,rowpts(k)+1))\squeeze(Hg(:,:,kj))),'fro');
                        end
                    end
                end
            end
            w(k) = norm(DeltaC(:));            
        end
    end
end

t = 0*p;
for i = 1:length(t)
    t(p(i)) = i;
end

mesh.BILU0ordered2unordered = p; % if p[i] = k then face k is pivoted at iteration i
mesh.BILU0unordered2ordered = t; % if t[k] = i then face k is pivoted at iteration i
% mesh.BILU0ordered2unordered = 1:mesh.nf;
% mesh.BILU0unordered2ordered = 1:mesh.nf;

% function i = indexmap(rowpts,k,j)    
% i = rowpts(k)+j;
% n = (rowpts(k)+1):rowpts(k+1); 
% i = n(j);




