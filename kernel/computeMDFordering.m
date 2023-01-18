function mesh = computeMDFordering(mesh, Hg)

nrows = mesh.cbsr_nrows;
nblks = mesh.cbsr_nblks;
maxcols = mesh.maxBlocksPerRow;
DeltaC = zeros(maxcols,maxcols);
rowpts    = mesh.cbsr_rowpts;
colind    = mesh.cbsr_colind;
w     = zeros(nrows,1);
%bsz  = size(Hg,1);
C  = zeros(nblks,1);
D  = zeros(nblks,1);
for i = 1:nrows
    if(rem(i,1000) == 0); disp(i); end
    for n = (rowpts(i)+1):(rowpts(i+1))
        tm = squeeze(Hg(:,:,rowpts(i)+1))\squeeze(Hg(:,:,n));
        C(n) = norm(tm(:));
        D(n) = norm(squeeze(Hg(:,:,n)),'fro');
    end
end
%max(abs(C1-C))

for k = 1:nrows
    if(rem(k,1000) == 0); disp(k); end
    n = (rowpts(k)+1):rowpts(k+1); % pointers to row k 
    c = colind(n);             % connectivity indices for row k    
    DeltaC = 0*DeltaC;
    for i = 2:length(n)
        ri  = c(i);              % row ri is a neighbor of row k  
        m = (rowpts(ri)+1):rowpts(ri+1); % pointers to row ri         
        a = find(colind(m)==k);      % find the position of row k in the neighboring list of row ri        
        ik = rowpts(ri)+a;  % find the connectivity index ik
        Cik = C(ik);             % get Cik
        Dik = D(ik);
        for j = 2:length(n)
            rj  = c(j);   % rj is a neighbor of row k  
            if j~=i                    % for each neigboring row rj that is different from row ri                
                kj = rowpts(k)+j; % find the connectivity index kj                                                             
                Ckj = C(kj);           % get Cik
                if isempty(find(colind(m)==rj, 1)) %If j does not neighbor i, then increase the discarded fill in
                    DeltaC(i,j) = Cik*Ckj;
                end
            end
        end
    end
    w(k) = norm(DeltaC(:));    
% % % % % %     w(k) = norm(inv(Hg(:,:,rowpts(k)+1)),'fro');
end

%max(abs(w1-w))

bign = 1e100;
p = zeros(nrows,1);
for l = 1:nrows
    if(rem(l,1000) == 0); disp(l); end
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
                Cik = C(ik);             % get Cik
                Dik = D(ik);
                for j = 2:length(ck)
                    rj  = ck(j); % rj is a neighbor of row k  
                    if i~=j              % for each neigboring row rj that is different from row ri                                
                        kj = rowpts(k)+nk(j); 
                        Ckj = C(kj);
                        if isempty(find(colind(m)==rj, 1)) % If j does not neighbor i, then increase the discarded fill in
                            DeltaC(i,j) = Cik*Ckj;
                        end
                    end
                end
            end
            w(k) = norm(DeltaC(:));            
        end
    end
end

% bign = 1e100;
% p = zeros(nrows,1);
% for rl = 1:nrows
%     rl
%     [~,pl] = min(w);  
%     p(rl) = pl;
%     w(pl) = bign;                % do not choose pl again        
%     nrl = rowpts(p(rl)+1)-rowpts(p(rl));    
%     for k = 2:nrl          % for each neighbor of row pl
%         plk = rowpts(rl)+k;
%         rk = colind(plk);                % row k is a neighbor to row pl          
%         if w(rk) < bign      
%             w(rk) = 0;
%             nrk = rowpts(rk+1)-rowpts(k); % pointers to row k             
%             for i = 2:nrk
%                 pki = rowpts(rk)+i;
%                 ri = colind(pki);
%                 if w(ri)<bign
%                     nri = rowpts(ri+1) - rowpts(ri);
%                     for j=2:nri
%                         pij = rowpts(ri)+j;
%                         if (colind(pij)==rk)
%                             break;
%                         end
%                     end                    
%                     ik = pij;
%                     Cik = C(ik);
%                     for j=2:nrk
%                         kj = rowpts(rk)+j;
%                         rj = colind(kj);
%                         if (j ~= i) && (w(rj)<bign)
%                             Ckj = C(kj);
%                             shared = 0;
%                             for m=2:nri
%                                 pim = rowpts(ri)+m;
%                                 if (colind(pim)==rj)
%                                     shared = 1;
%                                     break;
%                                 end
%                             end
%                             if shared==0
%                                 w(rk) = w(rk) + Cik*Cik*Ckj*Ckj;
%                             end
%                         end
%                     end
%                 end
%             end
%             w(rk) = sqrt(w(rk));            
%         end
%     end
% end

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




