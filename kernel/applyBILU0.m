function v = applyBILU0(mesh, Mg, x)

nrows = mesh.cbsr_nrows;
%nblks = mesh.cbsr_nblks;
rowpts    = mesh.cbsr_rowpts;
colind    = mesh.cbsr_colind;
ordered2unordered = mesh.BILU0ordered2unordered;
unordered2ordered = mesh.BILU0unordered2ordered;

y = x;
% solve L y = x
for j = 1:nrows
    rj = ordered2unordered(j);  % row rj    
    nj = (rowpts(rj)+1):rowpts(rj+1); % pointers to row rj 
    cj = colind(nj);               % connectivity indices for row rj            
            
    % Compute y_j <- x_j - L_ji*y_i for all i<j neighboring j 
    y(:,rj) = x(:,rj);
    for i = 2:length(nj)  % loop over each neighbor of row rj
        ri = cj(i);           % row ri is a neighbor of row rj
        if unordered2ordered(ri) < j
            ji = rowpts(rj)+i;
            y(:,rj) = y(:,rj) - Mg(:,:,ji)*y(:,ri);            
        end
    end                
end

v = y;
% solve U v = y
for j = 1:nrows
    rj = ordered2unordered(nrows-j+1);  % row rj    
    nj = (rowpts(rj)+1):rowpts(rj+1); % pointers to row rj 
    cj = colind(nj);               % connectivity indices for row rj            
    jj = rowpts(rj)+1;
    
    % Compute v_j <- y_j - U_ji*y_i for all i>j neighboring j 
    v(:,rj) = y(:,rj);
    for i = 2:length(nj)  % loop over each neighbor of row rj
        ri = cj(i);           % row ri is a neighbor of row rj
        if unordered2ordered(ri) > unordered2ordered(rj)
            ji = rowpts(rj)+i;
            v(:,rj) = v(:,rj) - Mg(:,:,ji)*v(:,ri);            
        end
    end    
    
    % Compute v_j <- U_jj \ v_j 
    v(:,rj) = Mg(:,:,jj)*v(:,rj);
end




