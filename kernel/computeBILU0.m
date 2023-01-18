function Mg = computeBILU0(mesh, Hg)

nrows = mesh.cbsr_nrows;
%nblks = mesh.cbsr_nblks;
rowpts    = mesh.cbsr_rowpts;
colind    = mesh.cbsr_colind;
ordered2unordered = mesh.BILU0ordered2unordered;
unordered2ordered = mesh.BILU0unordered2ordered;

Mg = Hg;
for j = 1:nrows
    rj = ordered2unordered(j);  % row rj    
    nj = (rowpts(rj)+1):rowpts(rj+1); % pointers to row rj 
    cj = colind(nj);               % connectivity indices for row rj        
    jj = rowpts(rj)+1;
    % inv(Ujj)
    Mg(:,:,jj) = inv(Mg(:,:,jj));
    for i = 2:length(nj)     % loop over each neighbor of row rj
        ri = cj(i);           % row ri is a neighbor of row rj 
        if unordered2ordered(ri) > j
           ni = (rowpts(ri)+1):rowpts(ri+1);
           ci = colind(ni);     % connectivity indices for row ri              
           aj = find(ci==rj);      % find the position of row rj in the neighboring list of row ri
           ij = rowpts(ri)+aj;
           % Lij = Uij*inv(Ujj)
           Mg(:,:,ij) = Mg(:,:,ij)*Mg(:,:,jj);
           ii = rowpts(ri)+1;
           ji = rowpts(rj)+i;
           % Uii = Uii - Lij*Uji
           Mg(:,:,ii) = Mg(:,:,ii) - Mg(:,:,ij)*Mg(:,:,ji);           
           for l = 2:length(ni)
               rij = ci(l); % row rij is a neighbor of both ri and rj
               if unordered2ordered(rij) > j
                   for m = 2:length(nj)
                       if rij == cj(m)       
                           il = rowpts(ri)+l;
                           jm = rowpts(rj)+m;
                           %  U_il <- U_il - L_ij*U_jm 
                           Mg(:,:,il) = Mg(:,:,il) - Mg(:,:,ij)*Mg(:,:,jm);
                       end
                   end
               end
           end
        end                    
    end
end

