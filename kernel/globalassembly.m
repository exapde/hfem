function [Hg, Rg] = globalassembly(mesh, He, Re)

ne    = mesh.ne;
nfe   = size(mesh.perm,2); % number of faces per element    
%nf    = mesh.nf;
nblks = mesh.cbsr_nblks;
nrows = mesh.cbsr_nrows;
rowpts    = mesh.cbsr_rowpts;
colind    = mesh.cbsr_colind;

if strcmp(mesh.hybrid,'hdg') 
    ncu = size(Re,1);
    npf = size(Re,2)/nfe;
    bsz = ncu*npf;        
    He = reshape(He,[ncu npf nfe ncu npf nfe ne]);
    Re = reshape(Re,[ncu npf nfe ne]);    
    elcon = reshape(mesh.elcon,[npf nfe ne]);
    
    Hg = zeros(bsz,bsz,nblks);
    Rg = zeros(bsz,nrows);
    for e = 1:ne
        t2f = mesh.t2f(e,:);                
        for a = 1:nfe
            i = t2f(a);  % block row i 
            ind = elcon(:,a,e) - (i-1)*npf;             
            tm = Re(:,ind,a,e);
            Rg(:,i) = Rg(:,i) + tm(:);            
            for b = 1:nfe            
                j = t2f(b); % block column j                           
                for k = (rowpts(i)+1):rowpts(i+1)
                    if colind(k) == j, break; end
                end
                jnd = elcon(:,b,e) - (j-1)*npf;  
                Hg(:,:,k) = Hg(:,:,k) + reshape(He(:,ind,a,:,jnd,b,e),[bsz bsz]);                
            end
        end
    end
    
%     ind = zeros(npf,1);
%     jnd = ind;    
%     Hg = zeros(bsz,bsz,nblks);
%     Rg = zeros(bsz,nrows);    
%     bsz = ncu*npf;
%     sz0 = npf*nfe;
%     sz1 = ncu*npf*nfe;
%     sz2 = ncu*bsz;
%     sz3 = bsz*bsz;
%     sz4 = ncu*sz1;
%     sz5 = npf*ncu*sz1;
%     sz6 = sz1*sz1;
%     for e = 1:ne                    
%         for a = 1:nfe
%             i = mesh.t2f((a-1)*ne+e);                
%             for m1 = 1:npf
%                 ind(m1) = elcon((e-1)*sz0+(a-1)*npf+m1) - (i-1)*npf;  
%             end
%             for m1 = 1:npf
%                 for n1 = 1:ncu
%                     Rg((i-1)*bsz+(m1-1)*ncu+n1) = Rg((i-1)*bsz+(m1-1)*ncu+n1) + Re((e-1)*sz1+(a-1)*bsz+(ind(m1)-1)*ncu+n1);                                                               
%                 end
%             end                        
%             for b = 1:nfe             
%                 j = mesh.t2f((b-1)*ne+e);      
%                 for k = (rowpts(i)+1):rowpts(i+1)
%                     if (colind(k) == j), break; end    
%                 end                                                                          
%                 for m2 = 1:npf 
%                     jnd(m2) = elcon((e-1)*sz0+(b-1)*npf+m2) - (j-1)*npf;
%                 end
%                 for m1 = 1:npf
%                     for n1 = 1:ncu
%                         for m2 = 1:npf
%                             for n2 = 1:ncu
%                                 t1 = (k-1)*sz3+(m1-1)*sz2+(n1-1)*bsz+(m2-1)*ncu+n2;
%                                 t2 = (e-1)*sz6+(b-1)*sz5+(jnd(m1)-1)*sz4+(n1-1)*sz1+(a-1)*bsz+(ind(m2)-1)*ncu+n2;
%                                 Hg(t1) = Hg(t1) + He(t2);                                 
%                             end
%                         end
%                     end
%                 end                
%             end                    
%         end               
%     end
elseif strcmp(mesh.hybrid,'edg') 
    ncu = size(Re,1);
    bsz = ncu;    
    nn = size(mesh.elcon,1);
    He = reshape(He,[bsz nn bsz nn ne]);
    Re = reshape(Re,[bsz nn ne]);
    Hg = zeros(bsz,bsz,nblks);
    Rg = zeros(bsz,nrows);    
    
    for e = 1:ne
        elc = mesh.elcon(:,e);
        for a = 1:nn
            i = elc(a);  % block row i
            Rg(:,i) = Rg(:,i) + Re(:,a,e);
            for b = 1:nn            
                j = elc(b); % block column j
                for k = (rowpts(i)+1):rowpts(i+1)
                    if colind(k) == j, break; end
                end                
                Hg(:,:,k) = Hg(:,:,k) + reshape(He(:,a,:,b,e),[bsz bsz]);                
            end
        end
    end    
    
%     sz0 = bsz*bsz;
%     sz1 = bsz*nn;
%     sz2 = bsz*nn*bsz;
%     sz3 = bsz*bsz*nn*nn;
%     for e = 1:ne
%         elc = mesh.elcon(:,e);
%         for a = 1:nn
%             i = elc(a);  % block row i
%             for m = 1:bsz                
%                 Rg((i-1)*bsz+m) = Rg((i-1)*bsz+m) + Re((e-1)*bsz*nn+(a-1)*bsz+m);
%             end
%             for b = 1:nn            
%                 j = elc(b); % block column j
%                 for k = (rowpts(i)+1):rowpts(i+1)
%                     if colind(k) == j, break; end
%                 end               
%                 for m = 1:bsz
%                     for n = 1:bsz                        
%                         t1 = (k-1)*sz0+(m-1)*bsz+n;
%                         t2 = (e-1)*sz3+(b-1)*sz2+(m-1)*sz1+(a-1)*bsz+n;
%                         Hg(t1) = Hg(t1) + He(t2);
%                     end
%                 end
%             end
%         end
%     end    
    
end




