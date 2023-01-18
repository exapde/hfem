function [extintelem,extintelempts,extintface,extintfacepts,outelem] = hdgpartition(t,t2f,elem2cpu,face2cpu,overlappinglevel,nproc)

extintelem = cell(nproc,1);
outelem = cell(nproc,1);
extintface = cell(nproc,1);
extintelempts = zeros(nproc,2);
extintfacepts = zeros(nproc,2);

ne = size(t,2);
[re,ce] = mkv2t(t+1,ne);                         
for i = 1:nproc
    intelem = find(elem2cpu==(i-1));     

    intface = find(face2cpu==(i-1));     

    elem = intelem;
    for j=1:overlappinglevel            
        elem = node2elem(t(:,elem)+1,re,ce);
    end
    extelem = setdiff(elem,intelem);                        
    
    extintelem{i} = [intelem; extelem];
    extintelempts(i,1:2) = [length(intelem) length(extelem)];        
    
    elem = intelem;
    for j=1:overlappinglevel+1            
        elem = node2elem(t(:,elem)+1,re,ce);
    end
    extelem = setdiff(elem,intelem);     
    
    outelem{i} = setdiff([intelem; extelem],extintelem{i}(1:sum(extintelempts(i,1:2))));
            
    face = t2f(:,extintelem{i}(1:sum(extintelempts(i,1:2))))+1;  
    face = unique(face(:));
    extface = setdiff(face,intface);                        
    extintface{i} = [intface; extface];
    extintfacepts(i,:)=[length(intface) length(extface)];    
    
    extintelem{i} = extintelem{i}-1;
    outelem{i} = outelem{i}-1;
    extintface{i}  = extintface{i}-1;  
end

