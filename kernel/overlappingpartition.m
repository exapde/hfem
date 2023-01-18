function [elem2cpu,extintelem,extintelempts,face2cpu,extintface,extintfacepts,ent2cpu,extintent,extintentpts,outelem] =...
         overlappingpartition(t2f,t,elcon,facecon,hybrid,overlappinglevel,nproc,check) 
     
if nargin < 8; check = 0; end

ne = size(t,2);
nf = size(facecon,2);
npf = size(facecon,1);
nfe = size(t2f,1);

% nonoverlapping element and face partition
[elem2cpu, face2cpu] = t2fpartition(t2f+1,ne,nproc);    

% overlapping HDG patition
[extintelem,extintelempts,extintface,extintfacepts,outelem] = hdgpartitionc(t,t2f,elem2cpu,face2cpu,[size(t) nf nfe overlappinglevel nproc]);                         

% overlapping EDG patition
if strcmp(hybrid,'hdg')==0                                 
    ndh = max(facecon(:))+1;
    [ent2cpu,extintent,extintentpts] = edgpartitionc(extintelem, elcon, facecon, t2f, face2cpu, [ne nf ndh nfe npf nproc]);    
else
    ent2cpu = []; extintent = []; extintentpts = [];
end

% TODO: Reorder elements and entities

for i = 1:nproc    
    if isempty(outelem{i})
        error('HDG element partition is incorrect.');
    end
    if ~isempty(intersect(extintelem{i},outelem{i}))
        error('HDG element partition is incorrect.');
    end
    extintelem{i} = [extintelem{i}; outelem{i}];
    extintelempts(i,3) = length(outelem{i});
end

if check==1
    % check METIS decomposition
    for i = 1:nproc                           
        intelem = find(elem2cpu==(i-1));     
        
        intface = find(face2cpu==(i-1));     
        
        face = t2f(:,intelem)+1;  
        face = unique(face(:));
        if ~isempty(setdiff(intface,face)) 
            error('Some faces in the processor are not connected to any element in the processor.');
        end            
        
%         n = size(f,1);
%         extintelem = f((n-1):n,intface)+1;  
%         extintelem = unique(extintelem(:));
%         extintelem = extintelem(extintelem>0);
%         if ~isempty(setdiff(intelem,extintelem)) 
%             error('Some elements in the processor are not connected to any face in the processor.');
%         end                            
    end
    
    % check HDG patition
    [re,ce] = mkv2t(t+1,ne);                         
    for i = 1:nproc
        intelem = find(elem2cpu==(i-1));     

        intface = find(face2cpu==(i-1));     
        
        elem = intelem;
        for j=1:overlappinglevel            
            elem = node2elem(t(:,elem)+1,re,ce);
        end
        extelem = setdiff(elem,intelem);                        
        if max(abs(extintelem{i}(1:sum(extintelempts(i,1:2)))+1-[intelem; extelem]))~=0
            error('HDG element partition is incorrect.');
        end
        if max(abs(extintelempts(i,1:2)-[length(intelem) length(extelem)]))~=0
            error('HDG element partition is incorrect.');
        end
        elem = intelem;
        for j=1:overlappinglevel+1            
            elem = node2elem(t(:,elem)+1,re,ce);
        end
        extelem = setdiff(elem,intelem);     
        if max(abs(outelem{i}+1-setdiff([intelem; extelem],extintelem{i}(1:sum(extintelempts(i,1:2)))+1)))~=0            
            error('HDG element partition is incorrect.');
        end
        
        face = t2f(:,extintelem{i}(1:sum(extintelempts(i,1:2)))+1)+1;  
        face = unique(face(:));
        extface = setdiff(face,intface);                        
        if max(abs(extintface{i}+1-[intface; extface]))~=0
            error('HDG face partition is incorrect.');
        end
        if max(abs(extintfacepts(i,:)-[length(intface) length(extface)]))~=0
            error('HDG face partition is incorrect.');
        end
    end

    if strcmp(hybrid,'hdg')==0                                 
        % check EDG patition
        ent2cpum = edgpartition(facecon, t, t2f, re, ce, elem2cpu, face2cpu, nproc, overlappinglevel);        
        if max(abs(ent2cpu-ent2cpum))~=0
            error('EDG entity partition is incorrect.');
        end
        for i = 1:nproc
            intent = find(ent2cpum==(i-1));            
            ent = elcon(:,extintelem{i}(1:sum(extintelempts(i,1:2)))+1)+1;
            ent = unique(ent(:));                                
            ent = ent(ent>0);                      
            extent = setdiff(ent,intent);      
            if max(abs(extintent{i}(1:extintentpts(i,1))+1-intent))~=0
                error('EDG entity partition is incorrect.');
            end
            if max(abs(extintent{i}+1-[intent; extent]))~=0
                error('EDG entity partition is incorrect.');
            end
            if max(abs(extintentpts(i,:)-[length(intent) length(extent)]))~=0
                error('EDG entity partition is incorrect.');
            end            
        end
    end
end

