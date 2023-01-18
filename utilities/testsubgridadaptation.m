function [udg,uhat,newsg] = testsubgridadaptation(mesh,mastersubgrid,udg,uh,c1,nlayer)

newsg = newssubgrid(mesh,mastersubgrid,udg,uh,c1,nlayer);

ns = length(mastersubgrid);
nc = size(udg{1},2);
ncu = size(uh,1);

oldsg = mesh.subgrids;
ds = oldsg - newsg;
if max(abs(ds))==0
    return;
else
    %es = find(ds~=0); % list of modified elements    
    for i=1:ns
        indi{i} = find(newsg==i);
        indj{i} = find(oldsg==i);
    end    
    for i=1:ns
        if isempty(indi{i})==0
            ni  = length(indi{i});
            ngi = size(mastersubgrid{i}.elemnodes,3); 
            npi = size(mastersubgrid{i}.shapvl,1); 
            udgt{i} = zeros(npi,nc,ngi*ni);
            jj = setdiff((1:ns),i);
            j1 = [];
            j2 = [];
            for k=1:ni
                e = indi{i}(k); % element e                           
                if ds(e)==0
                    m = find(indj{i}==e);
                    udgt{i}(:,:,(k-1)*ngi+1:k*ngi) = udg{i}(:,:,(m-1)*ngi+1:m*ngi);                                        
                else
                    j   = oldsg(e);                    
                    m   = find(indj{j}==e);  % position of element e in the list sj                                        
                    if j==jj(1)
                        j1 = [j1; e k m];                        
                    else
                        j2 = [j2; e k m];
                    end
                end
            end
            if isempty(j1)==0                    
                ngj = size(mastersubgrid{jj(1)}.elemnodes,3); % number of sub-elements for subgrid type j                                
                npj = size(mastersubgrid{jj(1)}.shapvl,1); 
                nj1 = size(j1,1);
                k1  = repmat((j1(:,2)'-1)*ngi,[ngi 1]) + repmat((1:ngi)',[1 nj1]);                
                k2  = repmat((j1(:,3)'-1)*ngj,[ngj 1]) + repmat((1:ngj)',[1 nj1]);                
                tm  = subgridprojection(mastersubgrid{jj(1)},mastersubgrid{i},mesh.dgnodes(:,:,j1(:,1)),reshape(udg{jj(1)}(:,:,k2),[npj nc ngj nj1]));                                                                                                                                
                udgt{i}(:,:,k1) = reshape(tm,[npi nc ngi*nj1]);
            end            
            if isempty(j2)==0                    
                ngj = size(mastersubgrid{jj(2)}.elemnodes,3); % number of sub-elements for subgrid type j                                
                npj = size(mastersubgrid{jj(2)}.shapvl,1); 
                nj2 = size(j2,1);
                k1  = repmat((j2(:,2)'-1)*ngi,[ngi 1]) + repmat((1:ngi)',[1 nj2]);                
                k2  = repmat((j2(:,3)'-1)*ngj,[ngj 1]) + repmat((1:ngj)',[1 nj2]);                
                tm  = subgridprojection(mastersubgrid{jj(2)},mastersubgrid{i},mesh.dgnodes(:,:,j2(:,1)),reshape(udg{jj(2)}(:,:,k2),[npj nc ngj nj2]));                                                                                                                                                
                udgt{i}(:,:,k1) = reshape(tm,[npi nc ngi*nj2]); 
            end            
        end
    end    
end
udg = udgt;
mesh.subgrids = newsg;

for i=1:ns
    npf  = size(mastersubgrid{i}.perm,1);        
    nf   = max(mastersubgrid{i}.elcon(:))/(npf);
    nes  = size(mastersubgrid{i}.t,1);    
    ind  = find(mesh.subgrids==i);
    ne   = length(ind);    
    uhat{i} = zeros(ncu,npf*nf,ne);        
    for k = 1:ne
        uhat{i}(:,:,k) = inituhat(mastersubgrid{i},mastersubgrid{i}.elcon,udg{i}(:,:,(k-1)*nes+1:k*nes),ncu);  
    end    
end




