function [udg,uhat,mesh] = subgridadaptation(mesh,mastersubgrid,udg,uhat,uh,c1)

ns = length(mastersubgrid);
if ns==1
    return;
end
if ns>3
    error('Subgrid types can not exceed 3');
end

u = udgproj(mesh,mastersubgrid,udg);
nc = size(u,2);
if exist('uh','var')
    ncu = size(uh,1);    
    q = getq(mesh,mastersubgrid{1},u,reshape(uh(:,mesh.elcon),ncu,[],size(u,3)));
    s = sensor(cat(2,u,q));    
else
    s = sensor(u);    
end
a = mean(squeeze(s));

ne    = length(mesh.subgrids);
newsg = ones(ne,1);

% update new subgrids
%c1 = -0.075; 
c2 = c1; 
c3 = 10;
ind0 = find(a <= c1); % p0 elements
if ns==2
    newsg(ind0) = 2; 
else
    ind2 = find(c2<a & a<=c3);               % high-order elements
    ind1 = setdiff((1:ne)',union(ind0,ind2)); % p1 elements    
    newsg(ind0) = 3;
    newsg(ind1) = 2;
end

% x = mesh.dgnodes(:,1,:);
% y = mesh.dgnodes(:,2,:);
% % dist = min(x.^2+(y).^2);
% % in = find(abs(dist-1)<1e-6);
% % newsg(in)=2;
% dist = min(x.^2+(y-1).^2);
% in = find(abs(dist)<1e-6);
% newsg(in)=3;
% dist = min(x.^2+(y+1).^2);
% in = find(abs(dist)<1e-6);
% newsg(in)=3;

for jj=1:0
    t2t = mkt2t(mesh.t,mesh.elemtype);
    ind = find(newsg==ns);
    for i=1:length(ind)    
        t = t2t(ind(i),:);
        k = t > 0;    
        q = t(k);
        j = ismember(q,ind);
        newsg(q(j==0))=ns;    
    end
end

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
%                     ngj = size(mastersubgrid{j}.elemnodes,3); % number of sub-elements for subgrid type j                                
%                     npj = size(mastersubgrid{j}.shapvl,1); 
%                     tm  = subgridprojection(mastersubgrid{j},mastersubgrid{i},mesh.dgnodes(:,:,e),reshape(udg{j}(:,:,(m-1)*ngj+1:m*ngj),[npj nc ngj 1]));                                                                                                
%                     udgt{i}(:,:,(k-1)*ngi+1:k*ngi) = tm;
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

% for high-order elements
% in = find(mesh.subgrids==1);
% if isempty(in)==0    
%     npf  = size(mastersubgrid{1}.perm,1);
%     nfe  = size(mastersubgrid{1}.perm,2);
%     UH = reshape(uh(:,mesh.elcon),ncu,npf,nfe,ne);
%     app.tdep = 0;
%     app.arg{end}=1;
%     app.itmax=10;
%     %[dudg,udg,dudg_duh,Run,BDn,Rln] = elementsolve(master{s0},app,mesh.dgnodes(:,:,ie),SH{s0}(:,:,id),UH(:,:,:,ie),UDG{s0}(:,:,id));
%     [~,tm] = elementsolve(mastersubgrid{1},app,mesh.dgnodes(:,:,in),0*udg{1},UH(:,:,:,in),udg{1});
%     udg{1} = tm;
% end

for i=1:ns
    npf  = size(mastersubgrid{i}.perm,1);        
    nf   = max(mastersubgrid{i}.elcon(:))/(npf);
    nes  = size(mastersubgrid{i}.t,1);    
    ind  = find(mesh.subgrids==i);
    ne   = length(ind);    
    uhat{i} = zeros(ncu,npf*nf,ne);        
    for k = 1:ne
        uhat{i}(:,:,k) = inituhat2(mastersubgrid{i},udg{i}(:,:,(k-1)*nes+1:k*nes),ncu);  
    end    
end

% for k = 1:length(newsg)
%     if a(k) <= c1
%         newsg(k) = porders;
%     elseif b(k) >= c2
%         newsg(k) = porders(end-1);
%     end
% end

% ind = find(a<=-2);      
% for i=1:length(mastersubgrid)        
%     if mastersubgrid{i}.porder==0                        
%         newsg(ind) = i;
%         break;
%     end
% end

% ind0 = find(a<=-2);       % elements
% st   = mesh.subgrids(ind0);  % subgrid types
% if isempty(ind0)==0          
%     for i=1:length(mastersubgrid)        
%         if mastersubgrid{i}.porder==0                        
%             for k=1:length(st)                                
%                 j = st(k);                
%                 if j~=i   
%                     e   = ind0(k); % element e
%                     sj  = find(mesh.subgrids == j); % list of elements with subgrid type j 
%                     m   = find(sj == e);  % position of element e in the list sj
%                     ngi = size(mastersubgrid{i}.elemnodes,3); % number of sub-elements for subgrid type i                          
%                     ngj = size(mastersubgrid{j}.elemnodes,3); % number of sub-elements for subgrid type j                                
%                     npj = size(udg{j},1);
%                     nc  = size(udg{j},2);                              
%                     % projection from subgrid type j to subgrid type i
%                     tm = subgridprojection(mastersubgrid{j},mastersubgrid{i},mesh.dgnodes(:,:,e),reshape(udg{j}(:,:,(m-1)*ngj+1:m*ngj),[npj nc ngj 1]));                                                                    
%                     % remove element e from subgrid type j 
%                     udg{j}(:,:,(m-1)*ngj+1:m*ngj) = [];
%                     % add element e to subgrid type i 
%                     mesh.subgrids(e) = i;
%                     si  = find(mesh.subgrids == i);
%                     n   = find(si == e);                    
%                     if n==length(si)
%                         udg{i}(:,:,(n-1)*ngi+1:n*ngi) = tm;                                        
%                         uhat{i}(:,:,n) = inituhat2(mastersubgrid{i},tm,ncu); 
%                     else
%                         udg{i}(:,:,n*ngi+1:end+ngi) = udg{i}(:,:,(n-1)*ngi+1:end);
%                         udg{i}(:,:,(n-1)*ngi+1:n*ngi) = tm;       
%                         uhat{i}(:,:,n+1:end+1) = uhat{i}(:,:,n:end); 
%                         uhat{i}(:,:,n) = inituhat2(mastersubgrid{i},tm,ncu); 
%                     end                     
%                 end                                
%             end            
%             break;
%         end
%     end
% end
