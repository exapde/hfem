function [dmd,meshp] = domaindecomposition(mesh,nproc)

if strcmp(mesh.hybrid,'hdg')
    dmd = domaindecomposition_hdg(mesh,nproc);
elseif strcmp(mesh.hybrid,'edg')    
    dmd = domaindecomposition_edg(mesh,nproc);    
end

[npf,nfe] = size(mesh.perm);
for i=1:nproc % loop over each subdomain
    meshp{i}.dgnodes = mesh.dgnodes(:,:,dmd{i}.elempart);
    meshp{i}.bf = mesh.bf(:,dmd{i}.elempart);
    % compute local element-to-entity connectivities
    elcon = reshape(mesh.elcon,[npf nfe mesh.ne]);
    meshp{i}.elcon = elcon(:,:,dmd{i}.elempart);
    ne = length(dmd{i}.elempart);        
    %reshape(meshp{i}.elcon,[npf*nfe ne])
    if strcmp(mesh.hybrid,'hdg')
        t2f = mesh.t2f(dmd{i}.elempart,:);    
        for j = 1:ne
            for k = 1:nfe
                t2f(j,k) = find(t2f(j,k) == dmd{i}.entpart);
            end
        end
        for j = 1:ne
            for k = 1:nfe
                ii = elcon(:,k,dmd{i}.elempart(j)) - (mesh.t2f(dmd{i}.elempart(j),k)-1)*npf;                
                ind = ((t2f(j,k)-1)*npf+1):(t2f(j,k)*npf);
                meshp{i}.elcon(:,k,j) = ind(ii);
            end
        end
        meshp{i}.t2f = t2f;
    elseif strcmp(mesh.hybrid,'edg')        
        for j = 1:ne
            for k = 1:nfe        
                for m = 1:npf                    
                    meshp{i}.elcon(m,k,j) = find(meshp{i}.elcon(m,k,j) == dmd{i}.entpart);
                end
            end
        end
        [~,meshp{i}.t2f] = mkt2f(mesh.t(dmd{i}.elempart,:),mesh.elemtype);
    end        
    meshp{i}.elcon = reshape(meshp{i}.elcon,[npf*nfe ne]);   
    meshp{i} = block_crs(meshp{i},mesh.hybrid);
    
    nownent = sum(dmd{i}.entnum(1:2));
    meshp{i}.BJ_rowpts = zeros(nownent+1,1);
    meshp{i}.BJ_colind = [];
    meshp{i}.BK_rowpts = 0;
    meshp{i}.BK_colind = [];
    meshp{i}.BK_rowind = 0;
    m = 1;
    for j = 1:nownent
        nj = (meshp{i}.cbsr_rowpts(j)+1):meshp{i}.cbsr_rowpts(j+1);
        cj = meshp{i}.cbsr_colind(nj);
        i1 = find(cj<=nownent);
        c1 = cj(i1);
        meshp{i}.BJ_rowpts(j+1) = meshp{i}.BJ_rowpts(j)+length(c1);
        meshp{i}.BJ_colind = [meshp{i}.BJ_colind; c1];
        i2 = find(cj>nownent);
        c2 = cj(i2);
        if ~isempty(c2)
            meshp{i}.BK_rowpts(m+1) = meshp{i}.BK_rowpts(m)+length(c2);
            meshp{i}.BK_colind = [meshp{i}.BK_colind; c2];
            meshp{i}.BK_rowind(m)  = j;
            m = m + 1;
        end
    end
    
%     dmd{i}.entpart'
%     dmd{i}.elempart'
%     meshp{i}.elcon 
%     pause
end

%[ae{i},fe{i},dudg{i},dudg_duh{i}] = hdg_elemental(masterp{i},appp{i},dgnodes{i},bf{i},udg{i},uh{i},sh{i});
    
function dmd = domaindecomposition_hdg(mesh,nproc)

[elempart,entpart] = meshpart(mesh.t2f,nproc);

% plot domain partition
plotpart(mesh, nproc, elempart, entpart);

elempart = elempart+1;
entpart = entpart+1;

t2t = mkt2t(mesh.t,mesh.elemtype);

% loop over each subdomain
for i=1:nproc
    
    % list of elements in subdomain i                       
    indelem{i} = find(elempart==i); 
    
    % list of faces in subdomain i                       
    indent{i} = find(entpart==i); 
end

% loop over each subdomain
for i=1:nproc    
        
    % numbering elements
    elemnumber = zeros(length(indelem{i}),1); 
    for j = 1:length(indelem{i})
        tj = t2t(indelem{i}(j),:);
        tj = tj(tj>0);
        if all(ismember(tj,indelem{i}))
            % elements inside the subdomain i
            elemnumber(j) = 2;
        else
            % elements on the interface between two subdomains
            elemnumber(j) = 1;
        end
    end
    tib = find(elemnumber==1);
    tin = find(elemnumber==2);
    
    % [interface elements, interior elements]   
    dmd{i}.elempart = [indelem{i}(tib); indelem{i}(tin)];
        
    % numbering faces 
    facenumber = zeros(length(indent{i}),1); 
    for j = 1:length(indent{i})
        fj = mesh.f(indent{i}(j),end-1:end);
        if fj(2) <= 0
            % faces on the domain boundary
            facenumber(j) = 2; 
        elseif all(ismember(fj,indelem{i})) % (ismember(fj(1),indelem{i}) && ismember(fj(2),indelem{i}))
            % faces inside the subdomain i
            facenumber(j) = 2;
        else
            % faces on the interface between two subdomains
            facenumber(j) = 1;            
        end
    end
    fib = find(facenumber==1);
    %fdb = find(facenumber==2);
    fin = find(facenumber==2);
        
    % [interior faces, own interface faces]   
    %dmd{i}.entpart = [indent{i}(fin); indent{i}(fib)];               
    
    elemrecv = zeros(length(fib),5);
    dmd{i}.entsend = zeros(length(fib),3);
    dmd{i}.entrecv = zeros(length(fib),3*(size(mesh.t2f,2)-1));
    for j = 1:length(fib)    
        %fj = dmd{i}.entpart(j);
        fj = indent{i}(fib(j));
        f2e = mesh.f(fj,end-1:end); % face-to-element 
        if ismember(f2e(1),indelem{i}) 
            % element f2e(1) belongs to subdomain i                                            
            % element f2e(2) belongs to subdomain k                     
            k = elempart(f2e(2));
            elemrecv(j,1) = k;                    % subdomain k is a neighbor of subdomain i
            elemrecv(j,2) = find(indelem{i}==f2e(1)); % local location of element f2e(1)   
            elemrecv(j,3) = f2e(1);                % global location of element f2e(1)   
            elemrecv(j,4) = find(indelem{k}==f2e(2)); % local location of element f2e(2)   
            elemrecv(j,5) = f2e(2);                % global location of element f2e(2)               
        else
            % element f2e(1) belongs to subdomain k         
            % element f2e(2) belongs to subdomain i                     
            k = elempart(f2e(1));
            elemrecv(j,1) = k;
            elemrecv(j,2) = find(indelem{i}==f2e(2)); % local location of element f2e(2)               
            elemrecv(j,3) = f2e(2);                % global location of element f2e(2)    
            elemrecv(j,4) = find(indelem{k}==f2e(1)); % local location of element f2e(1)   
            elemrecv(j,5) = f2e(1);                % global location of element f2e(1)         
        end
                
        dmd{i}.entsend(j,1) = k;
        dmd{i}.entsend(j,2) = find(indent{i}==fj); % local location of face fj
        dmd{i}.entsend(j,3) = fj;                % global location of face fj
                                
        f2f = mesh.f2f(:,fj); % face-to-face                      
        m = 1;
        for l=2:length(f2f) % for each neighboring face
            % determine if subdomain k contains the neighboring face
            for k = 1:nproc
                if ismember(f2f(l),indent{k}) && (i ~= k)                     
                    % face f2f(l) belongs to subdomain k                     
                    dmd{i}.entrecv(j,m) = k;
                    dmd{i}.entrecv(j,m+1) = find(indent{k}==f2f(l));
                    dmd{i}.entrecv(j,m+2) = f2f(l);
                    m = m + 3;
                end
            end
        end
    end
    
    % a list of neighboring subdomains
%     dmd{i}.nbsd = unique(elemrecv(:,1))';
% %     i1 = nbsd>i;    
% %     i2 = nbsd<i;
% %     dmd{i}.nbsd = [sort(nbsd(i1)); sort(nbsd(i2))];
%     dmd{i}.nbsd
%     elemrecv
%     pause
    
    % reordering 
    dmd{i}.nbsd = unique(elemrecv(:,1))';
    tm = [];
    for j = 1:length(dmd{i}.nbsd)
        ij = find(elemrecv(:,1)==dmd{i}.nbsd(j));
        tm = [tm; ij];
    end    
    
    % store elements received from neighboring subdomains to assemble the linear system
    dmd{i}.elemrecv = elemrecv(tm,:);    
            
    % store faces sent to neighboring subdomains to perform the matrix-vector product
    dmd{i}.entsend = dmd{i}.entsend(tm,:);
    
    % store faces received from neighboring subdomains to perform the matrix-vector product
    dmd{i}.entrecv = dmd{i}.entrecv(tm,:);
    
    % find interface faces that do not belong to subdomain i
    otherfaces = [];
    for j=1:length(tib)
        ej = indelem{i}(tib(j));
        t2fj = mesh.t2f(ej,:);
        fj = ismember(t2fj,indent{i});
        otherfaces = [otherfaces t2fj(fj==0)];
    end            
        
    % find faces that are connected to own interface faces and do not belong to subdomain i
    otherelements = [];
    for j = 1:length(fib)    
        fj = indent{i}(fib(j));
        f2f = mesh.f2f(:,fj); % face-to-face                              
        for l=2:length(f2f) % for each neighboring face
            % determine if subdomain i contains the neighboring face
            if ismember(f2f(l),[indent{i};otherfaces'])==0
                otherfaces = [otherfaces f2f(l)];
            end
        end
        f2e = mesh.f(fj,end-1:end); % face-to-element 
        if ismember(f2e(1),indelem{i}) 
            otherelements = [otherelements f2e(2)];
        else
            otherelements = [otherelements f2e(1)];
        end
    end    
    
    % [interior faces, own interface faces, other faces]   
    dmd{i}.entpart = [indent{i}(fin); indent{i}(fib(tm)); otherfaces(:)];  
    dmd{i}.entnum = [length(fin) length(fib) length(otherfaces)];
    
    % [interface elements, interior elements, other elements]  
    dmd{i}.elempart = [dmd{i}.elempart; otherelements(:)];
    dmd{i}.elemnum = [length(tib) length(tin) length(otherelements)];    
end

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc
    dmd{k}.elemsend = [];
end
for i = 1:nproc    
    for j = 1:length(dmd{i}.nbsd)
        k = dmd{i}.nbsd(j);
        ii = find(dmd{i}.elemrecv(:,1)==k);
        tm = dmd{i}.elemrecv(ii,[1 4 5 2 3]);
        tm(:,1) = i;
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end
end

for i = 1:nproc 
    dmd{i}.elemsend = unique(dmd{i}.elemsend,'rows');
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');   
end

% update dmd{i}.entrecv
nfe = size(mesh.t2f,2);
for i = 1:nproc        
    tm = [];
    for j = 1:nfe-1
        tm = [tm; dmd{i}.entrecv(:,(3*(j-1)+1):3*j)];        
    end
    dmd{i}.entrecv = tm;        
    r{i} = dmd{i}.entrecv; 
end
for i = 1:nproc    
    for j = 1:length(dmd{i}.nbsd)
        k = dmd{i}.nbsd(j);
        ii = find(dmd{i}.entsend(:,1)==k); 
        tm = dmd{i}.entsend(ii,:);
        tm(:,1) = i;
        dmd{k}.entrecv = [dmd{k}.entrecv; tm];        
    end
end

% update dmd{i}.entsend
for i = 1:nproc    
    for j = 1:length(dmd{i}.nbsd)
        k = dmd{i}.nbsd(j);
        ii = find(r{i}(:,1)==k); 
        tm = r{i}(ii,:);
        tm(:,1) = i;        
        dmd{k}.entsend = [dmd{k}.entsend; tm];        
    end
end


for i = 1:nproc 
    dmd{i}.entsend = unique(dmd{i}.entsend,'rows');
    dmd{i}.entrecv = unique(dmd{i}.entrecv,'rows');       
    for j = 1:size(dmd{i}.entsend,1)
        fj = dmd{i}.entsend(j,3);
        k = find(dmd{i}.entpart==fj);
        dmd{i}.entsend(j,2) = k;
    end
    %dmd{i}.entrecv'
    for j = 1:size(dmd{i}.entrecv,1)
        fj = dmd{i}.entrecv(j,3);
        k = find(dmd{i}.entpart==fj);
        dmd{i}.entrecv(j,2) = k;
    end            
    dmd{i}.nbsd = unique(dmd{i}.entrecv(:,1))';
end

% for i = 1:nproc    
%     dmd{i}.elemsend
%     dmd{i}.elemrecv
% end
% 
% for i = 1:nproc    
%     dmd{i}.entsend
%     dmd{i}.entrecv
% end
    
function dmd = domaindecomposition_edg(mesh,nproc)

[elempart,entpart] = meshpart(mesh.elcon',nproc);

% plot domain partition
plotpart(mesh, nproc, elempart, entpart);

elempart = elempart+1;
entpart = entpart+1;

mesh.f(:,end+1) = 1;    
[elcon,~,edgnodes] = elconnectivities(mesh);    
mesh.f(:,end) = [];    

% loop over each subdomain
for i=1:nproc
    
    % list of elements in subdomain i                       
    indelem{i} = find(elempart==i); 
    
    % list of edg nodes in subdomain i                       
    indent{i} = find(entpart==i); 
end

tm = [];
% loop over each subdomain
for i=1:nproc
    
    m = 1;
    % numbering edgnodes
    edgnumber = zeros(length(indent{i}),1); 
    for j = 1:length(indent{i})
        [~,edgj] = find(mesh.elcon == indent{i}(j));                
        edgj = unique(edgj);
        if all(ismember(edgj,indelem{i}))
            % edgnode indent{i}(j) is inside the subdomain i
            edgnumber(j) = 2;
        else
            % edgnode indent{i}(j) is on the interface between two subdomains
            edgnumber(j) = 1;            
            
            % find subdomains connected to edgnode indent{i}(j)
            dmd{i}.intfent(m,1) = indent{i}(j);
            for k = 1:nproc
                if any(ismember(edgj,indelem{k}))
                    dmd{i}.intfent(m,k+1) = 1;
                else
                    dmd{i}.intfent(m,k+1) = 0;
                end
            end
            m = m + 1;
        end            
    end
    e1{i} = find(edgnumber==1);
    e2{i} = find(edgnumber==2);
        
    % [interior edgnodes; own interface edgnodes]   
    dmd{i}.entpart = [indent{i}(e2{i}); indent{i}(e1{i})];          
    
    tm = [tm; dmd{i}.intfent];
end

for i=1:nproc
    
    % find all edgnodes on the interface including edgnodes from the neighbors
    ii = tm(:,i+1)==1;
    dmd{i}.allintfent = tm(ii,:);         
    
    % find all elements connected to the interface edgnodes
    intfelem = [];
    for j = 1:size(dmd{i}.allintfent,1)
        indnj = dmd{i}.allintfent(j,1);
        [~,edgj] = find(mesh.elcon == indnj);    
        edgj = unique(edgj); % all elements connected to edgnode indnj                
        intfelem = [intfelem; edgj(ismember(edgj,indelem{i}))];                
    end
    dmd{i}.intfelem = unique(intfelem);
    
    % [interface elements, interior elements]   
    dmd{i}.elempart = [dmd{i}.intfelem; setdiff(indelem{i},dmd{i}.intfelem)];
    
    % find subdomains which neighbor subdomain i
    m = 1;
    for j=1:nproc        
        if (i~=j) && any(dmd{i}.allintfent(:,j+1))
            dmd{i}.nbsd(m) = j;
            m = m + 1;
        end
    end    
    
    % find interface edgnodes that do not belong to subdomain i
    othernodes = [];
    for j=1:length(dmd{i}.intfelem)
        ej = dmd{i}.intfelem(j);
        edg = mesh.elcon(:,ej);
        fj = ismember(edg,indent{i});
        othernodes = [othernodes; edg(fj==0)];
    end        
    othernodes = unique(othernodes);    
    
    % [interior edgnodes, own interface edgnodes, other interface edgnodes]   
    %dmd{i}.entpart = [dmd{i}.entpart; othernodes]; 
    
    % find other edgnodes that are connected to own interface edgnodes and do not belong to subdomain i    
    otherelements = [];
    for j = 1:length(e1{i})    
        indnj = indent{i}(e1{i}(j));        
        [~,edgj] = find(mesh.elcon == indnj);    
        edgj = unique(edgj); % all elements connected to edgnode indnj                
        edg = mesh.elcon(:,edgj);
        edg = unique(edg(:));
        ii = ismember(edg,[dmd{i}.entpart; othernodes]);
        othernodes = [othernodes; edg(ii==0)];
        ii = ismember(edgj,indelem{i});
        otherelements = [otherelements; edgj(ii==0)];
    end    
    otherelements = unique(otherelements);

    % [interior edgnodes, own interface edgnodes, other interface edgnodes]   
    dmd{i}.entpart = [indent{i}(e2{i}); indent{i}(e1{i}); othernodes];    
    dmd{i}.entnum = [length(e2{i}) length(e1{i}) length(othernodes)];
    
    % [interface elements, interior elements, other elements]   
    dmd{i}.elempart = [dmd{i}.elempart; otherelements];    
    dmd{i}.elemnum = [length(dmd{i}.intfelem) length(indelem{i})-length(dmd{i}.intfelem) length(otherelements)];
    
%     dmd{i}.intfent
%     dmd{i}.intfelem
%     dmd{i}.allintfent
%     dmd{i}.nbsd
%     pause
end

% compute dmd{i}.elemrecv
for i = 1:nproc
    n = 1;
    for j = 1:size(dmd{i}.intfent,1) % for each interface edgnode
        indnj = dmd{i}.intfent(j,1);
        [~,edgj] = find(mesh.elcon == indnj);    
        edgj = unique(edgj); % all elements connected to edgnode indnj                        
        for m = 1:length(edgj) % for each element connected to edgnode indnj        
            for t = 1:length(dmd{i}.nbsd) % for each neighboring subdomain
                k = dmd{i}.nbsd(t);                        
                if ismember(edgj(m),indelem{k}) % if this element belongs to subdomain k
                    dmd{i}.elemrecv(n,1) = k; % subdomain k is a neighbor of subdomain i
                    dmd{i}.elemrecv(n,2) = find(edgj(m)==indelem{k}); % local location of element edgj(m)
                    dmd{i}.elemrecv(n,3) = edgj(m);  % global location of element edgj(m)
                    dmd{i}.elemrecv(n,4) = find(indnj==indent{i}); % local location of node indnj
                    dmd{i}.elemrecv(n,5) = indnj; % global location of node indnj                        
                    n = n + 1;
                end
            end                                   
        end        
    end
end

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc
    dmd{k}.elemsend = [];
end
for i = 1:nproc    
    for j = 1:length(dmd{i}.nbsd)
        k = dmd{i}.nbsd(j);
        ii = find(dmd{i}.elemrecv(:,1)==k);
        tm = dmd{i}.elemrecv(ii,:);
        tm(:,1) = i;
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end
end

for i = 1:nproc 
    dmd{i}.elemsend(:,end-1:end)=[];
    dmd{i}.elemrecv(:,end-1:end)=[];
    dmd{i}.elemsend = unique(dmd{i}.elemsend,'rows');
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');   
%     dmd{i}.elemsend
%     dmd{i}.elemrecv
end

% compute dmd{i}.entsend amd dmd{i}.entrecv
for i = 1:nproc
    q = 1;
    n = 1;
    for j = 1:size(dmd{i}.intfent,1) % for each interface edgnode
        indnj = dmd{i}.intfent(j,1);
        [~,edgj] = find(mesh.elcon == indnj);    
        edgj = unique(edgj); % all elements connected to edgnode indnj    
        edg = mesh.elcon(:,edgj);
        edg = unique(edg(:)); % all edgnodes connected to edgnode indnj    
        for t = 1:length(dmd{i}.nbsd) % for each neighboring subdomain
            k = dmd{i}.nbsd(t);     
            if ismember(indnj,dmd{k}.allintfent) % if node indnj is also on the subdomain k
                dmd{i}.entsend(q,1) = k;
                dmd{i}.entsend(q,2) = find(indnj==indent{i});
                dmd{i}.entsend(q,3) = indnj;
                q = q + 1;
            end
            for m = 1:length(edg)
                if ismember(edg(m),indent{k}) % if this edgnode belongs to subdomain k
                    dmd{i}.entrecv(n,1) = k; % subdomain k is a neighbor of subdomain i
                    dmd{i}.entrecv(n,2) = find(edg(m)==indent{k}); % local location of node edg(m)
                    dmd{i}.entrecv(n,3) = edg(m);  % global location of node edg(m)
                    %dmd{i}.entrecv(n,4) = find(indnj==indent{i}); % local location of node indnj
                    %dmd{i}.entrecv(n,5) = indnj; % global location of node indnj                        
                    n = n + 1;
                end
            end                
        end        
    end
end

% update dmd{i}.entsend
for i = 1:nproc    
    for j = 1:length(dmd{i}.nbsd)
        k = dmd{i}.nbsd(j);
        ii = find(dmd{i}.entrecv(:,1)==k); 
        tm = dmd{i}.entrecv(ii,:);
        tm(:,1) = i;        
        dmd{k}.entsend = [dmd{k}.entsend; tm(:,1:3)];        
    end
end

% update dmd{i}.entrecv
for i = 1:nproc    
    for j = 1:length(dmd{i}.nbsd)
        k = dmd{i}.nbsd(j);
        ii = find(dmd{i}.entsend(:,1)==k); 
        tm = dmd{i}.entsend(ii,:);
        tm(:,1) = i;
        dmd{k}.entrecv = [dmd{k}.entrecv; tm];        
    end
end

for i = 1:nproc     
    dmd{i}.entsend = unique(dmd{i}.entsend,'rows');
    dmd{i}.entrecv = unique(dmd{i}.entrecv,'rows');
end

for i = 1:nproc 
    dmd{i}.elemsend 
    dmd{i}.elemrecv 
end

for i = 1:nproc     
    dmd{i}.entsend
    dmd{i}.entrecv
end


function plotpart(mesh, np, elempart, entpart)
% plot the partitioned mesh 

p = mesh.p;
t = mesh.t;

bcol = [1 1 0;... % yellow 
        1 0 1;... % magneta
        0 1 1;... % cyan
        1 0 0;... % red
        0 1 0;... % green
        0 0 1;... % blue
        1,0.4,0.6;...
        0.4,0.6,1;...
       ];
figure(2); clf;
hold on;        
for i=0:np-1
    ind = find(elempart==i);
    ti = t(ind,:);
    simpplot(p,ti,[],bcol(i+1,:));                       
end
for it=1:size(t,1)
    pmid=mean(p(t(it,:),:),1);
    txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
    text(pmid(1),pmid(2),num2str(it),txtpars{:});
end
if strcmp(mesh.hybrid,'hdg')
    for it=1:size(mesh.f,1)
        pmid=mean(p(mesh.f(it,1:2),:),1);
        i = entpart(it)+1;
        txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',bcol(7-i,:)};
        text(pmid(1),pmid(2),num2str(it),txtpars{:});
    end
elseif strcmp(mesh.hybrid,'edg')
    mesh.f(:,end+1) = 1;    
    [~,~,edg] = elconnectivities(mesh);    
    mesh.f(:,end) = [];    
    for it=1:size(edg,1)
        pmid=edg(it,:);
        i = entpart(it)+1;
        txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',bcol(7-i,:)};
        text(pmid(1),pmid(2),num2str(it),txtpars{:});
    end
end


hold off;
axis equal;      
axis tight;
axis on;  







