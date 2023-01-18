
function meshp = domaindecomposition(mesh,nproc,RASlevel)
% Domain decomposition for 0-, 1- or 2-element overlap preconditioner

if nproc < 2; error('Domain decomposition required at least 2 processors.'); end

if strcmp(mesh.hybrid,'hdg')
    meshp = domaindecomposition_hdg(mesh,nproc,RASlevel);
elseif strcmp(mesh.hybrid,'edg')    
    meshp = domaindecomposition_edg(mesh,nproc,RASlevel);    
end

[npf,nfe] = size(mesh.perm);
for i=1:nproc % loop over each subdomain
    meshp{i}.perm = mesh.perm;
    meshp{i}.permgeom = mesh.permgeom;
    meshp{i}.dgnodes = mesh.dgnodes(:,:,meshp{i}.elempart);
    meshp{i}.bf = mesh.bf(:,meshp{i}.elempart);
    % compute local element-to-entity connectivities
    elcon = reshape(mesh.elcon,[npf nfe mesh.ne]);
    meshp{i}.elcon = elcon(:,:,meshp{i}.elempart);
    ne = length(meshp{i}.elempart);        
    %reshape(meshp{i}.elcon,[npf*nfe ne])
    if strcmp(mesh.hybrid,'hdg')
        t2f = mesh.t2f(meshp{i}.elempart,:);    
        for j = 1:ne
            for k = 1:nfe
                t2f(j,k) = find(t2f(j,k) == meshp{i}.entpart);
            end
        end
        for j = 1:ne
            for k = 1:nfe
                ii = elcon(:,k,meshp{i}.elempart(j)) - (mesh.t2f(meshp{i}.elempart(j),k)-1)*npf;                
                ind = ((t2f(j,k)-1)*npf+1):(t2f(j,k)*npf);
                meshp{i}.elcon(:,k,j) = ind(ii);
            end
        end
        meshp{i}.t2f = t2f;
    elseif strcmp(mesh.hybrid,'edg')        
        for j = 1:ne
            for k = 1:nfe        
                for m = 1:npf               
                    meshp{i}.elcon(m,k,j) = find(meshp{i}.elcon(m,k,j) == meshp{i}.entpart);
                end
            end
        end
        [~,meshp{i}.t2f] = mkt2f(mesh.t(meshp{i}.elempart,:),mesh.elemtype);
    end        
    meshp{i}.elcon = reshape(meshp{i}.elcon,[npf*nfe ne]);   
    meshp{i} = block_crs(meshp{i},mesh.hybrid);
    
    nownent = sum(meshp{i}.entpartpts(1:2));
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
            
    meshp{i}.nd  = mesh.nd;
    meshp{i}.ne  = size(meshp{i}.dgnodes,3);
    meshp{i}.nf  = max(meshp{i}.t2f(:));
    meshp{i}.ncd = size(meshp{i}.dgnodes,2);
    meshp{i}.nfe = size(meshp{i}.t2f,2);
    meshp{i}.nve = size(mesh.t,2);
    meshp{i}.nvf = size(mesh.tlocfc,2);
    t = mesh.t(meshp{i}.elempart,:);
    p = mesh.p(unique(t(:)),:);
    meshp{i}.nv  = size(p,1);
    if strcmp(mesh.hybrid,'hdg')
        npf = size(mesh.plocfc,1);
        meshp{i}.ndh = length(meshp{i}.entpart)*npf;
    elseif strcmp(mesh.hybrid,'edg')
        meshp{i}.ndh = length(meshp{i}.entpart);
    end
    meshp{i}.blkSize = mesh.blkSize;
    meshp{i}.BJ_nrows = length(meshp{i}.BJ_rowpts)-1;
    meshp{i}.BJ_nblks = length(meshp{i}.BJ_colind);
    meshp{i}.BK_nrows = length(meshp{i}.BK_rowpts)-1;
    meshp{i}.BK_nblks = length(meshp{i}.BK_colind);

    for j = 1:length(meshp{i}.nbsd)
        meshp{i}.entsendpts(j) = length(find(meshp{i}.entsend(:,1)==meshp{i}.nbsd(j)));
        meshp{i}.entrecvpts(j) = length(find(meshp{i}.entrecv(:,1)==meshp{i}.nbsd(j)));
        meshp{i}.elemsendpts(j) = length(find(meshp{i}.elemsend(:,1)==meshp{i}.nbsd(j)));
        meshp{i}.elemrecvpts(j) = length(find(meshp{i}.elemrecv(:,1)==meshp{i}.nbsd(j)));
    end
    meshp{i}.entsend = meshp{i}.entsend(:,2);
    meshp{i}.entrecv = meshp{i}.entrecv(:,2);
    meshp{i}.elemsend = meshp{i}.elemsend(:,2);
    meshp{i}.elemrecv = meshp{i}.elemrecv(:,2);
end

%[ae{i},fe{i},dudg{i},dudg_duh{i}] = hdg_elemental(masterp{i},appp{i},dgnodes{i},bf{i},udg{i},uh{i},sh{i});
    
function dmd = domaindecomposition_hdg(mesh,nproc,RASlevel)

[elempart,entpart] = meshpart(mesh.t2f,nproc);

% plot domain partition
%plotpart(mesh, nproc, elempart, entpart);

elempart = elempart+1;
entpart = entpart+1;

% loop over each subdomain
for i=1:nproc
    
    % list of elements in subdomain i                       
    indelem{i} = find(elempart==i); 
    
    % list of faces in subdomain i                       
    indent{i} = find(entpart==i); 
end

t2t = mkt2t(mesh.t,mesh.elemtype);
% loop over each subdomain
for i=1:nproc    
    % find all interface faces
    f2f = mesh.f2f(:,indent{i});
    i1 = find(f2f(:)>0);
    in = ones(size(f2f));
    in(i1) = ismember(f2f(i1),indent{i});
    [~,j1] = find(in==0);
    j1 = unique(j1);

    % list of all interface faces
    intf = indent{i}(j1);
    
    % find all faces that are connected to interface faces
    f2f = f2f(:,j1);
    i1 = find(f2f(:)>0);
    in = ones(size(f2f));
    in(i1) = ismember(f2f(i1),indent{i});
    f0 = unique(f2f(in==0));  % faces do not belong to subdomain i
    %f1 = unique(f2f(in==1));
    %f1 = f1(f1>0);            % faces belong to subdomain i        
    
    if RASlevel==2
        % find additional faces for RAS-2 preconditioner
        f2f = mesh.f2f(:,f0);
        i1 = find(f2f(:)>0);
        in = ones(size(f2f));
        in(i1) = ismember(f2f(i1),indent{i});
        f02 = unique(f2f(in==0));  % faces do not belong to subdomain i
        f0 = unique([f0; f02]);
    end
    
    % [interior faces, own interface faces, other faces]   
    dmd{i}.entpart = [setdiff(indent{i},intf); intf; f0];  
    dmd{i}.entpartpts = [length(indent{i})-length(intf) length(intf) length(f0)];    
        
    % store faces received from neighboring subdomains to perform the matrix-vector product
    dmd{i}.entrecv = [0*f0 length(indent{i})+(1:1:length(f0))' f0];
    for k = 1:nproc
        if k ~= i
            in = ismember(f0,indent{k});
            dmd{i}.entrecv(in,1) = k;
        end
    end
    dmd{i}.entrecv = unique(dmd{i}.entrecv,'rows');    
    dmd{i}.nbsd = unique(dmd{i}.entrecv(:,1))';        
    if dmd{i}.nbsd(1)==0
        error('something wrong');
    end
    
    % find all elements that are connected to interface faces
    elem = mesh.f(intf,end-1:end);
    i1 = find(elem(:)>0);
    in = ones(size(elem));
    in(i1) = ismember(elem(i1),indelem{i});
    elem0 = unique(elem(in==0)); % elements do not belong to subdomain i
    elem1 = unique(elem(in==1));
    elem1 = elem1(elem1>0);      % elements belong to subdomain i
    
    if RASlevel==2
        % find additional elements for RAS-2 preconditioner
        elem2 = t2t([elem1; elem0],:);    
        elem2 = unique(elem2(:));
        i1 = elem2==0;
        elem2(i1) = [];
        in = ismember(elem2,indelem{i});
        elem02 = unique(elem2(in==0)); % elements do not belong to subdomain i
        elem0 = unique([elem0; elem02]);
    end
    
    % [interior elements, interface elements, other elements]  
    dmd{i}.elempart = [setdiff(indelem{i},elem1); elem1; elem0];
    dmd{i}.elempartpts = [length(indelem{i})-length(elem1) length(elem1) length(elem0)];    
    
    % store elements received from neighboring subdomains to perform the matrix-vector assembly
    dmd{i}.elemrecv = [0*elem0 length(indelem{i})+(1:1:length(elem0))' elem0];
    for k = 1:nproc
        if k ~= i
            in = ismember(elem0,indelem{k});
            dmd{i}.elemrecv(in,1) = k;
        end
    end
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');        
end

% store faces sent to neighboring subdomains to perform the matrix-vector product
for k = 1:nproc
    dmd{k}.entsend = [];
end
for i = 1:nproc          
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.entrecv(:,1)==k;
        tm = dmd{i}.entrecv(ii,:);
        tm(:,1) = i;                
        for m = 1:size(tm,1)            
            tm(m,2) = find(dmd{k}.entpart==tm(m,3));
        end
        dmd{k}.entsend = [dmd{k}.entsend; tm];        
    end    
end

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc
    dmd{k}.elemsend = [];
end
for i = 1:nproc        
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.elemrecv(:,1)==k;
        tm = dmd{i}.elemrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)
            tm(m,2) = find(dmd{k}.elempart==tm(m,3));
        end
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end    
end

for i = 1:nproc 
    dmd{i}.elempart'
    dmd{i}.elemsend 
    dmd{i}.elemrecv 
    dmd{i}.entpart'
    dmd{i}.entsend
    dmd{i}.entrecv
end

    
function dmd = domaindecomposition_edg(mesh,nproc,RASlevel)

if RASlevel==2; error('Domain decomposition for EDG and 2+ element overlap not implemented yet.'); end

[elempart,entpart] = meshpart(mesh.elcon',nproc);

% plot domain partition
% plotpart(mesh, nproc, elempart, entpart);

elempart = elempart+1;
entpart = entpart+1;

% loop over each subdomain
for i=1:nproc
    
    % list of elements in subdomain i                       
    indelem{i} = find(elempart==i); 
    
    % list of edg nodes in subdomain i                       
    indent{i} = find(entpart==i); 
end

% loop over each subdomain
for i=1:nproc
    % find all interface edgnodes    
    edgnumber = zeros(length(indent{i}),1); 
    for j = 1:length(indent{i})
        [~,j1] = find(mesh.elcon==indent{i}(j));
        j1 = unique(j1);
        edg1 = mesh.elcon(:,j1);
        edg1 = unique(edg1(:));     % list of neirghboring edg nodes
        if all(ismember(edg1,indent{i}))
            % edgnode indent{i}(j) is fully inside the subdomain i (all
            % neighboring edg nodes are in the subdomain i)
            edgnumber(j) = 2;
        else
            % edgnode indent{i}(j) is on the interface between two subdomains
            edgnumber(j) = 1;            
        end
    end
        
    % list of all interface edgnodes in the subdomain
    intf = indent{i}(edgnumber==1);
    
    % find all neighboring edgnodes that are connected to interface edgnodes 
    edg0 = [];
    for j = 1:length(intf)
        [~,j1]=find(mesh.elcon==intf(j));
        j1 = unique(j1);
        edg1 = mesh.elcon(:,j1);
        edg1 = unique(edg1(:));
        in = ismember(edg1,indent{i});
        edg0 = [edg0; edg1(in==0)];
    end
    edg0 = unique(edg0);        % list of all edg nodes that are not in the subdomain but are connected to edg nodes in the subdomain
    
    % [interior edgnodes, own interface edgnodes, other edgnodes (neighbors in other processors)]   
    dmd{i}.entpart = [setdiff(indent{i},intf); intf; edg0];  
    dmd{i}.entpartpts = [length(indent{i})-length(intf) length(intf) length(edg0)];    
    
    % store edgnodes received from neighboring subdomains to perform the matrix-vector product
    dmd{i}.entrecv = [0*edg0 length(indent{i})+(1:1:length(edg0))' edg0];
    for k = 1:nproc
        if k ~= i
            in = ismember(edg0,indent{k});
            dmd{i}.entrecv(in,1) = k;
        end
    end
    dmd{i}.entrecv = unique(dmd{i}.entrecv,'rows');    
    dmd{i}.nbsd = unique(dmd{i}.entrecv(:,1))';     % list of subdomains from which info is required for matrix-vector product
    
    % find all elements that are connected to interface edgnodes
    elem0 = [];     % elements not in the process
    elem1 = [];     % elements in the process that neighbor edg nodes NOT in the process
    for j = 1:length(intf)
        [~,j1]=find(mesh.elcon==intf(j));
        j1 = unique(j1);        
        in = ismember(j1,indelem{i});
        elem0 = [elem0; j1(in==0)];
        elem1 = [elem1; j1(in==1)];
    end
    elem0 = unique(elem0);
    elem1 = unique(elem1);
    in = zeros(length(elem1),1);        
    for j = 1:length(elem1)        
        if all(ismember(mesh.elcon(:,elem1(j)),indent{i}))            
            in(j) = 1;
        end
    end  
    elem1 = elem1(in==0);    
    
    % [interior elements, interface elements, other elements (neighbors in other processors)]  
    dmd{i}.elempart = [setdiff(indelem{i},elem1); elem1; elem0];
    dmd{i}.elempartpts = [length(indelem{i})-length(elem1) length(elem1) length(elem0)];    
    
    % store elements received from neighboring subdomains to assemble linear system
    dmd{i}.elemrecv = [0*elem0 length(indelem{i})+(1:1:length(elem0))' elem0];
    for k = 1:nproc
        if k ~= i
            in = ismember(elem0,indelem{k});
            dmd{i}.elemrecv(in,1) = k;
        end
    end
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');        
end

% store edgnodes sent to neighboring subdomains to perform the matrix-vector product
for k = 1:nproc
    dmd{k}.entsend = [];
end
for i = 1:nproc
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.entrecv(:,1)==k;
        tm = dmd{i}.entrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)
            tm(m,2) = find(dmd{k}.entpart==tm(m,3));
        end
        dmd{k}.entsend = [dmd{k}.entsend; tm];        
    end    
end

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc
    dmd{k}.elemsend = [];
end
for i = 1:nproc        
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.elemrecv(:,1)==k;
        tm = dmd{i}.elemrecv(ii,:);
        tm(:,1) = i;        
        for m = 1:size(tm,1)            
            tm(m,2) = find(dmd{k}.elempart==tm(m,3));
        end
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end    
end

for i = 1:nproc 
    dmd{i}.elempart'
    dmd{i}.elemsend 
    dmd{i}.elemrecv 
    dmd{i}.entpart'
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
