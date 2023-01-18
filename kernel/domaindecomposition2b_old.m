
function [mesh,meshp] = domaindecomposition2b(mesh,nproc,DDobjectiveFlag,ent2entWeight)
% Domain decomposition for 0- or 1-element overlap preconditioner

if nargin<3; DDobjectiveFlag = 0; end
% DDflag:
%       0: Minimize edge cut
%       1: Minimized weighted edge-cut

if strcmp(mesh.hybrid,'hdg')
    if nargin < 4
        [mesh,meshp] = domaindecomposition_hdg(mesh,nproc,DDobjectiveFlag);
    else
        [mesh,meshp] = domaindecomposition_hdg(mesh,nproc,DDobjectiveFlag,ent2entWeight);
    end
elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
    if nargin < 4
        [mesh,meshp] = domaindecomposition_edg(mesh,nproc,DDobjectiveFlag);
    else
        [mesh,meshp] = domaindecomposition_edg(mesh,nproc,DDobjectiveFlag,ent2entWeight);
    end 
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
    elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
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
        if i1(end) ~= length(i1)
            error('ent2ent has neighboring entities in the middle of a row, instead of at the end');
        end
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
    elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
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
    meshp{i}.ent2entWeightLen = 2*(meshp{i}.BJ_nrows + meshp{i}.BJ_nblks + meshp{i}.BK_nblks);
end

%[ae{i},fe{i},dudg{i},dudg_duh{i}] = hdg_elemental(masterp{i},appp{i},dgnodes{i},bf{i},udg{i},uh{i},sh{i});

if strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
    for i=1:nproc % loop over each subdomain
        meshp{i}.matrecv = [];
        meshp{i}.matsend = [];
    end    
    for i=1:nproc
        for j = 1:size(dmd{i}.entrecv,1)
            ncpu = dmd{i}.entrecv(j,1);
            edgj = dmd{i}.entrecv(j,2);
            jr = (meshp{i}.cbsr_rowpts(edgj)+1):meshp{i}.cbsr_rowpts(edgj+1);
            edgk = meshp{i}.cbsr_colind(jr);    
            edgjg = meshp{i}.entpart(edgj);
            edgkg = meshp{i}.entpart(edgk);
            [~,je] = find(mesh.elcon==edgjg);            
            for k = 1:length(edgkg)
                [~,ke] = find(mesh.elcon==edgkg(k));
                jke = intersect(je,ke);
                if all(ismember(jke,meshp{i}.elempart))==0
                    meshp{i}.matrecv = [meshp{i}.matrecv; ncpu jr(k)]; 
                    m = find(dmd{ncpu}.entsend(:,3) == edgjg);                    
                    if dmd{ncpu}.entsend(m,1)~=i
                        error('something wrong');                        
                    end
                    nedgj = dmd{ncpu}.entsend(m,2);                                
                    njr = (meshp{ncpu}.cbsr_rowpts(nedgj)+1):meshp{ncpu}.cbsr_rowpts(nedgj+1);
                    nedgk = meshp{ncpu}.cbsr_colind(njr);    
                    nedgkg = meshp{ncpu}.entpart(nedgk);
                    nk = find(nedgkg==edgkg(k));
                    meshp{ncpu}.matsend = [meshp{ncpu}.matsend; i njr(nk)]; 
                end
            end
        end        
    end
end
for i=1:nproc
    for j = 1:length(meshp{i}.nbsd)
        meshp{i}.matsendpts(j) = length(find(meshp{i}.matsend(:,1)==meshp{i}.nbsd(j)));
        meshp{i}.matrecvpts(j) = length(find(meshp{i}.matrecv(:,1)==meshp{i}.nbsd(j)));
    end
    meshp{i}.matsend = meshp{i}.matsend(:,2);
    meshp{i}.matrecv = meshp{i}.matrecv(:,2);    
end

function [mesh,dmd] = domaindecomposition_hdg(mesh,nproc,DDobjectiveFlag,ent2entWeight)

if nargin < 3; DDobjectiveFlag = 0; end

if DDobjectiveFlag == 0
    [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart(mesh,nproc);
elseif DDobjectiveFlag == 1
    if nargin < 4;
        [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart_v2(mesh,nproc);
    else
        [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart_v2(mesh,nproc,ent2entWeight);
    end
else
    error('DDobjectiveFlag has invalid value');
end

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

% Check that all elements in the processor have at least one entity in the
% processor
for i=1:nproc
    elementsToRemove = [];
    for j=1:length(indelem{i})
        if length(setdiff(indent{i},mesh.t2f(indelem{i}(j),:))) == length(indent{i})
            entitiesInElement = mesh.t2f(indelem{i}(j),:);
            candidateProcessors = entpart(entitiesInElement);
            k = mode(candidateProcessors);
            indelem{k} = sort([indelem{k};indelem{i}(j)]);
            elementsToRemove = [elementsToRemove;indelem{i}(j)];
            warning('An element in the processor had did not have any entity in the processors. The issue has been fixed successfully.');
        end
    end
    indelem{i} = setdiff(indelem{i},elementsToRemove);
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
        
    % find all elements that are connected to interface faces
    elem = mesh.f(intf,end-1:end);
    i1 = find(elem(:)>0);
    in = ones(size(elem));
    in(i1) = ismember(elem(i1),indelem{i});
    elem0 = unique(elem(in==0)); % elements do not belong to subdomain i
    elem1 = unique(elem(in==1));
    elem1 = elem1(elem1>0);      % elements belong to subdomain i
        
    % find additional elements for RAS-1 preconditioner
    elem2 = t2t(elem1,:);    
    elem2 = unique(elem2(:));
    i1 = elem2==0;
    elem2(i1) = [];
    in = ismember(elem2,indelem{i});
    elem02 = unique(elem2(in==0)); % elements do not belong to subdomain i
    elem0 = unique([elem0; elem02]);    
    
    % [interior elements, interface elements, other elements]  
    elemi  = setdiff(indelem{i},elem1);
    elemj = t2t(elem1,:);    
    elemj = unique(elemj(:));
    i1 = elemj==0;
    elemj(i1) = [];
    in = ismember(elemj,elemi);
    elem2 = unique(elemj(in==1));     
    
    dmd{i}.elempart = [setdiff(elemi,elem2); [elem2; elem1]; elem0];
    dmd{i}.elempartpts = [length(elemi)-length(elem2) length(elem2)+length(elem1) length(elem0)];        
%     dmd{i}.elempart = [setdiff(indelem{i},elem1); elem1; elem0];
%     dmd{i}.elempartpts = [length(indelem{i})-length(elem1) length(elem1) length(elem0)];    
    
    % store elements received from neighboring subdomains to perform the matrix-vector assembly
    dmd{i}.elemrecv = [0*elem0 length(indelem{i})+(1:1:length(elem0))' elem0];
    for k = 1:nproc
        if k ~= i
            in = ismember(elem0,indelem{k});
            dmd{i}.elemrecv(in,1) = k;
        end
    end
    dmd{i}.elemrecv = unique(dmd{i}.elemrecv,'rows');            
    
    % find additional faces for RAS-1 preconditioner
    f2 = mesh.t2f(elem0,:);
    in = ismember(f2,indent{i});
    f02 = unique(f2(in==0));  % faces do not belong to subdomain i
    f0 = unique([f0; f02]);
    
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

% store faces sent to neigboring cpus to assemble the preconditioner
for k = 1:nproc
    dmd{k}.matrecv = dmd{k}.entrecv;
    on = [];
    for i = 1:size(dmd{k}.matrecv,1)
        fi = dmd{k}.matrecv(i,3);
        ei = mesh.f(fi,end-1:end);
        if (ei(2)<=0) || (all(ismember(ei,dmd{k}.elempart))==1)
            on = [on i];                    
        end
    end    
    dmd{k}.matrecv(on,:) = [];
end
for k = 1:nproc
    dmd{k}.matsend = [];
end
for i = 1:nproc          
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.matrecv(:,1)==k;
        tm = dmd{i}.matrecv(ii,:);
        tm(:,1) = i;                
        for m = 1:size(tm,1)            
            tm(m,2) = find(dmd{k}.entpart==tm(m,3));
        end
        dmd{k}.matsend = [dmd{k}.matsend; tm];        
    end    
end

% for i = 1:nproc 
%     dmd{i}.elempart'
%     dmd{i}.elemsend 
%     dmd{i}.elemrecv 
%     dmd{i}.entpart'
%     dmd{i}.entsend
%     dmd{i}.entrecv
% end

    
function [mesh,dmd] = domaindecomposition_edg(mesh,nproc,DDobjectiveFlag,ent2entWeight)

if nargin < 3; DDobjectiveFlag = 0; end

if DDobjectiveFlag == 0
    [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart(mesh,nproc);
elseif DDobjectiveFlag == 1
    if nargin < 4;
        [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart_v2(mesh,nproc);
    else
        [elempart, entpart, mesh.globalEnt2ent, mesh.globalEnt2entStart] = meshpart_v2(mesh,nproc,ent2entWeight);
    end
else
    error('DDobjectiveFlag has invalid value');
end

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

% Check that all elements in the processor have at least one entity in the
% processor
for i=1:nproc
    elementsToRemove = [];
    for j=1:length(indelem{i})
        if length(setdiff(indent{i},mesh.elcon(:,indelem{i}(j)))) == length(indent{i})
            entitiesInElement = mesh.elcon(:,indelem{i}(j));
            candidateProcessors = entpart(entitiesInElement);
            k = mode(candidateProcessors);
            indelem{k} = sort([indelem{k};indelem{i}(j)]);
            elementsToRemove = [elementsToRemove;indelem{i}(j)];
            warning('An element in the processor had did not have any entity in the processors. The issue has been fixed successfully.');
        end
    end
    indelem{i} = setdiff(indelem{i},elementsToRemove);
end

t2t = mkt2t(mesh.t,mesh.elemtype);

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
    elemi  = setdiff(indelem{i},elem1);
    elemj = t2t(elem1,:);    
    elemj = unique(elemj(:));
    i1 = elemj==0;
    elemj(i1) = [];
    in = ismember(elemj,elemi);
    elem2 = unique(elemj(in==1));         
    dmd{i}.elempart = [setdiff(elemi,elem2); [elem2; elem1]; elem0];
    dmd{i}.elempartpts = [length(elemi)-length(elem2) length(elem2)+length(elem1) length(elem0)];            
%     dmd{i}.elempart = [setdiff(indelem{i},elem1); elem1; elem0];
%     dmd{i}.elempartpts = [length(indelem{i})-length(elem1) length(elem1) length(elem0)];    
    
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

% for i = 1:nproc
%     dmd{i}.elempart'
%     dmd{i}.elemsend
%     dmd{i}.elemrecv 
%     dmd{i}.entpart'
%     dmd{i}.entsend
%     dmd{i}.entrecv
% end


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
        1 1 0;... % yellow 
        1 0 1;... % magneta
        0 1 1;... % cyan
        1 0 0;... % red
        0 1 0;... % green
        0 0 1;... % blue
        1,0.4,0.6;...
        0.4,0.6,1;...
        1 1 0;... % yellow 
        1 0 1;... % magneta
        0 1 1;... % cyan
        1 0 0;... % red
        0 1 0;... % green
        0 0 1;... % blue
        1,0.4,0.6;...
        0.4,0.6,1;...
        1 1 0;... % yellow 
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
elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
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
