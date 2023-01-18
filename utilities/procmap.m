function meshp = procmap(mesh,epart)
% MESHP = PROCMAP(MESH,EPART) 
% This function generates mesh structure MESHP for each subdomain according to the domain
% partition EPART. MESHP also includes subdomain-to-subdomain connectivities.
%
% MESH  : mesh structure for the whole domain
% EPART : domain partition from metis 
% MESHP : mesh structure for each subdomain

t2t = mkt2t(mesh.t,mesh.elemtype); % element-to-element connectivity for the whole domain

epart = epart+1; % domain partition from metis 
np = max(epart); % number of subdomains     
npf = size(mesh.perm,1); % number of points on one face
nfe = size(mesh.perm,2); % number of faces on one element

part = cell(np,1); % reordered partition for each subdomain
fpart = cell(np,1);
meshp = cell(np,1);
f = cell(np,1);
t2f = cell(np,1);
for i=1:np
    ind = find(epart==i); % list of elements in subdomain i               
    ti = mesh.t(ind,:);   % element connectivity of subdomain i
    t2ti = mkt2t(ti,mesh.elemtype); % element-to-element connectivity of subdomain i
    tj = []; % find boundary elements
    for j=1:size(t2ti,2)
        tj = [tj; find(t2ti(:,j)==0)];
    end
    tb = unique(tj); % list of boundary elements
    tk = []; % list of interface elements
    for j=1:length(tb)
        e = ind(tb(j)); % element e
        en = t2t(e,:);  % neighboring elements of element e
        ek = en>0;      % remove zero from en
        em = en(ek);    % neighboring elements of element e
        ei = ismember(em,ind);
        if all(ei)==0 % if any neighboring elements does not belong to subdomain i
            tk = [tk; tb(j)]; % then add this interface element to the list
        end
    end    
    tb = [tk; setdiff(tb,tk)]; % reorder tb to put interface elements first 
    tm = setdiff(1:length(ind),tb); % list of interior elements        
    part{i} = [ind(tb); ind(tm)]; % boundary elements then interior elements     
    
    % generate face-to-element and element-to-face connectivities for subdomain i 
    [f{i},t2f{i}] = mkt2f(mesh.t(part{i},:),mesh.elemtype);     
    
    nf = size(f{i},1);    
    % compute domain partitions for faces    
    fpart{i} = zeros(nf,1); % face partioning
    for n=1:nf
        p = sort(f{i}(n,1:end-2)); % nodes of face n
        [a,b] = ismember(p,sort(mesh.f(:,1:end-2),2),'rows');    
        if a==1
            fpart{i}(n) = b;
        else
            error('something wrong');
        end
    end
    
    % fixing f
    for n=1:nf
        p = f{i}(n,1:end-2);
        a = ismember(p,mesh.f(:,1:end-2),'rows');            
        if a==0             
            b = fpart{i}(n);
            f{i}(n,1:end-2) = mesh.f(b,1:end-2);
            e1 = part{i}(f{i}(n,end-1));            
            e2 = mesh.f(b,end-1);
            if (e1~=e2) && f{i}(n,end)>0
                f{i}(n,end-1:end) = f{i}(n,end:-1:end-1);
            end 
%             [i n e1 e2]
%             pause
%             [i n f{i}(n,:)]
%             [i n mesh.f(b,:)]
%             pause
        end
    end
    
    % check fpart and t2f
    e = fpart{i}(t2f{i})-mesh.t2f(part{i},:);
    if max(abs(e(:)))>0
        error('something wrong')
    end
    
    % generate mesh structure for subdomain i
    meshp{i}.nd = mesh.nd;
    meshp{i}.elemtype = mesh.elemtype;
    meshp{i}.nodetype = mesh.nodetype;
    meshp{i}.porder = mesh.porder;
    meshp{i}.perm = mesh.perm;
    meshp{i}.p = mesh.p;
    meshp{i}.t = mesh.t(part{i},:); 
    meshp{i}.f = [f{i} 0*f{i}(:,1)];
    meshp{i}.t2f = t2f{i};       
    [meshp{i}.f2f, meshp{i}.fe1, meshp{i}.fe2] = mkf2f(f{i}, t2f{i});
    meshp{i}.elcon = elconnectivities(meshp{i});    
    meshp{i}.elcon = reshape(meshp{i}.elcon,[npf nfe length(ind)]);
    meshp{i}.dgnodes = mesh.dgnodes(:,:,part{i});
    meshp{i}.bf = mesh.bf(:,part{i});
    meshp{i}.eint = length(tk); % number of interface elements for subdomain i
    meshp{i}.epart = part{i};   % global indices of elements of subdomain i
    meshp{i}.plocal = mesh.plocal;
    meshp{i}.tlocal = mesh.tlocal;
    meshp{i}.f = f{i};    
    meshp{i}.fpart = fpart{i};
%     for n=1:length(ind)
%         for j=1:nfe            
%             meshp{i}.elcon(:,j,n) = meshp{i}.elcon(:,j,n) - (meshp{i}.t2f(n,j)-1)*npf + (j-1)*npf;        
%         end
%     end
end

 
nf = size(mesh.f,1);
pmap = zeros(np,np);
emap = cell(np,np);
fmap = cell(np,np);
gmap = cell(np,np);
lmap = cell(np,np);
for i=1:np
    for j=1:np
        emap{i}{j} = [];
        fmap{i}{j} = [];
        gmap{i}{j} = [];
        lmap{i}{j} = [];
    end
end
for n=1:nf % for each face n
    fn = mesh.f(n,end-1:end); % two elements sharing face n
    if fn(2)>0 % if face n is interior, then do
        for i=1:np % for each subdomain
            i1 = find(part{i}==fn(1)); % find subdomain i that contains element fn(1)
            if isempty(i1)==0
                break;
            end
        end
        for j=1:np % for each subdomain           
            i2 = find(part{j}==fn(2)); % find subdomain j that contains element fn(2)
            if isempty(i2)==0
                break;
            end
        end        
        if i~=j % if subdomain i differs from subdomain j, then it is an interface face
            pmap(i,j) = 1; % connect subdomain i to subdomain j
            pmap(j,i) = 1; % connect subdomain j to subdomain i
            p = sort(mesh.f(n,1:end-2)); % nodes of face n
            [~,mi] = ismember(p,sort(f{i}(:,1:end-2),2),'rows'); % local position of face n in f{i}
            [~,mj] = ismember(p,sort(f{j}(:,1:end-2),2),'rows'); % local position of face n in f{j}            
            ia = find(t2f{i}(f{i}(mi,end-1),:)==mi);
            ii = [ia setdiff(1:nfe,ia)];            
            ja = find(t2f{j}(f{j}(mj,end-1),:)==mj);                        
            jj = [ja setdiff(1:nfe,ja)];
            emap{i}{j} = [emap{i}{j}; [i1 jj]];
            emap{j}{i} = [emap{j}{i}; [i2 ii]];
            fmap{i}{j} = [fmap{i}{j}; t2f{i}(f{i}(mi,end-1),ii)];
            fmap{j}{i} = [fmap{j}{i}; t2f{j}(f{j}(mj,end-1),jj)];                          
            gmap{i}{j} = [gmap{i}{j}; [fn(1) fn(2)]];
            gmap{j}{i} = [gmap{j}{i}; [fn(2) fn(1)]];            
            lmap{i}{j} = [lmap{i}{j}; [i1 i2]];
            lmap{j}{i} = [lmap{j}{i}; [i2 i1]];            
% %             fmap{i}{j} = [fmap{i}{j}; n];
% %             fmap{j}{i} = [fmap{j}{i}; n];  
%             p = sort(mesh.f(n,1:end-2));
%             [~,mi] = ismember(p,sort(f{i}(:,1:end-2),2),'rows');
%             [~,mj] = ismember(p,sort(f{j}(:,1:end-2),2),'rows');            
%             fmap{i}{j} = [fmap{i}{j}; mi];
%             fmap{j}{i} = [fmap{j}{i}; mj];              
%             [i1 mi t2f{i}(f{i}(mi,end-1),:)]
%             [i2 mj t2f{j}(f{j}(mj,end-1),:)]
%             pause
%             f{j}(mj,:)
%             [mi mj]
%             [i, j]
%             pause
        end
    end
end

% subdomain-to-subdomain connectivities
for i=1:np    
    meshp{i}.fint = 0; % number of faces on the interface for each subdomain   
    for j=1:np
        meshp{i}.pmap = find(pmap(i,:)==1);
        meshp{i}.emap{j} = emap{i}{j};
        meshp{i}.fmap{j} = fmap{i}{j};
        if pmap(i,j)==1
            meshp{i}.fint = meshp{i}.fint + size(fmap{i}{j},1);
        end
    end
end

% plot the mesh structures 
bcol = [1 1 0;... % yellow 
        1 0 1;... % magenta
        0 1 1;... % cyan
        1 0 0;... % red
        0 1 0;... % green 
        0 0 1;... % blue
        1,0.4,0.6;...
        0.4,0.6,1;...
       ];
   
figure(1); clf;
hold on;            
for i=1:np
    simpplot(mesh.p,mesh.t(part{i},:),[],bcol(i,:));                           
end
for i=1:np-1
    for j=i+1:np
        e = gmap{i}{j}(:);
        for it=1:length(e)
            pmid=mean(mesh.p(mesh.t(e(it),:),:),1);
            txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
            text(pmid(1),pmid(2),num2str(e(it)),txtpars{:});
        end        
    end
end
hold off;
axis equal;      
axis tight;
axis off;  

figure(2); clf;
hold on;            
for i=1:np
    simpplot(mesh.p,mesh.t(part{i},:),[],bcol(i,:));                           
end
for i=1:np-1
    for j=i+1:np
        e = gmap{i}{j}(:);
        l = lmap{i}{j}(:);
        for it=1:length(e)
            pmid=mean(mesh.p(mesh.t(e(it),:),:),1);
            txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
            text(pmid(1),pmid(2),num2str(l(it)),txtpars{:});
        end        
    end
end
hold off;
axis equal;      
axis tight;
axis off;  



% function [pmap,gmap,lmap,ind] = procmap(mesh,epart)
% 
% epart = epart+1;
% np = max(epart);
% 
% bcol = [1 1 0;...  
%         1 0 1;...
%         0 1 1;...
%         1 0 0;...
%         0 1 0;...
%         0 0 1;...
%         1,0.4,0.6;...
%         0.4,0.6,1;...
%        ];
% 
% figure(1); clf;
% hold on;            
% for i=1:np
%     ind{i} = find(epart==i);        
% %     [f{i},t2f{i}] = mkt2f(ti,mesh.elemtype);
% %     ane{i} = find(f{i}(:,end)==0);
% %     bne{i} = f{i}(ane{i},end-1); % local index
% %     cne{i} = ind{i}(bne{i}); % global index of the elements on the subdomain boundary
% %     tme = setdiff(1:length(ind{i}),bne{i});
% %     ind{i} = [cne{i}; ind{i}(tme)];
%     simpplot(mesh.p,mesh.t(ind{i},:),[],bcol(i,:));                           
%     ti = mesh.t(ind{i},:);
%     t2t = mkt2t(ti,mesh.elemtype);
%     tj = [];
%     for j=1:size(t2t,2)
%         tj = [tj; find(t2t(:,j)==0)];
%     end
%     tj = unique(tj); % list of boundary elements
%     tm = setdiff(1:length(ind{i}),tj);
%     e = ind{i}(tj);
%     ind{i} = [e; ind{i}(tm)];
% %     for it=1:length(e)
% %         pmid=mean(mesh.p(mesh.t(e(it),:),:),1);
% %         txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
% %         text(pmid(1),pmid(2),num2str(e(it)),txtpars{:});
% %     end            
% end
% hold off;
% axis equal;      
% axis tight;
% axis off;  
% 
% nf = size(mesh.f,1);
% pmap = zeros(np,np);
% for i=1:np
%     for j=1:np
%         gmap{i}{j} = [];
%         lmap{i}{j} = [];
%     end
% end
% for n=1:nf
%     fn = mesh.f(n,end-1:end);
%     if fn(2)>0
%         for i=1:np
%             i1 = find(ind{i}==fn(1));
%             if isempty(i1)==0
%                 break;
%             end
%         end
%         for j=1:np            
%             i2 = find(ind{j}==fn(2));
%             if isempty(i2)==0
%                 break;
%             end
%         end
%         if i~=j
%             pmap(i,j) = 1;
%             pmap(j,i) = 1;
%             gmap{i}{j} = [gmap{i}{j}; [fn(1) fn(2)]];
%             gmap{j}{i} = [gmap{j}{i}; [fn(2) fn(1)]];            
%             lmap{i}{j} = [lmap{i}{j}; [i1 i2]];
%             lmap{j}{i} = [lmap{j}{i}; [i2 i1]];
%         end
%     end
% end
% 
% % for it=1:size(mesh.t,1)
% %     pmid=mean(mesh.p(mesh.t(it,:),:),1);
% %     txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
% %     text(pmid(1),pmid(2),num2str(it),txtpars{:});
% % end
% 
% for i=1:np-1
%     for j=i+1:np
%         e = gmap{i}{j}(:);
%         for it=1:length(e)
%             pmid=mean(mesh.p(mesh.t(e(it),:),:),1);
%             txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
%             text(pmid(1),pmid(2),num2str(e(it)),txtpars{:});
%         end        
%     end
% end
% 
% % for i=1:np
% %     for j=1:length(cne{i})
% %         bj = bne{i}(j); % local index of the element j
% %         cj = cne{i}(j); % global index of the element j
% %         aj = ane{i}(j); % local index of the boundary face  
% %         dj = find(t2f{i}(bj,:) == aj); 
% %         fj = mesh.t2f(cj,dj); % global index of the boundary face
% %         ej = mesh.f(fj,end-1:end);
% %         if ej(2)>0 
% %             
% %         end
% %     end
% % end