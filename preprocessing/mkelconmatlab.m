function [elcon,edg,ndof] = mkelconmatlab(p,t,f,t2f,elemtype,isEDGface,hybrid,porder,dim)
%MKELCON compute element connectivities 
%   [ELCON,NDOF] = GETELCON(MESH)
%
%      MESH:      MESH STRUCTURE
%      ELCON:     Element connectivies
%      NDOF:      Number of degrees of freedom

if dim == 1
    elcon = t2f';
    if min(elcon(:)) == 0; error('Something wrong. Try increasing snap.'); end
    ndof = max(elcon(:));
    edg = [];
    return;
end

if porder==0
    npf   = 1;
    ne    = size(t,1);
    nfv   = size(t2f,2); % number of faces per element
    nf    = size(f,1);
    ncf = dim;               % number of corners of a face
    if dim == 3 && elemtype  % hex element
        ncf = dim+1;
    end    
    elcon = zeros(npf,nfv,ne);    
    for i = 1:nf        
        fe = f(i,ncf+1:ncf+2); % neighboring elements of face i        
        if1 = (t2f(fe(1),:)==i);           % location of face i on element fe(1)
        elcon(:,if1,fe(1)) = i;   % assign dof numbering of face i to elcon from element fe(1) 
        if fe(2)>0        
            if2 = (t2f(fe(2),:)==i);       % location of face i on element fe(2)
            elcon(:,if2,fe(2)) = i; % assign dof numbering of face i to elcon from element fe(2) 
        end                        
    end    
    elcon = reshape(elcon,[npf*nfv ne]);
    ndof = nf;
    edg =p(t',:);    
    edg=permute(reshape(edg,[size(t,2) ne dim]),[1 3 2]);
    if min(elcon(:)) == 0; error('Something wrong. Try increasing snap.'); end
    return;
end

if strcmp(hybrid,'hdg') 
    [elcon,ndof,edg] = hdgelcon(p,t,f,t2f,elemtype,porder,dim);    
elseif strcmp(hybrid,'edg') 
    [elcon,ndof,edg] = edgelcon(p,t,f,t2f,elemtype,porder,dim);
elseif strcmp(hybrid,'iedg')    
    [elcon,ndof,edg] = iedgelcon(p,t,f,t2f,elemtype,porder,dim);    
elseif strcmp(hybrid,'hedg')    
    [elcon,ndof,edg] = hedgelcon(p,t,f,t2f,elemtype,isEDGface,porder,dim);    
end

if min(elcon(:)) == 0; error('Something wrong. Try increasing snap.'); end

function in = xiny(x,y)
% Determine if each row of x is a member of y
% If row j of x is a member of y and x(j,:) = y(k,:) then in(j) = k
% Else in(j) = 0

[m,dim] = size(x);

in = zeros(m,1);
if dim==2
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2;
        [md,id] = min(d2);
        if md<1e-12, in(j)=id; end
    end
else
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2 + (y(:,3)-x(j,3)).^2;
        [md,id] = min(d2);
        if md<1e-12, in(j)=id; end
    end
end

function [elcon,ndof,hdg] = hdgelcon(p,t,f,t2f,elemtype,porder,dim)

ne       = size(t,1);
nfe      = size(t2f,2); % number of faces per element
nf       = size(f,1); % total number of faces
ncf = dim;               % number of corners of a face
if dim == 3 && elemtype  % hex element
    ncf = dim+1;
end

%[philocfc,philocvl] = localbasis(porder,dim,elemtype);
[philocvl,philocfc,~,~,perm] = localbasis(porder,dim,elemtype);
philocfc = philocfc{1};
perm = cell2mat(perm);

[npv,nnv] = size(philocvl);
npf       = size(philocfc,1);

% get dgnodes
dgnodes = zeros(npv,dim,ne);
for d=1:dim
  for n=1:nnv
    dp=philocvl(:,n)*p(t(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end
snap = 1e-8;
dgnodes = round(dgnodes/snap)*snap;
dgnodes = reshape(permute(dgnodes(perm,1:dim,:),[1,3,2]),[npf*nfe*ne,dim]);
dgnodes = reshape(dgnodes,[npf nfe ne dim]);

hdg = zeros(npf*nf,dim);
elcon = zeros(npf,nfe,ne);
ndof = 0;
for i = 1:nf
    fn = f(i,1:ncf);       % corners of face i 
    fe = f(i,ncf+1:ncf+2); % neighboring elements of face i            
    pf = philocfc*p(fn,:); % nodal points on face i    
    pf = round(pf/snap)*snap;
    ind = ((i-1)*npf+1):(i*npf);    
    hdg(ind,:) = pf;           % hdg nodes on faces
    
    % hdg face
    tof = ndof + (1:npf);   % dof numbering on face i
    ndof = ndof + npf;      % number of degrees of freedom     
    
    if1 = find(t2f(fe(1),:)==i);           % location of face i on element fe(1)
    dg1 = reshape(dgnodes(:,if1,fe(1),:),[npf dim]); 
    in = xiny(dg1,pf);              % match pf to dg1          
    elcon(:,if1,fe(1)) = tof(in);   % assign dof numbering of face i to elcon from element fe(1)             
    if fe(2)>0
        if2 = find(t2f(fe(2),:)==i);           % location of face i on element fe(2)
        dg2 = reshape(dgnodes(:,if2,fe(2),:),[npf dim]); 
        in = xiny(dg2,pf);              % match pf to dg1          
        elcon(:,if2,fe(2)) = tof(in);   % assign dof numbering of face i to elcon from element fe(1)                     
    end    
end
hdg = round(hdg/snap)*snap;
elcon = reshape(elcon,[npf nfe ne]);

function [elcon,ndof,edg] = edgelcon(p,t,f,t2f,elemtype,porder,dim)

[philocvl,philocfc,~,~,perm] = localbasis(porder,dim,elemtype);
philocfc = philocfc{1};
perm = cell2mat(perm);

nfe      = size(t2f,2); % number of faces per element
ne = size(t,1);
nf = size(f,1);
perm = perm(:);
bn = unique(perm);
nbn = length(bn);       % Number of cg nodes in the element faces

[npv,nnv] = size(philocvl);
[npf,nnf] = size(philocfc);

% get dgnodes on elements
dgnodes = zeros(npv,dim,ne);
for d=1:dim
  for n=1:nnv
    dp=philocvl(:,n)*p(t(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end
snap = 1e-8;
dgnodes = round(dgnodes/snap)*snap;

dgnodes = reshape(permute(dgnodes(bn,1:dim,:),[1,3,2]),[nbn*ne,dim]);
%snap = 1e-8;
%dgnodes = round(dgnodes/snap)*snap;

edg = flipdim(dgnodes,1);
[~,I] = unique(edg,'rows'); 
edg = edg(sort(I),:);
edg = flipdim(edg,1);
ndof = size(edg,1);

[~,b] = ismember(dgnodes,edg,'rows');
elcon = reshape(b,[nbn ne]);

bp = (1:1:nbn)';
pm = perm;
for i = 1:length(perm)
    j = bn == perm(i);
    pm(i) = bp(j);
end
elcon = elcon(pm,:);
elcon = reshape(elcon,[npf nfe ne]);

% % get dgnodes on faces
% dgnodes = zeros(npf,dim,nf);
% for d=1:dim
%   for n=1:nnf
%     dp=philocfc(:,n)*p(f(:,n),d)';
%     dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
%   end
% end
% snap = 1e-8;
% dgnodes = round(dgnodes/snap)*snap;
% 
% dgnodes = reshape(permute(dgnodes,[1,3,2]),[npf*nf,dim]);
% [~,b] = ismember(dgnodes,edg,'rows');
% facecon = reshape(b,[npf nf]);

function [elcon,ndof,edg] = iedgelcon(p,t,f,t2f,elemtype,porder,dim)

[philocvl,philocfc,~,~,perm] = localbasis(porder,dim,elemtype);
philocfc = philocfc{1};
perm = cell2mat(perm);

perm = perm(:);
ne       = size(t,1);
nfe      = size(t2f,2); % number of faces per element
[nf,nf2] = size(f); % total number of faces
ncf = dim;               % number of corners of a face
if dim == 3 && elemtype  % hex element
    ncf = dim+1;
end

[npv,nnv] = size(philocvl);
[npf,nnf] = size(philocfc);

% get dgnodes on faces
dgnodes = zeros(npf,dim,nf);
for d=1:dim
  for n=1:nnf
    dp=philocfc(:,n)*p(f(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end
snap = 1e-8;
dgnodes = round(dgnodes/snap)*snap;

% get edg nodes
iedg = find(f(:,ncf+2)>0);
iedg = sort(iedg);
edg = permute(dgnodes(:,:,iedg),[1 3 2]);
edg = reshape(edg,npf*length(iedg),dim);
edg = flipdim(edg,1);
[~,I] = unique(edg,'rows');
edg = edg(sort(I),:);
edg = flipdim(edg,1);

% get dgnodes
dgnodes = zeros(npv,dim,ne);
for d=1:dim
  for n=1:nnv
    dp=philocvl(:,n)*p(t(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end
snap = 1e-8;
dgnodes = round(dgnodes/snap)*snap;

dgnodes = reshape(permute(dgnodes(perm,1:dim,:),[1,3,2]),[npf*nfe*ne,dim]);
[~,b] = ismember(dgnodes,edg,'rows');
elcon = reshape(b,[npf nfe ne]);

dgnodes = reshape(dgnodes,[npf nfe ne dim]);
%dgnodes = permute(dgnodes,[1 2 4 3]);
ihdg = find(f(:,ncf+2)<=0);
ihdg = sort(ihdg);
nbf = length(ihdg);
ndof = size(edg,1);
pdg = zeros(npf,dim,nbf);
for j = 1:nbf
    i = ihdg(j);
    fn = f(i,1:ncf);       % corners of face i 
    fe = f(i,ncf+1:ncf+2); % neighboring elements of face i
    ff = f(i,nf2);         % 0 -> hdg face or 1 -> edg face        
    pf = philocfc*p(fn,:); % nodal points on face i
    pf = round(pf/snap)*snap;
    pdg(:,:,j) = pf;
    
    if (fe(2)>0) || (ff==1)
        error('something wrong');
    end                
    
    % hdg face
    tof = ndof + (1:npf);   % dof numbering on face i
    ndof = ndof + npf;      % number of degrees of freedom     
    
    if1 = find(t2f(fe(1),:)==i);           % location of face i on element fe(1)
     % dg nodes on face i from element fe(1)  
    %dgnodes(perm(:,if1),fe(1),:)
    dg1 = reshape(dgnodes(:,if1,fe(1),:),[npf dim]); 
    in = xiny(dg1,pf);              % match pf to dg1      
    elcon(:,if1,fe(1)) = tof(in);   % assign dof numbering of face i to elcon from element fe(1)             
end
pdg = permute(pdg,[1 3 2]);
edg = [edg; reshape(pdg,[npf*nbf dim])];

elcon = reshape(elcon,[npf nfe ne]);


function [elcon,ndof,edg] = hedgelcon(p,t,f,t2f,elemtype,isEDGface,porder,dim)

[philocvl,philocfc,~,~,perm] = localbasis(porder,dim,elemtype);
philocfc = philocfc{1};
perm = cell2mat(perm);

perm = perm(:);
ne       = size(t,1);
nfe      = size(t2f,2); % number of faces per element
[nf,nf2] = size(f); % total number of faces
ncf = dim;               % number of corners of a face
if dim == 3 && elemtype  % hex element
    ncf = dim+1;
end

[npv,nnv] = size(philocvl);
[npf,nnf] = size(philocfc);


% get dgnodes on faces
dgnodes = zeros(npf,dim,nf);
for d=1:dim
  for n=1:nnf
    dp=philocfc(:,n)*p(f(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end
snap = 1e-8;
dgnodes = round(dgnodes/snap)*snap;

% get edg nodes
iedg = find(isEDGface==1);
iedg = sort(iedg);
edg = permute(dgnodes(:,:,iedg),[1 3 2]);
edg = reshape(edg,npf*length(iedg),dim);
edg = flipdim(edg,1);
[~,I] = unique(edg,'rows');
edg = edg(sort(I),:);
edg = flipdim(edg,1);

% get dgnodes
dgnodes = zeros(npv,dim,ne);
for d=1:dim
  for n=1:nnv
    dp=philocvl(:,n)*p(t(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end
snap = 1e-8;
dgnodes = round(dgnodes/snap)*snap;

dgnodes = reshape(permute(dgnodes(perm,1:dim,:),[1,3,2]),[npf*nfe*ne,dim]);
[~,b] = ismember(dgnodes,edg,'rows');
elcon = reshape(b,[npf nfe ne]);

dgnodes = reshape(dgnodes,[npf nfe ne dim]);
%dgnodes = permute(dgnodes,[1 2 4 3]);
ihdg = find(isEDGface~=1);
ihdg = sort(ihdg);
if (length(iedg)+length(ihdg)) ~= nf
    error('Number of HDG and EDG does not add the total number of faces.');
end
nbf = length(ihdg);
ndof = size(edg,1);
pdg = zeros(npf,dim,nbf);
for j = 1:nbf
    i = ihdg(j);
    fn = f(i,1:ncf);       % corners of face i 
    fe = f(i,ncf+1:ncf+2); % neighboring elements of face i    
    pf = philocfc*p(fn,:); % nodal points on face i
    pf = round(pf/snap)*snap;           
    pdg(:,:,j) = pf;               
    
    if isEDGface(i)
        error('EDG face was detected when looping over non-EDG faces.');
    end
    
    % hdg face
    tof = ndof + (1:npf);   % dof numbering on face i
    ndof = ndof + npf;      % number of degrees of freedom     
    
    if1 = find(t2f(fe(1),:)==i);           % location of face i on element fe(1)
    dg1 = reshape(dgnodes(:,if1,fe(1),:),[npf dim]); 
    in = xiny(dg1,pf);              % match pf to dg1      
    elcon(:,if1,fe(1)) = tof(in);   % assign dof numbering of face i to elcon from element fe(1)     
    
    if fe(2)>0
        if2 = find(t2f(fe(2),:)==i);           % location of face i on element fe(2)
        dg2 = reshape(dgnodes(:,if2,fe(2),:),[npf dim]); 
        in = xiny(dg2,pf);              % match pf to dg2    
        elcon(:,if2,fe(2)) = tof(in);   % assign dof numbering of face i to elcon from element fe(2) 
    end
end
pdg = permute(pdg,[1 3 2]);
edg = [edg; reshape(pdg,[npf*nbf dim])];

elcon = reshape(elcon,[npf nfe ne]);


% function [elcon,ndof,edg] = hedgelcon(mesh)
% 
% dim      = nd;
% porder   = porder;
% elemtype = elemtype;
% perm = perm(:);
% 
% ne       = size(t,1);
% nfe      = size(t2f,2); % number of faces per element
% nf       = size(f,1); % total number of faces
% ncf = dim;               % number of corners of a face
% if dim == 3 && elemtype  % hex element
%     ncf = dim+1;
% end
% 
% [philocfc,philocvl] = localbasis(porder,dim,elemtype);
% [npv,nnv] = size(philocvl);
% npf       = size(philocfc,1);
% 
% % get dgnodes
% dgnodes = zeros(npv,dim,ne);
% for d=1:dim
%   for n=1:nnv
%     dp=philocvl(:,n)*p(t(:,n),d)';
%     dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
%   end
% end
% snap = 1e-8;
% dgnodes = round(dgnodes/snap)*snap;
% 
% hdg = zeros(npf*nf,dim);
% for i = 1:nf
%     fn = f(i,1:ncf);       % corners of face i     
%     pf = philocfc*p(fn,:); % nodal points on face i    
%     ind = ((i-1)*npf+1):(i*npf);
%     hdg(ind,:) = pf;           % hdg nodes on faces
% end
% hdg = reshape(hdg,[npf nf dim]);
% hdg = round(hdg/snap)*snap;
% 
% iedg = find(f(:,end)==1);
% iedg = sort(iedg);
% edg = hdg(:,iedg,:);
% edg = reshape(edg,npf*length(iedg),dim);
% edg = flipdim(edg,1);
% [~,I] = unique(edg,'rows'); 
% edg = edg(sort(I),:);
% edg = flipdim(edg,1);
% 
% dgnodes = reshape(permute(dgnodes(perm,1:dim,:),[1,3,2]),[npf*nfe*ne,dim]);
% [~,b] = ismember(dgnodes,edg,'rows');
% elcon = reshape(b,[npf nfe ne]);
% 
% dgnodes = reshape(dgnodes,[npf nfe ne dim]);
% %dgnodes = permute(dgnodes,[1 2 4 3]);
% ihdg = find(f(:,end)==0);
% if (length(iedg)+length(ihdg)) ~= nf
%     error('something wrong in the last column of f');
% end
% ihdg = sort(ihdg);
% nbf = length(ihdg);
% ndof = size(edg,1);
% pdg = zeros(npf,dim,nbf);
% for j = 1:nbf
%     i = ihdg(j);
%     fn = f(i,1:ncf);       % corners of face i 
%     fe = f(i,ncf+1:ncf+2); % neighboring elements of face i    
%     pf = philocfc*p(fn,:); % nodal points on face i
%     pf = round(pf/snap)*snap;           
%     pdg(:,:,j) = pf;
%     
%     % hdg face
%     tof = ndof + (1:npf);   % dof numbering on face i
%     ndof = ndof + npf;      % number of degrees of freedom     
%     
%     if1 = find(t2f(fe(1),:)==i);           % location of face i on element fe(1)
%     dg1 = reshape(dgnodes(:,if1,fe(1),:),[npf dim]); 
%     in = xiny(dg1,pf);              % match pf to dg1      
%     elcon(:,if1,fe(1)) = tof(in);   % assign dof numbering of face i to elcon from element fe(1)             
%     
%     if fe(2)>0
%         if2 = find(t2f(fe(2),:)==i);           % location of face i on element fe(2)
%         dg2 = reshape(dgnodes(:,if2,fe(2),:),[npf dim]); 
%         in = xiny(dg2,pf);              % match pf to dg1      
%         elcon(:,if2,fe(2)) = tof(in);   % assign dof numbering of face i to elcon from element fe(1)                     
%     end
% end
% pdg = permute(pdg,[1 3 2]);
% edg = [edg; reshape(pdg,[npf*nbf dim])];
% 
% elcon = reshape(elcon,[npf*nfe ne]);


% ne       = size(t,1);
% nfv      = size(t2f,2); % number of faces per element
% [nf,nf2] = size(f); % total number of faces
% ncf = dim;               % number of corners of a face
% 
% [philocfc,philocvl] = localbasis(porder,dim,elemtype);
% [npv,nnv] = size(philocvl);
% npf       = size(philocfc,1);
% 
% % get dgnodes
% % dgnodes = zeros(npv,dim,ne);
% % for d=1:dim
% %   for n=1:nnv
% %     dp=philocvl(:,n)*p(t(:,n),d)';
% %     dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
% %   end
% % end
% 
% elcon = zeros(npf,nfv,ne);
% facecon = zeros(npf,nf);
% edg   = [];
% ndof  = 0;
% for i = 1:nf
%     nvf = f(i,1);
%     fn = f(i,2:nvf+1);      % corners of face i 
%     e1 = f(i,nvf+2);
%     e2 = f(i,nvf+3);
%     ei = elementtype(e1);
%     ii = find(elem==ei);
%     npe = npes(ii);             
%     nfe = nfes(ii); 
%     nve = nves(ii);     
%     
%     if1 = find(t2f(e1,:)==i);
%     npf = npfs(ii,if1);        
%     pf = p(fn,:); % nodal points on face i
%     facenodes = philocfc{ii}{if1}*pf; % nodal points on face i
%         
%     if (isEDGface(i)==0)                  % hdg face         
%         tof = ndof + (1:npf);   % dof numbering on face i
%         ndof = ndof + npf;      % number of degrees of freedom 
%     elseif (isEDGface(i)==1)              % edg face 
%         if isempty(edg)
%             edg = facenodes;           % edg nodes on faces
%             tof = ndof + (1:npf);            
%             ieg = tof;          % dof numbering of edg nodes
%             ndof = ndof + npf; 
%         else                        
%             in = xiny(facenodes,edg);     % find which rows of pf are in edg 
%             jn = (in==0);          % new edg nodes have in==0
%             kn = (in>0);           % old edg nodes have in>0  
%             neg = sum(jn);         % number of new edg nodes                        
%             
%             edg = [edg; facenodes(jn,:)]; % update edg with new edg nodes            
%            
%             tof(jn) = ndof+(1:neg); % tof for new edge nodes            
%             tof(kn) = ieg(in(kn));  % tof for old edge nodes                                                      
%             
%             ieg = [ieg ndof+(1:neg)]; % update ieg with dof numbering of new edge nodes 
%             ndof = ndof + neg;        % update ndof     
%         end
%     end
% 
%     pv = p(t(e1,1:nve),:);
%     dgnodes = philocal{ii}*pv;        
%     dg1 = dgnodes(perm{ii}{if1},:);  % dg nodes on face i from element fe(1)
%     
%     in = xiny(dg1,pf);              % match pf to dg1      
%     elcon(1:npf,if1,fe(1)) = tof(in);   % assign dof numbering of face i to elcon from element fe(1) 
%     facecon(1:npf,i) = tof;    
%     if fe(2)>0
%         if2 = (t2f(fe(2),:)==i);       % location of face i on element fe(2)
%         dg2 = squeeze(dgnodes(perm{ii}(:,if2),:,fe(2))); % dg nodes on face i from element fe(2)
%         in = xiny(dg2,pf);            % match pf to dg2  
%         if all(in==0)   % case for faces on periodic boundaries (Hemant Chaurasia)
%             in1 = xiny(dg1,pf);
%             in = in1(size(dg1,1):-1:1);  % reverse of in1
%         end
%         elcon(:,if2,fe(2)) = tof(in); % assign dof numbering of face i to elcon from element fe(2) 
%     end                
% end
% %elcon = reshape(elcon,[npf*nfv ne]);
% 
% 
% function [philocfc,philocvl] = localbasis(porder,dim,elemtype) 
% 
% [plocvl,tlocal,plocfc] = mkmasternodes(porder,dim,elemtype,0);
% plocfc = plocfc{1};
% 
% if dim==2 && elemtype==0      % tri
%     xi  = plocfc(:,1);
%     philocfc(:,1) = 1 - xi;
%     philocfc(:,2) = xi;
%     xi  = plocvl(:,1);
%     eta = plocvl(:,2);    
%     philocvl(:,1) = 1 - xi - eta;
%     philocvl(:,2) = xi;
%     philocvl(:,3) = eta;
% elseif dim==2 && elemtype==1  % quad
%     xi  = plocfc(:,1);
%     philocfc(:,1) = 1 - xi;
%     philocfc(:,2) = xi;
%     xi  = plocvl(:,1);
%     eta = plocvl(:,2);    
%     philocvl(:,1) = (1-xi).*(1-eta);
%     philocvl(:,2) = xi.*(1-eta);
%     philocvl(:,3) = xi.*eta;
%     philocvl(:,4) = (1-xi).*eta;
% elseif dim==3 && elemtype==0  % tet
%     xi  = plocfc(:,1);
%     eta = plocfc(:,2);    
%     philocfc(:,1) = 1 - xi - eta;
%     philocfc(:,2) = xi;
%     philocfc(:,3) = eta;
%     xi   = plocvl(:,1);
%     eta  = plocvl(:,2);
%     zeta = plocvl(:,3);
%     philocvl(:,1) = 1 - xi - eta - zeta;
%     philocvl(:,2) = xi;
%     philocvl(:,3) = eta;
%     philocvl(:,4) = zeta;
% elseif dim==3 && elemtype==1   % hex
%     xi  = plocfc(:,1);
%     eta = plocfc(:,2);
%     philocfc(:,1) = (1-xi).*(1-eta);
%     philocfc(:,2) = xi.*(1-eta);
%     philocfc(:,3) = xi.*eta;
%     philocfc(:,4) = (1-xi).*eta;
%     xi   = plocvl(:,1);
%     eta  = plocvl(:,2);
%     zeta = plocvl(:,3);
%     philocvl(:,1) = (1-xi).*(1-eta).*(1-zeta);
%     philocvl(:,2) = xi.*(1-eta).*(1-zeta);
%     philocvl(:,3) = xi.*eta.*(1-zeta);
%     philocvl(:,4) = (1-xi).*eta.*(1-zeta);    
%     philocvl(:,5) = (1-xi).*(1-eta).*(zeta);
%     philocvl(:,6) = xi.*(1-eta).*(zeta);
%     philocvl(:,7) = xi.*eta.*(zeta);
%     philocvl(:,8) = (1-xi).*eta.*(zeta);        
% end
% 
% 
% function in = xiny(x,y)
% % Determine if each row of x is a member of y
% % If row j of x is a member of y and x(j,:) = y(k,:) then in(j) = k
% % Else in(j) = 0
% 
% [m,dim] = size(x);
% 
% in = zeros(m,1);
% % if dim==2
% %     for j=1:m
% %         d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2;
% %         [md,id] = min(d2);
% %         if md<1e-12, in(j)=id; end
% %     end
% % else
% %     for j=1:m
% %         d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2 + (y(:,3)-x(j,3)).^2;
% %         [md,id] = min(d2);
% %         if md<1e-12, in(j)=id; end
% %     end
% % end
% 
% if dim == 2
%     y1Rep = repmat(y(:,1),[1,m]);
%     y2Rep = repmat(y(:,2),[1,m]);
%     x1Rep = repmat(x(:,1)',[size(y,1),1]);
%     x2Rep = repmat(x(:,2)',[size(y,1),1]);
% 
%     dis = (y1Rep-x1Rep).^2 + (y2Rep-x2Rep).^2;
%     [md,id] = min(dis,[],1);
%     aux = md < (1e-10)^2;
%     in(aux) = id(aux);
% else
%     y1Rep = repmat(y(:,1),[1,m]);
%     y2Rep = repmat(y(:,2),[1,m]);
%     y3Rep = repmat(y(:,3),[1,m]);
%     x1Rep = repmat(x(:,1)',[size(y,1),1]);
%     x2Rep = repmat(x(:,2)',[size(y,1),1]);
%     x3Rep = repmat(x(:,3)',[size(y,1),1]);
% 
%     dis = (y1Rep-x1Rep).^2 + (y2Rep-x2Rep).^2 + (y3Rep-x3Rep).^2;
%     [md,id] = min(dis,[],1);
%     aux = md < (1e-10)^2;
%     in(aux) = id(aux);
% end
