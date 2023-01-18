
function [epart, npart, ent2ent, ent2entStart] = meshpart(mesh,np,p0,t0,plot)

if nargin < 5; plot = 0; end

% Compute connectivities "t" and number of entities in mesh
if strcmp(mesh.hybrid,'hdg')
    t = mesh.t2f;
    numTotalEntities = max(mesh.t2f(:));
    numRealEntities = length(unique(mesh.t2f(:)));
elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg') || strcmp(mesh.hybrid,'hedg')
    t = mesh.elcon';
    numTotalEntities = max(mesh.elcon(:));
    numRealEntities = length(unique(mesh.elcon(:)));
end
if numRealEntities ~= numTotalEntities
    error('There are ghost (missing) entities in the mesh.');
end

% Compute entity-to-entity connectivities of entire mesh:
% [ent2ent, ent2entStart] = computeEnt2ent(mesh);
% ent2entStart = ent2entStart + 1;        % Make indices start at 1 (instead of 0, as in C)
ent2ent = [];
ent2entStart = [];

% current directory
cdir = pwd;

% move to directory that contains metis programs
if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end
ps=strcat(pwd,sslash);
is=find(ps==sslash);
up=0;
while ~(strcmp(ps(is(end-up-1)+1:is(end-up)-1),'hdgv1.0') || strcmp(ps(is(end-up-1)+1:is(end-up)-1),'HDGv1.0'))
    up = up+1;
end
cd(strcat(ps(1:is(end-up)),'metis'));

% number of elements
nt = size(t,1);

% generate a temporary file to be used in metis
dlmwrite('temp.txt', nt, 'delimiter', ' ','precision',10);
dlmwrite('temp.txt', t, '-append', 'delimiter', ' ','precision',10);

% call mpmetis
str = ['!./mpmetis temp.txt ' num2str(np)];
eval(str);

% get mesh partitioning data
str = ['temp.txt.epart.' num2str(np)];
epart = textread(str,'%d');

% get node partitioning data
str = ['temp.txt.npart.' num2str(np)];
npart = textread(str,'%d');

% remove files
%delete('temp.txt');
str = ['temp.txt.epart.' num2str(np)];
delete(str);
str = ['temp.txt.npart.' num2str(np)];
delete(str);

% move back to current directory
cd(cdir);

nelem = zeros(1,np);
nent = nelem;
for i = 0:np-1
    nelem(i+1) = length(find(epart==i));
    nent(i+1) = length(find(npart==i));
end
% nelem
% nent
% [min(nelem) max(nelem)];
% [min(nent) max(nent)];
% nelem;
% nent;

% plot the partitioned mesh 
if plot
    p = p0;
    t = t0;
    bcol = [1 1 0;... % yellow 
            1 0 1;... % magneta
            0 1 1;... % cyan
            1 0 0;... % red
            0 1 0;... % green
            0 0 1;... % blue
            1,0.4,0.6;...
            0.4,0.6,1;...
           ];
    figure(1); clf;
    hold on;        
    for i=0:np-1
        ind = find(epart==i);
        ti = t(ind,:);
        simpplot(p,ti,[],bcol(i+1,:));                       
    end
    for it=1:size(t,1)
        pmid=mean(p(t(it,:),:),1);
        txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
        text(pmid(1),pmid(2),num2str(it),txtpars{:});
    end

    hold off;
    axis equal;      
    axis tight;
    axis on;  
    
    if strcmp(mesh.hybrid,'hdg')
        figure(2); clf;
        simpplot(p,t);
        hold on;
        for i=1:size(mesh.f,1)
            j = npart(i);        
            txtpars={'fontname','times','fontsize',20,'fontweight','bold', ...
                'horizontala','center','col','k', 'BackgroundColor',bcol(j+1,:)};            
            fi = mesh.f(i,1:2);
            pm = (mesh.p(fi(1),:)+mesh.p(fi(2),:))/2;
            text(pm(1),pm(2),num2str(i),txtpars{:});
        end    
    end
    
    
end

end


% function [ent2ent, ent2entStart] = computeEnt2ent(mesh)
% 
% if strcmp(mesh.hybrid,'hdg')
%     
%     ne    = size(mesh.t2f,1);
%     nfe   = size(mesh.t2f,2); 
%     il = zeros(nfe,nfe,ne);
%     jl = zeros(nfe,nfe,ne);
%     for i=1:ne    
%         con = mesh.t2f(i,:);
%         il(:,:,i) = repmat(con' , 1, nfe);
%         jl(:,:,i) = repmat(con , nfe, 1);        
%     end
%     il(:,:,1);
%     jl(:,:,1);
%     sl = ones(nfe,nfe,ne);
% 
% elseif strcmp(mesh.hybrid,'edg') || strcmp(mesh.hybrid,'iedg')
%     
%     [nn,ne] = size(mesh.elcon);
%     %ne = mesh.ne;
%     il = zeros(nn,nn,ne);
%     jl = zeros(nn,nn,ne);
%     for i=1:ne
%         con = mesh.elcon(:,i);
%         il(:,:,i) = repmat(con ,1,nn);
%         jl(:,:,i) = repmat(con',nn,1);
%     end    
%     sl = ones(nn,nn,ne);
%     
% end
% [rp, cj] = sparse_to_csr(sparse(il(:),jl(:),sl(:)));
% 
% ct = cj;
% for i = 1:length(rp)-1    
%     k = rp(i)+1:rp(i+1);    
%     ct(k) = [i; sort(setdiff(ct(k),i))]; % the first index is the block itself   
% end
% 
% if length(ct) ~= length(cj)
%     error('Matrix must have nonzero diagonal entries!');
% else
%     cj = ct;
% end
% 
% ent2entStart = rp;
% ent2ent = cj;
% 
% end
