function [elcon,t2f,bf,f,t,isEDGface] = periodicmatlab(pp,elcon,t2f,bf,f,t,isEDGface,porder,elemtype,perm,bndexpr,hybrid)
% PERIODIC modify elcon, f, t2f, nf, fcurved, bf to account for periodic boundary conditions
%
%    mesh = periodic(mesh,bndexpr)
%
%
%    bndexpr:  describes which boundaries are periodic:
%              bndexpr = {bnd1,expr1,bnd2,expr2;
%                   ....,.....,....,.....}
%          will make bnd1 and bnd2 periodic, where elements
%          are matched based on the expressions in expr1
%          (expressions depending on nodes p)
%
%          Ex: bndexpr = {1,'p(:,1),3,'p(:,1)'} will make bnds 1,3
%
% NOTE: Only valid for HDG faces on the periodic boundaries. At least the following
% three lines are based on that hypothesis:
% in  = xiny(q2,q1);
% tm1 = elcon(:,lf1,ebnd1(j));
% elcon(:,lf2,ebnd2(j)) = tm1(in);

tol = 1e-9; % Tolerance to match periodic boundaries.
nd  = size(pp,2);
nfv = size(t2f,2);    % number of faces per element
ne = size(t,1); % number of elements
npf = size(elcon,1)/nfv; % number of points per face

[~,philocvl,nvf] = localbasis(porder,nd,elemtype);
[npv,nnv] = size(philocvl);

% get dgnodes
dgnodes = zeros(npv,nd,ne);
for d=1:nd
  for n=1:nnv
    dp=philocvl(:,n)*pp(t(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end

nperiodic = size(bndexpr,1);
pbnd1 = cell(nperiodic,1);
fbnd1 = cell(nperiodic,1);
ebnd1 = cell(nperiodic,1);
vbnd1 = cell(nperiodic,1);
pbnd2 = cell(nperiodic,1);
fbnd2 = cell(nperiodic,1);
ebnd2 = cell(nperiodic,1);
vbnd2 = cell(nperiodic,1);
for i = 1:nperiodic

    if size(f,2)==nvf+3
        % get faces on the 1st boundary
        [pbnd1{i},fbnd1{i},ebnd1{i},vbnd1{i}] = getface(f(:,2:end),pp,bndexpr{i,1},bndexpr{i,2});
        % get faces on the 2nd boundary
        [pbnd2{i},fbnd2{i},ebnd2{i},vbnd2{i}] = getface(f(:,2:end),pp,bndexpr{i,3},bndexpr{i,4});
    elseif size(f,2)==nvf+2
        [pbnd1{i},fbnd1{i},ebnd1{i},vbnd1{i}] = getface(f,pp,bndexpr{i,1},bndexpr{i,2});
        [pbnd2{i},fbnd2{i},ebnd2{i},vbnd2{i}] = getface(f,pp,bndexpr{i,3},bndexpr{i,4});
    end

    % sort pbnd2 to match it to pbnd1
    ind = xiny(pbnd1{i}, pbnd2{i});
    pbnd2{i} = pbnd2{i}(ind,:);
    fbnd2{i} = fbnd2{i}(ind,:);
    ebnd2{i} = ebnd2{i}(ind,:);
    vbnd2{i} = vbnd2{i}(ind,:);

    if max(abs(pbnd2{i}(:)-pbnd1{i}(:))) > tol
        error('two periodic bundaries are not matched');
    end
end

t2fp = t2f;
tp = t;
elcon = reshape(elcon,[npf nfv ne]);
for i = 1:size(bndexpr,1)
    f(fbnd1{i},end) = f(fbnd2{i},end-1);
    if any(isEDGface(fbnd1{i}) ~= isEDGface(fbnd2{i}))
        error('Periodic faces must be either both HDG or both EDG.')
    end
    [isEDGface(fbnd1{i}), whichHDG] = min([isEDGface(fbnd1{i}) , isEDGface(fbnd2{i})], [], 2);
    isEDGface(fbnd2{i}) = isEDGface(fbnd1{i});

    nbf = length(fbnd1{i});
    for j=1:nbf
        lf1 = (t2fp(ebnd1{i}(j),:)==fbnd1{i}(j));  % location of face fbnd1(j) on element ebnd1(j)
        lf2 = (t2fp(ebnd2{i}(j),:)==fbnd2{i}(j));  % location of face fbnd2(j) on element ebnd2(j)
        t2f(ebnd2{i}(j),lf2) = t2f(ebnd1{i}(j),lf1);
        bf(lf1,ebnd1{i}(j)) = ebnd2{i}(j);
        bf(lf2,ebnd2{i}(j)) = ebnd1{i}(j);

        p   = squeeze(dgnodes(perm(:,lf1),:,ebnd1{i}(j)));  % dg nodes on face fbnd1(j) from element ebnd1(j)
        q1  = eval(bndexpr{i,2});        % evaluate periodic boundary expression
        p   = squeeze(dgnodes(perm(:,lf2),:,ebnd2{i}(j)));  % dg nodes on face fbnd2(j) from element ebnd2(j)
        q2  = eval(bndexpr{i,4});        % evaluate periodic boundary expression
        in  = xiny(q2,q1);               % match q2 to q1
%         in_  = xiny(q1,q2);               % match q1 to q2

%         if strcmp(hybrid,'edg')
        if isEDGface(fbnd2{i}(j))       % Both periodic faces are EDG
            tm1 = elcon(:,lf1,ebnd1{i}(j));     % dofs of face fbnd1(j)
            tm2 = elcon(:,:,ebnd2{i}(j));
            for k = 1:length(tm1)
                entk = tm2(k,lf2);
                tm2(tm2==entk) = tm1(in(k));
            end
            elcon(:,:,ebnd2{i}(j)) = tm2;
%         elseif strcmp(hybrid,'hdg') || strcmp(hybrid,'iedg')
        elseif ~isEDGface(fbnd2{i}(j))      % At least one peridoic face is HDG
            if whichHDG(j) == 1
                tm1 = elcon(:,lf1,ebnd1{i}(j));     % dofs of face fbnd1(j)
                elcon(:,lf2,ebnd2{i}(j)) = tm1(in);
            elseif whichHDG(j) == 2
                tm2 = elcon(:,lf2,ebnd2{i}(j));     % dofs of face fbnd1(j)
                elcon(:,lf1,ebnd1{i}(j)) = tm2(in_);
            else
                error('Something wrong.')
            end
        else
            error('hybrid flag not recognized.');
        end

        [~,vp1] = ismember(vbnd1{i}(j,:),tp(ebnd1{i}(j),:)); % location of vbnd1(j,:) on element ebnd1(j)
        [~,vp2] = ismember(vbnd2{i}(j,:),tp(ebnd2{i}(j),:)); % location of vbnd2(j,:) on element ebnd2(j)

        p = pp(tp(ebnd1{i}(j),vp1),:);
        q1  = eval(bndexpr{i,2});
        p = pp(tp(ebnd2{i}(j),vp2),:);
        q2  = eval(bndexpr{i,4});
        in = xiny(q2, q1);

        % Update tUpdated:
        tm1 = t(ebnd1{i}(j),vp1);
        t(ebnd2{i}(j),vp2) = tm1(in);
    end
end
f(cell2mat(fbnd2),:) = [];
isEDGface(cell2mat(fbnd2)) = [];

% Fix elcon (reorder nodes to eliminate holes)
mapping = zeros(max(elcon(:)),1);
uniqueNodes = unique(elcon(:));
mapping(uniqueNodes) = 1:length(uniqueNodes);
elcon = mapping(elcon);
elcon = reshape(elcon,[npf*nfv ne]);

mapping = zeros(max(t2f(:)),1);
uniqueNodes = unique(t2f(:));
mapping(uniqueNodes) = 1:length(uniqueNodes);
t2f = mapping(t2f);

mapping = zeros(max(t(:)),1);
uniqueNodes = unique(t(:));
mapping(uniqueNodes) = 1:length(uniqueNodes);
t = mapping(t);

function [pbnd,fbnd,ebnd,vbnd] = getface(f,p0,bn,expr)

ind = find(f(:,end)==0, 1);
if isempty(ind)==0
    bn = bn-1;
end

% get faces on the boundary bn
fbnd = find(f(:,end)==-bn);
ft   = f(fbnd,:);
% number of faces on the boundary bn
n    = length(fbnd);

p = p0(ft(1,1:end-2),:);
q = eval(expr);
[k, m] = size(q); % m = nd, k = number of vertices of a face

% get elements which contain faces on the boundary bn
ebnd = ft(:,end-1);

% get coordinate values of faces on the boundary bn
pbnd = zeros(k,m,n);
vbnd = zeros(n,k);
for i = 1:n % for each face i on the boundary bn
    p = p0(ft(i,1:end-2),:); % vertices of the face i
    q = eval(expr); % evaluate periodic expression of this face
    pbnd(:,:,i) = sortrows(snap(q));
    vbnd(i,:) = ft(i,1:end-2);
end
pbnd = reshape(permute(pbnd,[3,1,2]),n,[]);

function x=snap(x)

tol=sqrt(eps);
x=tol*round(x/tol);


function in = xiny(x,y)
% Determine if each row of x is a member of y
% If row j of x is a member of y and x(j,:) = y(k,:) then in(j) = k
% Else in(j) = 0

[m,dim] = size(x);
in = zeros(m,1);
if dim==1
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2;
        [md,id] = min(d2);
        if md<1e-12, in(j)=id; end
    end
elseif dim==2
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2;
        [md,id] = min(d2);
        if md<1e-12, in(j)=id; end
    end
elseif dim==3
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2 + (y(:,3)-x(j,3)).^2;
        [md,id] = min(d2);
        if md<1e-12, in(j)=id; end
    end
else
    n = size(y,1);
    for j=1:m
        d2 = sum((y - repmat(x(j,:),[n 1])).^2,2);
        [md,id] = min(d2);
        if md<1e-12, in(j)=id; end
    end
end

function [philocfc,philocvl,nvf] = localbasis(porder,dim,elemtype)

[plocvl,~,plocfc] = masternodes(porder,dim,elemtype,0);

if dim==2 && elemtype==0      % tri
    xi  = plocfc(:,1);
    philocfc(:,1) = 1 - xi;
    philocfc(:,2) = xi;
    xi  = plocvl(:,1);
    eta = plocvl(:,2);
    philocvl(:,1) = 1 - xi - eta;
    philocvl(:,2) = xi;
    philocvl(:,3) = eta;
    nvf = 2;
elseif dim==2 && elemtype==1  % quad
    xi  = plocfc(:,1);
    philocfc(:,1) = 1 - xi;
    philocfc(:,2) = xi;
    xi  = plocvl(:,1);
    eta = plocvl(:,2);
    philocvl(:,1) = (1-xi).*(1-eta);
    philocvl(:,2) = xi.*(1-eta);
    philocvl(:,3) = xi.*eta;
    philocvl(:,4) = (1-xi).*eta;
    nvf = 2;
elseif dim==3 && elemtype==0  % tet
    xi  = plocfc(:,1);
    eta = plocfc(:,2);
    philocfc(:,1) = 1 - xi - eta;
    philocfc(:,2) = xi;
    philocfc(:,3) = eta;
    xi   = plocvl(:,1);
    eta  = plocvl(:,2);
    zeta = plocvl(:,3);
    philocvl(:,1) = 1 - xi - eta - zeta;
    philocvl(:,2) = xi;
    philocvl(:,3) = eta;
    philocvl(:,4) = zeta;
    nvf = 3;
elseif dim==3 && elemtype==1   % hex
    xi  = plocfc(:,1);
    eta = plocfc(:,2);
    philocfc(:,1) = (1-xi).*(1-eta);
    philocfc(:,2) = xi.*(1-eta);
    philocfc(:,3) = xi.*eta;
    philocfc(:,4) = (1-xi).*eta;
    xi   = plocvl(:,1);
    eta  = plocvl(:,2);
    zeta = plocvl(:,3);
    philocvl(:,1) = (1-xi).*(1-eta).*(1-zeta);
    philocvl(:,2) = xi.*(1-eta).*(1-zeta);
    philocvl(:,3) = xi.*eta.*(1-zeta);
    philocvl(:,4) = (1-xi).*eta.*(1-zeta);
    philocvl(:,5) = (1-xi).*(1-eta).*(zeta);
    philocvl(:,6) = xi.*(1-eta).*(zeta);
    philocvl(:,7) = xi.*eta.*(zeta);
    philocvl(:,8) = (1-xi).*eta.*(zeta);
    nvf = 4;
end
