function mesh = mkmeshmixed(p,t,porder,bndexpr,elemtype,nodetype,fb,fd,varargin)
%MKMESH Creates mesh data structure 
%   MESH=MKMESH(p,t,porder,bndexpr,elemtype,nodetype,fd,fdparams)
%
%      MESH:      Mesh structure
%      P:         Node positions
%      T:         Node connectivities
%      PORDER:    Polynomial Order of Approximation (default=1)
%      BNDEXPR:   Cell Array of boundary expressions. The 
%                 number of elements in BNDEXPR determines 
%                 the number of different boundaries
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%      NODETYPE:  Flag determining node distribution 
%                 Flag = 0 uniform distribution (default)
%                 Flag = 1 nonuniform distribution
%      FB:        Index for the curved boundary
%      FD:        Distance Function d(x,y,z)
%      FDPARAMS:  Additional parameters passed to FD
%
%   See also: MKT2F, SETBNDNBRS, MASTERNODES, CREATENODES
%

if nargin<4
    error('Require at least four input aguments');
end
if nargin<5, elemtype=0; end
if nargin<6, nodetype=0; end
if nargin<7, fb=[]; end
if nargin<8, fd=[]; end

%[p,t]=fixmesh(p,t);
dim = size(p,2);

mesh.porder = porder; % polynomial degree
mesh.p = p;  % node positions
mesh.t = t;  % element-to-node connectivities 
mesh.elemtype = elemtype;

% compute face-to-node connectivities and element-to-face connectivities 
%[mesh.f,mesh.t2f] = mkt2f(mesh.t,elemtype);
[mesh.t2f,mesh.t2t,mesh.f,mesh.ne,mesh.nf,nfe,nve,nvf] = mkt2fmixed(mesh.t,mesh.elemtype,dim);
mesh.nfe = nfe; 
mesh.nve = nve;
mesh.nvf = nvf;

% make corrections for faces on the boundaries
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

% compute flags for curved faces and elements
if isempty(fb)
    mesh.fcurved = zeros(size(mesh.f,1),1); 
    mesh.tcurved = zeros(size(mesh.t,1),1);
else
    mesh.fcurved = (mesh.f(:,end)==-fb);
    ic = mesh.fcurved;
    mesh.tcurved = false(size(mesh.t,1),1);
    mesh.tcurved(mesh.f(ic,end-1)) = true;
end


% generate DG nodes
mesh.dgnodes = createnodes(mesh,fd,varargin);

mesh.np = size(mesh.p,1); % number of nodes 
mesh.nd = dim;            % problem dimension
mesh.nodetype = nodetype; % 0 -> uniform distribution or 1 -> nonuniform distribution 
mesh.bndexpr = bndexpr;
mesh.fb = fb;


function dgnodes=createnodes(mesh,fd,varargin)

if nargin < 2, fd=[]; end

nd = size(mesh.p,2);
nt = size(mesh.t,1);
elemtype = unique(mesh.elemtype);
nelem = length(elemtype);
npl = zeros(nelem,1);
for i = 1:nelem
    philocal{i} = localbasis(mesh.porder,nd,elemtype(i));
    npl(i) = size(philocal{i},1);
end
nplmax = max(npl);

% Allocate nodes
dgnodes=zeros(nplmax,nd,nt);
for i = 1:nelem
    if nelem>1
        ei = find(mesh.elemtype==elemtype(i));
    else
        ei = 1:1:nt;
    end
    npli = npl(i);
    npv = size(philocal{i},2);
    for dim=1:nd
      for node=1:npv
        dp=philocal{i}(:,node)*mesh.p(mesh.t(ei,node),dim)';
        dgnodes(1:npli,dim,ei)=dgnodes(1:npli,dim,ei)+permute(dp,[1,3,2]);
      end
    end
end

% Project nodes on the curved boundary
if ~isempty(fd) && mesh.porder>1
  tc=find(mesh.tcurved);
  for it=tc'
    p = dgnodes(:,:,it);
    deps=sqrt(eps)*max(max(p)-min(p));
    ed = find(mesh.f(abs(mesh.t2f(it,:)),4)<0);
    for id=ed'        
        e = find(philocal(:,id) < 1.e-6);        
        %d=feval(fd,p(e,:),varargin{:});        
        %dgradx=(feval(fd,[p(e,1)+deps,p(e,2)],varargin{:})-d)/deps;
        %dgrady=(feval(fd,[p(e,1),p(e,2)+deps],varargin{:})-d)/deps;
        d=feval(fd,p(e,:));
        dgradx=(feval(fd,[p(e,1)+deps,p(e,2)])-d)/deps;
        dgrady=(feval(fd,[p(e,1),p(e,2)+deps])-d)/deps;
        dgrad2=dgradx.^2+dgrady.^2;
        dgrad2(dgrad2==0)=1;
        p(e,:)=p(e,:)-[d.*dgradx./dgrad2,d.*dgrady./dgrad2];
    end
    dgnodes(:,:,it) = p;
  end
end

