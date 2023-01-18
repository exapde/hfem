function hh=meshplot(mesh)
%MESHPLOT  Plot Mesh Structure (with straight edges)
%    MESHPLOT(MESH,[OPTS])
%
%    MESH:       Mesh structure
%    OPTS:       (logical)
%      OPTS(1):  Plot dgnodes (default=0)
%      OPTS(2):  Plot triangle numbers (default=0)
%      OPTS(3):  Plot node numbers (default=0)


if nargin<2 | isempty(opts), opts=[0]; end
if length(opts)<2, opts=[opts,0]; end
if length(opts)<3, opts=[opts,0]; end

p=mesh.p;
t=mesh.t;

hh=[];
pars={'facecolor',[.8,1,.8],'edgecolor','k','Linew',1};
clf,hh=[hh;patch('faces',t,'vertices',p,pars{:})];
view(2),axis equal


if nargout<1, clear hh; end
