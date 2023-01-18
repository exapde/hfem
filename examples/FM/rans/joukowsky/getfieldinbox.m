function [udg,pdg,pm,xm,inde] = getfieldinbox(mesh,UDG,box,x,dmax,nref)

nd = mesh.nd;
[ne,nv] = size(mesh.t);
p = reshape(mesh.p(mesh.t',:),[nv ne nd]);
xm = reshape(mean(p,1),[ne nd]);

if nd==3
xmin = box(1,1); xmax = box(1,2);
ymin = box(2,1); ymax = box(2,2);
zmin = box(3,1); zmax = box(3,2);
ind1 = find( (xmin < xm(:,1)) & (xm(:,1)<xmax) & ...
             (ymin < xm(:,2)) & (xm(:,2)<ymax) & ...
             (zmin < xm(:,3)) & (xm(:,3)<zmax));         
elseif nd==2
xmin = box(1,1); xmax = box(1,2);
ymin = box(1,3); ymax = box(1,4);
%[xmin xmax ymin ymax]
ind1 = find( (xmin < xm(:,1)) & (xm(:,1)<xmax) & (ymin < xm(:,2)) & (xm(:,2)<ymax));    
end

xm = xm(ind1,:);         
if nd == 3
    distance = sqrt((x(1)-xm(:,1)).^2 + (x(2)-xm(:,2)).^2 + (x(3)-xm(:,3)).^2);    
else
    distance = sqrt((x(1)-xm(:,1)).^2 + (x(2)-xm(:,2)).^2);    
end
inde = distance<dmax;
inde = ind1(inde);

% get the velocity field and coordinate on those elements
udg = permute(UDG(:,:,inde),[1 3 2]);
pdg = permute(mesh.dgnodes(:,1:nd,inde),[1 3 2]);

% subdivision  
[pdg,udg] = scalarrefine(mesh,pdg,udg,nref);              
udg = round(udg/1e-10)*1e-10;
pdg = round(pdg/1e-10)*1e-10;
pm = reshape(mean(pdg,1),[size(pdg,2) nd]); 

function [pref,uref] = scalarrefine(mesh,p,u,nref)

[npl, nt, nd] = size(p);
porder=mesh.porder;
plocal=mesh.plocal;
tlocal=mesh.tlocal;

if isempty(nref), nref=ceil(log2(max(porder,1))); end
if mesh.elemtype==0  
    A0=koornwinder(plocal(:,1:nd),porder);
    [plocal,tlocal]=uniref3d(plocal,tlocal,nref);    
    A=koornwinder(plocal(:,1:nd),porder)/A0;
else
    A0=tensorproduct(plocal(:,1:nd),porder);
    m = porder*(nref+1)+1;     
    if nd==2
        [plocal,tlocal]=squaremesh(m,m,1);
    else
        [plocal,tlocal]=cubemesh(m,m,m,1);
    end
    A=tensorproduct(plocal(:,1:nd),porder)/A0;  
end

npln=size(plocal,1);
t = kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);
ne = size(t,1);
np = npln*nt;

nc = size(u,3);
uref = reshape(A*reshape(u,npl,nt*nc),[np nc]);
uref = reshape(uref(t',:), [size(t,2) ne nc]);

pref = reshape(A*reshape(p,npl,nt*nd),[np,nd]);
pref = reshape(pref(t',:),[size(t,2) ne nd]);

%function plotTetrahedron(x, p)






