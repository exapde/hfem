function V = getfieldatx(mesh,UDG,x,box)

nd = mesh.nd;
[ne,nv] = size(mesh.t);
p = reshape(mesh.p(mesh.t',:),[nv ne nd]);
xm = reshape(mean(p,1),[ne nd]);

if nargin > 3
    xmin = box(1,1); xmax = box(1,2);
    ymin = box(2,1); ymax = box(2,2);
    zmin = box(3,1); zmax = box(3,2);
    ind1 = find( (xmin < xm(:,1)) & (xm(:,1)<xmax) & ...
                 (ymin < xm(:,2)) & (xm(:,2)<ymax) & ...
                 (zmin < xm(:,3)) & (xm(:,3)<zmax));         
    ne = length(ind1);         
    xm = xm(ind1,:);         
end

dist = zeros(ne,1);
for i = 1:ne        
    s = sqrt((x(:,1)-xm(i,1)).^2 + (x(:,2)-xm(i,2)).^2 + (x(:,3)-xm(i,3)).^2);
    dist(i) = min(s);      
end

dmax = 0.5;
inde = find(dist<dmax);
if nargin > 3
    inde = ind1(inde);
end

% get the velocity field and coordinate on those elements
udg = permute(UDG(:,:,inde),[1 3 2]);
pdg = permute(mesh.dgnodes(:,1:nd,inde),[1 3 2]);

% subdivision
nref=1;
if mesh.porder>1         
    [pdg,udg] = scalarrefine(mesh,pdg,udg,nref);        
end                
udg = round(udg/1e-10)*1e-10;
pdg = round(pdg/1e-10)*1e-10;
pm = reshape(mean(pdg,1),[size(pdg,2) nd]); 

np = size(udg,1);
nc = size(udg,3);
P = reshape(pdg,[],nd);
U = reshape(udg,[],nc);
N = size(x,1);
V = zeros(N,nc);
for j = 2:1:N              
   d = sqrt((x(j,1)-P(:,1)).^2+(x(j,2)-P(:,2)).^2+(x(j,3)-P(:,3)).^2);
   if min(d)<1e-12
       V(j,:) = mean(U(d<1e-12,:),1);      
   else               
       d = (x(j,1)-pm(:,1)).^2+(x(j,2)-pm(:,2)).^2+(x(j,3)-pm(:,3)).^2;
       d = sqrt(d);
       [~,ks] = sort(d);          
       for k = 1:length(ks)               
           if k>1000
               V(j,:) = [0 0 0];
               v1 = reshape(pdg(1,ks(1),:),[1 3]);
               v2 = reshape(pdg(2,ks(1),:),[1 3]);
               v3 = reshape(pdg(3,ks(1),:),[1 3]);
               v4 = reshape(pdg(4,ks(1),:),[1 3]);
               v = [v1;v2;v3;v4]
               x(j,:)           
                figure(1);clf;
                plot3(v(:,1),v(:,2),v(:,3),'-');
                hold on;
                v = [v1;v3]; plot3(v(:,1),v(:,2),v(:,3),'-');
                v = [v1;v4]; plot3(v(:,1),v(:,2),v(:,3),'-');
                v = [v2;v4]; plot3(v(:,1),v(:,2),v(:,3),'-');
                plot3(x(j,1),x(j,2),x(j,3),'*');
                plot3(pm(ks(1),1),pm(ks(1),2),pm(ks(1),3),'o');           
               error('k>1000: something wrong');
               %in = PointInTetrahedron(v1, v2, v3, v4, x(j,:));
               break;
           end
           v1 = reshape(pdg(1,ks(k),:),[1 3]);
           v2 = reshape(pdg(2,ks(k),:),[1 3]);
           v3 = reshape(pdg(3,ks(k),:),[1 3]);
           v4 = reshape(pdg(4,ks(k),:),[1 3]);
           in = PointInTetrahedron(v1, v2, v3, v4, x(j,:));
           if in==1
               a = [v1 1; v2 1; v3 1; v4 1]\reshape(udg(:,ks(k),:),[np nc]);
               V(j,:)  = [x(j,:) 1]*a;
               break;
           end
       end
   end
    [j x(j,:) V(j,:)]
end

figure(1);clf;
plot3(x(:,1),x(:,2),x(:,3),'-b','LineWidth',2);
hold on;
plot3(x(:,1)+V(:,1),x(:,2)+V(:,2),x(:,3)+V(:,3),'-r','LineWidth',2);
for i = 1:N
    plot3([x(i,1) x(i,1)+V(i,1)],[x(i,2) x(i,2)+V(i,2)],[x(i,3) x(i,3)+V(i,3)],'--k','LineWidth',1);    
end
axis tight;

% figure(1);clf;
% plot3(x(:,1),x(:,2),x(:,3),'-');
% hold on;
% x = [v1;v3]; plot3(x(:,1),x(:,2),x(:,3),'-');
% x = [v1;v4]; plot3(x(:,1),x(:,2),x(:,3),'-');
% x = [v2;v4]; plot3(x(:,1),x(:,2),x(:,3),'-');
% plot3(p(j,1),p(j,2),p(j,3),'*');
% plot3(pm(ks(k),1),pm(ks(k),2),pm(ks(k),3),'o');


function sameside = IsSameSide(v1, v2, v3, v4, p)

normal = cross(v2 - v1, v3 - v1);
dotV4 = dot(normal, v4 - v1);
dotP = dot(normal, p - v1);

if abs(dotP)<=1e-14 
    sameside = 1;
elseif sign(dotP) == sign(dotV4)
    sameside = 1;
else
    sameside = 0;
end

function in = PointInTetrahedron(v1, v2, v3, v4, p)

if norm(p-v1)<=1e-12 || norm(p-v2)<=1e-12 || norm(p-v3)<=1e-12 || norm(p-v4)<=1e-12
    in = 1;
    return;
end

in = IsSameSide(v1, v2, v3, v4, p) & ...
     IsSameSide(v2, v3, v4, v1, p) & ...
     IsSameSide(v3, v4, v1, v2, p) & ...
     IsSameSide(v4, v1, v2, v3, p);               


function [pref,uref] = scalarrefine(mesh,p,u,nref)

[npl, nt, nd] = size(p);
porder=mesh.porder;
plocal=mesh.plocal;
tlocal=mesh.tlocal;

if isempty(nref), nref=ceil(log2(max(porder,1))); end
if size(tlocal,2)==4  
    A0=koornwinder(plocal(:,1:nd),porder);
    [plocal,tlocal]=uniref3d(plocal,tlocal,nref);    
    A=koornwinder(plocal(:,1:nd),porder)/A0;
else
    A0=tensorproduct(plocal(:,1:nd),porder);
    m = porder*(nref+1)+1;     
    [plocal,tlocal]=cubemesh(m,m,m,1);
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






