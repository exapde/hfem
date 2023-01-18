function mesh = mkmesh_circleinrect5(porder)

r = 1.0*1e-2;
R = 100;
xmin = -R;
xmax = R;
ymin = -15;
ylin   = -15;
ymax = R;
h = 13;

% n = 81;
% phi=linspace(0,2*pi,n);
% phi=phi(2:end);
% 
% pv1=[xmin xmax xmax xmin; ylin ylin ymax ymax]';
% pv2=r*[cos(phi); sin(phi)]';
% [p,t]=polymesh({pv1,pv2},[1,1],[1,0;0,1],[h,1.4],@href,h/5.5);

h = 10;
n = 81;
phi=linspace(0,2*pi,n);
phi=phi(2:end);

pv1=[xmin xmax xmax xmin; ylin ylin ymax ymax]';
pv2=r*[cos(phi); sin(phi)]';
[p,t]=polymesh({pv1,pv2},[1,1],[1,0;0,1],[h,1.35],@href,h/5);

% figure(1); clf; simpplot(p1,t1); axis on;
% figure(2); clf; simpplot(p2,t2); axis on;
% figure(3); clf; simpplot(p,t); axis on;
%[p,t] = fixmesh(p,t);

s1 = ['all(p(:,2)<' num2str(ymin) '+1e-7)'];
s2 = ['all(p(:,1)>' num2str(xmax) '-1e-7)'];
s3 = ['all(p(:,2)>' num2str(ymax) '-1e-7)'];
s4 = ['all(p(:,1)<' num2str(xmin) '+1e-7)'];
s5 = ['all(sum(p.^2,2)<' num2str(2*r) '^2)'];
bndexpr = {s1,s2,s3,s4,s5};   

fd=@(p) sqrt(sum(p.^2,2))-r;
mesh = mkmesh(p,t,porder,bndexpr,0,1,5,fd);


function t = tri(m,n,parity)

if parity==0
  t=[1,2,m+2; 1,m+2,m+1];
else
  t=[1,2,m+1; 2,m+2,m+1];
end

t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);

% Reorder triangles in Cartesian order

ix=[];
for i=1:n-1
  ix1=i+(n-1)*(0:m-2);
  ix2=ix1+(n-1)*(m-1);

  if parity==0
    ix12=[ix2;ix1];
  else
    ix12=[ix1;ix2];
  end

  ix=[ix,reshape(ix12,1,[])];
end

t=t(ix,:);

function h=href(p,hs)
b = 10;
c = 7;
pv1=[100,b;0,b;0,-c;100,-c;100,b];
inside1=inpolygon(p(:,1),p(:,2),pv1(:,1),pv1(:,2));
h=inf*ones(size(p,1),1);
h(inside1)=hs;

