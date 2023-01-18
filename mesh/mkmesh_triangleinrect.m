function mesh = mkmesh_triangleinrect(porder)

% xmin = -10;
% xmax =  11;
% ymin = -10;
% ymax =  10;

xmin = -5;
xmax =  6;
ymin = -5;
ymax =  5;

% rectangle
pv1 = [xmin ymin; xmax ymin; xmax ymax; xmin ymax];

% triangle
N = 5;
v1 = [0 0];
v2 = [1 -sqrt(3)/3];
v3 = [1 sqrt(3)/3];
[x1,y1] = linenodes(v1(1),v1(2),v2(1),v2(2),N);
[x2,y2] = linenodes(v2(1),v2(2),v3(1),v3(2),N);
[x3,y3] = linenodes(v3(1),v3(2),v1(1),v1(2),N);
pv2 = [x1 y1; x2(2:end) y2(2:end); x3(2:end-1) y3(2:end-1)];

H = 0.8;
[p,t]=polymesh({pv1,pv2},[1,1],[1,0;0,1],[H,1.75]);
[p,t] = fixmesh(p,t);

s1 = strcat('all(p(:,2)<',int2str(eval('ymin')),'+1e-6)');
s2 = strcat('all(p(:,1)>',int2str(eval('xmax')),'-1e-6)');
s3 = strcat('all(p(:,2)>',int2str(eval('ymax')),'-1e-6)');
s4 = strcat('all(p(:,1)<',int2str(eval('xmin')),'+1e-6)');
bndexpr={s1,s2,s3,s4,'all(sum(p.^2,2)<3^2)'};   

mesh = mkmesh(p,t,porder,bndexpr,0,1);

function [x,y] = linenodes(x1,y1,x2,y2,N)

xi = loginc(linspace(0,0.5,N),3);
xi = [xi 1-xi(end-1:-1:1)]';
x = x1 + xi*(x2-x1);
y = y1 + xi*(y2-y1);

function y=loginc(x,alpha)

a=min(x(:));b=max(x(:));
y = a + (b-a)*(exp(alpha*(x-a)/(b-a))-1)/(exp(alpha)-1);
