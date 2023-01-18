function mesh = mkmesh_square_imex(porder)

% xmin = -10;
% xmax =  11;
% ymin = -10;
% ymax =  10;

xmin = 0;
xmax =  1;
ymin = 0;
ymax =  1;

% rectangle
%pv1 = [xmin ymin; xmax ymin; xmax ymax; xmin ymax];
m  = 3;
a  = 2;
s1 = loginc(linspace(0,0.5,m),a);
s2 = logdec(linspace(0.5,1,m),a);
s  = [s1 s2(2:end)]';
n  = length(s);
pv1 = [s zeros(n,1); ones(n-1,1) s(2:end); s(end-1:-1:1) ones(n-1,1); zeros(n-2,1) s(end-1:-1:2)];

H = 0.5; %dictates element size
[p,t]=polymesh({pv1},[1],[1,0],[H,1.75]);
[p,t] = fixmesh(p,t);

s1 = strcat('all(p(:,2)<',int2str(eval('ymin')),'+1e-6)');
s2 = strcat('all(p(:,1)>',int2str(eval('xmax')),'-1e-6)');
s3 = strcat('all(p(:,2)>',int2str(eval('ymax')),'-1e-6)');
s4 = strcat('all(p(:,1)<',int2str(eval('xmin')),'+1e-6)');
bndexpr={s1,s2,s3,s4};   

mesh = mkmesh(p,t,porder,bndexpr,0,1);

function [x,y] = linenodes(x1,y1,x2,y2,N)

xi = loginc(linspace(0,0.5,N),3);
xi = [xi 1-xi(end-1:-1:1)]';
x = x1 + xi*(x2-x1);
y = y1 + xi*(y2-y1);

function y=loginc(x,alpha)

a=min(x(:));b=max(x(:));
y = a + (b-a)*(exp(alpha*(x-a)/(b-a))-1)/(exp(alpha)-1);
